#!/usr/bin/env python3
"""pose별 에피토프 예측 품질 — 예측 항체가 '진짜 에피토프 잔기'를 얼마나 맞게/특이적으로 접촉하나.
DockQ(스칼라, fold+배치 혼합)와 달리 '위치'만 본다 → MSA 깊이에 따라 예측 에피토프가 이동하나 직접 측정(2a 위치편향).
- true epitope = native에서 항체 heavy-atom 5Å 내 항원 잔기 → 참조 항원서열 위치로 표현
- pred epitope = pose에서 동일. 항원/항체는 서열매칭 식별(모델별 사슬ID 달라도 robust), 참조-서열 위치로 native/pose 정렬.
- 지표(불충분한 recall 단독 회피): recall(best/mean/min) · precision · F1 · MCC(불균형 강건) · AUPRC(threshold-free, 잔기별 항체거리 점수).
  ※ recall만 = 표면 전체 덮으면 100%(특이성0) → precision/MCC 필요. accuracy는 불균형(에피토프=소수)이라 오도 → MCC가 정도. AUPRC>AUROC(희소양성).
경로: pose=$DATA/<model>/<t>/rung<r>/results/**.cif · native=targets/<t>/native.cif · chains.json=targets/<t>/chains.json
사용(biopython+scipy env): python epitope_recall.py [--models boltz protenix] [--rungs 12] [--cutoff 5.0] [--out results/epitope_recall.csv]
"""
import argparse, csv, glob, json, math, os
import numpy as np
from scipy.spatial import cKDTree
from Bio.PDB import MMCIFParser, PDBParser
from Bio.Align import PairwiseAligner
from Bio.Data.IUPACData import protein_letters_3to1

T3 = {k.upper(): v for k, v in protein_letters_3to1.items()}
_al = PairwiseAligner(); _al.mode = "global"
_al.match_score = 1; _al.mismatch_score = -1; _al.open_gap_score = -3; _al.extend_gap_score = -0.5

def load(p): return (MMCIFParser(QUIET=True) if p.endswith(".cif") else PDBParser(QUIET=True)).get_structure("x", p)[0]

def seqca(ch):
    s, rr = [], []
    for r in ch:
        if "CA" in r:
            aa = T3.get(r.get_resname().upper())
            if aa: s.append(aa); rr.append(r)
    return "".join(s), rr

def best(model, seq, exclude=()):
    b = (None, -1e9, [])
    for ch in model:
        if ch.id in exclude: continue
        s, rr = seqca(ch)
        if len(s) < 5: continue
        sc = _al.score(s, seq) / max(len(s), len(seq), 1)
        if sc > b[1]: b = (ch.id, sc, rr)
    return b

def posmap(seq, ref):
    if not seq or not ref: return {}
    try: aln = _al.align(seq, ref)[0]
    except Exception: return {}
    m = {}
    for (a0, a1), (b0, b1) in zip(*aln.aligned):
        for k in range(a1 - a0): m[a0 + k] = b0 + k
    return m

def scored_epitope(ag_chains, ab_residues, cutoff):
    """반환: (접촉잔기 집합, {ref키: 항체까지 최소거리}). ag_chains=[(chain_key, residues, ref_seq)]."""
    ab_coords = np.array([a.coord for r in ab_residues for a in r if a.element != "H"], dtype=float)
    if len(ab_coords) == 0: return set(), {}
    tree = cKDTree(ab_coords)
    contacts, dist = set(), {}
    for ckey, residues, ref in ag_chains:
        seq = "".join(T3.get(r.get_resname().upper(), "X") for r in residues)
        pm = posmap(seq, ref)
        for i, res in enumerate(residues):
            if i not in pm: continue
            coords = np.array([a.coord for a in res if a.element != "H"], dtype=float)
            if len(coords) == 0: continue
            md = float(tree.query(coords, k=1)[0].min())
            key = (ckey, pm[i]); dist[key] = md
            if md <= cutoff: contacts.add(key)
    return contacts, dist

def auprc(scores, labels):
    """step-wise average precision. scores 큰 게 양성 추정. labels 0/1."""
    P = sum(labels)
    if P == 0: return float("nan")
    order = sorted(range(len(scores)), key=lambda i: -scores[i])
    tp = fp = 0; ap = 0.0; prev_rec = 0.0
    for i in order:
        if labels[i]: tp += 1
        else: fp += 1
        rec = tp / P; pre = tp / (tp + fp)
        ap += pre * (rec - prev_rec); prev_rec = rec
    return ap

def metrics(pred, dist, true, n_ag):
    tp = len(pred & true); fp = len(pred - true); fn = len(true - pred); tn = max(0, n_ag - tp - fp - fn)
    rec = tp / len(true) if true else 0.0
    pre = tp / len(pred) if pred else 0.0
    f1 = 2 * rec * pre / (rec + pre) if (rec + pre) else 0.0
    den = math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    mcc = (tp * tn - fp * fn) / den if den else 0.0
    keys = list(dist.keys())
    ap = auprc([-dist[k] for k in keys], [1 if k in true else 0 for k in keys]) if keys else float("nan")
    return dict(recall=rec, prec=pre, f1=f1, mcc=mcc, auprc=ap)

def antigen_refs(cj):  return [c["seq"] for c in cj["chains"] if c["role"] == "antigen"]
def antibody_refs(cj): return [c["seq"] for c in cj["chains"] if c["role"] in ("heavy", "light")]

def native_true(cj, native_path, cutoff):
    if not os.path.exists(native_path): return None
    nm = load(native_path); src = cj.get("src_chains", {})
    ag_ids = [str(x) for x in src.get("antigen", [])]
    ag = []
    for i, ref in enumerate(antigen_refs(cj)):
        cid = ag_ids[i] if i < len(ag_ids) else None
        rr = seqca(nm[cid])[1] if (cid and cid in nm) else best(nm, ref)[2]
        if rr: ag.append((i, rr, ref))
    ab_ids = [str(x) for x in list(src.get("H", [])) + list(src.get("L", []))]
    ab = [r for cid in ab_ids if cid in nm for r in seqca(nm[cid])[1]]
    if not ab:
        used = set()
        for ref in antibody_refs(cj):
            cid, _, rr = best(nm, ref, exclude=used)
            if cid: used.add(cid); ab.extend(rr)
    if not ag or not ab: return None
    true, _ = scored_epitope(ag, ab, cutoff)
    n_ag = sum(len(ref) for _, _, ref in ag)          # 총 항원 잔기(참조 서열 길이 합)
    return true, n_ag

def pose_metrics(cj, pose_path, cutoff, true, n_ag):
    m = load(pose_path); used = set(); ag = []
    for i, ref in enumerate(antigen_refs(cj)):
        cid, _, rr = best(m, ref, exclude=used)
        if cid is None: continue
        used.add(cid); ag.append((i, rr, ref))
    ab = []
    for ref in antibody_refs(cj):
        cid, _, rr = best(m, ref, exclude=used)
        if cid: used.add(cid); ab.extend(rr)
    if not ag or not ab: return None
    pred, dist = scored_epitope(ag, ab, cutoff)
    return metrics(pred, dist, true, n_ag)

def neff_of(target, ladders):
    for tsv in sorted(glob.glob(os.path.join(ladders, target, "*", "neff.tsv"))):
        m = {}
        for line in open(tsv):
            if line.startswith("rung") or not line.strip(): continue
            p = line.split()
            if len(p) >= 3: m[int(p[0])] = float(p[2])
        return m
    return {}

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--list", default="sweep_targets.csv")
    ap.add_argument("--targets-dir", default="targets")
    ap.add_argument("--data", default=os.environ.get("DATA", "/mnt/data/admuser/msadepth"))
    ap.add_argument("--models", nargs="+", default=["boltz", "protenix"])
    ap.add_argument("--rungs", type=int, default=12)
    ap.add_argument("--cutoff", type=float, default=5.0)
    ap.add_argument("--out", default="results/epitope_recall.csv")
    a = ap.parse_args()
    rows = []
    for r in csv.DictReader(open(a.list)):
        tgt = r["target"]; cjp = os.path.join(a.targets_dir, tgt, "chains.json")
        native = os.path.join(a.targets_dir, tgt, "native.cif")
        if not os.path.exists(cjp): continue
        cj = json.load(open(cjp))
        tr = native_true(cj, native, a.cutoff)
        if tr is None: print(f"{tgt}: native epitope 실패"); continue
        true, n_ag = tr
        nmap = neff_of(tgt, os.path.join(a.data, "ladders"))
        for model in a.models:
            for rung in range(a.rungs):
                poses = glob.glob(os.path.join(a.data, model, tgt, f"rung{rung}", "results", "**", "*.cif"), recursive=True)
                if not poses: continue
                ms = []
                for pose in poses:
                    try:
                        mm = pose_metrics(cj, pose, a.cutoff, true, n_ag)
                        if mm: ms.append(mm)
                    except Exception: continue
                if not ms: continue
                bp = max(ms, key=lambda x: x["recall"])            # best-recall pose
                recs = [x["recall"] for x in ms]
                aps = [x["auprc"] for x in ms if not math.isnan(x["auprc"])]
                rows.append(dict(target=tgt, group=r["group"], ab=r["ab"], model=model, rung=rung,
                                 neff80=nmap.get(rung, ""), n_true=len(true), n_ag=n_ag, n_pose=len(ms),
                                 best_recall=round(max(recs), 3), mean_recall=round(sum(recs) / len(recs), 3),
                                 min_recall=round(min(recs), 3),
                                 best_prec=round(bp["prec"], 3), best_f1=round(bp["f1"], 3), best_mcc=round(bp["mcc"], 3),
                                 best_auprc=round(bp["auprc"], 3) if not math.isnan(bp["auprc"]) else "",
                                 mean_auprc=round(sum(aps) / len(aps), 3) if aps else ""))
                print(f"  {tgt:14} {model:8} r{rung:<2} Neff80={nmap.get(rung,'?')} recall best/mean/min={max(recs):.2f}/{sum(recs)/len(recs):.2f}/{min(recs):.2f} F1={bp['f1']:.2f} MCC={bp['mcc']:.2f} AUPRC={bp['auprc']:.2f} (n_true={len(true)})")
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    cols = ["target", "group", "ab", "model", "rung", "neff80", "n_true", "n_ag", "n_pose",
            "best_recall", "mean_recall", "min_recall", "best_prec", "best_f1", "best_mcc", "best_auprc", "mean_auprc"]
    with open(a.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols); w.writeheader(); w.writerows(rows)
    print(f"\n→ {a.out} ({len(rows)}행). 지표: recall(best/mean/min)·precision·F1·MCC·AUPRC. 다음: 지표 vs depth 곡선(DockQ와 나란히).")

if __name__ == "__main__":
    main()
