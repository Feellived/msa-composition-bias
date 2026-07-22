#!/usr/bin/env python3
"""pose별 에피토프 recall — 예측 항체가 '진짜 에피토프 잔기'를 얼마나 맞게 접촉하나.
DockQ(스칼라, fold+배치 혼합)와 달리 '위치'만 본다 → MSA 깊이에 따라 예측 에피토프가 이동하나 직접 측정(2a 위치편향).
- true epitope = native에서 항체 heavy-atom 5Å 내 항원 잔기 → 참조 항원서열 위치로 표현
- pred epitope = pose에서 동일 계산(항원/항체는 서열매칭으로 식별, 모델별 사슬ID 달라도 robust)
- 참조-서열 위치 공간(항원 seq에 정렬)으로 native·pose 번호 불일치 해결. 다사슬 항원은 (src사슬, 위치)로 네임스페이스.
- recall=|pred∩true|/|true|, precision=|∩|/|pred|, F1. (target,rung)별 best-recall pose + mean.
경로: pose=$DATA/<model>/<t>/rung<r>/results/**.cif · native=targets/<t>/native.cif · chains.json=targets/<t>/chains.json
사용(biopython env): python epitope_recall.py [--models boltz protenix] [--rungs 12] [--cutoff 5.0] [--out results/epitope_recall.csv]
⚠️ 서버서 1복합체 스모크 먼저(ONLY 없음 → --list로 한 줄 csv 주거나 첫 타깃 로그 확인).
"""
import argparse, csv, glob, json, os
from Bio.PDB import MMCIFParser, PDBParser, NeighborSearch
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
    """seq 잔기 index → ref 서열 위치 (정렬)."""
    if not seq or not ref: return {}
    try: aln = _al.align(seq, ref)[0]
    except Exception: return {}
    m = {}
    for (a0, a1), (b0, b1) in zip(*aln.aligned):
        for k in range(a1 - a0): m[a0 + k] = b0 + k
    return m

def epitope(ag_chains, ab_residues, cutoff):
    """ag_chains: [(chain_key, residues, ref_seq)]. 반환 = 접촉 항원잔기의 (chain_key, ref위치) 집합."""
    ab_atoms = [a for r in ab_residues for a in r if a.element != "H"]
    if not ab_atoms: return set()
    ns = NeighborSearch(ab_atoms)
    epi = set()
    for ckey, residues, ref in ag_chains:
        seq = "".join(T3.get(r.get_resname().upper(), "X") for r in residues)
        pm = posmap(seq, ref)
        for i, res in enumerate(residues):
            if i not in pm: continue
            for atom in res:
                if atom.element == "H": continue
                if ns.search(atom.coord, cutoff): epi.add((ckey, pm[i])); break
    return epi

def native_true(cj, native_path, cutoff):
    if not os.path.exists(native_path): return None
    nm = load(native_path)
    ag = []
    for c in cj["chains"]:
        if c["role"] != "antigen": continue
        src = c.get("src")
        if src in nm: ag.append((src, seqca(nm[src])[1], c["seq"]))
    ab_src = [c.get("src") for c in cj["chains"] if c["role"] in ("heavy", "light")]
    ab = [r for cid in ab_src if cid in nm for r in seqca(nm[cid])[1]]
    if not ag or not ab: return None
    return epitope(ag, ab, cutoff), sum(1 for c in cj["chains"] if c["role"] == "antigen")

def pose_pred(cj, pose_path, cutoff):
    m = load(pose_path); used = set(); ag = []
    for c in cj["chains"]:
        if c["role"] != "antigen": continue
        cid, _, rr = best(m, c["seq"], exclude=used)
        if cid is None: continue
        used.add(cid); ag.append((c.get("src"), rr, c["seq"]))
    ab = []
    for c in cj["chains"]:
        if c["role"] in ("heavy", "light"):
            cid, _, rr = best(m, c["seq"], exclude=used)
            if cid: used.add(cid); ab.extend(rr)
    if not ag or not ab: return None
    return epitope(ag, ab, cutoff)

def prf(pred, true):
    inter = len(pred & true)
    rec = inter / len(true) if true else 0.0
    pre = inter / len(pred) if pred else 0.0
    f1 = 2 * rec * pre / (rec + pre) if (rec + pre) else 0.0
    return rec, pre, f1

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
        true_epi, n_agch = tr
        nmap = neff_of(tgt, os.path.join(a.data, "ladders"))
        for model in a.models:
            for rung in range(a.rungs):
                poses = glob.glob(os.path.join(a.data, model, tgt, f"rung{rung}", "results", "**", "*.cif"), recursive=True)
                if not poses: continue
                recs = []
                for pose in poses:
                    try:
                        pred = pose_pred(cj, pose, a.cutoff)
                        if pred is None: continue
                        recs.append(prf(pred, true_epi))
                    except Exception: continue
                if not recs: continue
                bestp = max(recs, key=lambda x: x[0])
                mean_rec = sum(x[0] for x in recs) / len(recs)
                rows.append(dict(target=tgt, group=r["group"], ab=r["ab"], model=model, rung=rung,
                                 neff80=nmap.get(rung, ""), n_true=len(true_epi), n_pose=len(recs),
                                 best_recall=round(bestp[0], 3), best_prec=round(bestp[1], 3), best_f1=round(bestp[2], 3),
                                 mean_recall=round(mean_rec, 3)))
                print(f"  {tgt:14} {model:8} rung{rung} Neff80={nmap.get(rung,'?')} recall(best)={bestp[0]:.3f} mean={mean_rec:.3f} (n_true={len(true_epi)}, {len(recs)} pose)")
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    with open(a.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["target", "group", "ab", "model", "rung", "neff80", "n_true", "n_pose",
                                          "best_recall", "best_prec", "best_f1", "mean_recall"])
        w.writeheader(); w.writerows(rows)
    print(f"\n→ {a.out} ({len(rows)}행). 다음: recall vs depth 곡선(DockQ와 나란히) — 깊이 줄일 때 에피토프 이동 여부.")

if __name__ == "__main__":
    main()
