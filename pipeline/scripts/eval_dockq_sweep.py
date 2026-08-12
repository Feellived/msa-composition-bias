#!/usr/bin/env python3
"""depth-sweep pose DockQ 채점 (pose vs native). 각 (target, model, rung)의 best DockQ + Neff80.
- 사슬은 chains.json 서열로 매칭(모델별 사슬ID 달라도 robust): 항원 사슬(들) + 항체 H/L.
- 병합: 항원 → chain A, 항체 H+L → chain B (양쪽 동일), DockQ CLI로 채점.
- RBD는 native를 chains.json의 crop 범위로 잘라 pose(크롭됨)와 맞춤.
경로: pose=$DATA/<model>/<target>/rung<r>/results/**.cif · native=targets/<target>/native.cif · Neff80=$DATA/ladders/<target>/<chain>/neff.tsv
사용(DockQ+biopython env): python eval_dockq_sweep.py [--models boltz protenix] [--out results/dockq_sweep.csv]
  ⭐ 이어달리기: --out CSV에 이미 있는 (target,model,rung)은 skip(재실행 값쌈). 강제 재채점=--rescore.
  ⭐ 생성과 동시 관찰: --watch 300  → 5분마다 새로 나온 pose만 채점·append (Ctrl-C 종료). protenix(GPU)와 병행(CPU) 안전.
"""
import argparse, csv, glob, json, os, re, subprocess, tempfile, time
from Bio.PDB import MMCIFParser, PDBParser, PDBIO
from Bio.PDB.Structure import Structure
from Bio.PDB.Model import Model as BModel
from Bio.PDB.Chain import Chain as BChain
from Bio.Align import PairwiseAligner
from Bio.Data.IUPACData import protein_letters_3to1

T3 = {k.upper(): v for k, v in protein_letters_3to1.items()}
_al = PairwiseAligner(); _al.mode = "global"
_al.match_score = 1; _al.mismatch_score = -1; _al.open_gap_score = -3; _al.extend_gap_score = -0.5
COLS = ["target", "group", "ab", "model", "rung", "neff80", "best_dockq", "n_pose"]

def load(p): return (MMCIFParser(QUIET=True) if p.endswith(".cif") else PDBParser(QUIET=True)).get_structure("x", p)[0]

def seqca(ch):
    s, rr = [], []
    for r in ch:
        if "CA" in r:
            aa = T3.get(r.get_resname().upper())
            if aa: s.append(aa); rr.append(r)
    return "".join(s), rr

def best(model, tgt_seq, exclude=()):
    """model 사슬 중 tgt_seq와 가장 잘 맞는 것 → (chain_id, score, residues)."""
    b = (None, -1e9, [])
    for ch in model:
        if ch.id in exclude: continue
        s, rr = seqca(ch)
        if len(s) < 5: continue
        sc = _al.score(s, tgt_seq) / max(len(s), len(tgt_seq), 1)
        if sc > b[1]: b = (ch.id, sc, rr)
    return b

def merged_pdb(ag_groups, ab_groups, out):
    s = Structure("m"); m = BModel(0); s.add(m); A = BChain("A"); Bc = BChain("B")
    i = 1
    for grp in ag_groups:
        for r in grp:
            if "CA" not in r: continue
            rr = r.copy(); rr.id = (" ", i, " "); A.add(rr); i += 1
    j = 1
    for grp in ab_groups:
        for r in grp:
            if "CA" not in r: continue
            rr = r.copy(); rr.id = (" ", j, " "); Bc.add(rr); j += 1
    m.add(A); m.add(Bc); io = PDBIO(); io.set_structure(s); io.save(out)

def dockq(model_pdb, native_pdb):
    try:
        r = subprocess.run(["DockQ", model_pdb, native_pdb], capture_output=True, text=True, timeout=300)
    except Exception:
        return None
    out = (r.stdout or "") + (r.stderr or "")
    vals = re.findall(r"\bDockQ\b[:\s]+([0-9]*\.[0-9]+)", out)
    return max(float(v) for v in vals) if vals else None

def _extract(model, cj, td, name, use_crop):
    ag_seqs = [(c["seq"], c.get("crop")) for c in cj["chains"] if c["role"] == "antigen"]
    ab_seqs = [c["seq"] for c in cj["chains"] if c["role"] in ("heavy", "light")]
    used = set(); ag = []; ab = []
    for seq, crop in ag_seqs:
        cid, _, rr = best(model, seq, exclude=used)
        if cid is None: return None
        used.add(cid)
        if use_crop and crop: rr = [r for r in rr if crop[0] <= r.id[1] <= crop[1]]
        ag.append(rr)
    for seq in ab_seqs:
        cid, _, rr = best(model, seq, exclude=used)
        if cid is None: continue
        used.add(cid); ab.append(rr)
    if not ab or not any(ag): return None
    out = os.path.join(td, f"{name}.pdb"); merged_pdb(ag, ab, out); return out

def _chain_res(model, cid):
    return seqca(model[cid])[1] if cid in model else None

def native_merged(cj, native_path, td):
    """native = prep이 기록한 정확한 원본 사슬(src_chains)로 추출 → 다중 copy에서 오짝(항원 copy1 + 항체 copy2) 방지."""
    if not os.path.exists(native_path): return None
    nm = load(native_path); src = cj.get("src_chains", {})
    crops = {c.get("src"): c.get("crop") for c in cj["chains"] if c["role"] == "antigen"}
    ag = []
    for cid in src.get("antigen", []):
        rr = _chain_res(nm, cid)
        if rr is None:                                   # 원본 사슬 못 찾으면 서열매칭 fallback
            seq = next((c["seq"] for c in cj["chains"] if c.get("src") == cid), None)
            if seq: rr = best(nm, seq)[2]
        if not rr: continue
        crop = crops.get(cid)
        if crop: rr = [r for r in rr if crop[0] <= r.id[1] <= crop[1]]
        ag.append(rr)
    ab = []
    for cid in list(src.get("H", [])) + list(src.get("L", [])):
        rr = _chain_res(nm, cid)
        if rr: ab.append(rr)
    if not any(ag) or not ab: return None
    out = os.path.join(td, "nat.pdb"); merged_pdb(ag, ab, out); return out

def pose_merged(cj, pose_path, td):
    return _extract(load(pose_path), cj, td, "pose", use_crop=False)    # pose는 서열매칭(모델별 사슬ID 다름) + 이미 crop된 서열

def neff_of(target, ladders):
    for tsv in sorted(glob.glob(os.path.join(ladders, target, "*", "neff.tsv"))):
        m = {}
        for line in open(tsv):
            if line.startswith("rung") or not line.strip(): continue
            p = line.split()
            if len(p) >= 3: m[int(p[0])] = float(p[2])
        return m
    return {}

def write_csv(out, rows_by_key):
    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLS); w.writeheader()
        for k in sorted(rows_by_key, key=lambda x: (x[1], x[0], x[2])):   # model,target,rung 순
            w.writerow(rows_by_key[k])

def score_pass(a, rows_by_key):
    """poses 있는데 아직 안 채점된 (target,model,rung)만 채점 → rows_by_key에 추가. 신규 건수 반환."""
    new = 0
    for r in csv.DictReader(open(a.list)):
        tgt = r["target"]; cjp = os.path.join(a.targets_dir, tgt, "chains.json")
        if not os.path.exists(cjp): continue
        pend = [(model, rung) for model in a.models for rung in range(a.rungs)
                if (tgt, model, rung) not in rows_by_key]
        if not pend: continue
        cj = json.load(open(cjp)); native = os.path.join(a.targets_dir, tgt, "native.cif")
        nmap = neff_of(tgt, os.path.join(a.data, "ladders"))
        with tempfile.TemporaryDirectory() as td:
            natm = native_merged(cj, native, td)
            if natm is None: print(f"{tgt}: native merge 실패"); continue
            for model, rung in pend:
                poses = glob.glob(os.path.join(a.data, model, tgt, f"rung{rung}", "results", "**", "*.cif"), recursive=True)
                if not poses: continue
                bq = None
                for pose in poses:
                    try:
                        pm = pose_merged(cj, pose, td)
                        if pm is None: continue
                        q = dockq(pm, natm)
                        if q is not None and (bq is None or q > bq): bq = q
                    except Exception: continue
                if bq is None: continue
                rows_by_key[(tgt, model, rung)] = dict(target=tgt, group=r["group"], ab=r["ab"], model=model,
                    rung=rung, neff80=nmap.get(rung, ""), best_dockq=round(bq, 3), n_pose=len(poses))
                new += 1
                print(f"  {tgt:14} {model:8} rung{rung} Neff80={nmap.get(rung,'?')} DockQ={bq:.3f} ({len(poses)} pose)")
        write_csv(a.out, rows_by_key)   # 타깃마다 flush → 중간에 끊겨도 보존, 실시간 관찰
    return new

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--list", default="sweep_targets.csv")
    ap.add_argument("--targets-dir", default="targets")
    ap.add_argument("--data", default=os.environ.get("DATA", "/mnt/data/admuser/msadepth"))
    ap.add_argument("--models", nargs="+", default=["boltz", "protenix"])
    ap.add_argument("--rungs", type=int, default=12)
    ap.add_argument("--out", default="results/dockq_sweep.csv")
    ap.add_argument("--watch", type=int, default=0, help="초. >0이면 그 간격으로 반복 채점(생성과 동시 관찰). Ctrl-C 종료")
    ap.add_argument("--rescore", action="store_true", help="--out 무시하고 전부 재채점")
    a = ap.parse_args()
    if subprocess.run(["which", "DockQ"], capture_output=True).returncode != 0:
        raise SystemExit("!! DockQ 없음 — DockQ 있는 env에서 실행 (pip install DockQ)")

    rows_by_key = {}
    if os.path.exists(a.out) and not a.rescore:
        for row in csv.DictReader(open(a.out)):
            rows_by_key[(row["target"], row["model"], int(row["rung"]))] = row
        print(f"기존 {len(rows_by_key)}건 로드 → 이미 채점된 건 skip (강제=--rescore)")

    while True:
        n = score_pass(a, rows_by_key)
        print(f"\n→ {a.out} (누적 {len(rows_by_key)}행, 이번 패스 신규 {n}). 성공기준 DockQ 0.23/0.49/0.80.")
        if a.watch <= 0: break
        print(f"[watch] {a.watch}s 대기 후 새 pose 확인... (Ctrl-C로 종료)")
        try: time.sleep(a.watch)
        except KeyboardInterrupt: print("\n중단."); break

if __name__ == "__main__":
    main()
