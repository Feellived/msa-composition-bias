#!/usr/bin/env python3
"""guided tFold(runs_guided_tfold, HADDOCK ICF) DockQ vs blind tFold(full-MSA, ungued). 항체 H+L 병합.
guided_best = patch별 best 중 최고(= 클러스터 오라클 상한). tFold 출력 = *.pdb.
blind = runs_rbd/<t>/out_tfold/*.pdb (full server MSA, ungued cofolder).
사용(DockQ 설치·torch+biopython env): python scripts/dockq_tfold_guided.py --targets "8XSI 9SBB 8SIS 9ML8 8SDF 8SIT 9ML9" """
import argparse, csv, glob, os, tempfile
import build_msafree_summary as B
from haddock_dockq import write_merged, dockq, native_merged


def best_dockq(poses, ni, natm, td):
    bq = None; n = 0
    for pose in poses:
        try:
            pm = B.load(pose)
            pid, _, ps, pca, _ = B.best(pm, ni["ag_seq"])
            abids = [ch.id for ch in pm if ch.id != pid]
            if not abids: continue
            mp = os.path.join(td, "m.pdb"); write_merged(pm, pid, abids, mp)
            q = dockq(mp, natm); n += 1
            if q is not None and (bq is None or q > bq): bq = q
        except Exception: continue
    return bq, n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--targets", default="8XSI 9SBB 8SIS 9ML8 8SDF 8SIT 9ML9")
    ap.add_argument("--guided-dir", default="runs_guided_tfold")
    ap.add_argument("--blind-dir", default="runs_rbd", help="ungued full-MSA tFold(out_tfold) 위치")
    ap.add_argument("--out", default="results/dockq_tfold_guided.csv")
    a = ap.parse_args()
    print(f"{'target':7}{'class':11}{'blind_full':11}{'guided_best':12}{'diff':>8}  #patch")
    print("-" * 58)
    rows = []
    for t in a.targets.split():
        t = t.upper(); ni = B.native_info(t)
        if ni is None: print(f"{t:7} native 없음"); continue
        with tempfile.TemporaryDirectory() as td:
            natm = native_merged(t, td)
            if natm is None: print(f"{t:7} native merge 실패"); continue
            blind = sorted(glob.glob(f"{a.blind_dir}/{t.lower()}/out_tfold/*.pdb"))
            pdirs = sorted(glob.glob(f"{a.guided_dir}/{t.lower()}/out_tfold_patch*"))
            bq_b, _ = best_dockq(blind, ni, natm, td)
            gbest = None; npatch = 0
            for pd in pdirs:
                gp = sorted(glob.glob(os.path.join(pd, "*.pdb")))
                if not gp: continue
                npatch += 1
                q, _ = best_dockq(gp, ni, natm, td)
                if q is not None and (gbest is None or q > gbest): gbest = q
            bs = f"{bq_b:.3f}" if bq_b is not None else "-"
            gs = f"{gbest:.3f}" if gbest is not None else "-"
            dd = f"{gbest - bq_b:+.3f}" if (gbest is not None and bq_b is not None) else "-"
            print(f"{t:7}{B.EPICLASS.get(t,''):11}{bs:>10} {gs:>11} {dd:>7}  {npatch}")
            rows.append({"target": t, "class": B.EPICLASS.get(t, ""), "blind_full": bs,
                         "guided_best": gs, "diff": dd, "n_patch": npatch})
    os.makedirs("results", exist_ok=True)
    with open(a.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["target", "class", "blind_full", "guided_best", "diff", "n_patch"])
        w.writeheader(); w.writerows(rows)
    print(f"\n→ {a.out}")
    print("해석: guided_best(=오라클 over 패치) > blind_full → 'HADDOCK ICF가 tFold를 near-native로 정교화'.")
    print("      ⚠️ tFold는 항체 MSA-free지만 항원 depth 편향 있음 — guided ICF가 그 편향을 계면 지정으로 눌러주나 확인.")


if __name__ == "__main__": main()
