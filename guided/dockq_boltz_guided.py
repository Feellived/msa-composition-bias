#!/usr/bin/env python3
"""guided Boltz(runs_guided, HADDOCK patch pocket 제약) DockQ vs blind Boltz(runs_msad_143). 항체 H+L 병합.
   guided_best = 패치별 best 중 최고(=oracle over 클러스터 상한). '올바른 패치를 주면 co-folder가 정교화하나' 테스트.
   사용(boltz env, DockQ 설치): python scripts/dockq_boltz_guided.py --targets "8XSI 9SBB 8SIS 9ML8 8SDF 8SIT 9ML9" """
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
    ap.add_argument("--guided-dir", default="runs_guided", help="guided 결과 폴더(변형 실험이면 예: runs_guided_top20u)")
    ap.add_argument("--out", default="results/dockq_boltz_guided.csv")
    a = ap.parse_args()
    print(f"{'target':7}{'class':11}{'blind_d143':11}{'guided_best':12}{'diff':>8}  #patch")
    print("-" * 56)
    rows = []
    for t in a.targets.split():
        t = t.upper(); ni = B.native_info(t)
        if ni is None: print(f"{t:7} native 없음"); continue
        with tempfile.TemporaryDirectory() as td:
            natm = native_merged(t, td)
            if natm is None: print(f"{t:7} native merge 실패"); continue
            blind = sorted(glob.glob(f"runs_msad_143/{t.lower()}/results_boltz/**/*_model_*.cif", recursive=True))
            pdirs = sorted(glob.glob(f"{a.guided_dir}/{t.lower()}/results_patch*"))
            bq_b, _ = best_dockq(blind, ni, natm, td)
            gbest = None; npatch = 0
            for pd in pdirs:
                gp = sorted(glob.glob(os.path.join(pd, "**", "*_model_*.cif"), recursive=True))
                if not gp: continue
                npatch += 1
                q, _ = best_dockq(gp, ni, natm, td)
                if q is not None and (gbest is None or q > gbest): gbest = q
            bs = f"{bq_b:.3f}" if bq_b is not None else "-"
            gs = f"{gbest:.3f}" if gbest is not None else "-"
            dd = f"{gbest - bq_b:+.3f}" if (gbest is not None and bq_b is not None) else "-"
            print(f"{t:7}{B.EPICLASS.get(t,''):11}{bs:>10} {gs:>11} {dd:>7}  {npatch}")
            rows.append({"target": t, "class": B.EPICLASS.get(t, ""), "blind_d143": bs,
                         "guided_best": gs, "diff": dd, "n_patch": npatch})
    os.makedirs("results", exist_ok=True)
    with open(a.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["target", "class", "blind_d143", "guided_best", "diff", "n_patch"])
        w.writeheader(); w.writerows(rows)
    print(f"\n→ {a.out}")
    print("해석: guided_best(=oracle over 클러스터) > blind_d143 이면 'HADDOCK 패치 힌트가 co-folder를 near-native로 정교화'.")
    print("      실제 배포는 '어느 패치인지' 자동선택 문제가 남음(재랭커 몫).")

if __name__ == "__main__": main()
