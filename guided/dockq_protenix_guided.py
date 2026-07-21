#!/usr/bin/env python3
"""guided Protenix(runs_guided_prot, HADDOCK 패치 pocket 제약) DockQ vs blind Protenix d143. 항체 H+L 병합.
   guided_best = 패치별 best 중 최고(= 클러스터 오라클 상한). Protenix 출력 = *_sample_*.cif.
   ⚠️ blind은 Protenix d143(runs_msad_143/<t>/out_protenix)이라 degraded — 진짜 배포비교는 server full(summary.csv)와.
   사용(DockQ 설치): python scripts/dockq_protenix_guided.py --targets "8XSI 9SBB 8SIS 9ML8 8SDF 8SIT 9ML9" """
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
    ap.add_argument("--guided-dir", default="runs_guided_prot")
    ap.add_argument("--blind-dir", default="runs_msad_143", help="Protenix blind(d143) 위치. 없으면 blind='-'")
    ap.add_argument("--out", default="results/dockq_protenix_guided.csv")
    a = ap.parse_args()
    print(f"{'target':7}{'class':11}{'blind_pd143':12}{'guided_best':12}{'diff':>8}  #patch")
    print("-" * 58)
    rows = []
    for t in a.targets.split():
        t = t.upper(); ni = B.native_info(t)
        if ni is None: print(f"{t:7} native 없음"); continue
        with tempfile.TemporaryDirectory() as td:
            natm = native_merged(t, td)
            if natm is None: print(f"{t:7} native merge 실패"); continue
            blind = sorted(glob.glob(f"{a.blind_dir}/{t.lower()}/out_protenix/**/*sample*.cif", recursive=True))
            pdirs = sorted(glob.glob(f"{a.guided_dir}/{t.lower()}/results_patch*"))
            bq_b, _ = best_dockq(blind, ni, natm, td)
            gbest = None; npatch = 0
            for pd in pdirs:
                gp = sorted(glob.glob(os.path.join(pd, "**", "*sample*.cif"), recursive=True))
                if not gp: continue
                npatch += 1
                q, _ = best_dockq(gp, ni, natm, td)
                if q is not None and (gbest is None or q > gbest): gbest = q
            bs = f"{bq_b:.3f}" if bq_b is not None else "-"
            gs = f"{gbest:.3f}" if gbest is not None else "-"
            dd = f"{gbest - bq_b:+.3f}" if (gbest is not None and bq_b is not None) else "-"
            print(f"{t:7}{B.EPICLASS.get(t,''):11}{bs:>10} {gs:>11} {dd:>7}  {npatch}")
            rows.append({"target": t, "class": B.EPICLASS.get(t, ""), "blind_pd143": bs,
                         "guided_best": gs, "diff": dd, "n_patch": npatch})
    os.makedirs("results", exist_ok=True)
    with open(a.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["target", "class", "blind_pd143", "guided_best", "diff", "n_patch"])
        w.writeheader(); w.writerows(rows)
    print(f"\n→ {a.out}")
    print("해석: guided_best(=오라클 over 패치) > blind → 'HADDOCK 패치가 Protenix를 near-native로 정교화'.")
    print("      ⚠️ blind=Protenix d143(degraded). 진짜 배포가치는 server full-MSA Protenix(summary.csv)와 비교.")
    print("      Protenix pocket은 soft라 Boltz(force:true)보다 steering 약할 수 있음.")


if __name__ == "__main__": main()
