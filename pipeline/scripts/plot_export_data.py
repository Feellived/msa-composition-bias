#!/usr/bin/env python3
"""[그림 데이터 추출] 노션 인수인계서 그림(A3·B1·B2·B3)용 compact 데이터를 stdout으로 뽑는다.

⚠️ 재랭커 아님 — 그림 재료만 만든다. stdlib only. 서버 results/ 에서 실행하고
   찍힌 출력을 그대로(4블록 전부) 붙여주면 로컬에서 그림을 그려 노션에 올린다.

  cd ~/projects/msa-composition-bias/pipeline && git pull
  python plot_export_data.py                 # results/pose_features.csv + results/crosscheck_depth.csv
  python plot_export_data.py --pf 다른.csv --cc 다른.csv

출력 4블록(각 [TAG]로 시작):
  [A3-HIST]   model,bin_lo,count                     # 모델별 DockQ 분포(0.1 구간)
  [A3-ORACLE] target,oracle,best_boltz,best_protenix # 타깃별 최고값: 전체(oracle) vs 모델별
  [B1]        crosscheck responsive 행(깊이-반응 요약)
  [B2B3]      target,model,group,rung,neff80,best_dq,max_recall,min_overrep,mean_recall,n
              (후보 복합체만 = crosscheck responsive ∪ 앵커 목록)
"""
import argparse, csv, math
from collections import defaultdict

CAND = ["8txu_HL", "9y0a_AB", "9azr_HL", "8wpy_AB", "8k3k_D", "8k46_I", "8y6a_CD", "8ulr_HL"]
B1COLS = ["target", "model", "dq_full", "dq_peak", "dq_peak_rung",
          "rec_full", "rec_peak", "rec_at_dqpk", "coincide", "responsive"]


def f(x):
    try:
        v = float(x)
        return None if math.isnan(v) else v
    except Exception:
        return None


def is_true(x):
    return str(x).strip().lower() in ("true", "1", "yes", "y", "t")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pf", default="results/pose_features.csv")
    ap.add_argument("--cc", default="results/crosscheck_depth.csv")
    a = ap.parse_args()
    pf = list(csv.DictReader(open(a.pf)))

    # ── [A3-HIST] 모델별 DockQ 히스토그램(0.1 구간) ──
    hist = defaultdict(lambda: [0] * 10)
    for r in pf:
        d = f(r.get("dockq"))
        if d is None:
            continue
        hist[r["model"]][min(9, max(0, int(d * 10)))] += 1
    print("[A3-HIST] model,bin_lo,count")
    for m in sorted(hist):
        for i, c in enumerate(hist[m]):
            print(f"{m},{i/10:.1f},{c}")

    # ── [A3-ORACLE] 타깃별 최고 DockQ: 전체 vs 모델별 ──
    bymt = defaultdict(list)
    byt = defaultdict(list)
    for r in pf:
        d = f(r.get("dockq"))
        if d is None:
            continue
        bymt[(r["target"], r["model"])].append(d)
        byt[r["target"]].append(d)
    print("[A3-ORACLE] target,oracle,best_boltz,best_protenix")
    for t in sorted(byt):
        orc = max(byt[t])
        bb = bymt.get((t, "boltz"))
        bp = bymt.get((t, "protenix"))
        print(f"{t},{orc:.3f},{max(bb):.3f}" if bb else f"{t},{orc:.3f},",
              end="")
        print(f",{max(bp):.3f}" if bp else ",")

    # ── [B1] crosscheck responsive 행 그대로 ──
    print("[B1] " + ",".join(B1COLS))
    cand = set(CAND)
    try:
        for r in csv.DictReader(open(a.cc)):
            if is_true(r.get("responsive", "")):
                cand.add(r["target"])
                print(",".join(str(r.get(k, "")) for k in B1COLS))
    except FileNotFoundError:
        print("(crosscheck_depth.csv 없음 — python analyze_crosscheck_depth.py 먼저 실행)")

    # ── [B2B3] 후보 복합체 rung 곡선(깊이별 best DockQ·recall·over-rep) ──
    grp = defaultdict(list)
    for r in pf:
        if r["target"] in cand:
            grp[(r["target"], r["model"], r.get("group", ""), int(float(r["rung"])))].append(r)
    print("[B2B3] target,model,group,rung,neff80,best_dq,max_recall,min_overrep,mean_recall,n")
    for (t, m, g, rg), ps in sorted(grp.items(), key=lambda x: (x[0][0], x[0][1], x[0][3])):
        dqs = [f(p.get("dockq")) for p in ps]; dqs = [x for x in dqs if x is not None]
        recs = [f(p.get("recall")) for p in ps]; recs = [x for x in recs if x is not None]
        ovs = [f(p.get("overrep")) for p in ps]; ovs = [x for x in ovs if x is not None]
        neff = f(ps[0].get("neff80"))
        bd = f"{max(dqs):.3f}" if dqs else ""
        mr = f"{max(recs):.3f}" if recs else ""
        mo = f"{min(ovs):.3f}" if ovs else ""
        mnr = f"{sum(recs)/len(recs):.3f}" if recs else ""
        nf = f"{neff:.1f}" if neff is not None else ""
        print(f"{t},{m},{g},{rg},{nf},{bd},{mr},{mo},{mnr},{len(ps)}")


if __name__ == "__main__":
    main()
