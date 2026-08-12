#!/usr/bin/env python3
"""[재분석 #1] 순열 null — Boltz 'rescue'가 진짜 깊이효과인가 best-of-N 착시인가.
   문헌(Ojala&Garriga 2010 라벨셔플 null; Porter et al. 2025 matched-depth 대조).
가설: rescue = "full-MSA서 실패 → 깊이 줄이면 회복". null = "깊이 무효 → full은 그냥 랜덤 rung".
방법(코스: rung별 best_recall만 사용; per-pose는 Step0 후 정밀판): 타깃 안 rung 값을 exchangeable로 보고,
   'full(rung0)'을 랜덤 rung으로 바꿔 같은 rescue/gain 재계산 10,000×. 관측이 null 밴드 밖이면 = full이 특별히 나쁨(진짜).
   ⚠️ rmax(=max over rungs)는 라벨불변 → 이 test는 'full이 특별히 낮은가'를 직접 검정.
사용: python analyze_perm_null.py [--recall results/epitope_recall.csv] [--perms 10000] [--fail 0.3] [--succ 0.5]
"""
import argparse, csv, math
from collections import defaultdict
import numpy as np

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--recall", default="epitope_recall.csv")
    ap.add_argument("--perms", type=int, default=10000)
    ap.add_argument("--fail", type=float, default=0.3, help="full 실패 기준")
    ap.add_argument("--succ", type=float, default=0.5, help="회복 성공 기준")
    a = ap.parse_args()
    rng = np.random.default_rng(0)

    rec = defaultdict(list)
    for r in csv.DictReader(open(a.recall)):
        rec[r["target"]].append(r)
    T = []
    for t, rows in rec.items():
        rows = sorted(rows, key=lambda z: int(z["rung"]))
        neffs = set(round(float(z["neff80"]), 1) for z in rows)
        if len(neffs) < 4:            # 깊이 range 없는 타깃 제외(작은 MSA)
            continue
        T.append(dict(t=t, grp=rows[0]["group"],
                      br=np.array([float(z["best_recall"]) for z in rows]),
                      mr=np.array([float(z["mean_recall"]) for z in rows])))
    print(f"깊이-range 타깃 {len(T)}개 (rung당 best_recall = pose 5개 max)")

    def observed(key):
        r0 = np.array([d[key][0] for d in T]); rmax = np.array([d[key].max() for d in T])
        gain = rmax - r0
        rescue = ((r0 < a.fail) & (rmax >= a.succ)).sum()
        rank0 = np.array([(d[key] < d[key][0]).sum() / (len(d[key]) - 1) for d in T])  # full보다 낮은 rung 비율(0=full이 최저)
        return gain.mean(), int(rescue), rank0.mean()

    def perm_once(key):
        r0p = np.array([d[key][rng.integers(len(d[key]))] for d in T])
        rmax = np.array([d[key].max() for d in T])
        gain = (rmax - r0p).mean()
        rescue = ((r0p < a.fail) & (rmax >= a.succ)).sum()
        return gain, rescue

    for key, name in [("br", "best_recall (max-of-5)"), ("mr", "mean_recall (표집통제)")]:
        og, orc, orank = observed(key)
        ng = np.empty(a.perms); nr = np.empty(a.perms)
        for i in range(a.perms):
            ng[i], nr[i] = perm_once(key)
        p_gain = (ng >= og).mean(); p_res = (nr >= orc).mean()
        print(f"\n=== {name} ===")
        print(f"  관측 mean-gain(rmax−full) = {og:.3f}   null 평균 {ng.mean():.3f} (95% {np.percentile(ng,2.5):.3f}~{np.percentile(ng,97.5):.3f})  → p={p_gain:.3f}")
        print(f"  관측 rescue수(full<{a.fail} & rmax≥{a.succ}) = {orc}   null 평균 {nr.mean():.1f} (95% {np.percentile(nr,2.5):.0f}~{np.percentile(nr,97.5):.0f})  → p={p_res:.3f}")
        print(f"  full이 특별히 낮은가: full보다 낮은 rung 평균비율 {orank:.2f} (0.5=랜덤, <0.5=full이 낮은쪽)")
    print("\n해석: best_recall p<0.05 = full이 특별히 나쁨(rescue 진짜 방향) / p>0.05 = full은 랜덤 rung = best-of-N 착시.")
    print("      mean_recall이 best보다 p 커지면 = 신호가 max-표집에 의존(착시 강함).")

if __name__ == "__main__":
    main()
