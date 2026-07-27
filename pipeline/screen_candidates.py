#!/usr/bin/env python3
"""[후보 재선별] "full MSA는 실패하는데 얕은 깊이 여러 칸에서 연속 성공"하는 복합체 찾기.

왜 기준을 바꿨나(2026-07-27): 예전 선별은 사다리 12칸 중 **가장 잘 나온 한 칸**을 골랐다.
그런데 한 칸은 추첨 1회이고, 실행 잡음이 커서(같은 입력 재실행 0.588→0.011) 한 칸만 튄 것은
운일 수 있다. 실제로:
  · 8ulr = rung2·3·4 **3칸 연속** 성공 → 재현됨
  · 9azr = rung4 **1칸만** 성공      → 재현 실패(0/40)
→ **연속 칸 수**를 기준으로 삼으면 운으로 튄 것을 거를 수 있다.

⚠️ boltz 데이터는 a3m 오염으로 무효 → 기본값이 protenix만 본다(--model 로 변경 가능).

사용(stdlib only):
  python screen_candidates.py                          # results/pose_features.csv
  python screen_candidates.py --succ 0.49 --fail 0.23  # 성공/‘full 실패’ 문턱
  python screen_candidates.py --model boltz            # boltz 재실행 후에 쓸 것
"""
import argparse, csv, math
from collections import defaultdict


def f(x):
    try:
        v = float(x)
        return None if math.isnan(v) else v
    except Exception:
        return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", default="results/pose_features.csv")
    ap.add_argument("--model", default="protenix")
    ap.add_argument("--succ", type=float, default=0.49, help="성공 문턱")
    ap.add_argument("--fail", type=float, default=0.23, help="full MSA가 '실패'라고 볼 문턱")
    ap.add_argument("--out", default="results/screen_candidates.csv")
    a = ap.parse_args()

    best = defaultdict(dict)   # target -> rung -> best dockq
    neff = defaultdict(dict)
    for r in csv.DictReader(open(a.csv)):
        if r["model"] != a.model:
            continue
        d = f(r.get("dockq"))
        if d is None:
            continue
        t, rg = r["target"], int(float(r["rung"]))
        if d > best[t].get(rg, -1):
            best[t][rg] = d
        n = f(r.get("neff80"))
        if n is not None:
            neff[t][rg] = n

    rows = []
    for t, bk in best.items():
        if 0 not in bk:
            continue
        full = bk[0]
        rgs = sorted(bk)
        # 얕은 칸(rung>0) 중 성공한 칸들의 '연속 구간' 최대 길이
        runlen = cur = 0
        streak, curstreak = [], []
        for rg in rgs:
            if rg == 0:
                continue
            if bk[rg] >= a.succ:
                cur += 1; curstreak.append(rg)
                if cur > runlen:
                    runlen, streak = cur, list(curstreak)
            else:
                cur = 0; curstreak = []
        n_succ = sum(1 for rg in rgs if rg > 0 and bk[rg] >= a.succ)
        peak_rg = max((rg for rg in rgs if rg > 0), key=lambda x: bk[x], default=None)
        rows.append(dict(target=t, full=round(full, 3),
                         full_fails=int(full < a.fail),
                         n_succ_rungs=n_succ, max_streak=runlen,
                         streak_rungs="|".join(map(str, streak)),
                         peak_rung=peak_rg, peak=round(bk[peak_rg], 3) if peak_rg is not None else None,
                         peak_neff=round(neff[t].get(peak_rg, float("nan")), 1) if peak_rg is not None else None))

    rows.sort(key=lambda r: (-r["max_streak"], -r["n_succ_rungs"], -(r["peak"] or 0)))
    print(f"모델 {a.model} · 성공≥{a.succ} · full 실패 기준 <{a.fail}\n")
    print(f"{'target':12}{'full':>7}{'실패?':>6}{'성공칸':>7}{'최대연속':>9}{'연속칸':>12}"
          f"{'peak':>7}{'@rung':>7}{'neff':>8}")
    print("-" * 82)
    strong = []
    for r in rows:
        if r["n_succ_rungs"] == 0 and r["max_streak"] == 0:
            continue
        mark = ""
        if r["full_fails"] and r["max_streak"] >= 2:
            mark = "  ★강후보"; strong.append(r)
        elif r["full_fails"] and r["n_succ_rungs"] >= 1:
            mark = "  (약함: 한 칸만)"
        print(f"{r['target']:12}{r['full']:>7.3f}{('예' if r['full_fails'] else '아니오'):>6}"
              f"{r['n_succ_rungs']:>7}{r['max_streak']:>9}{r['streak_rungs'] or '-':>12}"
              f"{(r['peak'] or 0):>7.3f}{r['peak_rung']:>7}{(r['peak_neff'] or 0):>8.1f}{mark}")

    print(f"\n검사한 타깃 {len(rows)}개 · ★강후보(full 실패 + 연속 2칸 이상) = {len(strong)}개")
    if strong:
        print("   " + ", ".join(r["target"] for r in strong))
        print("   → 이들에 comp_x_reps.sh (조성×반복)를 돌리면 8ulr과 같은 검정을 할 수 있다.")
    else:
        print("   (없음 — 8ulr 외에 연속 성공하는 타깃이 없다는 뜻)")
    print("\n⚠️ '한 칸만' 성공은 실행 운일 가능성이 큼(9azr가 그 경우였고 재현 실패).")

    import os
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    if rows:
        with open(a.out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            w.writeheader(); w.writerows(rows)
        print(f"→ {a.out}")


if __name__ == "__main__":
    main()
