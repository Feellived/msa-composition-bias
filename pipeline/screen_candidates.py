#!/usr/bin/env python3
"""[후보 선별 v2] "원래 MSA보다 얕은 깊이에서 얼마나 나아지나"로 전 복합체를 순위매김.

v1의 문제(2026-07-27 사용자 지적):
  ① 성공 문턱(0.49)을 넘어야만 후보로 셌다 → full 0.01 → 얕은깊이 0.40 인 복합체가 통째로 탈락.
  ② 칸마다 **1번씩만** 돌린 데이터로 걸렀다. 성공률이 30%인 복합체는 11칸 중 0~1칸만 성공해
     보이므로 놓친다(민감도 낮음).
→ v2는 **문턱 대신 '개선폭(gain)'으로 전 복합체를 줄세운다.** 검정은 나중에 반복 실행으로 하고,
  여기서는 **놓치지 않는 것**이 목적이다(관대하게 뽑고 실험으로 거른다).

등급: A = 개선폭 0.30↑ 이고 성공선(0.49) 넘김 / B = 개선폭 0.20↑ / C = 0.10↑ / D = 그 미만
  ⚠️ 등급은 '검정 결과'가 아니라 **실험 우선순위**다. B·C도 반복 실행하면 유의할 수 있다.

사용:
  python screen_candidates.py                       # protenix, DockQ 기준
  python screen_candidates.py --label recall        # 결합자리 회복률 기준
  python screen_candidates.py --model boltz         # boltz 재실행 후
  python screen_candidates.py --top 15
"""
import argparse, csv, math, os
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
    ap.add_argument("--label", default="dockq", choices=["dockq", "recall"])
    ap.add_argument("--succ", type=float, default=0.49)
    ap.add_argument("--top", type=int, default=0, help="상위 N개만 출력(0=전부)")
    ap.add_argument("--targets-dir", default="targets")
    ap.add_argument("--offsite-first", action="store_true",
                    help="B군(드문 자리에 붙는 항체)을 위로 — '편향 이탈' 이야기를 하려면 B군이 필요")
    ap.add_argument("--out", default="results/screen_candidates.csv")
    a = ap.parse_args()

    best, neff = defaultdict(dict), defaultdict(dict)
    for r in csv.DictReader(open(a.csv)):
        if r["model"] != a.model:
            continue
        v = f(r.get(a.label))
        if v is None:
            continue
        t, rg = r["target"], int(float(r["rung"]))
        if v > best[t].get(rg, -1):
            best[t][rg] = v
        n = f(r.get("neff80"))
        if n is not None:
            neff[t][rg] = n

    def grp_of(t):
        """chains.json의 AB(A=흔한 자리 on-site / B=드문 자리 off-site)와 결합자리 이름."""
        fp = os.path.join(a.targets_dir, t, "chains.json")
        try:
            cj = json.load(open(fp))
            return str(cj.get("AB", "?")), str(cj.get("label", ""))
        except Exception:
            return "?", ""

    rows = []
    for t, bk in best.items():
        if 0 not in bk:
            continue
        full = bk[0]
        red = {rg: v for rg, v in bk.items() if rg > 0}
        if not red:
            continue
        pk = max(red, key=lambda x: red[x])
        gain = red[pk] - full
        n_better = sum(1 for v in red.values() if v > full + 0.15)   # full 대비 뚜렷이 나은 칸 수
        n_succ = sum(1 for v in red.values() if v >= a.succ)
        if gain >= 0.30 and red[pk] >= a.succ:
            tier = "A"
        elif gain >= 0.20:
            tier = "B"
        elif gain >= 0.10:
            tier = "C"
        else:
            tier = "D"
        rows.append(dict(target=t, tier=tier, full=round(full, 3), peak=round(red[pk], 3),
                         gain=round(gain, 3), peak_rung=pk,
                         peak_neff=round(neff[t].get(pk, float("nan")), 1),
                         n_rungs_better=n_better, n_rungs_succ=n_succ))

    rows.sort(key=lambda r: -r["gain"])
    shown = rows[:a.top] if a.top else rows
    print(f"모델 {a.model} · 지표 {a.label} · 성공선 {a.succ}")
    print("개선폭 = (얕은 깊이 최고) − (원래 MSA).  등급은 검정 결과가 아니라 실험 우선순위.\n")
    print(f"{'등급':>4}{'target':>12}{'원래':>8}{'최고':>8}{'개선폭':>8}{'@rung':>7}{'neff':>9}"
          f"{'나은칸':>7}{'성공칸':>7}")
    print("-" * 74)
    for r in shown:
        print(f"{r['tier']:>4}{r['target']:>12}{r['full']:>8.3f}{r['peak']:>8.3f}{r['gain']:>+8.3f}"
              f"{r['peak_rung']:>7}{r['peak_neff']:>9.1f}{r['n_rungs_better']:>7}{r['n_rungs_succ']:>7}")

    from collections import Counter
    c = Counter(r["tier"] for r in rows)
    print(f"\n총 {len(rows)}개 — A {c['A']} · B {c['B']} · C {c['C']} · D {c['D']}")
    todo = [r for r in rows if r["tier"] in ("A", "B", "C")]
    if todo:
        print(f"\n▶ 반복 실행으로 검정할 가치가 있는 후보 = A+B+C {len(todo)}개")
        for r in todo:
            print(f"    {r['target']:12} 등급{r['tier']}  개선폭 {r['gain']:+.3f}  "
                  f"→  RUNG={r['peak_rung']} TARGET={r['target']} bash comp_x_reps.sh")
        print(f"\n  예비 검정(복합체당 얕은5회+원래5회 ≈ 16분): 총 약 {len(todo)*16}분")
    print("\n⚠️ 이 표는 칸마다 1회씩 돌린 데이터라 **놓친 것이 있을 수 있다**(민감도 낮음).")
    print("   따라서 등급 C까지 넉넉히 뽑아 실험으로 거르는 것이 맞다.")

    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    if rows:
        with open(a.out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            w.writeheader(); w.writerows(rows)
        print(f"\n→ {a.out}")


if __name__ == "__main__":
    main()
