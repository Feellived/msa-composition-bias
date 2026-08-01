#!/usr/bin/env python3
"""[게이트 탐색] 정답을 안 보고 "이 복합체에 이 방법이 통할까"를 미리 판별할 수 있나.

■ 왜
  59종 중 30종만 본 검정에 들어왔다. 거르는 규칙이 **진짜 결합자리 덮음(recall) ≥ 0.4**,
  즉 **정답을 쓴다.** 실전에서는 못 쓰는 규칙이므로 정직한 끝값은 17/59(29%)이고,
  17/30(57%)은 "어디에 통할지 미리 안다면"이라는 조건이 붙는다.
  → 정답 없는 신호로 그 판별이 되면 29%가 57%가 된다. 이 스크립트가 그 신호를 찾는다.

■ 어떻게
  pose_features.csv 의 **정답이 안 들어간 열만** 골라, 타깃마다 평균·최대·최소·표준편차로
  요약한 뒤, 적용 가능(본 검정 대상) 대 제외를 얼마나 가르는지 AUC 로 줄세운다.
  AUC 0.5 = 무작위, 1.0 = 완전 분리, 0.0 = 완전 반대(뒤집으면 완전 분리).

  ⚠️ n=59 에 양성 30 정도라 **학습 모델은 과적합한다.** 그래서 단일 피처 + 문턱만 본다.
     그리고 그 문턱조차 같은 자료에서 고르면 낙관적이므로 **한 타깃 빼기(LOO)** 로 다시 잰다.
  ⚠️ 아무것도 안 갈리면 그 자체가 결과다 — "현재 신뢰도 지표로는 적용 가능성을 예측할 수 없다".

사용 (stdlib only, CPU):
  python gate_scan.py
  python gate_scan.py --features results/pose_features.csv --maintest maintest.csv
  python gate_scan.py --top 15 --out results/gate_scan.csv
"""
import argparse
import csv
import math
import os
import statistics as st
from collections import defaultdict

# ⚠️ 정답이 들어간 열은 반드시 뺀다. 이름에 아래가 들어가면 제외.
BANNED = ("dockq", "recall", "true", "irms", "lrms", "fnat", "rmsd", "native",
          "cover", "precision", "hit", "succ", "label", "ok")


def is_banned(c):
    lc = c.lower()
    return any(b in lc for b in BANNED)


def auc(pos, neg):
    """순위 기반 AUC(= Mann-Whitney U / (n1*n2)). 동점은 절반으로 센다."""
    if not pos or not neg:
        return float("nan")
    allv = sorted(pos + neg)
    rank = {}
    i = 0
    while i < len(allv):
        j = i
        while j + 1 < len(allv) and allv[j + 1] == allv[i]:
            j += 1
        r = (i + j) / 2 + 1
        rank[allv[i]] = r
        i = j + 1
    s = sum(rank[v] for v in pos)
    n1, n2 = len(pos), len(neg)
    return (s - n1 * (n1 + 1) / 2) / (n1 * n2)


def best_split(pos, neg):
    """정확도가 가장 높은 문턱 하나. 반환 (정확도, 문턱, 방향) — 방향 +1이면 '크면 양성'."""
    vals = sorted(set(pos + neg))
    best = (0.0, float("nan"), 1)
    for k in range(len(vals)):
        t = vals[k]
        for d in (1, -1):
            ok = sum(1 for v in pos if (v >= t) == (d > 0)) + \
                 sum(1 for v in neg if (v >= t) != (d > 0))
            acc = ok / (len(pos) + len(neg))
            if acc > best[0]:
                best = (acc, t, d)
    return best


def loo_accuracy(pos, neg):
    """한 타깃 빼고 문턱을 고른 뒤 뺀 타깃에서 맞히나 — 문턱 고르기의 낙관을 걷어낸다."""
    items = [(v, 1) for v in pos] + [(v, 0) for v in neg]
    ok = 0
    for i in range(len(items)):
        tr = items[:i] + items[i + 1:]
        p = [v for v, y in tr if y]
        n = [v for v, y in tr if not y]
        if not p or not n:
            continue
        _, t, d = best_split(p, n)
        v, y = items[i]
        pred = 1 if ((v >= t) == (d > 0)) else 0
        ok += pred == y
    return ok / len(items)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--features", default="results/pose_features.csv")
    ap.add_argument("--maintest", default="maintest.csv",
                    help="본 검정 대상 명단(=적용 가능). 여기 있으면 양성, 없으면 음성")
    ap.add_argument("--targets-col", default="target")
    ap.add_argument("--top", type=int, default=12)
    ap.add_argument("--out", default="")
    a = ap.parse_args()

    rows = list(csv.DictReader(open(a.features)))
    if not rows:
        raise SystemExit(f"!! {a.features} 가 비었다")
    tcol = a.targets_col if a.targets_col in rows[0] else next(
        (c for c in rows[0] if "target" in c.lower()), None)
    if not tcol:
        raise SystemExit(f"!! 타깃 열을 못 찾음. 열 = {list(rows[0])[:12]}")

    pos_names = {r[tcol] for r in csv.DictReader(open(a.maintest))
                 if tcol in r} if os.path.exists(a.maintest) else set()
    if not pos_names:
        raise SystemExit(f"!! {a.maintest} 에서 적용 가능 명단을 못 읽었다 "
                         f"(pick_maintest_depth.py 를 --only 없이 다시 돌려 복구할 것)")

    num = [c for c in rows[0]
           if not is_banned(c) and c != tcol and
           sum(1 for r in rows[:200] if _f(r.get(c)) == _f(r.get(c))) > 100]
    if not num:
        raise SystemExit(f"!! 쓸 수 있는 숫자 열이 없다. 열 = {list(rows[0])}")

    by = defaultdict(lambda: defaultdict(list))
    for r in rows:
        for c in num:
            v = _f(r.get(c))
            if v == v:
                by[r[tcol]][c].append(v)

    tg = sorted(by)
    P = [t for t in tg if t in pos_names]
    N = [t for t in tg if t not in pos_names]
    print(f"■ 타깃 {len(tg)}종 (적용 가능 {len(P)} · 제외 {len(N)}) · 후보 피처 {len(num)}개")
    print(f"  ⚠️ 정답이 들어간 열은 이름으로 걸러냈다: {', '.join(BANNED)}\n")
    if not P or not N:
        raise SystemExit("!! 한쪽 무리가 비었다 — 명단을 확인할 것")

    out = []
    for c in num:
        for agg, fn in (("평균", st.mean), ("최대", max), ("최소", min),
                        ("표준편차", lambda v: st.pstdev(v) if len(v) > 1 else 0.0)):
            p = [fn(by[t][c]) for t in P if by[t].get(c)]
            n = [fn(by[t][c]) for t in N if by[t].get(c)]
            if len(p) < len(P) * 0.8 or len(n) < len(N) * 0.8:
                continue
            A = auc(p, n)
            if A != A:
                continue
            acc, thr, d = best_split(p, n)
            out.append(dict(feature=c, agg=agg, auc=round(A, 3),
                            sep=round(abs(A - 0.5) * 2, 3), acc=round(acc, 3),
                            thr=round(thr, 4), dir=("클수록 적용가능" if d > 0 else "작을수록 적용가능")))
    out.sort(key=lambda r: -r["sep"])

    print(f"{'피처':<24}{'집계':<6}{'AUC':>7}{'분리력':>7}{'문턱정확도':>10}   방향")
    print("-" * 82)
    for r in out[:a.top]:
        print(f"{r['feature'][:23]:<24}{r['agg']:<6}{r['auc']:>7.3f}{r['sep']:>7.3f}"
              f"{r['acc']:>10.3f}   {r['dir']}")

    print("\n■ 상위 3개를 한 타깃 빼기로 다시 재기 (문턱 고르기의 낙관을 걷어냄)")
    base = max(len(P), len(N)) / (len(P) + len(N))
    print(f"  기준선(많은 쪽으로 다 찍기) = {base:.3f}")
    for r in out[:3]:
        fn = {"평균": st.mean, "최대": max, "최소": min,
              "표준편차": lambda v: st.pstdev(v) if len(v) > 1 else 0.0}[r["agg"]]
        p = [fn(by[t][r["feature"]]) for t in P if by[t].get(r["feature"])]
        n = [fn(by[t][r["feature"]]) for t in N if by[t].get(r["feature"])]
        lo = loo_accuracy(p, n)
        mark = "✅ 기준선보다 나음" if lo > base + 0.05 else "✗ 기준선과 다를 바 없음"
        print(f"  {r['feature'][:23]:<24}{r['agg']:<6} 표 안 {r['acc']:.3f} → LOO {lo:.3f}   {mark}")

    print("\n⚠️ LOO 가 기준선을 못 넘으면 그것이 결과다 — 현재 신호로는 적용 가능성을 못 가른다.")
    if a.out and out:
        os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
        with open(a.out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(out[0])); w.writeheader(); w.writerows(out)
        print(f"→ {a.out}")


def _f(x):
    try:
        return float(x)
    except Exception:
        return float("nan")


if __name__ == "__main__":
    main()
