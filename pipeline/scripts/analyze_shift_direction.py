#!/usr/bin/env python3
"""[방향 검정] MSA를 줄이면 예측 결합 위치가 '흔한 자리에서 멀어지는' 쪽으로 일관되게 움직이나.

문제의식: "MSA를 바꾸면 위치가 변한다"는 그 자체로는 사소하다(입력을 바꾸면 출력이 바뀐다).
안 사소한 것은 **방향**이다. 그리고 방향이 있다면 다음 거울상 예측이 성립해야 한다.

  A군(진짜 자리 = 흔한 자리)  : 흔한 자리에서 멀어짐 → 정답에서 **멀어짐**(recall 하락)
  B군(진짜 자리 ≠ 흔한 자리)  : 흔한 자리에서 멀어짐 → 정답에 **가까워짐**(recall 상승)

두 군이 반대로 움직이면, 그건 "MSA가 예측을 자주 관측되는 부위로 끌어당긴다"는 기제의 증거다.
한쪽만 움직이거나 둘 다 같은 방향이면 그 해석은 성립하지 않는다.

⚠️ 선택 편향 방지: 각 rung의 대표값은 자세 5개의 **평균**, 타깃의 대표값은 rung>0 전체의 **평균**을
   쓴다(가장 좋은 칸을 고르지 않는다). 그래야 "좋아진 것만 골랐다"는 반박을 받지 않는다.

지표(둘 다 '예측 접촉면 중 그 영역인 비율'):
  recall  = 진짜 결합자리와 겹침   /   overrep = 흔한 자리와 겹침

사용(stdlib only):
  python analyze_shift_direction.py                    # protenix
  python analyze_shift_direction.py --model boltz
"""
import argparse, csv, json, math, os
import statistics as st
from collections import defaultdict


def f(x):
    try:
        v = float(x)
        return None if math.isnan(v) else v
    except Exception:
        return None


def sign_test(vals):
    """부호검정 — 0보다 큰 값이 우연보다 많은가(양측). 0은 제외."""
    pos = sum(1 for v in vals if v > 0); neg = sum(1 for v in vals if v < 0)
    n = pos + neg
    if n == 0:
        return float("nan"), 0, 0
    k = max(pos, neg)
    p = 2 * sum(math.comb(n, i) for i in range(k, n + 1)) / (2 ** n)
    return min(p, 1.0), pos, neg


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", default="results/pose_features.csv")
    ap.add_argument("--model", default="protenix")
    ap.add_argument("--targets-dir", default="targets")
    ap.add_argument("--out", default="results/shift_direction.csv")
    a = ap.parse_args()

    # (target, rung) -> 자세들의 recall/overrep 평균
    acc = defaultdict(lambda: defaultdict(lambda: ([], [])))
    for r in csv.DictReader(open(a.csv)):
        if r["model"] != a.model:
            continue
        rc, ov = f(r.get("recall")), f(r.get("overrep"))
        if rc is None or ov is None:
            continue
        t, rg = r["target"], int(float(r["rung"]))
        acc[t][rg][0].append(rc); acc[t][rg][1].append(ov)

    rows = []
    for t, bk in acc.items():
        if 0 not in bk or len(bk) < 3:
            continue
        r0, o0 = st.mean(bk[0][0]), st.mean(bk[0][1])
        red = [rg for rg in bk if rg > 0]
        rr = st.mean([st.mean(bk[rg][0]) for rg in red])
        oo = st.mean([st.mean(bk[rg][1]) for rg in red])
        try:
            cj = json.load(open(os.path.join(a.targets_dir, t, "chains.json")))
            g, lab = str(cj.get("AB", "?")), str(cj.get("label", ""))
        except Exception:
            g, lab = "?", ""
        rows.append(dict(target=t, grp=g, site=lab, n_rung=len(red),
                         full_recall=round(r0, 3), red_recall=round(rr, 3), d_recall=round(rr - r0, 3),
                         full_overrep=round(o0, 3), red_overrep=round(oo, 3), d_overrep=round(oo - o0, 3)))

    print(f"모델 {a.model} · MSA를 줄였을 때 예측 위치가 어디로 움직이나")
    print("  (각 값은 자세 평균 → 얕은 칸 전체 평균. 좋은 칸을 고르지 않음 = 선택 편향 없음)\n")
    for G, desc in [("A", "진짜 자리 = 흔한 자리"), ("B", "진짜 자리 ≠ 흔한 자리")]:
        sub = [r for r in rows if r["grp"] == G]
        if not sub:
            continue
        dr = [r["d_recall"] for r in sub]; do_ = [r["d_overrep"] for r in sub]
        pr, rp, rn = sign_test(dr); po, op, on = sign_test(do_)
        print(f"■ {G}군 ({desc}) — 복합체 {len(sub)}개")
        print(f"   진짜 자리 겹침 변화  평균 {st.mean(dr):+.3f}  ·  오른 것 {rp} / 내린 것 {rn}  ·  부호검정 p={pr:.4f}")
        print(f"   흔한 자리 겹침 변화  평균 {st.mean(do_):+.3f}  ·  오른 것 {op} / 내린 것 {on}  ·  부호검정 p={po:.4f}")
        print()
    A = [r["d_recall"] for r in rows if r["grp"] == "A"]
    B = [r["d_recall"] for r in rows if r["grp"] == "B"]
    if A and B:
        print(f"[거울상 확인] 진짜 자리 겹침 변화: A군 {st.mean(A):+.3f} vs B군 {st.mean(B):+.3f}")
        if st.mean(A) < 0 < st.mean(B):
            print("   → 두 군이 반대 방향. 'MSA가 흔한 자리로 끌어당긴다'는 기제와 일치.")
        else:
            print("   → 거울상 아님. 그 기제로는 설명되지 않는다(다른 설명 필요).")
    print("\n[타깃별]")
    print(f"  {'target':11}{'군':>3}{'결합자리':>11}{'진짜자리 변화':>14}{'흔한자리 변화':>14}")
    for r in sorted(rows, key=lambda x: (x["grp"], -x["d_recall"])):
        print(f"  {r['target']:11}{r['grp']:>3}{r['site'][:10]:>11}{r['d_recall']:>+14.3f}{r['d_overrep']:>+14.3f}")

    if rows:
        os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
        with open(a.out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            w.writeheader(); w.writerows(rows)
        print(f"\n→ {a.out}")


if __name__ == "__main__":
    main()
