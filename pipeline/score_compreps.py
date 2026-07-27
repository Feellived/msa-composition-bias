#!/usr/bin/env python3
"""[채점] 반복 실행 결과에서 '원래 MSA' vs '얕은 깊이'의 성공률·점수를 비교한다.

단위 = **실행 1회**(자세 5개 중 최고). 한 실행 안의 자세들은 서로 상관되어 있어 독립 표본이 아니다.

두 가지로 비교한다:
  ① 성공 횟수 — DockQ 0.49(중간 등급)·0.23(대략 맞음) 두 문턱에서 Fisher 정확검정(단측)
  ② 점수 순위 — 문턱을 안 쓰고 값 자체를 비교(Mann–Whitney, 단측). 성공선을 못 넘는
     복합체(예: 0.03 → 0.34)도 이걸로는 잡힌다.

입력 = dump_seedrep_full.py 가 만든 자세 단위 CSV.
  python dump_seedrep_full.py --data $DATA/compreps --only 8ulr_HL --csv-out results/compreps_8ulr_HL.csv
  python score_compreps.py --csv results/compreps_*.csv
"""
import argparse, csv, glob, math, os
import statistics as st
from collections import defaultdict, Counter
from itertools import combinations


def fisher(a, b, c, d):
    n1, n2, k = a + b, c + d, a + c
    if n1 == 0 or n2 == 0 or math.comb(n1 + n2, k) == 0:
        return float("nan")
    return sum(math.comb(n1, x) * math.comb(n2, k - x)
               for x in range(a, min(n1, k) + 1)) / math.comb(n1 + n2, k)


def ranks_of(vals):
    idx = sorted(range(len(vals)), key=lambda i: vals[i])
    r = [0.0] * len(vals)
    i = 0
    while i < len(idx):
        j = i
        while j + 1 < len(idx) and vals[idx[j + 1]] == vals[idx[i]]:
            j += 1
        avg = (i + j) / 2.0 + 1
        for k in range(i, j + 1):
            r[idx[k]] = avg
        i = j + 1
    return r


def mannwhitney(x, y):
    """단측 p — x가 y보다 큰 쪽. 작은 표본은 정확계산, 크면 정규근사."""
    n1, n2 = len(x), len(y)
    if n1 == 0 or n2 == 0:
        return float("nan")
    allv = list(x) + list(y)
    rk = ranks_of(allv)
    R1 = sum(rk[:n1])
    U1 = R1 - n1 * (n1 + 1) / 2
    if math.comb(n1 + n2, n1) <= 60000:                 # 정확계산
        N = n1 + n2
        cnt = 0; tot = 0
        for comb in combinations(range(N), n1):
            tot += 1
            r = sum(rk[i] for i in comb) - n1 * (n1 + 1) / 2
            if r >= U1:
                cnt += 1
        return cnt / tot
    mu = n1 * n2 / 2                                     # 정규근사(동점 보정)
    N = n1 + n2
    ties = sum(t ** 3 - t for t in Counter(allv).values())
    var = n1 * n2 / 12 * ((N + 1) - ties / (N * (N - 1))) if N > 1 else 0
    if var <= 0:
        return float("nan")
    z = (U1 - mu - 0.5) / math.sqrt(var)
    return 0.5 * math.erfc(z / math.sqrt(2))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", nargs="+", default=["results/compreps_*.csv"])
    ap.add_argument("--label", default="dockq", choices=["dockq", "recall"])
    ap.add_argument("--out", default="results/compreps_summary.csv")
    a = ap.parse_args()
    files = sorted({p for pat in a.csv for p in glob.glob(pat)})
    if not files:
        raise SystemExit("!! CSV를 못 찾음. dump_seedrep_full.py를 먼저 실행할 것.")

    print(f"단위 = 실행 1회(자세 중 최고).  지표 = {a.label}\n")
    print(f"{'target':11}{'조건':>7}{'실행':>5}{'중앙값':>8}{'최소~최대':>14}"
          f"{'≥0.49':>7}{'≥0.23':>7}")
    print("-" * 62)
    rows = []
    for fp in files:
        runs = defaultdict(list)
        tgt = ""
        for r in csv.DictReader(open(fp)):
            tgt = r.get("target", "") or tgt
            try:
                runs[r["seed"]].append(float(r[a.label]))
            except Exception:
                pass
        if not runs:
            print(f"{os.path.basename(fp):11}  (자료 없음)"); continue
        full = [max(v) for k, v in runs.items() if k.startswith("seedfull")]
        red = [max(v) for k, v in runs.items() if not k.startswith("seedfull")]
        if not full or not red:
            print(f"{tgt:11}  한쪽 조건이 없음 (원래 {len(full)}회 · 얕은 {len(red)}회) — 건너뜀")
            continue
        for nm, v in (("원래", full), ("얕은", red)):
            print(f"{tgt if nm=='원래' else '':11}{nm:>7}{len(v):>5}{st.median(v):>8.3f}"
                  f"{f'{min(v):.2f}~{max(v):.2f}':>14}"
                  f"{sum(1 for x in v if x>=0.49):>7}{sum(1 for x in v if x>=0.23):>7}")
        p49 = fisher(sum(1 for x in red if x >= .49), sum(1 for x in red if x < .49),
                     sum(1 for x in full if x >= .49), sum(1 for x in full if x < .49))
        p23 = fisher(sum(1 for x in red if x >= .23), sum(1 for x in red if x < .23),
                     sum(1 for x in full if x >= .23), sum(1 for x in full if x < .23))
        pmw = mannwhitney(red, full)
        verdict = ("✅ 얕은 쪽이 유의하게 높음" if min(p49, p23, pmw) < 0.05
                   else ("△ 방향은 맞으나 유의하지 않음" if st.median(red) > st.median(full)
                         else "✗ 차이 없음/반대"))
        print(f"{'':11}  → 성공수 검정 p(0.49)={p49:.3f} · p(0.23)={p23:.3f} · "
              f"점수순위 p={pmw:.4f}   {verdict}\n")
        rows.append(dict(target=tgt, n_full=len(full), n_red=len(red),
                         med_full=round(st.median(full), 3), med_red=round(st.median(red), 3),
                         succ49_full=sum(1 for x in full if x >= .49), succ49_red=sum(1 for x in red if x >= .49),
                         succ23_full=sum(1 for x in full if x >= .23), succ23_red=sum(1 for x in red if x >= .23),
                         p_fisher49=round(p49, 4), p_fisher23=round(p23, 4),
                         p_ranktest=round(pmw, 4), verdict=verdict))
    if rows:
        os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
        with open(a.out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            w.writeheader(); w.writerows(rows)
        print(f"→ {a.out}")
    print("※ 세 p 중 '점수순위'가 가장 민감하다(성공선을 못 넘는 복합체도 잡음). 문턱 검정은 해석이 쉬움.")


if __name__ == "__main__":
    main()
