#!/usr/bin/env python3
"""[채점] 반복 실행 결과에서 '원래 MSA' vs '얕은 깊이'의 성공률·점수를 비교한다.

단위 = **실행 1회**(자세 5개 중 최고). 한 실행 안의 자세들은 서로 상관되어 있어 독립 표본이 아니다.

두 가지로 비교한다:
  ① 성공 횟수 — DockQ 0.49(중간 등급)·0.23(대략 맞음) 두 문턱에서 Fisher 정확검정(단측)
  ② 점수 순위 — 문턱을 안 쓰고 값 자체를 비교(Mann–Whitney, 단측). 성공선을 못 넘는
     복합체(예: 0.03 → 0.34)도 이걸로는 잡힌다.

입력 = eval_dump_seedrep.py 가 만든 자세 단위 CSV.
  python eval_dump_seedrep.py --data $DATA/compreps --only 8ulr_HL --csv-out results/compreps_8ulr_HL.csv
  python eval_compreps.py --csv results/compreps_*.csv
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


def heterogeneity_p(counts, sizes, max_states=3_000_000):
    """'모든 조성의 성공률이 같다'는 귀무가설의 정확검정(다변량 초기하).

    통계량 = 조성별 성공수의 제곱합(클수록 조성 간에 갈림). 성공 총수를 실행들에
    무작위 배분했을 때 관측만큼 치우칠 확률을 낸다. ⚠️ 조성당 반복이 2회 이상이어야 의미 있다.
    """
    S = sum(counts); obs = sum(c * c for c in counts); n = len(sizes)
    states = 1
    for x in sizes:
        states *= (x + 1)
    if states > max_states:                      # 너무 크면 몬테카를로(재현 위해 고정 시드)
        import random
        rnd = random.Random(0); runs = []
        for i, sz in enumerate(sizes):
            runs += [i] * sz
        lab = [1] * S + [0] * (len(runs) - S)
        ge = tot = 0
        for _ in range(200_000):
            rnd.shuffle(lab)
            c = [0] * n
            for r, l in zip(runs, lab):
                c[r] += l
            tot += 1
            if sum(x * x for x in c) >= obs:
                ge += 1
        return ge / tot
    tot = [0]; ge = [0]
    def rec(i, rem, w, sq):
        if i == n:
            if rem == 0:
                tot[0] += w
                if sq >= obs:
                    ge[0] += w
            return
        for k in range(0, min(sizes[i], rem) + 1):
            rec(i + 1, rem - k, w * math.comb(sizes[i], k), sq + k * k)
    rec(0, S, 1, 0)
    return ge[0] / tot[0] if tot[0] else float("nan")


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
    ap.add_argument("--label", default="dockq", choices=["dockq", "recall", "overrep", "n_contact"],
                    help="overrep=예측 중 흔한 자리 비율(낮을수록 좋음, --lower-better와 함께)")
    ap.add_argument("--lower-better", action="store_true",
                    help="값이 낮을수록 좋은 지표(overrep 등)일 때. 부호를 뒤집어 검정")
    ap.add_argument("--succ-th", type=float, default=None,
                    help="성공 문턱(기본 dockq/recall=0.49). overrep이면 '이 값 미만'이 성공")
    ap.add_argument("--out", default="results/compreps_summary.csv")
    a = ap.parse_args()
    SGN = -1.0 if a.lower_better else 1.0          # 낮을수록 좋으면 부호 반전 후 동일 로직
    TH = a.succ_th if a.succ_th is not None else 0.49
    files = sorted({p for pat in a.csv for p in glob.glob(pat)})
    if not files:
        raise SystemExit("!! CSV를 못 찾음. eval_dump_seedrep.py를 먼저 실행할 것.")

    print(f"단위 = 실행 1회(자세 중 {'최저' if a.lower_better else '최고'}).  지표 = {a.label}"
          f"{'  (낮을수록 좋음)' if a.lower_better else ''}  성공 문턱 {'<' if a.lower_better else '≥'}{TH}\n")
    print(f"{'target':11}{'조건':>7}{'실행':>5}{'중앙값':>8}{'최소~최대':>14}"
          f"{'성공':>7}{'강한성공':>7}")
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
        pick = min if a.lower_better else max
        # 값이 하나도 없는 실행은 뺀다. 예전에는 빈 리스트에 max()를 걸어 ValueError로 죽었다 —
        # DockQ 환경 없이 돌리면 dockq 열이 통째로 비어서 실제로 그렇게 된다(2026-07-29).
        runs = {k: v for k, v in runs.items() if v}
        if not runs:
            print(f"{tgt or os.path.basename(fp):11}  '{a.label}' 값이 한 개도 없다 — 채점 환경 확인"
                  f" (dockq 지표면 conda activate DockQ)"); continue
        full = [pick(v) for k, v in runs.items() if k.startswith("seedfull")]
        red = [pick(v) for k, v in runs.items() if not k.startswith("seedfull")]
        if not full or not red:
            print(f"{tgt:11}  한쪽 조건이 없음 (원래 {len(full)}회 · 얕은 {len(red)}회) — 건너뜀")
            continue
        for nm, v in (("원래", full), ("얕은", red)):
            ok1 = sum(1 for x in v if (x < TH if a.lower_better else x >= TH))
            ok2 = sum(1 for x in v if (x < TH/2 if a.lower_better else x >= 0.23))
            print(f"{tgt if nm=='원래' else '':11}{nm:>7}{len(v):>5}{st.median(v):>8.3f}"
                  f"{f'{min(v):.2f}~{max(v):.2f}':>14}{ok1:>7}{ok2:>7}")
        # 조성별 성공률 — 반복이 2회 이상일 때만 이질성 검정이 의미 있다
        bycomp = defaultdict(list)
        for k, v in runs.items():
            if not k.startswith("seedfull"):
                bycomp[k.split("_r")[0]].append(max(v))
        reps = [len(v) for v in bycomp.values()]
        if bycomp and min(reps) >= 2:
            names = sorted(bycomp)
            cnt = [sum(1 for x in bycomp[c] if (x < TH if a.lower_better else x >= TH)) for c in names]
            ph = heterogeneity_p(cnt, [len(bycomp[c]) for c in names])
            det = " ".join(f"{c.replace('seed','')}:{k}/{len(bycomp[c])}" for c, k in zip(names, cnt))
            print(f"{'':11}  조성별 성공  {det}")
            print(f"{'':11}  → ⭐조성 간 이질성 정확검정 p = {ph:.4f}"
                  f"   {'(조성이 성공률을 좌우함)' if ph < 0.05 else '(조성 간 차이 불충분)'}")
        else:
            ph = float("nan")
            print(f"{'':11}  (조성당 반복 {min(reps) if reps else 0}회 — 이질성 검정 불가, 2회 이상 필요)")

        gp = (lambda x: x < TH) if a.lower_better else (lambda x: x >= TH)
        gp2 = (lambda x: x < TH/2) if a.lower_better else (lambda x: x >= 0.23)
        p49 = fisher(sum(1 for x in red if gp(x)), sum(1 for x in red if not gp(x)),
                     sum(1 for x in full if gp(x)), sum(1 for x in full if not gp(x)))
        p23 = fisher(sum(1 for x in red if gp2(x)), sum(1 for x in red if not gp2(x)),
                     sum(1 for x in full if gp2(x)), sum(1 for x in full if not gp2(x)))
        pmw = mannwhitney([SGN * x for x in red], [SGN * x for x in full])
        verdict = ("✅ 얕은 쪽이 유의하게 높음" if min(p49, p23, pmw) < 0.05
                   else ("△ 방향은 맞으나 유의하지 않음" if SGN * st.median(red) > SGN * st.median(full)
                         else "✗ 차이 없음/반대"))
        print(f"{'':11}  → 성공수 검정 p(0.49)={p49:.3f} · p(0.23)={p23:.3f} · "
              f"점수순위 p={pmw:.4f}   {verdict}\n")
        rows.append(dict(target=tgt, n_full=len(full), n_red=len(red),
                         med_full=round(st.median(full), 3), med_red=round(st.median(red), 3),
                         succ49_full=sum(1 for x in full if x >= .49), succ49_red=sum(1 for x in red if x >= .49),
                         succ23_full=sum(1 for x in full if x >= .23), succ23_red=sum(1 for x in red if x >= .23),
                         p_fisher49=round(p49, 4), p_fisher23=round(p23, 4),
                         p_ranktest=round(pmw, 4), p_heterogeneity=(round(ph,4) if ph==ph else ''), verdict=verdict))
    if rows:
        os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
        with open(a.out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            w.writeheader(); w.writerows(rows)
        print(f"→ {a.out}")
    print("※ 세 p 중 '점수순위'가 가장 민감하다(성공선을 못 넘는 복합체도 잡음). 문턱 검정은 해석이 쉬움.")


if __name__ == "__main__":
    main()
