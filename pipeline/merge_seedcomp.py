#!/usr/bin/env python3
"""[지표 통일 ②-완결] rerun_seedcomp_all.sh 가 타깃별로 낸 CSV를 합쳐
4.3절 문장("N종 중 M종에서 조성 쪽이 더 많은 자리를 찾았다")을 30종 기준으로 다시 낸다.

이항검정은 variance_all.py 의 binom_tail 과 같은 닫힌 형태를 그대로 쓴다
(scipy 의존 없이, 이 저장소의 기존 관례와 통일).

사용 (conda activate boltz · pipeline/ 에서, rerun_seedcomp_all.sh 다음에):
  python -u merge_seedcomp.py
"""
import csv
import glob
import math
import statistics as st


def binom_tail(k, n, p=0.5):
    return sum(math.comb(n, i) * p ** i * (1 - p) ** (n - i) for i in range(k, n + 1))


def main():
    files = sorted(glob.glob("results/seedcomp_by_target/*.csv"))
    if not files:
        raise SystemExit("!! results/seedcomp_by_target/*.csv 가 없다 — rerun_seedcomp_all.sh 먼저")

    d_site, d_reach, rows_out = [], [], []
    for f in files:
        tgt = f.split("/")[-1][:-4]
        rows = list(csv.DictReader(open(f)))
        if not rows:
            print(f"  ! {tgt}: 빈 결과 — 건너뜀"); continue
        nmax = max(int(r["n_run"]) for r in rows)
        last = [r for r in rows if int(r["n_run"]) == nmax][0]
        ds = float(last["comp_site"]) - float(last["seed_site"])
        dr = float(last["comp_reach"]) - float(last["seed_reach"])
        d_site.append(ds); d_reach.append(dr)
        rows_out.append((tgt, nmax, ds, dr))
        print(f"  {tgt:<10} n={nmax:<3} 자리차이 {ds:+.2f}  정답도달차이 {dr:+.2f}")

    n = len(d_site)
    win_site = sum(1 for v in d_site if v > 0)
    tie_site = sum(1 for v in d_site if v == 0)
    win_reach = sum(1 for v in d_reach if v > 0)

    print("\n" + "=" * 60)
    print(f"  읽은 타깃 {n}/30종 (건너뛴 것은 위에 ! 로 표시)")
    print(f"  자리 수: 조성군 승 {win_site}/{n}종 (동률 {tie_site})"
          f" · 차이 중앙값 {st.median(d_site):+.3f}")
    p_site = 2 * min(binom_tail(win_site, n - tie_site), 1 - binom_tail(win_site - 1, n - tie_site))
    print(f"    부호검정(동률 제외 n={n - tie_site}) 양측 p = {min(p_site, 1.0):.4f}")
    print(f"  정답 도달: 조성군 승 {win_reach}/{n}종 · 차이 중앙값 {st.median(d_reach):+.3f}")
    p_reach = 2 * min(binom_tail(win_reach, n), 1 - binom_tail(win_reach - 1, n))
    print(f"    부호검정 양측 p = {min(p_reach, 1.0):.4f}")

    with open("results/seedcomp_merged.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["target", "n_run", "site_diff", "reach_diff"])
        w.writerows(rows_out)
    print(f"\n→ results/seedcomp_merged.csv ({n}행)")


if __name__ == "__main__":
    main()
