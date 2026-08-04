#!/usr/bin/env python3
"""[지표 통일 ③] 4.6절 "정답이 존재하는 N종 중 M종을 식별했다" 통계를 새로 정의해 낸다.

이 통계를 냈던 원 코드는 저장소에 남아 있지 않다. 대신 이 저장소의 다른 모든 절과
같은 23종 기준(후보가 2개 이상이라 실제로 고를 게 있는 타깃 — eval_selectors.py 의
기본값 --min-cand 2, 4.4·4.5절이 이미 쓰는 것과 같은 집합)을 그대로 써서 처음부터
다시 정의한다. "27종"이라는 예전 분모는 쓰지 않는다.

정의 (본문 3.1절의 "절반 이상 덮는다"와 통일):
  식별 성공 = AbEpiScore-max 로 고른 후보의 결합 자리 겹침이 0.5 이상
  무작위 기대값 = 타깃마다 후보 하나를 무작위로 골랐을 때 성공할 확률(=성공 후보 수÷전체 후보 수)의 합
  유의성 = 순열검정(타깃마다 후보를 무작위로 다시 뽑기를 반복해 우연히 이 이상 성공하는 빈도를 잰다),
           site_reproducibility.py 의 뒤섞기 검정과 같은 방식.

사용 (conda activate boltz · pipeline/ 에서, rerun_sites_all.sh 다음에):
  python -u identify_stat.py --sites results \
      --abepi ../../bk21-antibody-ml/consensus_docking/results/abepiscore_all.csv
"""
import argparse
import csv
import glob
import json
import os
import random
import re
import statistics as st
from collections import defaultdict

import eval_selectors as ES


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sites", default="results")
    ap.add_argument("--abepi", required=True)
    ap.add_argument("--min-cand", type=int, default=2)
    ap.add_argument("--succ-th", type=float, default=0.5, help="절반 이상 덮으면 식별 성공")
    ap.add_argument("--nperm", type=int, default=10000)
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args()

    comp = defaultdict(lambda: defaultdict(list))
    for r in csv.DictReader(open(a.abepi)):
        try:
            comp[r["target"]][re.sub(r"_r\d+$", "", r["run"])].append(float(r["score"]))
        except Exception:
            pass

    T = []
    for p in sorted(glob.glob(os.path.join(a.sites, "sites_*.json"))):
        d = json.load(open(p))
        if len(d.get("candidates", [])) >= a.min_cand:
            T.append(d)
    if not T:
        raise SystemExit("!! 후보 2개 이상인 타깃이 없다 — rerun_sites_all.sh 먼저 돌렸는지 확인")

    n = len(T)
    print(f"타깃 {n}종 (후보 {a.min_cand}개 이상)\n")

    picked_ok = 0
    exp_random = 0.0
    per_target_success_rate = []
    rng = random.Random(a.seed)
    obs_hits = 0

    for d in T:
        cs = comp.get(d["target"], {})
        k = ES.abepi_pick(d, cs, "max")
        cand = d["candidates"]
        succ_flags = [1 if ES.recall(c) >= a.succ_th else 0 for c in cand]
        rate = sum(succ_flags) / len(cand)
        per_target_success_rate.append(rate)
        exp_random += rate
        picked = next((c for c in cand if c["cand"] == k), None) if k is not None else None
        hit = 1 if picked is not None and ES.recall(picked) >= a.succ_th else 0
        obs_hits += hit
        print(f"  {d['target']:<10} 후보 {len(cand):>2}개 · 성공후보비율 {rate:.2f} · "
              f"AbEpiScore 픽 {'성공' if hit else '실패'}")

    print(f"\nAbEpiScore-max 가 식별한 타깃: {obs_hits}/{n}")
    print(f"무작위로 골랐을 때 기대되는 성공 수: {exp_random:.1f}/{n}")

    # 순열검정 — 타깃마다 후보를 무작위로 하나 뽑는 것을 nperm 회 반복
    ge = 0
    for _ in range(a.nperm):
        h = sum(1 for d in T if rng.random() < (
            sum(1 for c in d["candidates"] if ES.recall(c) >= a.succ_th) / len(d["candidates"])))
        if h >= obs_hits:
            ge += 1
    p = (ge + 1) / (a.nperm + 1)
    print(f"순열검정(무작위 픽 {a.nperm}회 반복, 단측) p = {p:.4f}")


if __name__ == "__main__":
    main()
