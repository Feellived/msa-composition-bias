#!/usr/bin/env python3
"""[선택기 유의성 검정] 각 선택기가 '후보를 무작위로 하나 고르기'보다 나은지 검정한다.

select_eval_selectors.py 는 평균 F1 만 나란히 보여줄 뿐 유의성을 말하지 않는다. 그런데
4.6절은 "우연과 구별되지 않았다" 같은 판정을 하므로 검정이 필요하다. 선택기 정의는
eval_selectors 에서 그대로 가져와, 같은 판·같은 집계로 비교한다.

두 가지를 함께 낸다(둘 다 이 저장소의 기존 관례를 따른다).
  ① 순열검정 — 귀무가설 "타깃마다 후보를 균등하게 무작위로 하나 고른다".
     그 상태에서 평균 F1 을 nperm 번 만들어 관측값이 분포의 어디쯤인지 본다(단측).
     analyze_site_reproducibility.py 의 뒤섞기 검정과 같은 방식이다.
  ② 부호검정 — 타깃마다 (선택기가 고른 후보의 F1) 대 (그 타깃 후보들의 평균 F1)을
     짝지어 비교. 3.4절 규약대로 양측이며 동률은 제외한다.

⚠️ 선택기를 여러 개 동시에 검정하므로 다중비교 보정이 필요하다. 본 스크립트는 보정
   전 p 와 함께 Bonferroni 문턱(0.05/선택기 수)을 같이 찍는다.

사용 (conda activate boltz · pipeline/ 에서):
  python -u select_test_selectors.py --sites results \
      --abepi ../../epitope-guided-docking/pipeline/results/abepiscore_all.csv \
      --iptm results/iptm_all.csv
"""
import argparse
import csv
import glob
import json
import math
import os
import random
import re
import statistics as st
from collections import defaultdict

import eval_selectors as ES


def binom_tail(k, n, p=0.5):
    return sum(math.comb(n, i) * p ** i * (1 - p) ** (n - i) for i in range(k, n + 1))


def sign_test(wins, losses):
    """양측 부호검정. 동률은 이미 제외된 상태로 받는다."""
    n = wins + losses
    if n == 0:
        return float("nan")
    return min(1.0, 2 * min(binom_tail(wins, n), binom_tail(losses, n)))


def load_scores(path):
    m = defaultdict(lambda: defaultdict(list))
    for r in csv.DictReader(open(path)):
        try:
            m[r["target"]][re.sub(r"_r\d+$", "", r["run"])].append(float(r["score"]))
        except Exception:
            pass
    return m


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sites", default="results")
    ap.add_argument("--abepi", required=True)
    ap.add_argument("--iptm", default="")
    ap.add_argument("--min-cand", type=int, default=2)
    ap.add_argument("--nperm", type=int, default=20000)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", default="results/selector_tests.csv")
    a = ap.parse_args()

    comp = load_scores(a.abepi)
    iptm = load_scores(a.iptm) if a.iptm else None

    T = []
    for p in sorted(glob.glob(os.path.join(a.sites, "sites_*.json"))):
        d = json.load(open(p))
        if len(d.get("candidates", [])) >= a.min_cand:
            T.append(d)
    if not T:
        raise SystemExit("!! 후보가 2개 이상인 타깃이 없다")

    SEL = ES.build_selectors(lambda d: comp.get(d["target"], {}),
                             (lambda d: iptm.get(d["target"], {})) if iptm else None)

    # 타깃별 후보 F1 목록과 무작위 기대값(후보 평균)
    cand_f1 = [[ES.f1(c) for c in d["candidates"]] for d in T]
    rnd_exp = [st.mean(v) for v in cand_f1]
    n = len(T)
    print(f"타깃 {n}종 · 선택기 {len(SEL)}개 · 순열 {a.nperm}회\n")

    # 순열 귀무분포는 선택기와 무관하므로 한 번만 만든다
    rng = random.Random(a.seed)
    null = [st.mean([rng.choice(v) for v in cand_f1]) for _ in range(a.nperm)]
    null_mean = st.mean(null)

    rows = []
    print(f"{'선택기':<16}{'평균F1':>8}{'무작위':>8}{'순열p':>9}{'승':>4}{'패':>4}{'동률':>5}{'부호p':>9}")
    print("-" * 66)
    for nm, fn in SEL.items():
        vals = []
        for d in T:
            k = fn(d)
            vals.append(0.0 if k is None else
                        ES.f1(next(c for c in d["candidates"] if c["cand"] == k)))
        obs = st.mean(vals)
        p_perm = (sum(1 for x in null if x >= obs) + 1) / (a.nperm + 1)
        w = sum(1 for v, e in zip(vals, rnd_exp) if v > e)
        l = sum(1 for v, e in zip(vals, rnd_exp) if v < e)
        t = n - w - l
        p_sign = sign_test(w, l)
        print(f"{nm:<16}{obs:>8.3f}{null_mean:>8.3f}{p_perm:>9.4f}{w:>4}{l:>4}{t:>5}{p_sign:>9.4f}")
        rows.append(dict(selector=nm, mean_f1=round(obs, 4), random_mean=round(null_mean, 4),
                         perm_p=round(p_perm, 5), wins=w, losses=l, ties=t,
                         sign_p=round(p_sign, 5)))

    bonf = 0.05 / len(SEL)
    print(f"\n선택기 {len(SEL)}개를 동시에 검정했으므로 Bonferroni 문턱 = 0.05/{len(SEL)} = {bonf:.4f}")
    print("순열검정은 단측(선택기가 무작위보다 낫다), 부호검정은 양측이다.")

    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    with open(a.out, "w", newline="") as fh:
        w_ = csv.DictWriter(fh, fieldnames=list(rows[0])); w_.writeheader(); w_.writerows(rows)
    print(f"→ {a.out}")


if __name__ == "__main__":
    main()
