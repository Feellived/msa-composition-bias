#!/usr/bin/env python3
"""[C 유의성] 상위 k개를 내놓을 때, 무작위로 k개 뽑는 것보다 유의하게 쓸 만한가.

■ 왜
  select_eval_topk.py 가 낸 표는 아직 기술 통계다. 특히 정제된 후보(sites_refined)의 상위 3개에서
  iptm_max·iptm_mean·abepi_max·abepi_top3 가 천장(전체 후보 중 최선)과 거의 같은 값에
  도달했는데(F1 0.548 대 천장 0.549, 쓸만함 비율 0.70 = 천장과 동일), 이게 우연이 아닌지
  검정해야 발표에 쓸 수 있다.

■ 무엇을 재나 — '쓸 만한 후보'가 상위 k개 안에 있는가(이진)
  select_test_selectors.py 는 평균 F1 을 검정했지만, 여기서는 **쓸만함 비율**
  (진짜 결합자리를 analyze_site_reproducibility.py 와 같은 기준 0.5 이상 덮는 후보가 상위 k개 안에
  있는가)을 검정한다. select_eval_topk.py §6 의 실용 기준(≥0.70)이 재는 것과 같은 지표다.

  귀무가설 = "상위 k개를 무작위로 뽑아도 이만큼 쓸 만한 게 걸린다".
  타깃마다 무작위로 k개를 뽑는 걸 nperm 번 반복해 평균 쓸만함 비율의 귀무분포를 만들고,
  선택기의 관측값이 그 분포의 어디쯤인지 본다(단측, analyze_site_reproducibility.py 의 뒤섞기 검정과 같은 방식).

  ⚠️ 다중비교: 이 진단은 select_eval_topk.py 표에서 상위 3개 성능이 좋아 보인 것들만 검정한다.
     사후에 좋아 보이는 걸 고르는 것 자체가 다중비교이므로, 표에 나온 선택기 전부를
     Bonferroni 로 같이 보정한다(0.05 / 검정한 선택기 수).

사용 (conda activate boltz · pipeline/ 에서):
  python -u select_test_topk.py --sites results/sites_refined \
      --abepi ../../epitope-guided-docking/pipeline/results/abepiscore_all.csv \
      --iptm results/iptm_all.csv --k 3
"""
import argparse
import glob
import json
import os
import random
import statistics as st

import eval_selectors as ES
import eval_topk as ET


def usable_rate_of(T, scorer, k, cover_of, cover_ok):
    """타깃마다 상위 k개 안에 '쓸 만한 후보'가 있으면 1, 없으면 0."""
    out = []
    for d in T:
        sc = scorer(d)
        if not sc:
            out.append(0)
            continue
        order = sorted(sc, key=lambda c: sc[c], reverse=True)[:k]
        cs = [c for c in d["candidates"] if c["cand"] in order]
        out.append(int(any(cover_of(c) >= cover_ok for c in cs)))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sites", default="results")
    ap.add_argument("--abepi", required=True)
    ap.add_argument("--iptm", default="")
    ap.add_argument("--min-cand", type=int, default=2)
    ap.add_argument("--k", type=int, default=3)
    ap.add_argument("--cover-ok", type=float, default=ET.COVER_OK)
    ap.add_argument("--nperm", type=int, default=20000)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", default="results/topk_tests.csv")
    a = ap.parse_args()

    comp, iptm = ET.load_scores(a.abepi), ET.load_scores(a.iptm)
    T = []
    for p in sorted(glob.glob(os.path.join(a.sites, "sites_*.json"))):
        d = json.load(open(p))
        if len(d.get("candidates", [])) >= a.min_cand:
            T.append(d)
    if not T:
        raise SystemExit(f"!! 후보가 {a.min_cand}개 이상인 타깃이 없다")

    rng = random.Random(a.seed)
    comp_of = lambda d: comp.get(d["target"], {})
    iptm_of = (lambda d: iptm.get(d["target"], {})) if iptm is not None else None
    S = ET.build_scorers(comp_of, iptm_of, rng)
    cover_of = ES.recall

    n = len(T)
    print(f"타깃 {n}종 · 상위 {a.k}개 · '쓸 만함' = 덮음 ≥ {a.cover_ok} · 순열 {a.nperm}회\n")

    # 귀무분포 — 선택기와 무관하므로 한 번만 만든다(타깃마다 후보 목록은 같다)
    cand_ids = [[c["cand"] for c in d["candidates"]] for d in T]
    def rand_hit(d, ids):
        if len(ids) <= a.k:
            picked = ids
        else:
            picked = rng.sample(ids, a.k)
        cs = [c for c in d["candidates"] if c["cand"] in picked]
        return int(any(cover_of(c) >= a.cover_ok for c in cs))
    null = [st.mean(rand_hit(d, ids) for d, ids in zip(T, cand_ids)) for _ in range(a.nperm)]
    null_mean = st.mean(null)

    rows = []
    print(f"{'선택기':<14}{'쓸만함비율':>11}{'무작위':>9}{'순열p':>9}")
    print("-" * 45)
    for nm, fn in S.items():
        hits = usable_rate_of(T, fn, a.k, cover_of, a.cover_ok)
        obs = st.mean(hits)
        p = (sum(1 for x in null if x >= obs) + 1) / (a.nperm + 1)
        print(f"{nm:<14}{obs:>11.3f}{null_mean:>9.3f}{p:>9.4f}")
        rows.append(dict(selector=nm, k=a.k, usable_rate=round(obs, 4),
                         random_mean=round(null_mean, 4), perm_p=round(p, 5)))

    bonf = 0.05 / len(S)
    print(f"\nBonferroni 문턱 = 0.05/{len(S)} = {bonf:.4f}  (표에 나온 선택기 전부를 검정한 것으로 보정)")
    print("순열검정은 단측(선택기가 무작위 k개 뽑기보다 낫다).")

    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    with open(a.out, "w", newline="") as fh:
        import csv
        w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    print(f"→ {a.out}")


if __name__ == "__main__":
    main()
