#!/usr/bin/env python3
"""[선택기 평가] 후보 자리를 고르는 규칙들을 크기-정규화 지표로 나란히 비교한다.

⚠️ 왜 이 스크립트가 필요한가 (2026-08-01 발견):
  후보를 덮음(recall)만으로 평가하면 **넓게 던진 후보가 자동으로 이긴다.**
  실제로 정답이 26잔기인데 후보가 62잔기, 정답 16잔기인데 후보가 110잔기인 경우가 있고,
  덮음 기준으로 줄세우면 "가장 큰 후보를 고르기"가 모든 선택기를 이긴다(0.728).
  → 선택기 비교에는 **F1 / Jaccard / 정밀도** 같은 크기-정규화 지표를 써야 한다.

집계 방식이 결정적이다:
  AbEpiScore 자세 점수를 후보로 올릴 때 **최대값**을 쓰면 ncomp보다 낫고,
  **평균**을 쓰면 무작위보다 나쁘다. 한 조성의 자세 5개가 다 같은 자리로 가지 않기 때문이다.

입력 (둘 다 커밋돼 있음):
  results/sites_<타깃>.json                        (site_reproducibility.py --dump-sites)
  <consensus_docking>/results/abepiscore_all.csv    (score_abepitope.py)

사용:
  python eval_selectors.py --sites results --abepi ../bk21-antibody-ml/consensus_docking/results/abepiscore_all.csv
  python eval_selectors.py ... --pick ncomp_x_abemax --out results/selected_sites.csv
"""
import argparse
import csv
import glob
import json
import os
import re
import statistics as st
from collections import defaultdict


# ── 지표 ─────────────────────────────────────────────────────────────────────
def recall(c):
    return c.get("true_covered") or 0.0


def precision(c):
    return c.get("precision") or 0.0


def f1(c):
    p, r = precision(c), recall(c)
    return 0.0 if p + r == 0 else 2 * p * r / (p + r)


def jaccard(c):
    p, r = precision(c), recall(c)
    return 0.0 if p == 0 or r == 0 else 1 / (1 / p + 1 / r - 1)


METRICS = {"f1": f1, "jaccard": jaccard, "precision": precision, "recall": recall}


# ── 선택기 ───────────────────────────────────────────────────────────────────
def abepi_pick(d, comp_scores, how):
    """자세 점수를 후보로 올린다. how = mean | max | cmax(조성평균의 최대)."""
    sc = {}
    for c in d["candidates"]:
        flat = [v for cm in c["comps"] for v in comp_scores.get(cm, [])]
        per = [comp_scores[cm] for cm in c["comps"] if comp_scores.get(cm)]
        if not flat:
            continue
        sc[c["cand"]] = {"mean": st.mean(flat), "max": max(flat),
                         "cmax": max(st.mean(v) for v in per)}[how]
    return max(sc, key=sc.get) if sc else None


def ncomp_pick(d):
    """조성이 가장 많이 모인 후보. 동점이면 잔기 수가 적은 쪽."""
    return max(d["candidates"], key=lambda c: (c["n_comp"], -len(c["residues"])))["cand"]


def build_selectors(comp_scores_of, iptm_scores_of=None):
    sel = {
        "ncomp": lambda d: ncomp_pick(d),
        "largest": lambda d: max(d["candidates"], key=lambda c: len(c["residues"]))["cand"],
        "smallest": lambda d: min(d["candidates"], key=lambda c: len(c["residues"]))["cand"],
        "abepi_mean": lambda d: abepi_pick(d, comp_scores_of(d), "mean"),
        "abepi_max": lambda d: abepi_pick(d, comp_scores_of(d), "max"),
        "abepi_cmax": lambda d: abepi_pick(d, comp_scores_of(d), "cmax"),
        # 조성 수가 갈리면 그것을 믿고, 동점일 때만 AbEpiScore 로 가른다
        "ncomp_x_abemax": lambda d: (
            ncomp_pick(d) if len({c["n_comp"] for c in d["candidates"]}) > 1
            else abepi_pick(d, comp_scores_of(d), "max")),
    }
    if iptm_scores_of is not None:
        # 모델이 스스로 매기는 계면 신뢰도. 집계 방식은 AbEpiScore 와 동일하게 맞춘다.
        sel["iptm_max"] = lambda d: abepi_pick(d, iptm_scores_of(d), "max")
        sel["iptm_mean"] = lambda d: abepi_pick(d, iptm_scores_of(d), "mean")
    return sel


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sites", default="results", help="sites_<타깃>.json 이 있는 폴더")
    ap.add_argument("--abepi", required=True, help="abepiscore_all.csv")
    ap.add_argument("--min-cand", type=int, default=2, help="후보가 이보다 적으면 제외")
    ap.add_argument("--iptm", default="", help="collect_iptm.py 가 낸 CSV. 주면 ipTM 선택기도 함께 비교")
    ap.add_argument("--pick", default="", help="이 선택기로 고른 결과를 CSV 로 남긴다")
    ap.add_argument("--out", default="results/selected_sites.csv")
    ap.add_argument("--summary-out", default="results/eval_selectors_summary.csv",
                    help="화면에 찍는 지표×선택기 표를 CSV 로도 남긴다(로그 안 붙여도 되게)")
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
        raise SystemExit("!! 후보 2개 이상인 타깃이 없다")
    iptm = None
    if a.iptm:
        iptm = defaultdict(lambda: defaultdict(list))
        for r in csv.DictReader(open(a.iptm)):
            try:
                iptm[r["target"]][re.sub(r"_r\d+$", "", r["run"])].append(float(r["score"]))
            except Exception:
                pass
    SEL = build_selectors(lambda d: comp.get(d["target"], {}),
                          (lambda d: iptm.get(d["target"], {})) if iptm is not None else None)

    print(f"타깃 {len(T)}종 · 선택기 {len(SEL)}개\n")
    summary_rows = []
    for mname, m in METRICS.items():
        ceil = st.mean([max(m(c) for c in d["candidates"]) for d in T])
        rnd = st.mean([st.mean([m(c) for c in d["candidates"]]) for d in T])
        print(f"■ {mname}")
        print(f"   {'선택기':<16}{'평균':>8}{'천장대비':>9}")
        for nm, fn in SEL.items():
            v = []
            for d in T:
                k = fn(d)
                v.append(0.0 if k is None else
                         m(next(c for c in d["candidates"] if c["cand"] == k)))
            print(f"   {nm:<16}{st.mean(v):>8.3f}{st.mean(v)/ceil:>9.2f}")
            summary_rows.append([mname, nm, round(st.mean(v), 4), round(st.mean(v) / ceil, 4)])
        summary_rows.append([mname, "random", round(rnd, 4), round(rnd / ceil, 4)])
        summary_rows.append([mname, "oracle", round(ceil, 4), 1.0])
        print(f"   {'random':<16}{rnd:>8.3f}{rnd/ceil:>9.2f}")
        print(f"   {'oracle':<16}{ceil:>8.3f}{1.0:>9.2f}\n")

    # 데모 성립(고른 후보 ≠ 원래 MSA 자리) — 선택기 품질 지표가 아님에 주의
    print("■ 데모 성립 수 (고른 후보 ≠ 원래 MSA 자리)")
    print("   ⚠️ 이건 선택기가 기준선과 얼마나 다르게 고르나이지 품질이 아니다.")
    for nm, fn in SEL.items():
        n = 0
        for d in T:
            full = next((c for c in d["candidates"] if c.get("from_full_msa")), None)
            k = fn(d)
            if full and k is not None and k != full["cand"]:
                n += 1
        print(f"   {nm:<16}{n:>3}/{len(T)}")

    os.makedirs(os.path.dirname(a.summary_out) or ".", exist_ok=True)
    with open(a.summary_out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["n_target", len(T)])
        w.writerow(["metric", "selector", "mean", "ceil_ratio"])
        w.writerows(summary_rows)
    print(f"\n→ {a.summary_out}")

    if a.pick:
        if a.pick not in SEL:
            raise SystemExit(f"!! 모르는 선택기: {a.pick} (가능: {', '.join(SEL)})")
        os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
        with open(a.out, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["target", "selector", "cand", "n_comp", "n_res",
                        "true_covered", "precision", "f1", "from_full_msa"])
            for d in T:
                k = SEL[a.pick](d)
                if k is None:
                    continue
                c = next(x for x in d["candidates"] if x["cand"] == k)
                w.writerow([d["target"], a.pick, k, c["n_comp"], len(c["residues"]),
                            c.get("true_covered"), c.get("precision"), round(f1(c), 4),
                            bool(c.get("from_full_msa"))])
        print(f"\n→ {a.out}  ({a.pick} 로 고른 후보)")


if __name__ == "__main__":
    main()
