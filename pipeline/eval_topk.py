#!/usr/bin/env python3
"""[C] 후보를 '하나 고르기' 대신 '상위 k개 내놓기' + 자세 수 편향을 보정한 선택기.

두 가지를 한 스크립트에 넣는다. 둘 다 sites_*.json 재집계이고 서로 맞물린다.

■ (1) 상위 k개 — 산출물의 형태를 바꾼다
  선택기 9개가 전부 무작위와 구별되지 않았다(test_selectors.py). 그러면 하나 고르기를
  접고 **상위 2~3개를 다 내놓는 게** 실용적으로 낫다. 랩이 epitope binning 실험을 설계할 때
  "후보 자리 3개"는 충분히 쓸 만하다 — 변이체를 몇 개 만들지가 그걸로 정해진다.
  실패를 감추는 게 아니라 **고르기가 안 되니 안 골라도 되는 형태로 내는** 설계 변경이다.

  두 가지를 잰다.
    · 상위 k개 중 최고 F1 의 평균      (얼마나 좋은 걸 손에 쥐나)
    · 상위 k개 안에 '진짜 자리를 절반 이상 덮는 후보'가 있는 비율  (쓸 만한가)
  무작위로 k개 뽑기 · 천장(전체 중 최선)과 나란히 놓는다.

■ (2) 자세 수 편향 보정 — 지금 결론을 뒤집을 수 있는 유일한 재분석
  `abepi_max` 는 그 후보에 딸린 **모든 자세 점수의 최댓값**이다. 그런데 후보마다 딸린
  자세 수가 다르다. **자세가 많은 후보는 뽑기를 많이 했다는 이유만으로 최댓값이 커진다.**
  실제로 조성이 많이 모인 후보를 고르는 `ncomp` 는 0.284 로 무작위(0.326)보다 나쁘다.
  즉 "조성이 많이 모인 큰 후보"가 오히려 틀린 자리인데 `abepi_max` 가 그쪽으로 끌려가고
  있을 수 있다.

  보정 둘을 새 선택기로 넣는다.
    · abepi_top3   상위 3개 자세 점수의 평균 (개수에 덜 흔들린다)
    · abepi_maxeq  후보마다 같은 수의 자세를 무작위로 뽑은 뒤 최댓값, 여러 번 평균

  ⚠️ 이미 선택기 9개를 검정했다. **새로 검정할 것은 이 둘뿐**이라고 미리 정해 두었다
     (plan/PLAN_202608_final8days.md §6). Bonferroni 문턱 0.05/2 = 0.025.
     유의성 자체는 test_selectors.py 가 낸다. 여기서는 점수만 만든다.

■ 안전장치
  기존 `eval_selectors.py` 는 건드리지 않는다(그 숫자가 이미 보고서에 들어가 있다).
  대신 여기서 만든 점수의 최댓값이 eval_selectors 의 선택과 **같은지 자동으로 대조**하고,
  어긋나면 화면에 알린다. 어긋나면 둘 중 하나가 틀린 것이다.

사용 (conda activate boltz · pipeline/ 에서):
  python -u eval_topk.py --sites results \
      --abepi ../../bk21-antibody-ml/consensus_docking/results/abepiscore_all.csv \
      --iptm results/iptm_all.csv --out results/topk.csv
  python -u eval_topk.py --sites results/sites_refined ... --out results/topk_refined.csv
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

COVER_OK = 0.5          # '쓸 만한 후보' = 진짜 자리를 이만큼 덮는다 (site_reproducibility 와 같은 기준)


def load_scores(path):
    m = defaultdict(lambda: defaultdict(list))
    if not path:
        return None
    for r in csv.DictReader(open(path)):
        try:
            m[r["target"]][re.sub(r"_r\d+$", "", r["run"])].append(float(r["score"]))
        except Exception:
            pass
    return m


def pose_scores(c, sc):
    """후보에 딸린 자세 점수 전부 (조성별 목록을 펼친다)."""
    return [v for cm in c["comps"] for v in sc.get(cm, [])]


def build_scorers(comp_of, iptm_of=None, rng=None, ndraw=200):
    """이름 → (타깃 d → {후보id: 점수}). 최댓값을 고르면 eval_selectors 의 선택과 같아야 한다."""
    def agg(how, of):
        def f(d):
            sc, out = of(d), {}
            for c in d["candidates"]:
                flat = pose_scores(c, sc)
                if not flat:
                    continue
                if how == "mean":
                    out[c["cand"]] = st.mean(flat)
                elif how == "max":
                    out[c["cand"]] = max(flat)
                elif how == "cmax":
                    per = [sc[cm] for cm in c["comps"] if sc.get(cm)]
                    out[c["cand"]] = max(st.mean(v) for v in per) if per else None
                elif how == "top3":
                    out[c["cand"]] = st.mean(sorted(flat, reverse=True)[:3])
            return {k: v for k, v in out.items() if v is not None}
        return f

    def maxeq(of):
        """⭐ 후보마다 같은 수의 자세만 뽑아 최댓값 — 자세 수 편향을 없앤다."""
        def f(d):
            sc = of(d)
            per = {c["cand"]: pose_scores(c, sc) for c in d["candidates"]}
            per = {k: v for k, v in per.items() if v}
            if not per:
                return {}
            n = min(len(v) for v in per.values())
            if n <= 0:
                return {}
            r = rng or random.Random(0)
            out = {}
            for k, v in per.items():
                if len(v) == n:
                    out[k] = max(v)
                else:
                    out[k] = st.mean(max(r.sample(v, n)) for _ in range(ndraw))
            return out
        return f

    S = {
        # eval_selectors 와 같은 것들 — 대조용
        "ncomp":     lambda d: {c["cand"]: (c["n_comp"], -len(c["residues"])) for c in d["candidates"]},
        "largest":   lambda d: {c["cand"]: len(c["residues"]) for c in d["candidates"]},
        "smallest":  lambda d: {c["cand"]: -len(c["residues"]) for c in d["candidates"]},
        "abepi_mean": agg("mean", comp_of),
        "abepi_max":  agg("max", comp_of),
        "abepi_cmax": agg("cmax", comp_of),
        # ⭐ 새로 검정할 둘 (사전 지정)
        "abepi_top3":  agg("top3", comp_of),
        "abepi_maxeq": maxeq(comp_of),
    }
    if iptm_of is not None:
        S["iptm_max"] = agg("max", iptm_of)
        S["iptm_mean"] = agg("mean", iptm_of)
    return S


def topk_stats(T, scorer, k, f1_of, cover_of):
    """상위 k개 중 최고 F1 · 상위 k개 안에 쓸 만한 후보가 있는 비율."""
    best, hit, used = [], [], 0
    for d in T:
        sc = scorer(d)
        if not sc:
            best.append(0.0); hit.append(0)
            continue
        used += 1
        order = sorted(sc, key=lambda c: sc[c], reverse=True)[:k]
        cs = [c for c in d["candidates"] if c["cand"] in order]
        best.append(max((f1_of(c) for c in cs), default=0.0))
        hit.append(int(any(cover_of(c) >= COVER_OK for c in cs)))
    return st.mean(best), st.mean(hit), used


def baselines(T, k, f1_of, cover_of, rng, ndraw=2000):
    """무작위로 k개 뽑기(몬테카를로) · 천장(전체 중 최선)."""
    rb, rh = [], []
    for d in T:
        cs = d["candidates"]
        if len(cs) <= k:
            rb.append(max(f1_of(c) for c in cs))
            rh.append(int(any(cover_of(c) >= COVER_OK for c in cs)))
            continue
        b = h = 0.0
        for _ in range(ndraw):
            s = rng.sample(cs, k)
            b += max(f1_of(c) for c in s)
            h += any(cover_of(c) >= COVER_OK for c in s)
        rb.append(b / ndraw); rh.append(h / ndraw)
    ob = st.mean(max(f1_of(c) for c in d["candidates"]) for d in T)
    oh = st.mean(int(any(cover_of(c) >= COVER_OK for c in d["candidates"])) for d in T)
    return (st.mean(rb), st.mean(rh)), (ob, oh)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sites", default="results")
    ap.add_argument("--abepi", required=True)
    ap.add_argument("--iptm", default="")
    ap.add_argument("--min-cand", type=int, default=2)
    ap.add_argument("--ks", default="1 2 3")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", default="results/topk.csv")
    a = ap.parse_args()

    comp, iptm = load_scores(a.abepi), load_scores(a.iptm)
    T = []
    for p in sorted(glob.glob(os.path.join(a.sites, "sites_*.json"))):
        d = json.load(open(p))
        if len(d.get("candidates", [])) >= a.min_cand:
            T.append(d)
    if not T:
        raise SystemExit(f"!! 후보가 {a.min_cand}개 이상인 타깃이 없다 — --sites {a.sites}")

    rng = random.Random(a.seed)
    comp_of = lambda d: comp.get(d["target"], {})
    iptm_of = (lambda d: iptm.get(d["target"], {})) if iptm is not None else None
    S = build_scorers(comp_of, iptm_of, rng)
    f1_of, cover_of = ES.f1, ES.recall

    ncand = st.mean(len(d["candidates"]) for d in T)
    print(f"타깃 {len(T)}종 · 후보 평균 {ncand:.1f}개 · 선택기 {len(S)}개")
    print(f"'쓸 만한 후보' = 진짜 결합자리를 {COVER_OK:.0%} 이상 덮는 후보\n")

    # ── eval_selectors 와 대조 (여기 점수의 최댓값 = 거기 선택이어야 한다) ──
    OLD = ES.build_selectors(comp_of, iptm_of)
    bad = []
    for nm, fn in OLD.items():
        if nm not in S:
            continue
        for d in T:
            sc = S[nm](d)
            mine = max(sc, key=sc.get) if sc else None
            if fn(d) != mine:
                bad.append((nm, d["target"], fn(d), mine))
    if bad:
        print(f"⚠️ eval_selectors 의 선택과 어긋난 곳 {len(bad)}건 — 둘 중 하나가 틀렸다:")
        for nm, t, o, m in bad[:8]:
            print(f"   {nm:<12} {t:<11} eval_selectors={o}  eval_topk={m}")
        print()
    else:
        print("✅ eval_selectors 의 선택과 전부 일치 — 점수 정의가 같다\n")

    KS = [int(x) for x in a.ks.split()]
    rows = []
    for k in KS:
        (rb, rh), (ob, oh) = baselines(T, k, f1_of, cover_of, rng)
        print(f"■ 상위 {k}개를 내놓을 때")
        print(f"   {'선택기':<14}{'최고F1':>9}{'쓸만함비율':>12}{'천장대비':>10}")
        for nm, fn in S.items():
            b, h, used = topk_stats(T, fn, k, f1_of, cover_of)
            mark = "" if used == len(T) else f"  ({len(T)-used}종 점수없음)"
            print(f"   {nm:<14}{b:>9.3f}{h:>12.2f}{b/ob if ob else 0:>10.2f}{mark}")
            rows.append(dict(k=k, selector=nm, best_f1=round(b, 4),
                             usable_rate=round(h, 4), ceil_ratio=round(b / ob, 4) if ob else "",
                             n_target=len(T), n_scored=used))
        print(f"   {'무작위 k개':<14}{rb:>9.3f}{rh:>12.2f}{rb/ob if ob else 0:>10.2f}")
        print(f"   {'천장(전체)':<14}{ob:>9.3f}{oh:>12.2f}{1.0:>10.2f}\n")
        rows.append(dict(k=k, selector="random", best_f1=round(rb, 4),
                         usable_rate=round(rh, 4), ceil_ratio=round(rb / ob, 4) if ob else "",
                         n_target=len(T), n_scored=len(T)))
        rows.append(dict(k=k, selector="oracle", best_f1=round(ob, 4),
                         usable_rate=round(oh, 4), ceil_ratio=1.0,
                         n_target=len(T), n_scored=len(T)))

    print("⚠️ 상위 k개는 후보를 여러 개 내놓는 산출물이다. k 를 키우면 당연히 오르므로,")
    print("   무작위 k개와의 차이만이 선택기의 기여다.")
    print("   유의성은 test_selectors.py 가 따로 답한다(새로 검정할 선택기는 "
          "abepi_top3 · abepi_maxeq 둘뿐, Bonferroni 0.025).")

    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    print(f"\n→ {a.out}")


if __name__ == "__main__":
    main()
