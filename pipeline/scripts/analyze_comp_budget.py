#!/usr/bin/env python3
"""[E] 조성을 몇 개까지 줄여도 결과가 유지되나 — 랩이 실제로 돌릴 수 있는지의 근거.

■ 왜
  지금은 타깃마다 조성을 여러 번 재추첨해 예측한다. 이게 이 파이프라인에서 가장 비싼 부분이다.
  **조성을 절반으로 줄여도 후보 품질이 유지되면 랩에서 쓸 때 비용이 반이 된다.**
  발표에서 "이건 실제로 돌릴 수 있다"를 뒷받침하는 유일한 정량 근거다.

■ 새 예측이 필요 없다
  이미 돌린 조성들 중에서 **부분집합을 뽑아 다시 조립**할 뿐이다. 조성 n개를 무작위로 골라
  후보를 만들고 채점하기를 여러 번 반복해 평균을 낸다.

■ 무엇을 재나
  · F1 천장  — 만들어진 후보 중 최선. 조성이 줄면 정답을 담은 후보를 놓칠 수 있다.
  · 후보 개수 — 조성이 줄면 후보도 줄어 **고를 것이 없어진다**(후보 1개면 선택기가 죽는다).
  둘을 같이 봐야 한다. 천장만 보면 "조성 2개로도 충분"이라는 잘못된 결론이 난다.

■ 판정 기준 (plan/PLAN_202608_final8days.md §6, 결과 보기 전 고정)
  성공 = 조성 수를 **절반으로** 줄였을 때 F1 천장 하락 ≤ 0.03

  ⚠️ 이건 "조성을 줄여도 된다"는 뜻이지 "조성이 중요하지 않다"는 뜻이 아니다.
     조성이 자리를 정한다는 주장(4.2절)은 조성 **간 차이**에서 나온 것이고,
     여기서 재는 것은 몇 번 재추첨해야 그 차이를 다 보나이다.

사용 (conda activate boltz · pipeline/ 에서):
  python -u analyze_comp_budget.py --all --out results/comp_budget.csv
  python -u analyze_comp_budget.py --targets "8ulr_HL 8k3k_D" --reps 200
"""
import argparse
import csv
import glob
import os
import random
import statistics as st
from collections import defaultdict

import refine_sites as RS


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--targets", default="")
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--dir", default="results")
    ap.add_argument("--data", default=os.environ.get("DATA", "/mnt/data/msadepth") + "/compreps")
    ap.add_argument("--targets-dir", default="targets")
    ap.add_argument("--cutoff", type=float, default=5.0)
    ap.add_argument("--link", type=float, default=0.5)
    ap.add_argument("--cons-frac", type=float, default=0.5)
    ap.add_argument("--merge-frac", type=float, default=0.75,
                    help="analyze_site_reproducibility.py 의 기본값과 맞춰 둔다(합집합 아님)")
    ap.add_argument("--max-res", type=int, default=0)
    ap.add_argument("--reps", type=int, default=100, help="조성 수마다 무작위 부분집합을 몇 번 뽑나")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", default="results/comp_budget.csv")
    a = ap.parse_args()

    tg = a.targets.replace(",", " ").split()
    if a.all or not tg:
        SKIP = {"summary", "all"}
        tg = sorted(t for t in (os.path.basename(p)[9:-4]
                                for p in glob.glob(os.path.join(a.dir, "compreps_*.csv")))
                    if t not in SKIP)

    loaded = []
    for t in tg:
        p = os.path.join(a.dir, f"compreps_{t}.csv")
        if not os.path.exists(p):
            print(f"  ! {t}: {p} 없음 — 건너뜀"); continue
        d = RS.load_runs(t, p, a.data, a.targets_dir, a.cutoff)
        if d and d["groups"]:
            loaded.append(d)
            print(f"  {t:<11} 조성 {len(d['groups']):>2}가지 · 정답 {len(d['true']):>3}잔기")
        else:
            print(f"  ! {t}: 실행/자세를 못 읽음 — 건너뜀")
    if not loaded:
        raise SystemExit("!! 읽은 타깃이 없다")

    full = max(len(d["groups"]) for d in loaded)
    rng = random.Random(a.seed)
    print(f"\n타깃 {len(loaded)}종 · 조성 최대 {full}가지 · 부분집합 {a.reps}회씩")
    print(f"조립 설정 cons={a.cons_frac} merge={a.merge_frac} max_res={a.max_res or '-'}\n")
    print(f"{'조성수':>6}{'F1천장':>9}{'덮음':>8}{'정밀도':>8}{'후보수':>8}{'전량대비':>10}{'쓴타깃':>8}")
    print("-" * 57)

    rows, base = [], None
    for n in range(2, full + 1):
        f1s, recs, pres, ncs, used = [], [], [], [], 0
        for d in loaded:
            keys = list(d["groups"])
            if len(keys) < n:
                continue                       # 조성이 n개도 안 되는 타깃은 이 칸에서 뺀다
            used += 1
            vf, vr, vp, vn = [], [], [], []
            for _ in range(a.reps):
                sub = {k: d["groups"][k] for k in rng.sample(keys, n)}
                cands, ncons = RS.build(sub, a.cons_frac, a.merge_frac, a.max_res, a.link)
                s = RS.score(cands, d["true"], ncons, n)
                vf.append(s["f1"]); vr.append(s["rec"]); vp.append(s["pre"]); vn.append(s["n_cand"])
            f1s.append(st.mean(vf)); recs.append(st.mean(vr))
            pres.append(st.mean(vp)); ncs.append(st.mean(vn))
        if not f1s:
            continue
        m = st.mean(f1s)
        if n == full:
            base = m
        rows.append(dict(n_comp=n, f1=round(m, 4), recall=round(st.mean(recs), 4),
                         precision=round(st.mean(pres), 4), n_cand=round(st.mean(ncs), 2),
                         n_target=used, reps=a.reps))
        print(f"{n:>6}{m:>9.3f}{st.mean(recs):>8.3f}{st.mean(pres):>8.3f}"
              f"{st.mean(ncs):>8.2f}{'':>10}{used:>8}")

    if base:
        for r in rows:
            r["vs_full"] = round(r["f1"] - base, 4)
        print(f"\n■ 전량({full}가지) 대비")
        for r in rows:
            warn = "  ⚠️ 후보가 말라 선택기가 죽는다" if r["n_cand"] < 1.8 else ""
            print(f"   조성 {r['n_comp']:>2}가지 → F1 천장 {r['f1']:.3f} ({r['vs_full']:+.3f})"
                  f" · 후보 {r['n_cand']:.2f}개{warn}")
        half = max(2, full // 2)
        hit = next((r for r in rows if r["n_comp"] == half), None)
        if hit:
            ok = hit["vs_full"] >= -0.03 and hit["n_cand"] >= 1.8
            print(f"\n■ 판정 (사전 기준: 조성 절반에서 하락 ≤ 0.03)")
            print(f"   조성 {full}가지 → {half}가지 : F1 천장 {hit['vs_full']:+.3f}"
                  f" · 후보 {hit['n_cand']:.2f}개")
            print("   " + ("✅ 성공 — 조성을 절반으로 줄여도 된다. 랩 적용 비용이 반으로 준다."
                           if ok else
                           "✗ 실패 — 조성을 줄이면 손해다. 현재 조성 수가 필요하다."))

    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    print(f"\n→ {a.out}")


if __name__ == "__main__":
    main()
