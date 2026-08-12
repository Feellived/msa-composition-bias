#!/usr/bin/env python3
"""[요약] analyze_site_reproducibility.py 결과를 타깃별 한 줄로 모아 세 질문에 답한다.

  ① 조성이 자리를 정하는가      뒤섞기 검정 p < 0.05 인 타깃 수
  ② 정답 후보가 만들어지는가    후보 중 하나라도 정답을 절반 이상 덮는 타깃 수
  ③ 그 후보를 고를 수 있는가    선택기(ncomp)가 고른 것이 정답인 타깃 수

⚠️ ②와 ③의 격차가 이 연구의 병목이다. ②는 되는데 ③이 안 되면
   "자리는 만들어지는데 고르지를 못한다" = 다음 과제는 후보 단위 재랭커.

옛 판(정답으로 자세를 고른 것)과 새 판(AbEpiScore)을 나란히 놓으려면 두 폴더를 다 준다:
  python select_summarize_sites.py --dir results/honest
  python select_summarize_sites.py --dir results/honest --vs results       # 옛 판과 비교
"""
import argparse
import csv
import glob
import os
from collections import defaultdict


def load(d):
    """폴더의 site_repro_*.csv → {타깃: {stats, cands}}"""
    out = {}
    for p in sorted(glob.glob(os.path.join(d, "site_repro_*.csv"))):
        rows = list(csv.DictReader(open(p)))
        if not rows:
            continue
        t = rows[0]["target"]
        f = lambda k, r=rows[0]: (float(r[k]) if r.get(k) not in (None, "", "nan") else float("nan"))
        cands = []
        for r in rows:
            try:
                cands.append(dict(cand=r["cand"], n_comp=int(r["n_comp"]), n_res=int(r["n_res"]),
                                  rec=float(r["true_covered"]), pre=float(r["precision"])))
            except Exception:
                pass
        if cands:
            out[t] = dict(within=f("within"), between=f("between"), ratio=f("ratio"),
                          p=f("perm_p"), cands=cands)
    return out


def ncomp_pick(cands):
    """조성이 가장 많이 모인 후보. 동점이면 잔기 수가 적은 쪽."""
    return max(cands, key=lambda c: (c["n_comp"], -c["n_res"]))


def score(D, hit=0.5):
    sig = [t for t, v in D.items() if v["p"] == v["p"] and v["p"] < 0.05]
    gen = [t for t, v in D.items() if max(c["rec"] for c in v["cands"]) >= hit]
    pick = [t for t, v in D.items() if ncomp_pick(v["cands"])["rec"] >= hit]
    return sig, gen, pick


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default="results/honest")
    ap.add_argument("--vs", default="", help="비교할 다른 폴더(예: 옛 판이 있는 results)")
    ap.add_argument("--hit", type=float, default=0.5, help="정답 후보로 칠 덮음 기준")
    ap.add_argument("--out", default="")
    a = ap.parse_args()

    D = load(a.dir)
    if not D:
        raise SystemExit(f"!! {a.dir} 에 site_repro_*.csv 가 없다")

    print(f"■ {a.dir}   타깃 {len(D)}종\n")
    print(f"{'타깃':<11}{'후보':>4}{'조성내':>7}{'조성간':>7}{'비율':>6}{'p':>8}"
          f"{'최고덮음':>9}{'ncomp덮음':>10}  판정")
    print("-" * 88)
    rows = []
    for t, v in sorted(D.items()):
        best = max(v["cands"], key=lambda c: c["rec"])
        pk = ncomp_pick(v["cands"])
        sig = v["p"] == v["p"] and v["p"] < 0.05
        mark = ("✅ 고름" if pk["rec"] >= a.hit else
                "⚠️ 있는데 못 고름" if best["rec"] >= a.hit else "✗ 후보에 없음")
        print(f"{t:<11}{len(v['cands']):>4}{v['within']:>7.3f}{v['between']:>7.3f}"
              f"{v['ratio']:>6.2f}{v['p']:>8.4f}{best['rec']:>9.2f}{pk['rec']:>10.2f}  "
              f"{'유의 ' if sig else '     '}{mark}")
        rows.append(dict(target=t, n_cand=len(v["cands"]), within=v["within"], between=v["between"],
                         ratio=v["ratio"], perm_p=v["p"], best_recall=best["rec"],
                         best_res=best["n_res"], ncomp_recall=pk["rec"], ncomp_res=pk["n_res"],
                         significant=int(sig), generated=int(best["rec"] >= a.hit),
                         picked=int(pk["rec"] >= a.hit)))

    sig, gen, pick = score(D, a.hit)
    n = len(D)
    print(f"\n■ 세 질문")
    print(f"  ① 조성이 자리를 정하는가 (뒤섞기 p<0.05)     {len(sig):>3}/{n}")
    print(f"  ② 정답 후보가 만들어지는가 (덮음≥{a.hit})       {len(gen):>3}/{n}")
    print(f"  ③ 그것을 고를 수 있는가 (ncomp)               {len(pick):>3}/{n}")
    print(f"\n  ⚠️ ②와 ③의 격차 = {len(gen) - len(pick)}종 — 후보에 정답이 있는데 못 골랐다.")
    if gen:
        print(f"     못 고른 타깃: {' '.join(t for t in gen if t not in pick) or '없음'}")

    if a.vs:
        O = load(a.vs)
        com = sorted(set(D) & set(O))
        if com:
            s2, g2, p2 = score({t: O[t] for t in com}, a.hit)
            s1, g1, p1 = score({t: D[t] for t in com}, a.hit)
            print(f"\n■ 두 판 비교 (공통 {len(com)}종)   {a.vs} → {a.dir}")
            print(f"  ① 유의        {len(s2):>3} → {len(s1):>3}")
            print(f"  ② 후보 생성    {len(g2):>3} → {len(g1):>3}")
            print(f"  ③ 선택 성공    {len(p2):>3} → {len(p1):>3}")
            print("  ⚠️ 옛 판은 자세를 정답(DockQ)으로 골랐다. 떨어지는 게 정상이고, 그 폭이 곧 대가다.")

    if a.out:
        os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
        with open(a.out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
        print(f"\n→ {a.out}")


if __name__ == "__main__":
    main()
