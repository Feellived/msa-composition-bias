#!/usr/bin/env python3
"""iDist 과대표집 점수 vs 깊이-취약성 — 가설 검정.
검정1: 과대표집(near-dup 이웃수)이 A(우세)>B(비우세) 인가 = A/B 라벨·'B=희귀에피토프' 검증.
검정2: 과대표집 클수록 full-MSA에서 잘 붙나(prior가 modal 자리 도움) + 낮을수록 깊이-취약한가(위치편향).
입력: results/overrep_idist.csv (idist_overrep.py) + epitope_recall.csv (깊이별 지표).
사용: python idist_analyze.py [--overrep results/overrep_idist.csv] [--recall epitope_recall.csv]
"""
import argparse
import numpy as np, pandas as pd
from scipy.stats import mannwhitneyu, spearmanr
import matplotlib; matplotlib.use("Agg")
import matplotlib.font_manager as fm, matplotlib.pyplot as plt
for p in ["/System/Library/Fonts/Supplemental/AppleGothic.ttf"]:
    try: fm.fontManager.addfont(p); matplotlib.rcParams["font.family"] = "AppleGothic"
    except Exception: pass
matplotlib.rcParams["axes.unicode_minus"] = False

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--overrep", default="results/overrep_idist.csv")
    ap.add_argument("--recall", default="epitope_recall.csv")
    ap.add_argument("--metric", default="n_0.15", help="과대표집 컬럼(n_0.04 near-dup / n_0.15 등)")
    ap.add_argument("--out", default="results/idist_overrep.png")
    a = ap.parse_args()

    ov = pd.read_csv(a.overrep)
    rc = pd.read_csv(a.recall)
    # 복합체별 깊이 지표: full(rung0) best_recall, 깊이-이동폭
    g = rc.sort_values("rung").groupby("target")
    dep = g.agg(rung0=("best_recall", "first"),
                rng=("best_recall", lambda x: x.max() - x.min()),
                ab=("ab", "first"), group=("group", "first")).reset_index()
    d = ov.merge(dep[["target", "rung0", "rng"]], on="target", how="left")
    d["logov"] = np.log1p(d[a.metric])

    print("=== 검정1: 과대표집 A vs B (항원군별) ===")
    for fam in sorted(d.family.unique()):
        sub = d[d.family == fam]
        A = sub[sub.ab == "A"][a.metric]; B = sub[sub.ab == "B"][a.metric]
        if len(A) and len(B):
            try: u, pv = mannwhitneyu(A, B, alternative="greater")
            except Exception: pv = float("nan")
            print(f"  [{fam}] A중앙 {A.median():.0f} (n{len(A)}) vs B중앙 {B.median():.0f} (n{len(B)})  MWU(A>B) p={pv:.3f}")
    A = d[d.ab == "A"][a.metric]; B = d[d.ab == "B"][a.metric]
    if len(A) and len(B):
        u, pv = mannwhitneyu(A, B, alternative="greater")
        print(f"  [전체] A중앙 {A.median():.0f} vs B중앙 {B.median():.0f}  MWU(A>B) p={pv:.4f}")

    print("\n=== 검정2: 과대표집 vs 깊이 지표 (Spearman) ===")
    for col, name in [("rung0", "full-MSA recall"), ("rng", "깊이-이동폭")]:
        m = d.dropna(subset=[col, a.metric])
        if len(m) > 3:
            r, pv = spearmanr(m[a.metric], m[col])
            print(f"  과대표집({a.metric}) vs {name}: rho={r:+.3f} p={pv:.3f} (n{len(m)})")

    # 그림
    col = {"A": "#2c7fb8", "B": "#d95f0e", "C": "#999999"}
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.4))
    fams = sorted(d.family.unique()); xt = []
    for i, fam in enumerate(fams):
        for k, ab in enumerate(["A", "B"]):
            v = d[(d.family == fam) & (d.ab == ab)][a.metric].values
            if not len(v): continue
            x = i * 2.5 + k
            ax[0].scatter(np.full(len(v), x) + np.random.RandomState(0).uniform(-.1, .1, len(v)), v,
                          color=col[ab], s=36, edgecolor="k", linewidth=.4, alpha=.85)
            ax[0].plot([x - .25, x + .25], [np.median(v)] * 2, "k", lw=2)
            xt.append((x, f"{fam}\n{ab}"))
    ax[0].set_xticks([t[0] for t in xt]); ax[0].set_xticklabels([t[1] for t in xt], fontsize=8)
    ax[0].set_ylabel(f"과대표집 (레퍼런스 이웃수, {a.metric})")
    ax[0].set_title("검정1: 우세(A) > 비우세(B) 인가", fontsize=10)
    for ab in ["A", "B", "C"]:
        s = d[d.ab == ab].dropna(subset=["rung0"])
        if len(s): ax[1].scatter(s[a.metric], s.rung0, color=col[ab], s=40, edgecolor="k", linewidth=.4, alpha=.85, label=ab)
    ax[1].set_xlabel(f"과대표집 ({a.metric})"); ax[1].set_ylabel("full-MSA recall (rung0)")
    ax[1].set_title("검정2: 과대표집 클수록 full서 잘 붙나", fontsize=10); ax[1].legend(fontsize=8)
    plt.tight_layout(); plt.savefig(a.out, dpi=140)
    print(f"\n→ {a.out}")

if __name__ == "__main__":
    main()