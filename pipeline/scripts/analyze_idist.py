#!/usr/bin/env python3
"""iDist 과대표집 점수 vs 깊이-취약성 — 가설 검정.
검정1: 과대표집이 A(우세)>B(비우세)인가 = A/B 라벨·'B=희귀 에피토프' 검증.
검정2: 과대표집 클수록 full-MSA에서 잘 붙나(prior가 modal 자리 도움) + 깊이 이동폭과 관계.
⚠️ 신호마다 '과대표집' 방향이 다름: min_dist·mean_knn = 작을수록 과대표집(-1), frac_ndup·n_0.04 = 클수록(+1).
   → 모두 sign을 곱해 '클수록 과대표집'으로 통일한 뒤 A>B(greater)로 검정.
입력: results/overrep_idist.csv (analyze_idist_overrep.py) + epitope_recall.csv (깊이별 지표).
사용: python analyze_idist.py [--overrep results/overrep_idist.csv] [--recall results/epitope_recall.csv]
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

# (컬럼, 방향[+1=클수록 과대표집 / -1=작을수록], 설명)
METRICS = [("min_dist",  -1, "가장 가까운 훈련 쌍둥이까지 거리(작을수록 과대표집)"),
           ("mean_knn",  -1, "최근접 5개 평균거리(작을수록 과대표집)"),
           ("frac_ndup", +1, "near-dup(≤0.04) 비율(클수록 과대표집)"),
           ("n_0.04",    +1, "near-dup(≤0.04) 개수(클수록 과대표집)")]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--overrep", default="results/overrep_idist.csv")
    ap.add_argument("--recall", default="results/epitope_recall.csv")
    ap.add_argument("--primary", default="min_dist", help="그림·주 신호로 쓸 컬럼")
    ap.add_argument("--out", default="results/idist_overrep.png")
    a = ap.parse_args()

    ov = pd.read_csv(a.overrep)
    rc = pd.read_csv(a.recall)
    metrics = [m for m in METRICS if m[0] in ov.columns]
    # 복합체별 깊이 지표: full(rung0) best_recall + 이동폭(max-min)
    g = rc.sort_values("rung").groupby("target")
    dep = g.agg(rung0=("best_recall", "first"),
                rng=("best_recall", lambda x: x.max() - x.min())).reset_index()
    d = ov.merge(dep, on="target", how="left")
    print(f"과대표집 {len(ov)}행 · 깊이지표 병합 {d['rung0'].notna().sum()}/{len(d)}행\n")

    print("=== 검정1: 과대표집 A(우세) > B(비우세)? (방향 통일 후 MWU greater) ===")
    for mcol, sgn, desc in metrics:
        line = f"  {mcol:9}"
        for fam in sorted(d.family.unique()) + ["전체"]:
            sub = d if fam == "전체" else d[d.family == fam]
            A = (sgn * sub[sub.ab == "A"][mcol]).dropna(); B = (sgn * sub[sub.ab == "B"][mcol]).dropna()
            if len(A) and len(B):
                try: _, pv = mannwhitneyu(A, B, alternative="greater")
                except Exception: pv = float("nan")
                mark = "*" if pv < 0.05 else " "
                line += f"  {fam}:A{sgn*A.median():+.3g}/B{sgn*B.median():+.3g} p={pv:.2f}{mark}"
        print(line)
    print("  (raw 중앙값 표기; A/B는 원래 단위. * = MWU p<0.05)")

    print("\n=== 검정2: 과대표집 vs 깊이 지표 (Spearman, +=과대표집↑일수록↑) ===")
    for mcol, sgn, _ in metrics:
        row = f"  {mcol:9}"
        for col, name in [("rung0", "full recall"), ("rng", "깊이 이동폭")]:
            m = d.dropna(subset=[col, mcol])
            if len(m) > 3:
                r, pv = spearmanr(sgn * m[mcol], m[col])
                row += f"  {name}: rho={r:+.2f} p={pv:.2f}"
        print(row)
    print("  기대: full recall과 +상관(과대표집→full서 이미 잘 붙음), 깊이 이동폭과 -상관(과대표집→깊이 둔감).")

    # 그림 (주 신호 = primary, raw 단위)
    pm = a.primary if a.primary in ov.columns else metrics[0][0]
    psgn = dict((m[0], m[1]) for m in METRICS)[pm]
    note = "낮을수록" if psgn < 0 else "높을수록"
    col = {"A": "#2c7fb8", "B": "#d95f0e", "C": "#999999"}
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.4))
    fams = sorted(d.family.unique()); xt = []; rng = np.random.RandomState(0)
    for i, fam in enumerate(fams):
        for k, ab in enumerate(["A", "B", "C"]):
            v = d[(d.family == fam) & (d.ab == ab)][pm].dropna().values
            if not len(v): continue
            x = i * 3.2 + k
            ax[0].scatter(np.full(len(v), x) + rng.uniform(-.1, .1, len(v)), v,
                          color=col[ab], s=36, edgecolor="k", linewidth=.4, alpha=.85)
            ax[0].plot([x - .28, x + .28], [np.median(v)] * 2, "k", lw=2)
            xt.append((x, f"{fam}\n{ab}"))
    ax[0].set_xticks([t[0] for t in xt]); ax[0].set_xticklabels([t[1] for t in xt], fontsize=8)
    ax[0].set_ylabel(f"과대표집 신호 ({pm})")
    ax[0].set_title(f"검정1: 우세(A) vs 비우세(B) — {note} 과대표집", fontsize=10)
    for ab in ["A", "B", "C"]:
        s = d[d.ab == ab].dropna(subset=["rung0", pm])
        if len(s): ax[1].scatter(s[pm], s.rung0, color=col[ab], s=40, edgecolor="k", linewidth=.4, alpha=.85, label=ab)
    ax[1].set_xlabel(f"과대표집 ({pm}, {note} 큼)"); ax[1].set_ylabel("full-MSA recall (rung0)")
    ax[1].set_title("검정2: 과대표집 vs full서 결합", fontsize=10); ax[1].legend(fontsize=8)
    plt.tight_layout(); plt.savefig(a.out, dpi=140)
    print(f"\n→ {a.out}")

if __name__ == "__main__":
    main()
