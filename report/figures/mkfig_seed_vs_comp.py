# -*- coding: utf-8 -*-
"""그림 7 — 같은 실행 예산에서 조성군이 시드군보다 더 많은 자리를 찾는다 (4.3절).

그림 2(자카드 재현성)와 한 쌍으로 읽히도록 같은 형식으로 그린다:
대각선 산점도 · 같은 색 · 같은 글꼴 크기 · 제목 없음 · 범례에 종 수 표기.
"""
import csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams["font.family"] = "AppleGothic"
plt.rcParams["axes.unicode_minus"] = False

BLUE = "#1F4E79"
GRAY = "#999999"

rows = list(csv.DictReader(open("/tmp/fig7/data.csv")))
up   = [(float(r["seed_site"]), float(r["comp_site"])) for r in rows
        if float(r["comp_site"]) > float(r["seed_site"])]
rest = [(float(r["seed_site"]), float(r["comp_site"])) for r in rows
        if float(r["comp_site"]) <= float(r["seed_site"])]

fig, ax = plt.subplots(figsize=(3.30, 3.10), dpi=400)

lo, hi = 0.7, 4.3
ax.plot([lo, hi], [lo, hi], ls="--", lw=0.9, color=GRAY, zorder=1)

ax.scatter([x for x, _ in up], [y for _, y in up], s=26, color=BLUE,
           zorder=3, label=f"조성 쪽이 많음 ({len(up)})")
ax.scatter([x for x, _ in rest], [y for _, y in rest], s=26, facecolors="white",
           edgecolors=BLUE, linewidths=1.1, zorder=3,
           label=f"같거나 시드 쪽이 많음 ({len(rest)})")

ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
ax.set_xticks([1, 2, 3, 4]); ax.set_yticks([1, 2, 3, 4])
ax.set_xlabel("시드만 바꿨을 때 찾은 자리 수", fontsize=9.5)
ax.set_ylabel("조성을 바꿨을 때 찾은 자리 수", fontsize=9.5)
ax.tick_params(labelsize=8.5, length=3, width=0.8)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
for s in ("left", "bottom"):
    ax.spines[s].set_linewidth(0.9)

ax.legend(fontsize=7.6, loc="lower right", frameon=False,
          handletextpad=0.4, borderpad=0.2, labelspacing=0.35)

fig.tight_layout(pad=0.4)
fig.savefig("/tmp/fig7/F7_seed_vs_comp.png", dpi=400)
print("저장 완료 · 조성 우세", len(up), "· 나머지", len(rest))
