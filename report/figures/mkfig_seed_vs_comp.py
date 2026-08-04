# -*- coding: utf-8 -*-
"""그림(4.3절) — 실행 수를 늘릴수록 조성군이 시드군보다 자리를 더 빨리 넓힌다.

그림 2(대각선 산점도)와 형태가 겹치지 않도록 성장 곡선으로 그린다.
두 팔이 n=1 에서 같은 점에서 출발해 갈라지는 것이 이 절의 주장 그대로다.
"""
import csv
import glob
import statistics as st
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams["font.family"] = "AppleGothic"
plt.rcParams["axes.unicode_minus"] = False

BLUE, GRAY = "#1F4E79", "#8C8C8C"
SRC = "/Users/zzuhyeong2/projects/bk21-msa-depth-bias/pipeline/results/seedcomp_by_target/*.csv"

byn = {n: [[], []] for n in (1, 2, 3, 4)}
for f in sorted(glob.glob(SRC)):
    for r in csv.DictReader(open(f)):
        n = int(r["n_run"])
        if n in byn:
            byn[n][0].append(float(r["comp_site"]))
            byn[n][1].append(float(r["seed_site"]))

ns = [1, 2, 3, 4]
comp = [st.mean(byn[n][0]) for n in ns]
seed = [st.mean(byn[n][1]) for n in ns]

fig, ax = plt.subplots(figsize=(3.40, 2.45), dpi=400)

ax.fill_between(ns, seed, comp, color=BLUE, alpha=0.12, lw=0, zorder=1)
ax.plot(ns, comp, "-o", color=BLUE, lw=1.6, ms=5, zorder=3,
        label="조성을 바꿈")
ax.plot(ns, seed, "--s", color=GRAY, lw=1.4, ms=4.2, zorder=3,
        markerfacecolor="white", markeredgewidth=1.1, label="시드만 바꿈")

ax.set_xticks(ns)
ax.set_xlim(0.85, 4.15)
ax.set_ylim(0.95, 2.42)
ax.set_yticks([1.0, 1.5, 2.0])
ax.set_xlabel("실행 수", fontsize=9.5)
ax.set_ylabel("서로 구별되는 결합 자리 수", fontsize=9.5)
ax.tick_params(labelsize=8.5, length=3, width=0.8)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
for s in ("left", "bottom"):
    ax.spines[s].set_linewidth(0.9)

ax.legend(fontsize=8, loc="upper left", frameon=False,
          handletextpad=0.5, borderpad=0.2, labelspacing=0.35)

fig.tight_layout(pad=0.4)
fig.savefig("/tmp/fig7/F7b_growth.png", dpi=400)
print("조성군", [round(v, 2) for v in comp])
print("시드군", [round(v, 2) for v in seed])
