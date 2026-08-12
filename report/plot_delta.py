#!/usr/bin/env python3
# 렌즈5 — full MSA 대비 DockQ 변화(Δ) 히트맵
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.font_manager as fm, matplotlib.pyplot as plt
fm.fontManager.addfont('/System/Library/Fonts/Supplemental/AppleGothic.ttf')
matplotlib.rcParams['font.family'] = 'AppleGothic'
matplotlib.rcParams['axes.unicode_minus'] = False
import pandas as pd, numpy as np
from matplotlib.gridspec import GridSpec

BASE = os.environ.get('ANALYSIS_DIR', 'analysis')   # 중간 산출 CSV 폴더
GC = {'A': '#2c7fb8', 'B': '#d95f0e', 'C': '#7a7a7a'}

d = pd.read_csv(f'{BASE}/dockq_sweep_boltz.csv')
full = d[d.rung == 0].set_index('target').best_dockq
d['delta'] = d.apply(lambda r: r.best_dockq - full[r.target], axis=1)

# rung 0~9 는 47 복합체 전부 존재 -> 결측 없는 완전 행렬 (rung10은 34개만 있어 제외)
sub = d[d.rung <= 9].copy()
piv = sub.pivot_table(index='target', columns='rung', values='delta')      # Δ
abmap = d.drop_duplicates('target').set_index('target').ab

# 정렬: ab 그룹(A->B->C) 블록, 각 블록 내부는 '깊이 축소로 얻는 최대 이득' 내림차순
red = piv.loc[:, 1:9]                       # 축소 rung(1~9)만
maxgain = red.max(axis=1)
argrung = red.idxmax(axis=1)                # 최대 이득이 나온 깊이
order = sorted(piv.index, key=lambda t: ({'A':0,'B':1,'C':2}[abmap[t]], -maxgain[t]))
piv = piv.loc[order]
M = piv.values                              # (47,10)
targets = piv.index.tolist()
abs_ = [abmap[t] for t in targets]

# 발산 컬러맵 대칭. vmax=0.30 -> 흔한 ±0.05~0.20 이 또렷, 극단셀(±0.8)은 포화
VMAX = 0.30
cmap = plt.get_cmap('RdBu_r').copy(); cmap.set_bad('#ffffff')

fig = plt.figure(figsize=(10.2, 13.2))
gs = GridSpec(2, 2, width_ratios=[40, 1.1], height_ratios=[47, 2.2],
              hspace=0.06, wspace=0.03, left=0.20, right=0.90, top=0.905, bottom=0.075)
ax = fig.add_subplot(gs[0, 0])
axm = fig.add_subplot(gs[1, 0], sharex=ax)     # 하단: 열(깊이)별 평균 Δ = '평평한 평균'
cax = fig.add_subplot(gs[0, 1])

im = ax.imshow(M, aspect='auto', cmap=cmap, vmin=-VMAX, vmax=VMAX,
               interpolation='nearest')

# 각 복합체가 '어느 깊이에서 최대 이득'인지 표시 (이득 > 0.05 인 행만)
for i, t in enumerate(targets):
    if maxgain[t] > 0.05:
        ax.plot(argrung[t], i, marker='o', ms=4.2, mfc='none',
                mec='black', mew=1.1, zorder=5)

# y축: 복합체명 + ab 그룹색
ax.set_yticks(range(len(targets)))
ax.set_yticklabels(targets, fontsize=5.4)
for tick, a in zip(ax.get_yticklabels(), abs_):
    tick.set_color(GC[a])
ax.set_xticks(range(10))
ax.set_xticklabels(['full\n0'] + [str(r) for r in range(1, 10)], fontsize=8)
ax.tick_params(length=0)
ax.set_ylabel('항체–항원 복합체 47개  (ab 그룹별: 파랑 A · 주황 B · 회색 C)', fontsize=9)

# ab 그룹 경계선 + 왼쪽 그룹 라벨
bounds = {}
for i, a in enumerate(abs_):
    bounds.setdefault(a, [i, i]); bounds[a][1] = i
for a, (lo, hi) in bounds.items():
    if hi + 1 < len(abs_):
        ax.axhline(hi + 0.5, color='white', lw=2.2)
        ax.axhline(hi + 0.5, color='0.35', lw=0.7)
    ax.text(-0.09, (lo + hi) / 2, f'그룹 {a}\n(n={hi-lo+1})', transform=ax.get_yaxis_transform(),
            ha='right', va='center', fontsize=8.5, color=GC[a], fontweight='bold')
ax.set_xlim(-0.5, 9.5)

# 하단 스트립: 열(깊이)별 평균 Δ  — 이게 '평평'해서 신호가 없어 보이는 함정
colmean = np.nanmean(M, axis=0).reshape(1, -1)
axm.imshow(colmean, aspect='auto', cmap=cmap, vmin=-VMAX, vmax=VMAX, interpolation='nearest')
axm.set_yticks([0]); axm.set_yticklabels(['47개 평균 Δ'], fontsize=7.5)
for j in range(10):
    axm.text(j, 0, f'{colmean[0, j]:+.02f}', ha='center', va='center',
             fontsize=6.2, color='0.15')
axm.set_xticks(range(10))
axm.set_xticklabels(['full\n(최심)'] + [str(r) for r in range(1, 9)] + ['9\n(얕음)'], fontsize=7.5)
axm.tick_params(length=0)
axm.set_xlabel('MSA 깊이 rung  (0 = full MSA 최심 → 오른쪽으로 갈수록 얕음)', fontsize=9.5)
for s in axm.spines.values(): s.set_visible(False)
for s in ax.spines.values(): s.set_visible(False)

# 컬러바
cb = fig.colorbar(im, cax=cax)
cb.set_label('Δ DockQ  =  DockQ(축소) - DockQ(full)', fontsize=8.5)
cb.ax.tick_params(labelsize=7)
cb.ax.text(0.5, 1.015, '빨강 = 축소로 이득', transform=cb.ax.transAxes, ha='center',
           va='bottom', fontsize=7, color='#b2182b')
cb.ax.text(0.5, -0.02, '파랑 = 손해', transform=cb.ax.transAxes, ha='center',
           va='top', fontsize=7, color='#2166ac')

fig.suptitle('렌즈5 — full MSA 대비 DockQ 변화(Δ) 히트맵: 평균이 지우는 복합체별 깊이 신호',
             fontsize=13, fontweight='bold', y=0.968)
fig.text(0.20, 0.928,
         '위=복합체별 셀(빨강 이득·파랑 손해, ○=최대 이득 깊이)  ·  아래=47개 평균(거의 백색=평평)'
         '   →   개별 이득/손해가 평균에서 상쇄됨',
         fontsize=8.6, color='0.25')

fig.savefig(f'{BASE}/viz_delta.png', dpi=170, facecolor='white')
print('saved viz_delta.png')

# ---- 핵심 수치 ----
gain = (maxgain > 0.02).sum()
big = (np.abs(red.values) >= 0.3).sum()
print('mean-per-rung range:', round(np.nanmin(colmean),4), 'to', round(np.nanmax(colmean),4))
print('cell delta range:', round(np.nanmin(red.values),3), 'to', round(np.nanmax(red.values),3))
print('complexes gaining >0.02:', gain, '/47')
print('complexes with big shift |Δ|>=0.3:', np.sum(np.abs(red).max(axis=1)>=0.3))
print('argmax rung spread (gain>0.05):', dict(argrung[maxgain>0.05].value_counts().sort_index()))