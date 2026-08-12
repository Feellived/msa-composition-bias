#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""렌즈3: 자기정규화 스파게티 (A/B/C 그룹 facet)
각 복합체의 DockQ를 자기 안에서 min-max 정규화(0~1)한 뒤 정규화 깊이에 대해 옅은 선으로.
평균이 상쇄되는 이유(반응 '모양'이 제각각: 상승/하강/피크형)를 직접 보여준다."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
fm.fontManager.addfont('/System/Library/Fonts/Supplemental/AppleGothic.ttf')
matplotlib.rcParams['font.family'] = 'AppleGothic'
matplotlib.rcParams['axes.unicode_minus'] = False
import pandas as pd
import numpy as np

BASE = '/private/tmp/claude-501/-Users-zzuhyeong2-Library-CloudStorage-GoogleDrive-a01056371120-gmail-com-My-Drive-SNU-BK-Summer-2026/be7059af-14d1-4717-a37d-c316fbae23f5/scratchpad/analysis'
df = pd.read_csv(f'{BASE}/dockq_sweep_boltz.csv')

COL = {'A': '#2c7fb8', 'B': '#d95f0e', 'C': '#7a7a7a'}
GRID = np.linspace(0, 1, 41)   # 공통 정규화-깊이 격자 (중앙값 계산용)

# 각 복합체(target,ab)를 자기 안에서 정규화
curves = {'A': [], 'B': [], 'C': []}   # (nd, ynorm, peak_nd)
peak_locs = {'A': [], 'B': [], 'C': []}
flat_n = 0
for (t, ab), g in df.groupby(['target', 'ab']):
    g = g.sort_values('rung')
    v = g.best_dockq.values.astype(float)
    r = g.rung.values.astype(float)
    nd = r / r.max()                      # 정규화 깊이: 0=full MSA(깊음) → 1=가장 얕음
    rng = v.max() - v.min()
    if rng < 1e-9:                        # 완전 평평 → 중앙 고정
        y = np.full_like(v, 0.5)
        flat_n += 1
    else:
        y = (v - v.min()) / rng
    peak_nd = nd[int(np.argmax(v))]
    curves[ab].append((nd, y))
    peak_locs[ab].append(peak_nd)

fig, axes = plt.subplots(1, 3, figsize=(15, 5.6), sharey=True)
titles = {'A': 'A 항체', 'B': 'B 항체', 'C': 'C 항체'}

for ax, ab in zip(axes, ['A', 'B', 'C']):
    c = COL[ab]
    n = len(curves[ab])
    # 옅은 스파게티 선 (복합체별)
    for nd, y in curves[ab]:
        ax.plot(nd, y, color=c, alpha=0.28, lw=1.1, zorder=2)
    # 각 복합체 최적점(peak): y=1 위치가 x축을 따라 흩어짐 = 깊이 신호의 다양성
    for nd, y in curves[ab]:
        pk = nd[int(np.argmax(y))]
        ax.plot(pk, 1.0, marker='v', color=c, alpha=0.55,
                ms=7, mec='white', mew=0.6, zorder=4)
    # 그룹 중앙값 (공통 격자에 보간 후 중앙값) — 굵게
    M = np.vstack([np.interp(GRID, nd, y) for nd, y in curves[ab]])
    med = np.nanmedian(M, axis=0)
    ax.plot(GRID, med, color='black', lw=3.6, zorder=6)
    ax.plot(GRID, med, color=c, lw=2.1, zorder=7)
    ax.text(GRID[-1], med[-1] + 0.03, '그룹 중앙값', color=c, fontsize=9.5,
            fontweight='bold', ha='right', va='bottom', zorder=8)

    # peak 위치 분포 요약
    pl = np.array(peak_locs[ab])
    f_full = (pl < 0.05).mean()
    f_int = ((pl >= 0.05) & (pl <= 0.95)).mean()
    f_shal = (pl > 0.95).mean()
    ax.text(0.5, 0.03,
            f'각 복합체 최적 깊이:  full MSA {f_full*100:.0f}%   ·   중간 {f_int*100:.0f}%   ·   가장 얕음 {f_shal*100:.0f}%',
            transform=ax.transAxes, ha='center', va='bottom', fontsize=9.4,
            color='#222', bbox=dict(boxstyle='round,pad=0.4', fc='white',
                                    ec=c, alpha=0.95, lw=1.3))

    ax.set_title(f'{titles[ab]}  (n={n} 복합체)', fontsize=13, color=c, fontweight='bold')
    ax.set_xlabel('정규화 MSA 깊이\n(0 = full MSA / 깊음  →  1 = 가장 얕음)', fontsize=10.5)
    ax.set_xlim(-0.03, 1.03)
    ax.set_ylim(-0.05, 1.12)
    ax.grid(True, alpha=0.25, lw=0.6)

axes[0].set_ylabel('자기정규화 DockQ\n(0 = 그 복합체의 최악  →  1 = 최선)', fontsize=10.5)

fig.suptitle('자기정규화 스파게티 — 복합체마다 최적 MSA 깊이가 달라 그룹 평균이 상쇄된다',
             fontsize=14.5, fontweight='bold', y=0.99)
fig.text(0.5, 0.005,
         '옅은 선 = 복합체 1개 (자기 안에서 0~1 정규화) · ▽ = 각 복합체의 최적 깊이 (x축에 흩어짐) · '
         '굵은 선 = 그룹 중앙값(거의 평평)',
         ha='center', fontsize=9.3, color='#444')

fig.tight_layout(rect=[0, 0.03, 1, 0.955])
OUT = f'{BASE}/viz_spaghetti.png'
fig.savefig(OUT, dpi=150, bbox_inches='tight')
print('saved', OUT)

# 콘솔 키수치
allpk = np.concatenate([np.array(peak_locs[a]) for a in 'ABC'])
print('total complexes', len(allpk), 'flat', flat_n)
print('peak full  %.0f%%' % ((allpk < 0.05).mean()*100))
print('peak inter %.0f%%' % (((allpk >= 0.05) & (allpk <= 0.95)).mean()*100))
print('peak shal  %.0f%%' % ((allpk > 0.95).mean()*100))