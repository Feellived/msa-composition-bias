#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
렌즈1: '오라클 이득' (적응적 깊이 선택의 가치)
복합체별 수평 덤벨: full(rung0) DockQ vs best-over-depth DockQ.
그룹평균 꺾은선이 숨기는 per-complex 깊이 신호를 드러낸다.
"""
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import pandas as pd
import numpy as np

fm.fontManager.addfont('/System/Library/Fonts/Supplemental/AppleGothic.ttf')
matplotlib.rcParams['font.family'] = 'AppleGothic'
matplotlib.rcParams['axes.unicode_minus'] = False

BASE = os.environ.get('ANALYSIS_DIR', 'analysis')   # 중간 산출 CSV 폴더
CSV = f'{BASE}/dockq_sweep_boltz.csv'
OUT = f'{BASE}/viz_oracle.png'

AB_COLOR = {'A': '#2c7fb8', 'B': '#d95f0e', 'C': '#7a7a7a'}
THRESH = [0.23, 0.49, 0.80]
GROUP_ORDER = ['RBD', 'HA', 'Env', 'C']
BAR_COLOR = '#3b6b7a'   # 그룹 막대 (ab 색과 구분되는 중립 슬레이트)
FULL_COLOR = '#b0b0b0'  # full(rung0) 마커 (연한 회색)

# ---------- 데이터 준비 ----------
df = pd.read_csv(CSV)
full = df[df.rung == 0].set_index('target')['best_dockq']
best = df.groupby('target')['best_dockq'].max()
meta = df.groupby('target').agg(group=('group', 'first'), ab=('ab', 'first'))
o = meta.copy()
o['full'] = full
o['best'] = best
o['gain'] = o['best'] - o['full']

# 요약 수치
grp_gain = o.groupby('group')['gain'].mean().reindex(GROUP_ORDER)
mean_gain = o['gain'].mean()
median_gain = o['gain'].median()
# 임계 통과 카운트
cross = {t: int(((o.best >= t) & (o.full < t)).sum()) for t in THRESH}
full_pass = {t: int((o.full >= t).sum()) for t in THRESH}
best_pass = {t: int((o.best >= t).sum()) for t in THRESH}

# ---------- 행 배치: 그룹 블록, 블록 내부 best 오름차순 ----------
rows = []          # (target, ab, full, best, gain, y)
group_spans = []   # (group, y_top, y_bot, y_center)
y = 0.0
GAP = 1.4          # 그룹 사이 간격
for g in GROUP_ORDER:
    sub = o[o.group == g].sort_values('best')  # best 기준 정렬(그룹 내부)
    y_start = y
    for tgt, r in sub.iterrows():
        rows.append((tgt, r.ab, r.full, r.best, r.gain, y))
        y += 1.0
    y_end = y - 1.0
    group_spans.append((g, y_start, y_end, (y_start + y_end) / 2.0))
    y += GAP
ymax = y

rows_df = pd.DataFrame(rows, columns=['target', 'ab', 'full', 'best', 'gain', 'y'])

# ---------- Figure 레이아웃 ----------
fig = plt.figure(figsize=(13.5, 12.2))
gs = fig.add_gridspec(
    2, 2, width_ratios=[2.35, 1.0], height_ratios=[1.0, 1.0],
    left=0.135, right=0.975, top=0.905, bottom=0.075,
    wspace=0.28, hspace=0.30)
ax = fig.add_subplot(gs[:, 0])       # 왼쪽 큰 덤벨 패널
ax_bar = fig.add_subplot(gs[0, 1])   # 우상: 그룹 평균 이득 막대
ax_hist = fig.add_subplot(gs[1, 1])  # 우하: best-full 분포

# ===== (1) 덤벨 패널 =====
# 그룹 블록 배경 음영(교대) + 그룹 라벨 + 경계선
for i, (g, y0, y1, yc) in enumerate(group_spans):
    if i % 2 == 0:
        ax.axhspan(y0 - 0.5, y1 + 0.5, color='#000000', alpha=0.035, zorder=0)
    # 그룹 라벨(왼쪽 여백)
    ax.text(-0.135, yc, f'{g}\n(n={len(o[o.group==g])})',
            transform=ax.get_yaxis_transform(),
            ha='center', va='center', fontsize=11, fontweight='bold',
            color='#333333', linespacing=1.3)
    # 그룹 경계선
    if i < len(group_spans) - 1:
        yb = y1 + 0.5 + GAP / 2.0
        ax.axhline(yb, color='#cccccc', lw=0.8, zorder=1)

# 임계 세로선 (라벨은 하단 여백에 배치 — 제목과 충돌 방지)
th_labels = {0.23: '0.23\nacceptable', 0.49: '0.49\nmedium', 0.80: '0.80\nhigh'}
y_thlab = -0.55
for t in THRESH:
    ax.axvline(t, color='#999999', ls='--', lw=1.0, zorder=1)
    ax.text(t, y_thlab, th_labels[t], ha='center', va='top',
            fontsize=8.5, color='#777777', linespacing=1.1)

# 덤벨: 연결선(=이득) + full(회색 hollow) + best(ab색 filled)
for _, r in rows_df.iterrows():
    c = AB_COLOR[r.ab]
    lo, hi = min(r.full, r.best), max(r.full, r.best)
    # 연결선 = 이득 폭 (ab 색, 반투명)
    ax.plot([lo, hi], [r.y, r.y], color=c, lw=2.4, alpha=0.45,
            solid_capstyle='round', zorder=2)
    # full 마커 (rung0)
    ax.plot(r.full, r.y, 'o', mfc='white', mec=FULL_COLOR, mew=1.6,
            ms=6.5, zorder=3)
    # best 마커 (깊이 최적)
    ax.plot(r.best, r.y, 'o', color=c, ms=8.0, zorder=4,
            markeredgecolor='white', markeredgewidth=0.6)
    # full<0.23 인데 best가 통과하면 별 표시(punchline)
    if r.best >= 0.23 and r.full < 0.23:
        ax.plot(r.best + 0.028, r.y, marker='*', color='#c0392b',
                ms=9.5, zorder=5)

ax.set_yticks(rows_df['y'])
ax.set_yticklabels(rows_df['target'], fontsize=7.2)
ax.tick_params(axis='y', length=0)
ax.set_ylim(-1.35, ymax - GAP + 0.3)
ax.set_xlim(-0.02, 0.92)
ax.set_xlabel('DockQ', fontsize=11)
ax.set_title('복합체별 full(rung0) → 깊이-최적(best) DockQ',
             fontsize=12.5, fontweight='bold', pad=8)
ax.grid(axis='x', color='#eeeeee', lw=0.6, zorder=0)
for s in ['top', 'right', 'left']:
    ax.spines[s].set_visible(False)

# 범례
legend_elems = [
    Line2D([0], [0], marker='o', mfc='white', mec=FULL_COLOR, mew=1.6,
           ms=7, ls='', label='full = rung0 (전체 MSA)'),
    Line2D([0], [0], marker='o', color='#555555', ms=8, ls='',
           label='best = 깊이-최적(오라클)'),
    Line2D([0], [0], marker='*', color='#c0392b', ms=10, ls='',
           label='full<0.23 → best≥0.23 (구제)'),
    Line2D([0], [0], marker='o', color=AB_COLOR['A'], ms=8, ls='', label='항체 A'),
    Line2D([0], [0], marker='o', color=AB_COLOR['B'], ms=8, ls='', label='항체 B'),
    Line2D([0], [0], marker='o', color=AB_COLOR['C'], ms=8, ls='', label='항체 C'),
]
ax.legend(handles=legend_elems, loc='lower right', fontsize=8.3,
          framealpha=0.92, edgecolor='#dddddd', ncol=1)

# ===== (2) 우상: 그룹 평균 이득 막대 =====
gx = np.arange(len(GROUP_ORDER))
bars = ax_bar.bar(gx, grp_gain.values, color=BAR_COLOR, width=0.62, zorder=3)
ax_bar.axhline(mean_gain, color='#c0392b', ls='--', lw=1.3, zorder=2,
               label=f'전체 평균 {mean_gain:.3f}')
for xi, v in zip(gx, grp_gain.values):
    ax_bar.text(xi, v + 0.004, f'{v:.3f}', ha='center', va='bottom',
                fontsize=9, color='#333333')
ax_bar.set_xticks(gx)
ax_bar.set_xticklabels([f'{g}' for g in GROUP_ORDER], fontsize=10)
ax_bar.set_ylabel('평균 (best - full)', fontsize=10)
ax_bar.set_title('그룹 평균 오라클 이득', fontsize=11.5, fontweight='bold', pad=6)
ax_bar.set_ylim(0, max(grp_gain.max(), mean_gain) * 1.28)
ax_bar.legend(fontsize=8.2, loc='upper right', framealpha=0.9,
              edgecolor='#dddddd')
ax_bar.grid(axis='y', color='#eeeeee', lw=0.6, zorder=0)
for s in ['top', 'right']:
    ax_bar.spines[s].set_visible(False)

# ===== (3) 우하: best-full 분포 (히스토그램 + 러그, ab 색 스트립) =====
gains = o['gain'].values
bins = np.linspace(0, max(0.8, gains.max()), 21)
ax_hist.hist(gains, bins=bins, color='#9fb9c2', edgecolor='white', zorder=2)
ax_hist.axvline(mean_gain, color='#c0392b', ls='--', lw=1.3, zorder=3,
                label=f'평균 {mean_gain:.3f}')
ax_hist.axvline(median_gain, color='#2c3e50', ls=':', lw=1.3, zorder=3,
                label=f'중앙값 {median_gain:.3f}')
# ab 색 러그(개별 복합체)
for _, r in o.iterrows():
    ax_hist.plot(r.gain, -0.5, marker='|', color=AB_COLOR[r.ab],
                 ms=9, mew=1.4, zorder=4, clip_on=False)
ax_hist.set_xlabel('오라클 이득  (best - full DockQ)', fontsize=10)
ax_hist.set_ylabel('복합체 수', fontsize=10)
ax_hist.set_title('오라클 이득 분포 (n=47)', fontsize=11.5,
                  fontweight='bold', pad=6)
ax_hist.set_xlim(-0.01, max(0.8, gains.max()) + 0.02)
ax_hist.legend(fontsize=8.2, loc='upper right', framealpha=0.9,
               edgecolor='#dddddd')
ax_hist.grid(axis='y', color='#eeeeee', lw=0.6, zorder=0)
for s in ['top', 'right']:
    ax_hist.spines[s].set_visible(False)

# ===== 상단 제목 + 하단 캡션 =====
fig.suptitle('오라클 이득: 깊이를 복합체마다 맞게 고르면 얻는 값 (Boltz)',
             fontsize=15.5, fontweight='bold', x=0.5, y=0.972)
fig.text(0.135, 0.937,
         '핵심: 복합체마다 최적 MSA 깊이가 달라서 그룹평균 꺾은선에서는 서로 '
         '상쇄돼 이 신호가 사라진다. 여기서는 상쇄를 우회해 per-complex 이득을 그대로 보인다.',
         fontsize=10, color='#444444', ha='left')

cap = (f'임계 통과 수 — DockQ≥0.23: full {full_pass[0.23]} → best {best_pass[0.23]} '
       f'(구제 +{cross[0.23]})   |   ≥0.49: {full_pass[0.49]} → {best_pass[0.49]} (+{cross[0.49]})   |   '
       f'≥0.80: {full_pass[0.80]} → {best_pass[0.80]} (+{cross[0.80]})       '
       f'★ = full로는 실패, 깊이 선택으로 acceptable 통과한 복합체')
fig.text(0.5, 0.022, cap, fontsize=9, color='#555555', ha='center')

fig.savefig(OUT, dpi=155, facecolor='white')
print('saved', OUT)
print('mean_gain', round(mean_gain, 3), 'median', round(median_gain, 3),
      'max', round(o.gain.max(), 3))
print('grp_gain', grp_gain.round(3).to_dict())
print('cross', cross, 'full_pass', full_pass, 'best_pass', best_pass)
print('gain>0', int((o.gain > 0.001).sum()), 'gain>=0.1', int((o.gain >= 0.1).sum()))