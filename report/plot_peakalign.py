#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
렌즈2 = '피크 정렬 평균' (peak-aligned average).

목적: 그룹평균 vs MSA깊이(rung) 꺾은선은 평평해서 "깊이 신호가 없다"처럼 보이지만,
실제로는 개별 복합체가 서로 다른 깊이(rung)에서 최적이라 평균에서 상쇄되는 것.
각 복합체의 최적 깊이(argmax DockQ rung)를 x=0으로 재정렬하면 뾰족한 피크가 살아남는다.
같은 그림에 '정렬 안 한 원평균'(고정깊이 정렬)도 겹쳐 대비.
"""
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import numpy as np, pandas as pd

fm.fontManager.addfont('/System/Library/Fonts/Supplemental/AppleGothic.ttf')
matplotlib.rcParams['font.family'] = 'AppleGothic'
matplotlib.rcParams['axes.unicode_minus'] = False

BASE = os.environ.get('ANALYSIS_DIR', 'analysis')   # 중간 산출 CSV 폴더
CSV = f'{BASE}/dockq_sweep_boltz.csv'
OUT = f'{BASE}/viz_peakalign.png'

COL = {'A': '#2c7fb8', 'B': '#d95f0e', 'C': '#7a7a7a', 'ALL': '#222222'}
RAW = '#b0b0b0'
CENTER_RUNG = 5          # 원평균(고정깊이 정렬)을 rung-5 기준으로 겹침 (argmax 중앙값=5)
MIN_FRAC = 0.5           # offset별 복합체 수가 이 비율 이상일 때만 정렬곡선 표시

df = pd.read_csv(CSV)
df = df[df.rung.between(0, 9)].copy()   # 0~9는 47복합체 전부 존재 (균일)

def per_complex(sub):
    """복합체별 rung->dockq 시리즈 + argmax rung."""
    out = {}
    for tgt, g in sub.groupby('target'):
        s = g.set_index('rung').best_dockq.reindex(range(10))
        amax = int(s.idxmax())
        out[tgt] = (s, amax)
    return out

def aligned_curve(comp):
    """offset별 평균/표준오차/개수 (피크 정렬)."""
    N = len(comp)
    rows = {}
    for off in range(-9, 10):
        vals = []
        for tgt, (s, amax) in comp.items():
            r = amax + off
            if 0 <= r <= 9 and not np.isnan(s.iloc[r]):
                vals.append(s.iloc[r])
        if vals:
            v = np.array(vals)
            rows[off] = (v.mean(), v.std(ddof=0) / np.sqrt(len(v)), len(v))
    offs = sorted(k for k, (_, _, n) in rows.items() if n >= max(4, MIN_FRAC * N))
    m = np.array([rows[o][0] for o in offs])
    se = np.array([rows[o][1] for o in offs])
    n = np.array([rows[o][2] for o in offs])
    return np.array(offs), m, se, n, N

def raw_curve(comp):
    """rung별 원평균(고정깊이 정렬) -> x = rung - CENTER_RUNG."""
    mat = np.vstack([s.values for s, _ in comp.values()])   # (Ncomplex, 10)
    mean = np.nanmean(mat, axis=0)
    x = np.arange(10) - CENTER_RUNG
    return x, mean

subsets = [('전체 (A+B+C, N=47)', df, 'ALL'),
           ('항체 A', df[df.ab == 'A'], 'A'),
           ('항체 B', df[df.ab == 'B'], 'B'),
           ('항체 C', df[df.ab == 'C'], 'C')]

fig, axes = plt.subplots(2, 2, figsize=(12.4, 9.4))
axes = axes.ravel()

summary = {}
for ax, (title, sub, key) in zip(axes, subsets):
    comp = per_complex(sub)
    offs, m, se, n, N = aligned_curve(comp)
    xr, mr = raw_curve(comp)
    c = COL[key]

    # 원평균 (평평) — 고정깊이 정렬
    ax.plot(xr, mr, color=RAW, lw=2.4, ls='--', marker='o', ms=4,
            zorder=2, label='정렬 안 한 원평균\n(고정깊이 rung-5 기준)')

    # 피크 정렬 평균 (뾰족)
    ax.fill_between(offs, m - se, m + se, color=c, alpha=0.18, zorder=3)
    ax.plot(offs, m, color=c, lw=2.8, marker='o', ms=5, zorder=4,
            label='피크 정렬 평균\n(각 복합체 최적깊이=0)')
    # 피크 지점 강조
    i0 = list(offs).index(0)
    ax.scatter([0], [m[i0]], s=120, facecolor=c, edgecolor='white',
               linewidth=1.6, zorder=5)

    flat = float(np.nanmean(mr))
    peak = float(m[i0])
    lift = peak - flat
    summary[key] = dict(N=N, peak=peak, flat=flat, lift=lift,
                        m=m, offs=offs,
                        drop1=peak - float(m[list(offs).index(1)]) if 1 in offs else np.nan,
                        drop_1=peak - float(m[list(offs).index(-1)]) if -1 in offs else np.nan)

    ax.axvline(0, color='#cccccc', lw=1, zorder=1)
    ax.set_title(title, fontsize=13, fontweight='bold', color=c if key != 'ALL' else '#111')
    ax.set_xlabel('깊이 오프셋  (rung - 복합체별 최적 rung)', fontsize=10.5)
    ax.set_ylabel('평균 최적 DockQ', fontsize=10.5)
    ax.set_xlim(-6.4, 6.4)
    ax.set_xticks(range(-6, 7))
    ax.grid(alpha=0.25, lw=0.6)

    # 피크 상승폭 주석
    ax.annotate(f'피크 상승폭\n+{lift:.3f}',
                xy=(0, peak), xytext=(1.7, peak + 0.006),
                fontsize=10, color=c, fontweight='bold',
                ha='left', va='bottom')
    ax.legend(fontsize=8.4, loc='lower center', framealpha=0.9,
              handlelength=1.8, borderpad=0.6)

    # 전체 패널: argmax rung 분포 인셋 (왜 원평균이 평평한지)
    if key == 'ALL':
        amax_list = [amax for _, amax in comp.values()]
        iax = inset_axes(ax, width='40%', height='30%', loc='upper left',
                         borderpad=1.4)
        cnt = pd.Series(amax_list).value_counts().reindex(range(10), fill_value=0)
        iax.bar(range(10), cnt.values, color='#555', width=0.8)
        iax.set_title('복합체별 최적 깊이(rung) 분포', fontsize=8)
        iax.set_xlabel('최적 rung (0=full MSA)', fontsize=7)
        iax.set_ylabel('복합체 수', fontsize=7)
        iax.set_xticks(range(0, 10, 2))
        iax.tick_params(labelsize=6.5)
        iax.grid(axis='y', alpha=0.25)

fig.suptitle('렌즈2 · 피크 정렬 평균 — MSA 깊이 반응은 실재하나 복합체마다 최적 깊이가 달라 원평균에선 상쇄된다',
             fontsize=13.5, fontweight='bold', y=0.985)
fig.text(0.5, 0.005,
         '읽는 법: 회색 점선(원평균)은 거의 평평 → "깊이 무관"처럼 보임.  '
         '색 실선(각 복합체 최적깊이를 0으로 재정렬한 평균)은 0에서 뾰족한 피크 → 복합체마다 깊이 반응이 실재.',
         ha='center', fontsize=9.5, color='#333')
fig.tight_layout(rect=[0, 0.02, 1, 0.965])
fig.savefig(OUT, dpi=150, bbox_inches='tight')
print('SAVED', OUT)

# 핵심 수치 출력
for k, v in summary.items():
    print(f"[{k}] N={v['N']} peak(off0)={v['peak']:.3f} flat(raw)={v['flat']:.3f} "
          f"lift=+{v['lift']:.3f} drop_+1={v['drop1']:.3f} drop_-1={v['drop_1']:.3f}")