#!/usr/bin/env python3
# 렌즈4 = '깊이효과 분포'
# 요지: 그룹평균 깊이곡선은 평평(≈0) 하지만, 개별 복합체의 변동폭은 크고
#       최적 깊이(argmax rung)가 full에 안 몰리고 흩어져 있다 = 신호가 분산돼 숨음.
import matplotlib
matplotlib.use('Agg')
import matplotlib.font_manager as fm, matplotlib.pyplot as plt
fm.fontManager.addfont('/System/Library/Fonts/Supplemental/AppleGothic.ttf')
matplotlib.rcParams['font.family'] = 'AppleGothic'
matplotlib.rcParams['axes.unicode_minus'] = False
import pandas as pd, numpy as np

BASE = '/private/tmp/claude-501/-Users-zzuhyeong2-Library-CloudStorage-GoogleDrive-a01056371120-gmail-com-My-Drive-SNU-BK-Summer-2026/be7059af-14d1-4717-a37d-c316fbae23f5/scratchpad/analysis'
OUT = BASE + '/viz_dist.png'
COL = {'A': '#2c7fb8', 'B': '#d95f0e', 'C': '#7a7a7a'}
ABS = ['A', 'B', 'C']

d = pd.read_csv(BASE + '/dockq_sweep_boltz.csv')
d['cx'] = d.target + '_' + d.ab

# ---- per-complex 요약 ----
def summ(x):
    full = x.loc[x.rung == 0, 'best_dockq']
    full = full.iloc[0] if len(full) else np.nan
    imax = x.best_dockq.idxmax()
    return pd.Series({'ab': x.ab.iloc[0], 'group': x.group.iloc[0],
                      'full': full, 'max': x.best_dockq.max(), 'min': x.best_dockq.min(),
                      'argmax_rung': int(x.loc[imax, 'rung'])})
s = d.groupby('cx').apply(summ, include_groups=False)
s['range'] = s['max'] - s['min']
s['gain'] = s['max'] - s['full']
N = len(s)

# ---- 그룹평균 깊이곡선 (평평) : 모든 47복합체가 존재하는 rung 0~9로 balanced 계산 ----
bal = d[d.rung <= 9]
mean_curve = bal.groupby('rung').best_dockq.mean()
amp = mean_curve.max() - mean_curve.min()   # 평균곡선 진폭 ≈ 0.04

# ---- figure ----
fig = plt.figure(figsize=(15.2, 5.3))
gs = fig.add_gridspec(1, 3, wspace=0.28, left=0.055, right=0.985, top=0.86, bottom=0.135)

# ========== (a) 복합체별 깊이-변화폭 (max-min) : ab별 바이올린 + 스트립 ==========
axa = fig.add_subplot(gs[0, 0])
rng = np.random.default_rng(7)
positions = [1, 2, 3]
data_by_ab = [s.loc[s.ab == ab, 'range'].values for ab in ABS]
vp = axa.violinplot(data_by_ab, positions=positions, widths=0.78,
                    showextrema=False, showmedians=False)
for body, ab in zip(vp['bodies'], ABS):
    body.set_facecolor(COL[ab]); body.set_alpha(0.22)
    body.set_edgecolor(COL[ab]); body.set_linewidth(1.3)
for pos, ab in zip(positions, ABS):
    vals = s.loc[s.ab == ab, 'range'].values
    jit = rng.uniform(-0.13, 0.13, len(vals))
    axa.scatter(pos + jit, vals, s=34, color=COL[ab], alpha=0.85,
                edgecolor='white', linewidth=0.5, zorder=3)
    med = np.median(vals)
    axa.plot([pos - 0.28, pos + 0.28], [med, med], color=COL[ab], lw=2.4, zorder=4)
    axa.text(pos - 0.33, med, f'중앙값 {med:.2f}', ha='right', va='center',
             fontsize=8.2, color=COL[ab], zorder=5)
# 평균곡선 진폭 참조선 : "평균만 보면 이만큼밖에 안 움직임"
axa.axhline(amp, ls='--', lw=1.5, color='#444')
axa.text(0.32, amp + 0.008, f'그룹평균 곡선 진폭 ≈ {amp:.2f}', ha='left', va='bottom',
         fontsize=8.5, color='#444', style='italic')
axa.set_xticks(positions); axa.set_xticklabels([f'항체 {a}\n(n={int((s.ab==a).sum())})' for a in ABS])
axa.set_ylabel('복합체별 깊이-변화폭  (max - min DockQ)', fontsize=10.5)
axa.set_title('(a) 개별 복합체는 깊이에 따라 크게 움직인다', fontsize=11.5, fontweight='bold', pad=8)
axa.set_ylim(-0.03, s['range'].max() * 1.08)
axa.set_xlim(0.25, 3.75)
axa.grid(axis='y', ls=':', alpha=0.4)
axa.spines[['top', 'right']].set_visible(False)

# 인셋 : 그룹평균 깊이곡선 (평평) — 상쇄되어 평평해 보이는 그 곡선
axin = axa.inset_axes([0.50, 0.60, 0.46, 0.34])
axin.plot(mean_curve.index, mean_curve.values, '-o', color='#444', ms=3, lw=1.4)
axin.set_ylim(0, 0.5)
axin.set_title('그룹평균 곡선 (평평)', fontsize=8, pad=2)
axin.set_xlabel('rung (0=full→얕음)', fontsize=7)
axin.tick_params(labelsize=6.5)
axin.set_yticks([0, 0.25, 0.5])
axin.grid(ls=':', alpha=0.35)

# ========== (b) 복합체별 최적깊이 argmax rung 히스토그램 (ab 누적) ==========
axb = fig.add_subplot(gs[0, 1])
bins = np.arange(-0.5, 11.5, 1)
stack = [s.loc[s.ab == ab, 'argmax_rung'].values for ab in ABS]
axb.hist(stack, bins=bins, stacked=True,
         color=[COL[a] for a in ABS], edgecolor='white', linewidth=0.6,
         label=[f'항체 {a}' for a in ABS])
n_full = int((s.argmax_rung == 0).sum())
axb.axvline(0, color='#444', ls='--', lw=1.4)
axb.set_ylim(0, 7.9)
axb.annotate(f'full(rung0)이 최적인 복합체\n= {n_full}/{N}  ({100*n_full/N:.0f}%)',
             xy=(0.05, 4.9), xytext=(2.6, 6.2), fontsize=9, color='#222',
             arrowprops=dict(arrowstyle='->', color='#444', lw=1.2))
axb.set_xlabel('최적 깊이  (argmax rung;  0 = full MSA  →  10 = 가장 얕음)', fontsize=10)
axb.set_ylabel('복합체 수', fontsize=10.5)
axb.set_title('(b) 최적 깊이는 full에 안 몰리고 흩어진다', fontsize=11.5, fontweight='bold', pad=8)
axb.set_xticks(range(0, 11))
axb.legend(fontsize=8.5, frameon=False, loc='upper right')
axb.grid(axis='y', ls=':', alpha=0.4)
axb.spines[['top', 'right']].set_visible(False)

# ========== (c) full 대비 최대이득 (best − full) 분포 (ab 누적) ==========
axc = fig.add_subplot(gs[0, 2])
bins2 = np.arange(0, s['gain'].max() + 0.06, 0.05)
stack2 = [s.loc[s.ab == ab, 'gain'].values for ab in ABS]
axc.hist(stack2, bins=bins2, stacked=True,
         color=[COL[a] for a in ABS], edgecolor='white', linewidth=0.6,
         label=[f'항체 {a}' for a in ABS])
mean_g, med_g = s.gain.mean(), s.gain.median()
axc.axvline(med_g, color='#333', ls='-', lw=1.8)
axc.axvline(mean_g, color='#333', ls='--', lw=1.6)
axc.text(med_g + 0.01, axc.get_ylim()[1]*0.92, f'중앙값 {med_g:.2f}', fontsize=8.5, color='#333')
axc.text(mean_g + 0.01, axc.get_ylim()[1]*0.78, f'평균 {mean_g:.2f}', fontsize=8.5, color='#333')
n_gain = int((s.gain > 0.05).sum())
axc.set_xlabel('full 대비 최대이득  (best - full DockQ)', fontsize=10)
axc.set_ylabel('복합체 수', fontsize=10.5)
axc.set_title('(c) 깊이를 바꾸면 얻는 이득의 분포', fontsize=11.5, fontweight='bold', pad=8)
axc.legend(fontsize=8.5, frameon=False, loc='upper right')
axc.grid(axis='y', ls=':', alpha=0.4)
axc.spines[['top', 'right']].set_visible(False)
axc.text(0.98, 0.60, f'이득 > 0.05 인 복합체\n= {n_gain}/{N}',
         transform=axc.transAxes, ha='right', fontsize=9, color='#222',
         bbox=dict(boxstyle='round', fc='#f4f4f4', ec='#ccc'))

fig.suptitle('렌즈4 · 깊이효과 분포 :  평균 변화 ≈ 0 이지만 개별 복합체의 변동폭은 크고 최적 깊이는 분산 → 신호가 숨는다',
             fontsize=13, fontweight='bold', y=0.975)
fig.savefig(OUT, dpi=150)
print('saved', OUT)

# ---- 핵심 수치 ----
print('N complexes =', N)
print('mean_curve amplitude (rung0-9) =', round(amp, 3))
print('range: median=%.3f mean=%.3f max=%.3f' % (s.range.median(), s.range.mean(), s.range.max()))
print('gain : median=%.3f mean=%.3f max=%.3f' % (s.gain.median(), s.gain.mean(), s.gain.max()))
print('argmax==0 (full best) = %d/%d (%.0f%%)' % (n_full, N, 100*n_full/N))
print('gain>0.05 = %d/%d' % (n_gain, N))
print('mean range by ab:', {a: round(s.loc[s.ab==a,'range'].mean(),3) for a in ABS})