import pandas as pd, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.font_manager as fm, matplotlib.pyplot as plt
fm.fontManager.addfont("/System/Library/Fonts/Supplemental/AppleGothic.ttf")
matplotlib.rcParams["font.family"]="AppleGothic"; matplotlib.rcParams["axes.unicode_minus"]=False
d=pd.read_csv("epitope_recall.csv")

# ---------- (1) 지표별 분리: 케이스 4종을 지표마다 별도 그림 ----------
cases=[("8q7s_O","A 우세·견고","#08519c"),("8y6a_CD","B 비우세·얕게 rescue","#d95f0e"),
       ("8wpy_AB","B 비우세·깊은MSA 필수(절벽)","#238b45"),("8t4d_OQ","A·중간깊이 최적","#6a51a3")]
metrics=[("best_recall","recall (best pose)","cases_recall.png"),
         ("best_auprc","AUPRC (best pose)","cases_auprc.png"),
         ("best_mcc","MCC (best pose)","cases_mcc.png")]
for col,ylab,fn in metrics:
    fig,ax=plt.subplots(figsize=(7.2,4.4))
    for t,lab,c in cases:
        g=d[d.target==t].sort_values("rung")
        ax.plot(g.rung, g[col], "-o", color=c, ms=5, lw=1.8, label=f"{t} · {lab}")
    if col=="best_mcc": ax.axhline(0,color="grey",lw=.6,ls="--")
    ax.set_xlabel("rung  (0 = full MSA → 오른쪽 = 얕음/단일서열)")
    ax.set_ylabel(ylab); ax.set_ylim(-0.25 if col=="best_mcc" else -0.02, 1.05)
    ax.set_title(f"{ylab} — MSA 깊이에 따른 반응", fontsize=12)
    ax.legend(fontsize=8, loc="lower center", ncol=1)
    ax.grid(alpha=.25)
    plt.tight_layout(); plt.savefig(fn,dpi=140); print("saved",fn)

# ---------- (2) 히트맵: 깊이축 정규화(결측 제거) + 결론 정렬 ----------
GRID=np.linspace(0,1,11)   # 0=full → 1=최대축소, 모든 행 동일 11칸
def resample(col):
    rows={}
    for t,g in d.groupby("target"):
        g=g.sort_values("rung"); n=len(g)
        if n<2: continue
        frac=g.rung.values/(n-1)                    # 0..1 축소진행
        rows[t]=np.interp(GRID, frac, g[col].values)
    return pd.DataFrame(rows, index=[f"{x:.1f}" for x in GRID]).T   # index=target, cols=frac

meta=d.groupby("target").agg(ab=("ab","first"),group=("group","first"),
                             r0=("best_recall","first")).reset_index()
# rung0 = full 값을 정확히(첫 rung)
r0=d.sort_values("rung").groupby("target").best_recall.first()
meta["r0"]=meta.target.map(r0)
# 정렬: A → B → C, 그룹 내 full-MSA(r0) 내림차순 (결론: A는 위·밝음, B/C 아래)
meta=meta.sort_values(["ab","r0"],ascending=[True,False])
order=meta.target.tolist()

for col,ttl,fn,cmap,vlim in [("best_recall","best_recall","heatmap_recall.png","viridis",(0,1)),
                             ("best_mcc","best_MCC (음수=진짜자리 회피)","heatmap_mcc.png","RdBu_r",(-0.3,1))]:
    P=resample(col).reindex(order)
    fig,ax=plt.subplots(figsize=(7.5,11))
    im=ax.imshow(P.values,aspect="auto",cmap=cmap,vmin=vlim[0],vmax=vlim[1])
    ax.set_xticks(range(11)); ax.set_xticklabels([f"{x:.1f}" for x in GRID],fontsize=7)
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels([f"{t} [{meta.set_index('target').ab[t]}]" for t in order],fontsize=5.8)
    ax.set_xlabel("MSA 축소 진행  (0 = full → 1 = 최대축소/거의 단일서열)")
    ax.set_title(f"{ttl}\n위=A(우세) · 아래로 B(비우세)/C(대조), 그룹내 full-MSA 내림차순",fontsize=10)
    abs_=meta.ab.tolist()
    for j in range(1,len(abs_)):
        if abs_[j]!=abs_[j-1]: ax.axhline(j-0.5,color="k",lw=1.5)
    fig.colorbar(im,ax=ax,fraction=0.04,pad=0.03)
    plt.tight_layout(); plt.savefig(fn,dpi=140); print("saved",fn)