import pandas as pd, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.font_manager as fm, matplotlib.pyplot as plt
fm.fontManager.addfont("/System/Library/Fonts/Supplemental/AppleGothic.ttf")
matplotlib.rcParams["font.family"]="AppleGothic"; matplotlib.rcParams["axes.unicode_minus"]=False
d=pd.read_csv("epitope_recall.csv")

# ---------- Fig1: 케이스별 지표 vs 깊이 (recall/AUPRC/MCC 일치 + 깊이반응) ----------
cases=[("8q7s_O","A 우세부위·견고"),("8y6a_CD","B 비우세·얕게 rescue"),
       ("8wpy_AB","B 비우세·깊은MSA 필수(절벽)"),("8t4d_OQ","A·중간깊이 최적")]
fig,ax=plt.subplots(2,2,figsize=(11,7)); ax=ax.ravel()
for i,(t,ttl) in enumerate(cases):
    g=d[d.target==t].sort_values("rung"); x=g.rung.values
    ax[i].plot(x,g.best_recall,"-o",color="#08519c",label="recall(best)",ms=4)
    ax[i].plot(x,g.best_auprc,"--s",color="#d95f0e",label="AUPRC(best)",ms=4)
    ax[i].plot(x,g.best_mcc,":^",color="#238b45",label="MCC(best)",ms=4)
    ax[i].plot(x,g.mean_recall,"-",color="#9ecae1",lw=1,alpha=.9,label="recall(mean=재현)")
    ax[i].axhline(0,color="grey",lw=.5); ax[i].set_ylim(-0.25,1.05)
    ax[i].set_title(f"{t}  ({ttl})",fontsize=10)
    ax[i].set_xlabel("rung  (0=full MSA → 오른쪽=얕음)"); ax[i].set_ylabel("지표값")
    n=g.neff80.values
    ax[i].set_xticks(x); ax[i].set_xticklabels([f"{r}\n{nf:g}" for r,nf in zip(x,n)],fontsize=6)
ax[0].legend(fontsize=7,loc="lower left",ncol=2)
fig.suptitle("케이스별: recall·AUPRC·MCC가 함께 움직이며 MSA 깊이에 반응 (x밑=Neff80)",fontsize=11)
plt.tight_layout(rect=[0,0,1,0.97]); plt.savefig("epitope_cases.png",dpi=140); print("saved epitope_cases.png")

# ---------- Fig2: 전체 히트맵 (best_recall & best_AUPRC) A/B/C 그룹 ----------
order=d.groupby("target").agg(ab=("ab","first"),group=("group","first"),lvl=("best_recall","mean")).reset_index()
order=order.sort_values(["ab","group","lvl"],ascending=[True,True,False])
tg=order.target.tolist()
def pivot(col):
    p=d.pivot_table(index="target",columns="rung",values=col).reindex(tg)
    return p.reindex(columns=range(0,11))
fig,axes=plt.subplots(1,2,figsize=(12,10))
for ax0,col,ttl in zip(axes,["best_recall","best_auprc"],["best_recall","best_AUPRC"]):
    P=pivot(col)
    im=ax0.imshow(P.values,aspect="auto",cmap="viridis",vmin=0,vmax=1)
    ax0.set_xticks(range(11)); ax0.set_xticklabels(range(11),fontsize=7)
    ax0.set_yticks(range(len(tg)))
    ax0.set_yticklabels([f"{t} [{a}]" for t,a in zip(tg,order.ab)],fontsize=5.5)
    ax0.set_xlabel("rung (0=full → 얕음)"); ax0.set_title(ttl,fontsize=11)
    # 그룹 경계선
    abs_=order.ab.tolist()
    for j in range(1,len(abs_)):
        if abs_[j]!=abs_[j-1]: ax0.axhline(j-0.5,color="red",lw=1.2)
    fig.colorbar(im,ax=ax0,fraction=0.046,pad=0.04)
fig.suptitle("전체 47복합체 × MSA깊이 — 위=A(우세) 중=B(비우세) 아래=C(대조), 빨간선=그룹경계",fontsize=11)
plt.tight_layout(rect=[0,0,1,0.98]); plt.savefig("epitope_heatmap.png",dpi=140); print("saved epitope_heatmap.png")