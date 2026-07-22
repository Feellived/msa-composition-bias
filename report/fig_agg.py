import pandas as pd, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.font_manager as fm, matplotlib.pyplot as plt
fm.fontManager.addfont("/System/Library/Fonts/Supplemental/AppleGothic.ttf")
matplotlib.rcParams["font.family"]="AppleGothic"; matplotlib.rcParams["axes.unicode_minus"]=False
d=pd.read_csv("epitope_recall.csv")
d=d[d.rung<=9]                                  # rung0~9 = 47복합체 전부 존재(결측 0)
col={"A":"#2c7fb8","B":"#d95f0e","C":"#7a7a7a"}
lab={"A":"A 우세부위(on-site)","B":"B 비우세(off-site)","C":"C 대조(비과대표집)"}

metrics=[("best_recall","recall (best pose)"),("best_auprc","AUPRC (best pose)"),("best_mcc","MCC (best pose)")]

# ---------- (1) 꺾은선: 그룹별 rung당 평균 ± SEM ----------
fig,axes=plt.subplots(1,3,figsize=(14,4.4))
for ax,(mc,ylab) in zip(axes,metrics):
    for ab in ["A","B","C"]:
        g=d[d.ab==ab].groupby("rung")[mc]
        mean=g.mean(); sem=g.std()/np.sqrt(g.count()); n=int(d[d.ab==ab].target.nunique())
        ax.plot(mean.index,mean.values,"-o",color=col[ab],ms=4,label=f"{lab[ab]} (n={n})")
        ax.fill_between(mean.index,mean-sem,mean+sem,color=col[ab],alpha=.18)
    if mc=="best_mcc": ax.axhline(0,color="grey",lw=.6,ls="--")
    ax.set_xlabel("rung  (0 = full MSA → 9 = 얕음)"); ax.set_ylabel(ylab)
    ax.set_title(ylab,fontsize=11); ax.grid(alpha=.25)
    ax.set_ylim(-0.15 if mc=="best_mcc" else 0,1.0)
axes[0].legend(fontsize=8,loc="upper right")
fig.suptitle("그룹별 rung당 평균 ± SEM — A/B/C 분리는 크나 깊이에 따른 그룹-평균 변화는 작음(개별 rescue는 서로 다른 깊이의 스파이크라 평균서 상쇄)",fontsize=10)
plt.tight_layout(rect=[0,0,1,0.95]); plt.savefig("agg_lines.png",dpi=140); print("saved agg_lines.png")

# ---------- (2) 막대: 그룹×깊이대역(full/중간/얕음) 평균 ----------
d["band"]=pd.cut(d.rung,[-1,2,6,9],labels=["full(r0-2)","중간(r3-6)","얕음(r7-9)"])
fig,axes=plt.subplots(1,3,figsize=(14,4.2))
bands=["full(r0-2)","중간(r3-6)","얕음(r7-9)"]; x=np.arange(3); w=0.26
for ax,(mc,ylab) in zip(axes,metrics):
    for k,ab in enumerate(["A","B","C"]):
        vals=[d[(d.ab==ab)&(d.band==b)][mc].mean() for b in bands]
        errs=[d[(d.ab==ab)&(d.band==b)][mc].sem() for b in bands]
        ax.bar(x+(k-1)*w,vals,w,yerr=errs,capsize=2,color=col[ab],label=lab[ab] if mc=="best_recall" else None)
    if mc=="best_mcc": ax.axhline(0,color="grey",lw=.6)
    ax.set_xticks(x); ax.set_xticklabels(bands,fontsize=8); ax.set_ylabel(ylab); ax.set_title(ylab,fontsize=11)
    ax.set_ylim(-0.1 if mc=="best_mcc" else 0,1.0)
axes[0].legend(fontsize=8)
fig.suptitle("그룹 × 깊이대역 평균 (막대) — 전 대역에서 A > B > C, 대역 간 변화는 작음",fontsize=10)
plt.tight_layout(rect=[0,0,1,0.95]); plt.savefig("agg_bars.png",dpi=140); print("saved agg_bars.png")

# 콘솔 요약
print("\n그룹별 전체 평균(rung0-9):")
for mc,yl in metrics:
    print(f"  {yl:20}",{ab:round(d[d.ab==ab][mc].mean(),3) for ab in ["A","B","C"]})