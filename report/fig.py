import pandas as pd, numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
fm.fontManager.addfont("/System/Library/Fonts/Supplemental/AppleGothic.ttf")
matplotlib.rcParams["font.family"]="AppleGothic"
matplotlib.rcParams["axes.unicode_minus"]=False
d = pd.read_csv("epitope_recall.csv")

def per(g):
    br=g.sort_values("rung"); 
    return pd.Series(dict(rung0=br.best_recall.values[0], ceil=g.best_recall.max(),
                          mean_ceil=g.mean_recall.max(), rng=g.best_recall.max()-g.best_recall.min(),
                          ab=g.ab.iloc[0]))
s=d.groupby("target").apply(per,include_groups=False).reset_index()
order=["A","B","C"]; lab={"A":"A 우세부위\n(on-site)","B":"B 비우세\n(off-site)","C":"C 대조\n(비과대표집)"}
col={"A":"#2c7fb8","B":"#d95f0e","C":"#999999"}

fig,ax=plt.subplots(1,3,figsize=(13,4.2))
# P1: full(rung0) vs 축소천장 vs 재현천장 — 그룹 평균
x=np.arange(3); w=0.26
for i,(k,c) in enumerate(zip(["rung0","ceil","mean_ceil"],["#bdd7e7","#6baed6","#08519c"])):
    vals=[s[s.ab==a][k].mean() for a in order]
    ax[0].bar(x+(i-1)*w, vals, w, color=c, label={"rung0":"full MSA(rung0)","ceil":"깊이축소 천장(best)","mean_ceil":"재현 천장(mean)"}[k])
ax[0].set_xticks(x); ax[0].set_xticklabels([lab[a] for a in order]); ax[0].set_ylabel("epitope best_recall")
ax[0].set_title("① full에선 off-site가 낮고,\n축소하면 천장 열림 (단 재현천장은 낮음)",fontsize=10)
ax[0].legend(fontsize=7,loc="upper right"); ax[0].set_ylim(0,1)
# P2: 복합체별 깊이-이동폭 strip
for a in order:
    v=s[s.ab==a].rng.values; xj=order.index(a)+np.random.RandomState(0).uniform(-.12,.12,len(v))
    ax[1].scatter(xj,v,color=col[a],s=42,alpha=.8,edgecolor="k",linewidth=.4)
    ax[1].plot([order.index(a)-.2,order.index(a)+.2],[v.mean()]*2,color="k",lw=2)
ax[1].set_xticks(range(3)); ax[1].set_xticklabels([lab[a] for a in order]); ax[1].set_ylabel("깊이-이동폭 (max-min best_recall)")
ax[1].set_title("② off-site가 깊이에 가장 크게 움직임",fontsize=10); ax[1].set_ylim(0,1)
# P3: best천장 vs mean천장 = 스파이크성(대각선 아래일수록 1-pose 운)
for a in order:
    sub=s[s.ab==a]; ax[2].scatter(sub.ceil,sub.mean_ceil,color=col[a],s=42,alpha=.8,edgecolor="k",linewidth=.4,label=lab[a].split("\n")[0])
ax[2].plot([0,1],[0,1],"k--",lw=.8); ax[2].fill_between([0,1],[0,0],[0,1],color="grey",alpha=.06)
ax[2].set_xlabel("best_recall 천장 (아무 pose나)"); ax[2].set_ylabel("mean_recall 천장 (재현)")
ax[2].set_title("③ 대각선서 멀수록 best-of-5 운\n(off-site rescue 대부분 여기)",fontsize=10)
ax[2].legend(fontsize=7,loc="lower right"); ax[2].set_xlim(0,1); ax[2].set_ylim(0,1)
plt.tight_layout(); plt.savefig("epitope_recall_ABC.png",dpi=140)
print("saved epitope_recall_ABC.png")