import pandas as pd, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.font_manager as fm, matplotlib.pyplot as plt
fm.fontManager.addfont("/System/Library/Fonts/Supplemental/AppleGothic.ttf")
matplotlib.rcParams["font.family"]="AppleGothic"; matplotlib.rcParams["axes.unicode_minus"]=False
rng=np.random.RandomState(42)
d=pd.read_csv("dockq_sweep_boltz.csv")
col={"A":"#2c7fb8","B":"#d95f0e","C":"#7a7a7a"}; lab={"A":"A 우세부위","B":"B 비우세","C":"C 대조"}

# ---------- (1) 그룹평균 꺾은선 (부트스트랩 95% CI) + DockQ 임계선 ----------
dd=d[d.rung<=9]
def boot(vals,B=2000):
    if len(vals)<2: return vals.mean(), vals.mean(), vals.mean()
    m=[rng.choice(vals,len(vals),replace=True).mean() for _ in range(B)]
    return vals.mean(), np.percentile(m,2.5), np.percentile(m,97.5)
fig,ax=plt.subplots(figsize=(8,5))
for ab in ["A","B","C"]:
    g=dd[dd.ab==ab]; n=g.target.nunique()
    xs=sorted(g.rung.unique()); mean=[];lo=[];hi=[]
    for r in xs:
        m,l,h=boot(g[g.rung==r].best_dockq.values); mean.append(m);lo.append(l);hi.append(h)
    ax.plot(xs,mean,"-o",color=col[ab],ms=5,label=f"{lab[ab]} (n={n})")
    ax.fill_between(xs,lo,hi,color=col[ab],alpha=.18)
for y,t in [(0.23,"acceptable 0.23"),(0.49,"medium 0.49"),(0.80,"high 0.80")]:
    ax.axhline(y,color="grey",lw=.7,ls="--"); ax.text(9.05,y,t,fontsize=7,va="center",color="grey")
ax.set_xlabel("rung  (0 = full MSA → 9 = 얕음)"); ax.set_ylabel("DockQ (best-of-5 pose)")
ax.set_title("Boltz MSA depth-sweep — 그룹별 DockQ 평균 (부트스트랩 95% CI)",fontsize=11)
ax.set_ylim(0,1); ax.legend(fontsize=9,loc="upper right"); ax.grid(alpha=.25)
plt.tight_layout(); plt.savefig("dockq_agg_lines.png",dpi=140); print("saved dockq_agg_lines.png")

# ---------- (2) 복합체별 요약 통계 ----------
rows=[]
for t,g in d.groupby("target"):
    g=g.sort_values("rung"); v=g.best_dockq.values
    rows.append(dict(target=t, 항원=g.group.iloc[0], ab=g.ab.iloc[0], rung수=len(v),
                     평균=round(v.mean(),3), 최소=round(v.min(),3), 최대=round(v.max(),3),
                     변화량=round(v.max()-v.min(),3), full_r0=round(v[0],3),
                     최고=round(v.max(),3), 최고rung=int(g.rung.values[v.argmax()])))
S=pd.DataFrame(rows).sort_values(["ab","항원","target"])
S.to_csv("dockq_summary.csv",index=False)
# 그룹 전체 평균
print("\n그룹 평균 DockQ:",{ab:round(d[d.ab==ab].best_dockq.mean(),3) for ab in ["A","B","C"]})
print("성공률(전 pose,rung 중 best):")
for th in [0.23,0.49,0.80]:
    print(f"  DockQ≥{th}: 복합체수 =",{ab:int((S[S.ab==ab].최대>=th).sum()) for ab in ["A","B","C"]},"/ (A20 B18 C9)")
print(f"\ndockq_summary.csv ({len(S)}행)  columns={list(S.columns)}")
# wide (raw matrix) 도 저장 → 첨부용
W=d.pivot_table(index=["target","group","ab"],columns="rung",values="best_dockq").reset_index()
W.to_csv("dockq_wide.csv",index=False); print(f"dockq_wide.csv ({len(W)}행 × rung {list(d.rung.unique())[:3]}...)")