import pandas as pd, numpy as np
d = pd.read_csv("epitope_recall.csv")
print("총 복합체:", d.target.nunique(), " 행:", len(d))

# 1) 지표 합의 (pooled Pearson)
print("\n[1] 지표 상관(pooled):")
cols = ["best_recall","best_f1","best_mcc","best_auprc","mean_recall"]
print(d[cols].corr().round(3).to_string())

# 2) 복합체별 요약: 깊이-민감도(rung간 best_recall 범위), 수준, 폭(폭=recall>0.5 rung수), 지속성(mean이 best 따라가나)
def summ(g):
    br = g.best_recall.values; mr = g.mean_recall.values
    r0 = g.sort_values("rung").best_recall.values[0]
    return pd.Series(dict(
        n_true=int(g.n_true.iloc[0]), n_ag=int(g.n_ag.iloc[0]),
        max_br=br.max(), min_br=br.min(), rng=br.max()-br.min(),
        rung0=r0, best_over_rungs=br.max(),
        n_hit_br=int((br>=0.5).sum()),           # best_recall≥0.5 인 rung 수 (폭)
        max_mr=mr.max(),                          # mean_recall 최고 (지속성=여러 pose가 맞나)
        sustained=int((mr>=0.5).sum()),           # mean_recall≥0.5 rung 수 (재현성 있는 폭)
        mean_mcc=g.best_mcc.mean().round(3)))
s = d.groupby(["target","group","ab"]).apply(summ, include_groups=False).reset_index()

# 3) A/B/C 대비
print("\n[2] on-site(A) vs off-site(B) vs 대조(C) — 복합체 평균:")
agg = s.groupby("ab").agg(
    n=("target","size"),
    평균_rung0=("rung0","mean"),
    평균_best천장=("best_over_rungs","mean"),
    평균_깊이범위=("rng","mean"),
    평균_지속폭_mr50=("sustained","mean"),
    평균_스파이크폭_br50=("n_hit_br","mean")).round(2)
print(agg.to_string())

# 4) reduced-better(축소가 더 좋음): 얕은쪽 최고가 rung0보다 유의미하게↑ 인 복합체
s["reduced_better"] = (s.best_over_rungs - s.rung0 >= 0.3)
print("\n[3] '축소가 rung0보다 +0.30 이상' 복합체(ab별):")
print(s[s.reduced_better].groupby("ab").target.apply(list).to_string())

# 5) deep-required(깊은쪽 필요): rung0 높은데 천장이 rung0 근처(즉 얕으면 떨어짐) — 8wpy류
s["deep_req"] = (s.rung0>=0.8) & (s.min_br<=0.1)
print("\n[4] 'deep 높다가 얕으면 붕괴(min_br≤0.1)' 복합체:")
print(s[s.deep_req][["target","ab","rung0","min_br"]].to_string(index=False))

# 6) 스파이크 vs 지속: best는 높은데 mean 낮은(=best-of-5 운) 복합체 판별
s["spikey"] = (s.best_over_rungs>=0.7) & (s.max_mr<0.4)
print("\n[5] 스파이크형(best천장≥0.7 이나 mean천장<0.4 = 1개 pose 운):")
print(s[s.spikey][["target","ab","best_over_rungs","max_mr"]].to_string(index=False))

# 7) n_true 작은 복합체(recall 신뢰 낮음)
print("\n[6] n_true<8 (recall=1이 무의미할 수 있음):", s[s.n_true<8][["target","ab","n_true"]].to_dict("records"))

# 8) 견고한 on-site: 전 rung 높은
s["robust"] = (s.min_br>=0.6)
print("\n[7] 견고형(min_br≥0.6, 전 깊이 잘 맞음):")
print(s[s.robust].groupby("ab").target.apply(list).to_string())