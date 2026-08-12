# 에피토프 예측 품질(위치) vs MSA 깊이 — Boltz (v1, 2026-07-22)

> DockQ(fold+배치 혼합 스칼라)와 달리 **"항체가 진짜 에피토프 자리에 붙었나"(위치)만** 재는 지표.
> pose를 native 에피토프(항체 5Å 내 항원 잔기, 참조 항원서열 위치)와 비교. 47 복합체 × 최대 11 rung × 5 pose.
> 데이터: `report/epitope_recall.csv` · 그림: `report/figures/epitope_recall_ABC.png` · 코드: `pipeline/lib_epitope_recall.py`.

## 지표 (recall 단독 불충분 → 세트로)
`recall(best/mean/min)` · `precision` · `F1` · **`MCC`**(불균형 강건 = 정직한 accuracy) · **`AUPRC`**(threshold-free, 잔기별 항체거리 점수).
근거: recall만 보면 표면 전체를 덮는 pose가 100%로 뜸(특이성 0). accuracy는 에피토프=소수라 항상 높게 나와 오도. ARID-sf류가 recall 대신 정확도+AUROC를 같이 보는 이유와 동일(우리는 희소양성이라 AUROC 대신 AUPRC).

## 결과 1 — 지표 4종이 거의 완전 일치
pooled 상관: best_recall↔F1 0.93, ↔MCC 0.95, ↔AUPRC 0.91 / F1↔MCC 0.99, F1↔AUPRC 0.98.
→ 대표값으로 **F1·MCC** 하나만 봐도 됨. 단 어긋나는 예외를 정확히 잡음:
- **8vye_AD (n_true=3)**: recall=1.0 인데 F1 0.26·MCC 0.37 → 큰 patch가 우연히 3잔기 포함(특이성 없음). recall만이면 "완벽"으로 오판. **recall 단독의 위험을 보여주는 대표 사례**(+ n_true=3은 데이터 품질도 재확인 대상).

## 결과 2 — on-site / off-site / 대조가 full MSA에서 갈린다 (핵심)
| 그룹 | full(rung0) best_recall | 깊이축소 천장(best) | 깊이-이동폭(max−min) | 재현폭(mean_recall≥0.5 rung수) |
|---|---|---|---|---|
| **A 우세부위(on-site)** (n=20) | **0.72** | 0.89 | 0.34 | 6.0 |
| **B 비우세(off-site)** (n=18) | **0.42** | 0.66 | **0.41** | 3.6 |
| **C 대조(비과대표집)** (n=9) | 0.30 | 0.64 | 0.38 | 2.2 |

- **on-site는 full MSA에서 이미 잘 붙지만(0.72), off-site는 full에서 못 붙고(0.42) 깊이를 낮춰야 천장(0.66)이 열림.** off-site가 깊이에 가장 크게 움직임(이동폭 0.41).
- **같은 항원(RBD) 안에서** A(8q7s_O·8kep·8k5g·8k5h = RBM) 전부 견고, B(8y6a·8wpy·8xsj·8k46·8k3k = off-RBM) 전부 취약 → **항원 난이도가 아니라 "에피토프 위치(우세/비우세)"가 깊이-취약성을 결정.**
- 해석: "깊은 MSA가 과대표집(우세) 자리로 끌고 비우세는 억눌린다"는 위치편향의 **A/B 대비 정량화** (DockQ 스칼라론 안 보였던 구조).
- 견고형(전 깊이 min_recall≥0.6): 거의 다 A(on-site). off-site·대조는 소수.

## 결과 3 — 냉정한 한계 (앞 DockQ 분석과 일치)
- **rescue 대부분 스파이크형**: best천장≥0.7 이나 mean천장<0.4 인 복합체 8개(8y6a·8q7s_H·8txu·9y0a·8t4d·8t4b·8kdm·8ume) = **5 pose 중 1개 운, 재현 안 됨.** (그림 ③ 대각선 아래 대량 분포.)
- **방향 양쪽**: 8wpy·8k46·8t49 = deep 필요(얕으면 붕괴), 8y6a·8xsj = shallow 필요. → **깊이(양)가 아니라 조성(어느 서열)** 이라는 결론 재확인.

## 함의 → 다음
1. **A/B 위치-취약성 대비**는 재현성 있는 그룹-수준 신호(스파이크와 별개) → 발표 headline 가능.
2. rescue의 스파이크성 → **seed-복제 통제실험**(같은 개수 다른 draw N회)로 "깊이 vs 조성" 확정 필요성 재확인.
3. Protenix/Chai 스윕 완료 후 같은 지표로 재계산 → 모델 반응성 랭킹(Exp0) → 통제실험 1모델 선택.

## 재현
```bash
python lib_epitope_recall.py --models boltz --rungs 12    # → results/epitope_recall.csv
# 분석·그림: report 폴더 analyze.py / fig.py 참조 (pandas+matplotlib)
```