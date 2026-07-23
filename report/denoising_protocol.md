# MSA-깊이 rescue — 표집잡음 제거 프로토콜 (문헌 차용, 인용) — 2026-07-23

> 문제: rescue/gain을 **best-of-N(pose 5 × rung 11 = 55개 중 max)**로 정의해서, 실패 복합체도 우연히 threshold를
> 넘음(P≈0.7). 평균으로 보면 신호 사라짐. + rung마다 **MSA 추첨 1번**이라 "깊이 vs 조성" 미분리.
> 아래는 5각도 문헌탐색(워크플로 wvf2r40hq)이 수렴한 **차용 방법 4종 + 5단계 프로토콜.** 방법 novelty=0, 전부 인용.

## 우리와 같은 현상을 이름 붙인 선행
- **AFsample2 (Kalakoti & Wallner 2025, Commun Biol)**: MSA 마스킹이 모델 1개면 오히려 나쁨(DockQ 0.138 vs vanilla 0.163), N 늘려야 이득 = "평균 사라지고 max만 산다".
- **Porter et al. 2025 (J Mol Biol, "Does sequence clustering confound AF2?")**: AF-cluster 효과가 **같은 깊이 랜덤 추첨**으로 재현 → 깊이 주장은 matched-depth 랜덤 대조를 이겨야 함(= 우리 MSA-추첨-1번 문제).
- **Guan & Keating 2025 (Protein Science 34:e70331)**: 타깃별 짝 Wilcoxon + per-pose 성공확률(이항 CI) + top-1 vs oracle regret = 베낄 템플릿.

## 핵심 추정량 = pass@k (Chen et al. 2021, Codex, arXiv:2107.03374)
`success@k = 1 − C(n−c, k)/C(n, k)`  (n=pose 수, c=성공 pose 수, 성공 = epitope recall > θ=0.3)
- **pass@1 = c/n** = 뽑은 pose 하나가 성공할 확률 = **고정 예산이라 best-of-N이 못 부풀림.** 주 효과.
- 진단: 깊이 낮출 때 **pass@1이 오르면 진짜 rescue** / pass@55만 높고 pass@1 평평하면 순전히 shots-on-goal.
- naive `1−(1−p̂)^k`는 편향(복원추출 가정) → 위 조합식이 무편향(without-replacement MVUE).
- 안정 pass@5엔 n≈10–20 표본 필요(관례) → 지금 5는 부족.

## null (go/no-go) = matched best-of-N 순열검정 (Ojala & Garriga 2010; Porter 2025)
타깃 안 55 표본(11 rung × 5 pose) rung 라벨만 셔플(개수 유지) → **똑같은 max-over-rung rescue 재계산** 10,000×.
관측 rescue가 null 밴드 안이면 = 깊이 무효. **null도 반드시 best-of-N이어야 P≈0.7이 양쪽에 떠 상쇄.**
보완: 순서추세(Jonckheere–Terpstra / per-target Spearman→Stouffer); best-rung 점추정은 sample-splitting(Andrews–Kitagawa–McCloskey 2024).

## 정규화 (3축)
- x축 = **log(Neff)**(원 서열수 X; RBD~30 vs Env~3000 비교 위해 항원별 full-Neff 분수로 rung 인덱싱). (Jumper 2021)
- y축 = **full 대비 짝지은 Δrecall** = recall(rung)−recall(full). (AFsample2 improvement-over-no-mask)
- 타깃 pooling = per-target z-score(CASP outlier 규칙) 또는 within-target 분위수.
- 흡수기 = 혼합모델 `recall ~ log(depth) + (1|target) [+ (1|target:draw)]`(⚠️도킹 선행 없음, 표준 분산분해로 제시).

## 표집 계획 (유일한 새 GPU) — Step 4
del Alamo 2022는 깊이당 ~수십 draw. 최소 신뢰 floor:
**5 복합체(easy→borderline→impossible) × 3 깊이(full/mid/single) × 5 MSA추첨 × 15 pose = 1,125 fold** (GPU 1~3일).
GPU 부족 시: 3×2×5×10 = 300 fold. → 분산분해(깊이/추첨/pose) + **재현성 필터(del Alamo): off-site 에피토프가 독립 추첨 ≥2번 재현될 때만 rescue 인정.**

## 5단계
0. **pose별 recall 5개 재생성**(지금 best/mean/min만) — 저장 pose 재채점, CPU. (Chen 2021 요건)
1. **pass@1 ± Wilson CI vs log(Neff)** 그림. pass@1 안 오르면 헤드라인 = "rescue=best-of-N 착시".
2. **순열 null**(rung 라벨 셔플, max-over-rung 재계산) + BH-FDR. **= go/no-go.**
3. **혼합모델 + 짝 Wilcoxon** (null이면 TOST 등가로 "±0.05 내 무효" 적극 진술; Lakens 2017).
4. **표적 multi-draw**(위) — 분산분해 + 재현성 필터 + 교차적합.
5. 최종 패널: rung별 **oracle − ipTM-top1 격차** → "낮은 깊이서 진짜 자리 만들지만 ipTM이 못 고름" = Phase-1 재랭커로 연결(AFsample; Hitawala–Gray 2025).

## 유일한 함정
best-of-N rescue를 **best-of-N이 아닌 null(고정 threshold·평균·per-sample)과 비교 금지** — 항상 rescue처럼 보임.
Step 2 null은 동일 max-over-rung을 셔플 데이터로 재계산해야 상쇄. 보고 효과는 pass@1(고정예산).

**출처:** 워크플로 wvf2r40hq(6에이전트, 5각도 수렴). Chen2021·Ojala2010·Porter2025·delAlamo2022·Wayment-Steele2024·Wallner2023·Kalakoti&Wallner2025·Guan&Keating2025·Jumper2021·Truchon&Bayly2007·Andrews–Kitagawa–McCloskey2024·Lakens2017.
