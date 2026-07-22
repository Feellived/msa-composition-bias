# MSA 조성 통제 실험 계획 — "깊이(개수) vs 특정 서열(조성)" 인과 분리 (v1, 2026-07-22)

> 목적: depth-response가 **MSA의 양(깊이/개수)** 때문인지 **특정 동족서열(조성)** 때문인지를 통제 실험으로 가른다.
> 후자면 "어떤 서열(그 서열의 공변 잔기)이 배치를 좌우한다"는 **mutation/coevolution insight** → DMS 대비 사전선별(랩 차별화)로 연결.
> ⚠️ 현재 EDA(`boltz_depth_sweep_analysis.md`)의 사다리는 **rung마다 독립 랜덤 subsample**이라 "개수 vs 조성"을 못 가른다 — 아래는 그걸 가르는 재설계.

---

## 0. 모델 선택 — 반응성 최고 1개만 (앞으로 통제 실험 전부 이 모델로)
- Boltz·Chai·**Protenix**(현재 스윕 중) 3모델 depth-sweep 완주 후, **가장 깊이-반응 좋은 모델 1개** 선택. 통제 실험은 **그 하나만** (compute 절약, 사용자 지시).
- **선택 기준(정량):** ① 지속형(sustained) 깊이-전이 수 ② 반응 복합체 평균 진폭 ③ 단발 스파이크 비율(낮을수록 좋음) ④ sweet-spot(reduced>full) 복합체 수.
- **예상:** Protenix가 MSA 강의존(consensus data: core full 0.9→depth 0.03 붕괴)이라 반응성 최고 가능성. Boltz 현재 = sweet-spot 후보 14 + deep-required 1(8wpy)이나 대부분 스파이크.

## 대상 복합체
선택 모델에서 **반응 복합체**(sweet-spot + deep-required) 중 가장 깨끗·큰 것 **3~5개**. (Boltz면 8wpy[deep-required 절벽], 8t4d[mid-peak], 8y6a/9y0a/8txu[reduced-better] 등.)

---

## 실험 순서

### Exp1 — Seed 복제 (깊이 vs 조성 1차 필터, 제일 중요)
각 대상 복합체의 **전이 깊이(개수)에서, 같은 서열 개수를 다른 seed로 N=10회** 독립 subsample → DockQ 분포.
- **분산 작음(다 비슷)** → **깊이(개수)**가 원인. 조성 무관 → 그 복합체는 조성 실험 불필요.
- **분산 큼(어떤 draw는 높고 어떤 건 낮음)** → **조성(어느 서열)**이 원인 → Exp2·3 진행.
- **이게 8wpy 절벽이 "18서열이라서"인지 "그 랜덤 draw가 핵심서열을 놓쳐서"인지 가르는 결정적 검정.**

### Exp2 — Nested 사다리 (전이=제거된 특정 서열)
사다리를 **nested**로 재빌드: rung_{k+1} ⊂ rung_k (이전 rung에서 서열을 **빼기만**, 재추첨 X). → rung→rung DockQ 변화 = **그 단계에서 빠진 특정 서열군**에 귀속. 전이가 나는 단계의 제거 서열 = 인과 후보.
- 구현: build_ladder에 `--nested` 모드(누적 제거) 추가.

### Exp3 — 인과 서열 식별 (leave-one-out / add-one-in, 클러스터 단위)
full MSA를 **동일성(예 62%) 또는 taxonomy로 M개 클러스터**로 묶고:
- **LOCO(leave-one-cluster-out):** full − cluster_i → DockQ. 큰 하락 = cluster_i **필요**(제거하면 무너짐).
- **AOCI(add-one-cluster-in):** single-seq + cluster_i → DockQ. 큰 rescue = cluster_i가 **정답 유도**.
- 교차 확인: "인과 클러스터" 제거가 깊이 전이를 재현하나? 그 추가가 single-seq를 rescue하나?
- 클러스터 단위라 compute 유한(서열별 아님).

### Exp4 — 인과 클러스터 특성화 → mutation insight
인과 cluster의: ① taxa/생물종 ② query와의 %동일성 ③ 커버하는 **MSA 열(column) = 에피토프/paratope 잔기인가** ④ 그 열의 **공변(APC-보정 MI)이 cluster 유무로 바뀌나**. → **어느 잔기의 공진화가 배치를 결정하나** = mutation/DMS 표적.

## Readout
**DockQ(배치 품질) + epitope recall(항체가 진짜 자리로 이동하나)** 이원. (2a의 "깊을 때 과대표집 자리로 끌림"도 여기서 pose 파싱으로 확인.)

## Compute (유한, 1모델)
- Exp1: 5타깃 × 10 seed × 5 sample = 250 run.
- Exp3: 5타깃 × ~15 cluster × 2(LOCO+AOCI) × 5 sample ≈ 750 run.
- 1 GPU 청크씩. 클러스터 단위 + 5타깃으로 유계.

## 정직한 한계
- 클러스터 granularity가 해상도 좌우(너무 크면 뭉뚱, 너무 잘면 compute↑).
- MSA subsample 확률성 → Exp1(seed 복제)을 반드시 먼저.
- 1모델 → 일반화 제한(명시). Protenix에서 되면 Boltz 교차확인.
- epitope recall = pose 파싱 필요(서버).

## 다음 착수
Protenix 스윕 완주 → **모델 반응성 랭킹(Exp0)** → 모델·타깃 확정 → **Exp1(seed 복제)** 부터. build_ladder `--nested`·LOCO/AOCI 서브샘플러·epitope recall 파서 = 구현 대상.
