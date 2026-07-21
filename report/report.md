# 보고서 — 항원 MSA 깊이 편향과 depth-response 기반 pose 선택

> 작업 중 초안. 계획서 = [`../plan/research_plan.md`](../plan/research_plan.md).

## 1. 배경
- (문제) 항체–항원 도킹의 병목 + 흔한 자리 편향.
- (동기) 신규·비면역우세 에피토프 항체가 더 틀린다.

## 2. 방법
- 데이터셋 3그룹(A/B/C), leakage-free 2023-06 이후, 49 복합체.
- 생성: Boltz-2·Chai-1·Protenix-base + HADDOCK·SnugDock × Neff80 사다리 × seed.
- 지표: merged DockQ(0.23/0.49/0.80) · epitope recall · over-rep overlap · 항체 내부 RMSD.

## 3. 결과

### 3.1 편향의 존재 (over-representation map)
- (예정) 항원별 잔기 접촉 빈도 지도.

### 3.2 depth-response — rescue 특이성 (E1)
- (예정) 그룹별 depth–DockQ 곡선. B rescue / A 악화 / C 무변.

### 3.3 sweet-spot (E2)
- (예정) 그룹·모델별 최적 깊이 분포.

### 3.4 커버리지 (E3)
- (예정) depth-sweep oracle vs single-MSA, B 집중.

### 3.5 선택 (E4)
- (예정) 기성 재랭커 × enriched vs single 풀.

### 3.6 편향 센서 (E5)
- (예정) 예측 자리 이탈 신호 ↔ 그룹 라벨 상관.

## 4. 논의
- 정직한 위치(방법 선점, 기여=결합·적용·wet).
- 한계(표본 ~49, case study).

## 5. 참고문헌
- 계획서 §참고문헌 참조.
