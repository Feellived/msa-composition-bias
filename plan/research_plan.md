# 항원 MSA 깊이 편향의 진단과 depth-response 기반 항체–항원 도킹 pose 선택

*SNU BK21 여름 학부연구 2026 · 연구 계획안 (v1, 2026-07-20)*

---

## 초록

딥러닝 기반 복합체 구조예측 모델(co-folder; AlphaFold3 및 그 오픈 재현판)은 항체–항원 도킹에서 뛰어난 성능을 보이나, 표적에 따라 정확도가 크게 흔들리고 특히 신규·비면역우세 에피토프에서 실패가 잦다. 본 연구는 이 실패의 상당 부분이 **모델이 학습 데이터에 과대표집된 에피토프로 예측을 쏠리게 하는 편향**에서 비롯되며, 그 편향이 **항원 다중서열정렬(MSA)의 깊이**를 통해 전달된다는 가설을 검증한다. 핵심 실험은 동일한 항체–항원 쌍에 대해 항원 MSA 깊이를 full에서 single-sequence까지 log 스케일로 낮추며(depth-sweep) 예측 pose의 궤적(depth-response)을 관측하는 것이다. 우리는 (1) 과대표집 에피토프가 정답이 아닌 항체에서 MSA 감소가 올바른 pose를 되살리고(rescue), (2) 과도한 감소는 공진화 신호 붕괴로 pose 품질을 떨어뜨려 최적 깊이(sweet-spot)가 존재하며, (3) depth-response 자체가 편향을 감지하고 near-native pose를 선택하는 데 유효한 신호임을 보인다. 결과물은 depth-sweep 데이터셋과, 그 위에서 기성 재랭커(DeepRank-Ab / ARID-sf)가 단일-MSA 풀보다 더 나은 pose를 고르게 만드는 pipeline이며, 이는 향후 항체–항원 재랭커·co-folder 학습에 넣을 수 있는 새로운 특징(feature)의 근거를 제공한다.

---

## 1. 배경

### 1.1 문제 정의: 도킹은 되지만 선택이 안 된다

In silico 항체 신약개발은 수백~수만 개의 항체 후보를 생성할 수 있으나, 모든 후보에 실제 결합 assay·epitope binning·구조 분석을 수행하는 것은 비현실적이다. 따라서 실험 진입 전, 항체가 항원의 **어느 에피토프에 어떤 방향·형태로 결합하는지**를 높은 정확도로 예측해 후보를 선별해야 한다. 신뢰할 수 있는 항체–항원 복합체 구조는 에피토프·선택성 분석, 상호작용 잔기 규명, epitope binning 가설, affinity maturation 후보 제안, escape mutation 평가의 구조적 근거가 되며, 초기 pose가 실제와 크게 다르면 이후 모든 계산이 잘못된 구조 위에서 이루어진다.

AlphaFold3와 그 오픈 재현판(Boltz, Chai, Protenix 등)은 항체–항원 복합체 예측을 크게 향상시켰으나 여전히 취약 범주에 속한다. 보고된 성공률은 seed 수·평가 기준에 크게 의존하며, single-seed 기준 acceptable(DockQ≥0.23) 성공률은 대략 30–35%, high-accuracy는 약 10% 수준으로 알려져 있다. 예비 분석에서 우리는 두 가지 실패 양식을 관찰했다. 첫째, 여러 co-folder의 자기신뢰 점수(ipTM)와 물리 도킹 점수(HADDOCK score)가 **자신이 생성한 near-native pose를 선택하지 못한다**(생성은 되는데 선택은 실패한다). 둘째, 실패의 상당수는 fold 붕괴가 아니라 **배치(placement) 오류** — 즉 항체가 잘못된 에피토프로 놓이는 것이다.

### 1.2 편향의 원인: 학습 데이터 과대표집과 MSA 채널

특정 항원에 방대한 데이터가 축적되어 있어도, 그 데이터는 흔히 일부 잘 알려진(면역우세) 에피토프에 편중된다. 예컨대 SARS-CoV-2 수용체 결합 도메인(RBD)은 공개 구조가 매우 많지만 항체 복합체는 수용체 결합 모티프(RBM, 잔기 437–508) 주변에 몰려 있다. 새로운 항체가 이 과대표집 부위 밖(cryptic·측면 보존 에피토프)에 결합하는 경우, 모델은 학습에서 자주 관찰된 결합 양식으로 치우친 pose를 제시하기 쉽다. 즉 데이터의 양이 많다고 새 에피토프 예측이 정확한 것은 아니며, 오히려 과대표집이 잘못된 방향으로 예측을 끌어당긴다.

이 편향이 전달되는 주요 통로가 **항원 MSA의 깊이와 조성**이다. 깊은 MSA는 강한 공진화 신호를 제공하지만, 그 신호는 학습 분포에 과대표집된 결합 양식의 통계를 함께 담고 있어 과대표집 부위로의 인력(attraction)을 강화한다. 반대로 MSA를 얕게 하면 그 인력이 약해지지만, 동시에 정확한 접힘·배치에 필요한 공진화 정보도 잃는다.

### 1.3 왜 도킹에서 이 편향이 특히 중요한가

MSA 깊이 조절은 기존에도 연구되었으나, 대개 **단백질의 대안 conformation을 샘플링**하기 위한 목적이었다(§2). 항체–항원 도킹은 결이 다르다. 여기서 잘못된 에피토프 예측은 단순한 정확도 저하를 넘어, 이후 binding energy 계산·interface 분석·mutation 제안·epitope binning 등 **downstream wet 실험 설계 전체를 오도**한다. 따라서 도킹 맥락에서 MSA 편향은 "다양성"의 문제가 아니라 "정확한 에피토프 탐색"의 문제이며, 이를 진단하고 완화하는 것이 실용적 가치를 가진다.

---

## 2. 선행 연구

### 2.1 co-folder의 항체–항원 예측 한계: 배치·선택의 문제

AlphaFold3(Abramson et al., 2024)를 비롯한 co-folder는 성능이 seed 수에 강하게 의존하며, 저자들은 항체–항원을 seed 예산에 특히 민감한 별도 범주로 취급한다. Hitawala & Gray(2025, mAbs)는 단일 seed 기준 high-accuracy 도킹 성공률이 항체 약 10%, 나노바디 약 13%에 불과하고 단일 seed 실패율이 약 60%임을 보고했으며, 이는 AF3 원 논문의 1000-seed 조건 ~60% 성공과 크게 대비된다. 성능차의 본질은 항체의 fold 품질이 아니다: 이들은 항체 사슬 내부 배향 confidence는 높으나(inter-chain ipTM ≈ 0.8) 계면 confidence는 낮고(평균 ipTM ≈ 0.3), 결정적으로 **처음에 잘못 도킹된 항체도 에피토프 위치를 제공하면 올바르게 도킹**됨을 보였다. 즉 실패의 상당 부분은 접힘이 아니라 **에피토프 오배치(placement)**이며, 이는 본 연구가 겨냥하는 대상과 정확히 일치한다. 다중 모델 벤치마크 FoldBench(2025, Nat Commun) 역시 항체–항원 계면을 벤치 내 가장 어려운 범주 중 하나로 규정하며, AF3가 DockQ≥0.23 성공률 45.4%로 최고이나 나머지 co-folder는 평균 약 30%, 다섯 모델 모두 실패율이 50%를 넘음을 정량화했다.

이 한계의 한 축은 **자기신뢰 점수(confidence)가 올바른 결합을 담보하지 못한다**는 데 있다. Yin & Pierce(2024, Protein Science)는 429개 복합체에서 AlphaFold의 confidence가 최상위 pose를 신뢰성 있게 고르지 못함을 보였고, 최근 연구(bioRxiv 2026; 106개 나노바디 복합체와 11,342개 비-동족 쌍에 대해 AF3·Boltz-2·Chai-1 평가)는 세 모델 모두 기하학적으로 그럴듯한 복합체를 자주 생성하나 ipTM이 동족(정답)과 비-동족(오답)을 구분하지 못하며, 샘플링을 늘려도 기하 품질만 개선될 뿐 판별력은 정체됨을 직접 입증했다. Ibex(Dreyer et al., 2025)는 seed를 1000까지 늘려도 항체 loop 예측 향상이 거의 없음을 보여, 단순한 샘플 증가로는 이 문제가 풀리지 않음을 시사한다.

### 2.2 MSA 깊이 조절: 대안 conformation을 위한 도구

MSA 깊이를 낮추면 co-folder가 단일 basin이 아니라 여러 구조를 생성한다는 사실은 잘 확립되어 있다. del Alamo et al.(2022, eLife)은 AF2의 내부 파라미터(max_seq·extra_seq)를 16에서 5120까지 조절해 막단백질·수용체의 대안 conformation을 샘플링했고, Monteiro da Silva et al.(2024, Nat Commun)은 얕은 subsampling으로 conformation 간 상대 population을 80% 이상 정확도로 재현했다. AF-Cluster(Wayment-Steele et al., 2024, Nature)는 MSA를 진화적 subfamily로 클러스터링하여(예: KaiB) metamorphic 단백질의 두 fold를 모두 재현했으며, AFsample(Wallner, 2023)은 MSA를 건드리지 않고 대량 sampling과 dropout으로, AFsample2(Kalakoti & Wallner, 2025)는 MSA column masking으로 공진화 신호를 희석해 다양성을 얻었다. **그러나 이 방법들의 목적은 예외 없이 단량체·다량체의 conformational diversity 확보이지, 항체–항원 에피토프 편향의 진단이 아니다.**

깊이–정확도 관계는 선형이 아니라 포화형이며, 유효 서열 수 Neff가 대략 30 미만으로 떨어지면 정확도가 급감한다(NEFFy; Rajabi et al., 2026). 따라서 depth를 원자료 서열 수가 아니라 Neff로 정규화하는 것이 표적 간 공정 비교에 필요하다. 결정적으로 McCoy et al.(2024, Protein Science)은 엄격히 필터링된 항체–항원 벤치마크에서 **Neff와 DockQ 사이에 유의한 상관이 없음**을 보고했다. 저자들은 항원이 항체와 공진화하지 않고 오히려 회피하도록 진화하기 때문이라 해석하는데, 이는 "깊은 항원 MSA가 항체–항원 도킹에 반드시 유익하지는 않다"는 본 연구 가설의 발판이 된다. 같은 연구는 정확한 예측일수록 PDB에 흔한 계면 기하를, 부정확한 예측일수록 드문 기하를 보이는 **체계적 편향**을 정량화하여, 모델이 과대표집된 흔한 계면으로 예측을 쏠리게 함을 시사했다.

### 2.3 항체–항원 pose 재랭킹과 모델 품질평가

"좋은 pose는 생성되지만 랭킹이 틀린다"는 문제의식은 다수의 전용 재랭커를 낳았다. DeepRank-Ab(2026, Comms Biol)는 등변 그래프 신경망으로 계면을 채점해 AF3의 Top-1 성공률을 35.5%p 끌어올렸고, ABAG-Rank(2026)는 기하 기술자와 AF confidence만으로 learning-to-rank를 수행하며, ARID-sf(2026)는 물리 정보를 결합한 스코어링으로 150만 개 이상의 도킹 모델을 학습했다. 나노바디 특화 NbX(Tam et al., 2021)와 통합 confidence 지표 AntiConf(2026, pDockQ2+pTM), 그리고 PAE 기반 reference-free 지표 pDockQ2(Zhu et al., 2023)·ipSAE(Dunbrack, 2025)가 이 계보를 이룬다. 더 넓게는 구조 모델 품질평가(EMA)가 단량체에서 복합체로 확장되어(MULTICOM_qa, CASP15), CASP16의 QMODE3는 대규모 모델 풀에서 최적 모델을 **선택**하는 능력을 별도로 평가하기에 이르렀다 — 이는 본 연구가 다루는 문제(큰 pose 풀에서 near-native 선택)와 정확히 일치한다.

이 방법들의 공통점이자 본 연구와의 대비축은 명확하다: **모두 co-folder가 내놓은 단일 정지 구조(또는 decoy)를 사후에 채점하거나 confidence를 조합할 뿐, 예측이 섭동에 어떻게 반응하는지는 사용하지 않는다.**

### 2.4 학습 데이터 과대표집과 에피토프 면역우세 편중

co-folder의 정확도가 물리적 이해보다 **학습 계면의 암기**에 상당 부분 의존한다는 증거가 축적되고 있다. Guan & Keating(2025, Protein Science)은 AF3·Boltz-1·Chai-1이 학습 컷오프 이후 공개된 구조에서 이전 구조보다 정확도가 유의하게 하락함을(원자 수준 정확도 pre-cutoff 22–38% vs post-cutoff 6–13%), 그리고 학습셋에 결합부위 매치가 충분하면(≥200) 89–100%가 성공하지만 매치가 전무하면 어떤 모델도 원자 정확 예측에 실패함을 보였다. Glukhov, Vajda & Kozakov(2026, Curr Opin Struct Biol)는 이를 "암기 → 일반화 실패 → 물리의 필요"로 프레이밍한다.

이 암기의 근원에는 구조 데이터의 편향이 있다. SARS-CoV-2 RBD를 인식하는 항체의 구조 분석은 접근 가능 표면의 상당 부분이 잠재적 에피토프임에도 인식 빈도가 소수의 면역우세 사이트에 강하게 쏠림을 보였으며, 애초에 구조가 풀린 항체는 전체 서열 항체 레퍼토리의 약 5%에 불과하다(23-epitope 분석, 2023). AsEP 벤치마크(Liu et al., 2024, NeurIPS)는 1,723개 항체–항원 복합체가 641개 항원·973개 에피토프 그룹으로 long-tail 재사용됨을 정량화하고(예: 한 항원에 수십 개 항체가 반복 결합), 테스트 에피토프를 학습에서 완전히 배제하는 epitope-group split으로 novel 결합부위 일반화를 평가하는 프로토콜을 제공한다. 데이터 누수가 지표를 부풀린다는 인접 분야의 경고(Nat Mach Intell, 2025) 역시 비중복·일반화 평가의 필요를 뒷받침한다.

### 2.5 본 연구의 위치와 차별점

종합하면 선행 연구는 세 가지 사실을 확립했다. 첫째, MSA 깊이가 co-folder의 출력 basin을 지배한다(§2.2). 둘째, 그러나 그 깊이 조절은 지금까지 conformational diversity를 위해서만 쓰였고, 항체–항원에서는 깊은 항원 MSA가 유익하지 않으며(McCoy et al., 2024) 오히려 학습 과대표집 편향의 통로가 될 수 있다(Guan & Keating, 2025). 셋째, 실패는 fold가 아니라 배치·선택의 문제이나(Hitawala & Gray, 2025), 이를 겨냥한 재랭커들은 모두 단일 정지 구조만 채점한다(§2.3).

본 연구는 이 세 흐름이 만나는 지점에 있다. 우리는 (a) 항체–항원을 대상으로, (b) 항원 MSA 깊이를 편향의 **개입 변수**로 삼아 full에서 single-sequence까지 스윕하고, (c) 그 섭동에 대한 예측의 반응(depth-response)을 과대표집 에피토프 편향의 진단 및 near-native pose 선택 신호로 사용한다. MSA 서브샘플링(§2.2), 재랭커(§2.3), 학습 편향 규명(§2.4)은 각각 선점되어 있으므로, 본 연구의 기여는 새 알고리즘이 아니라 **이 세 요소를 항체–항원 도킹의 진단·선택 문제에 결합·적용**하는 데 있다. 우리가 아는 한, 항원 MSA 깊이 스윕으로 에피토프 편향을 진단하고 depth-response를 pose 선택 신호로 사용한 선행 연구는 없다.

---

## 3. 가설

- **H1 (편향의 존재·인과).** 항원 MSA 깊이가 co-folder의 예측을 과대표집 에피토프로 편향시킨다.
- **H2 (rescue).** 과대표집 부위가 정답이 아닌 항체(off-site)에서, 항원 MSA를 낮추면 올바른(near-native) pose가 되살아난다.
- **H3 (sweet-spot).** MSA를 지나치게 낮추면 공진화 신호가 붕괴해 pose 품질(DockQ)이 떨어진다. 따라서 편향 완화와 pose 품질을 동시에 만족하는 **최적 깊이(sweet-spot)**가 존재하며, 그 위치는 표적별로 다르다.
- **H4 (depth-response = 신호).** 깊이 섭동에 대한 pose의 반응(이탈 여부·이탈 깊이·안정성)은 정답을 모르는 추론 시점에도 (a) 그 표적이 편향에 취약한지, (b) 어느 pose가 near-native인지를 예측하는 데 유효한 신호이다.

---

## 4. 방법

### 4.1 데이터셋: 3그룹 설계

편향의 인과를 특이성으로 검증하기 위해 항체–항원 복합체를 세 그룹으로 구성한다. "on-site/off-site"는 과대표집 부위를 기준으로 정의되므로, 비과대표집 항원은 이 축으로 나뉘지 않는다(따라서 2×2가 아니라 3그룹이다).

| 그룹 | 정의 | 예측(항원 MSA 감소 시) | 역할 |
|---|---|---|---|
| A | 과대표집 항원 + 정답 에피토프가 그 부위 **안** | pose 품질 저하 (맞는 신호 제거) | reverse control |
| B | 과대표집 항원 + 정답 에피토프가 그 부위 **밖** | near-native pose rescue | 주효과 |
| C | 비과대표집 항원 | 에피토프 이동 없음(flip 없음) | no-flip control |

- 규모: A·B 합쳐 약 15개(각 그룹 균형), C 약 5개, 총 약 20개. 소규모이므로 **통계적 일반화가 아니라 메커니즘 입증(case study)**으로 프레이밍한다.
- 분류의 객관성: 각 항원의 과대표집 점수를 공개 구조 DB에서 잔기별 항체 접촉 빈도로 정량(over-representation map)하고, 정답 에피토프 위치와 대조하여 A/B/C를 기계적으로 배정한다.
- **누수(leakage) 규율:** 사용 co-folder들의 학습 컷오프 이후 공개된 복합체만 사용한다. 본 실험은 "모델이 MSA에 의존"해야 편향 효과가 드러나므로, 컷오프 이후 복합체를 쓰는 것은 공정성뿐 아니라 **효과 측정의 전제**다(암기한 복합체는 MSA를 줄여도 구조를 recall해 편향 효과를 가린다).

### 4.2 모델(생성 포트폴리오)

MSA를 입력으로 사용하는 co-folder만이 depth-sweep의 대상이 된다. 물리·Rosetta 도킹은 MSA를 쓰지 않으므로 depth-sweep이 무의미하며, 서로 다른 축(에너지·기하)의 decorrelated 생성기로서 풀 다양성에만 기여한다.

| 역할 | 도구 | MSA 사용 | depth-sweep | 학습 컷오프 |
|---|---|---|---|---|
| co-folder | Boltz-2 | 예 | 예(궤적 축) | 2023-06-01 |
| co-folder | Chai-1 | 예(single-seq 모드 병행 가능) | 예 | 2021-01-12 |
| co-folder | Protenix-base | 예 | 예 | 2021-09-30 |
| (옵션) 물리 | HADDOCK3 | 아니오 | 아니오 | — |
| (옵션) Rosetta | SnugDock | 아니오 | 아니오 | — |

- Protenix는 **base(2021-09) 체크포인트**를 사용한다. 2025-06 컷오프 변이체는 테스트 복합체를 암기해 편향 효과를 오염시킬 수 있어 본 실험에 부적합하다.
- 세 co-folder 중 가장 늦은 컷오프가 Boltz-2(2023-06)이므로, 테스트 복합체는 **2023-06 이후 공개분**으로 한정한다.
- AlphaFold3는 가중치 승인·비상업 제약으로 직접 실행 대상에서 제외한다.

### 4.3 MSA 깊이 사다리

- **천장:** full MSA(100%) — 각 co-folder가 해당 항원에 대해 실제 소비하는 최대 깊이(모델별로 다를 수 있음, 명시).
- **축(보고 단위):** 원자료 정렬 수가 아니라 **Neff80**(80% 동일성 cutoff, 잔기별 median effective sequences — AlphaFold 자체가 쓰는 depth 지표)로 정규화하여 표적 간 비교를 보장한다.
- **간격:** full → single-sequence(Neff≈1)까지 **log 등간격(배증) 6–7단계.** 원자료 근사로는 대략 [full · 1/2 · 1/4 · … · single-seq].
- **해석 랜드마크:** Neff≈30 부근을 공진화 신호 붕괴 임계로 표시한다("여기서부터 편향이 풀리기 시작").
- **seed:** 각 (모델, 깊이)마다 3–5 seed. 같은 깊이에서도 diffusion 확률성으로 pose가 흔들리므로, 깊이 효과와 seed 노이즈를 분리하기 위한 필수 요소다.

### 4.4 지표

- **주지표:** merged-chain DockQ(항원=A, 항체 H+L 병합=B), CAPRI 3-tier(0.23 Acceptable / 0.49 Medium / 0.80 High).
- **보조:** epitope recall(예측 접촉 잔기와 정답 에피토프의 겹침), over-representation overlap(예측 접촉이 과대표집 부위와 겹치는 비율), 항체 내부 CA-RMSD(fold 품질 대 배치 오류 분리).
- **깊이 축:** Neff80.

### 4.5 depth-response 특징

각 co-folder pose에 대해 궤적 기반 특징을 계산한다.

- flip depth: 깊이를 낮추며 예측 에피토프가 과대표집 부위에서 이탈하는 지점.
- over-rep overlap의 깊이-기울기: 얕아질수록 과대표집 부위 겹침이 줄어드는 정도.
- pose 변위·앙상블 분산: 깊이 간 pose 이동량, seed 간 변동.
- 안정성: 깊이 전 구간에서 동일 에피토프를 유지하는지(안정 = 신뢰, 섭동에도 안 움직이면서 과대표집 부위와 일치하면 오히려 편향 의심).

### 4.6 선택 실험(off-the-shelf)

새 모델을 학습하지 않고, 기성 재랭커(DeepRank-Ab / ARID-sf, 보조로 pDockQ2·ipSAE)를 두 풀에 적용해 비교한다.

- **single 풀:** full-MSA에서 생성한 pose만.
- **enriched 풀:** depth-sweep 전 구간(+옵션 물리/Rosetta) pose.

각 풀에서 재랭커의 top-1 DockQ를 비교하고, 궤적 특징을 추가했을 때의 이득을 별도로 평가한다(§5 ablation).

---

## 5. 실험 설계

### 5.1 특이성 검증(H1·H2·H3)

핵심은 rescue의 **특이성**이다. 항원 MSA를 낮출 때, near-native rescue가 **B에서만** 나타나고 A에서는 오히려 악화, C에서는 이동이 없어야 한다. 만약 rescue가 A·C에서도 나타나면 이는 편향이 아니라 단순한 "MSA 감소 = 다양성 증가"이며, **B에 국한된 rescue만이 과대표집 편향을 원인으로 지목하는 인과 증거**가 된다. 따라서 A(reverse)와 C(no-flip)는 본 설계에 내장된 두 개의 negative control이다.

- (선택 보강) Neff-matched shuffling 대조: 동일 Neff에서 깊이 감소 대신 서열을 섞어, 효과가 "깊이/편향"인지 "무작위 노이즈"인지 분리한다. 이는 선점 문헌(MSA shuffling)과의 직접 비교를 제공한다.

### 5.2 baseline 및 ablation(H4)

- **baseline:** (1) full-MSA + ipTM 선택(표준), (2) single-MSA 풀 재랭킹, (3) best-single-model.
- **capability(주장):** enriched 풀 + 재랭커의 top-1 DockQ가 위 baseline을 이기는가, 그 이득이 B(off-site)에 집중되는가.
- **ablation(머니 실험):** 재랭커에 depth-response 특징을 넣은 경우와 뺀 경우를 비교하여, 궤적 특징이 단일 구조 재랭킹 대비 추가 이득을 주는지 판정한다. 이득이 있으면 궤적이 본 연구의 novelty이며, 없으면 정직하게 제거하고 커버리지(§6.3)로 후퇴한다.

---

## 6. 예상 결과

### 6.1 sweet-spot 분포
그룹별로 최적 깊이의 분포가 다르게 나타날 것으로 예상한다. A는 깊은 쪽(full 근처), B는 중간~얕은 쪽(Neff≈30 이하)으로 sweet-spot이 이동하고, C는 깊이에 대체로 둔감할 것이다.

### 6.2 rescue의 특이성(H1·H2·H3)
B에서 항원 MSA 감소가 과대표집 부위 이탈과 함께 DockQ를 유의하게 올리고(near-native rescue), A에서는 DockQ가 떨어지며, C에서는 에피토프 이동이 없을 것이다. 또한 B에서도 지나친 감소(Neff≈1 부근)에서는 pose 품질이 다시 떨어져, sweet-spot의 존재가 드러날 것이다.

### 6.3 커버리지(생성 단계 진단)
depth-sweep 풀의 oracle(달성 가능한 최고) DockQ가 단일-깊이 풀보다 높고, 그 격차가 B(off-site)에 집중될 것이다. 이는 full-MSA co-folder가 off-site near-native pose를 **애초에 생성하지 못하며, MSA 섭동이 그것을 풀에 넣어준다**는 것을 정량한다.

### 6.4 선택 capability
enriched 풀 + 기성 재랭커의 top-1 성공률이 baseline(full+ipTM, single-MSA 재랭킹)을 넘고, 이득이 B에 집중될 것이다. 즉 재랭커가 단일-MSA 풀에서는 고르지 못하던 near-native를, MSA-섭동으로 넓힌 풀에서는 골라낼 것이다.

### 6.5 depth-response = 편향 센서(H4)
flip 신호(이탈 여부·이탈 깊이)가 그룹 라벨(B vs A·C)과 강하게 상관하여, 정답을 모르는 상태에서도 "이 케이스가 편향에 취약한가"를 예측할 것이다. ablation에서 궤적 특징이 추가 이득을 주면, depth-response가 pose 선택의 유효 신호임을 확립한다.

---

## 7. 결론 및 기여

MSA는 과대표집을 통해 항체–항원 예측에 에피토프 편향을 부여한다. 기존 연구는 MSA 감소·shuffling을 주로 **대안 conformation 다양성** 확보에 사용했으나, 도킹에서는 이 편향이 **에피토프 정확도를 직접 저해하며 downstream wet 실험을 오도**한다. 본 연구는 (1) 이 편향의 인과를 3그룹 특이성 설계로 진단하고, (2) MSA 깊이 섭동으로 off-site near-native pose를 생성해 기성 재랭커가 이를 골라내게 하며, (3) depth-response가 편향 감지·pose 선택에 유효한 특징임을 보인다. 산출물인 **depth-sweep 데이터셋과 depth-response 특징**은 향후 항체–항원 재랭커·co-folder 학습에 넣을 수 있는 새로운 신호의 근거를 제공한다.

방법 각 요소는 선점되어 있으므로(§2), 기여는 새 알고리즘이 아니라 **"depth-response를 항체–항원 도킹의 선택·진단 신호로 결합·적용"**하는 데 있으며, 이를 정직하게 명시한다.

---

## 8. 한계 및 향후 연구

- **표본 규모.** 약 20개는 메커니즘 입증용이며 통계적 일반화의 검정력은 낮다. 학습형 선택기(depth-response 특징 학습)는 별도의 대규모 leakage-free SAbDab 세트가 필요하며(누적 컷오프 이후 비중복 복합체 규모는 수백 수준), 본 20개 세트로 모델을 학습했다고 주장하지 않는다.
- **추론 시 sweet-spot 자동 선택은 어려운 가설.** 메커니즘(§6.2·6.5)은 견고하나, 궤적 특징이 단일 구조 재랭킹 대비 추가 이득을 주는지는 ablation이 판정한다. 이득이 없어도 커버리지(§6.3)와 편향 진단은 성립한다(graceful degradation).
- **선점.** MSA 서브샘플링(AFsample2·AF-Cluster·subsampled-AF2), 재랭커(DeepRank-Ab·ARID-sf·ABAG-Rank·pDockQ2·ipSAE), 학습 편향(Guan & Keating 2025)은 모두 인용·활용한다.
- **wet 연결(전망).** 예측된 에피토프·pose는 alanine scanning 또는 epitope binning(BLI 경쟁 assay)으로 검증 가능하며, 특히 연구실의 novel 표적(예: Beacon 스크리닝 산출물)에 대한 prospective 예측→wet 검증이 향후 핵심 확장이다.

---

## 부록: 실행 순서(요약)

1. A/B/C 복합체 확정(RBD로 A·B, 비과대표집 clean 복합체로 C; 2023-06 이후).
2. Boltz-2·Chai-1·Protenix-base × depth 사다리(full→single-seq, Neff80 log) × seed 3–5 생성.
3. pose별 merged DockQ·epitope recall·over-rep overlap·Neff80 산출 → depth-sweep 데이터셋.
4. depth-response 특징 추출 → 3그룹 특이성·sweet-spot·커버리지 분석(§6.1–6.3, 6.5).
5. 기성 재랭커로 enriched vs single 비교 + 궤적 특징 ablation(§6.4).
6. 결과 정리 → 그림 5종(그룹별 depth-DockQ 곡선, rescue 특이성, oracle 커버리지, 재랭커 성공률, flip↔그룹 상관).

---

## 참고문헌

**co-folder 항체–항원 한계·선택 문제**

1. Abramson J. et al. (2024). Accurate structure prediction of biomolecular interactions with AlphaFold 3. *Nature* 630:493–500. https://www.nature.com/articles/s41586-024-07487-w
2. Hitawala F. N. & Gray J. J. (2025). What does AlphaFold3 learn about antibody and nanobody docking, and what remains unsolved? *mAbs* 17:2545601. https://www.tandfonline.com/doi/full/10.1080/19420862.2025.2545601
3. FoldBench: an all-atom benchmark for biomolecular structure prediction (2025). *Nature Communications*. https://www.nature.com/articles/s41467-025-67127-3
4. Structural Plausibility Without Binding Specificity: Limits of AI-Based Antibody–Antigen Structure Prediction Confidence Scores (2026). *bioRxiv* 2026.03.02.709004. https://www.biorxiv.org/content/10.64898/2026.03.02.709004v1
5. Yin R. & Pierce B. G. (2024). Evaluation of AlphaFold antibody–antigen modeling with implications for improving predictive accuracy. *Protein Science* 33:e4865. https://pmc.ncbi.nlm.nih.gov/articles/PMC10349958/
6. Yin R., Wang X. & Pierce B. G. (2025). Evaluating deep learning based structure prediction methods on antibody–antigen complexes. *bioRxiv* 2025.07.11.662141. https://www.biorxiv.org/content/10.1101/2025.07.11.662141v2.full
7. Dreyer F. A. et al. (2025). Ibex: Conformation-Aware Structure Prediction of Antigen-Recognizing Immune Proteins. *arXiv*:2507.09054. https://arxiv.org/abs/2507.09054
8. McCoy K. M. et al. (2024). A comparison of antibody–antigen complex sequence-to-structure prediction methods and their systematic biases. *Protein Science* 33:e5127. https://onlinelibrary.wiley.com/doi/10.1002/pro.5127

**MSA 깊이 조절·Neff**

9. del Alamo D., Sala D., Mchaourab H. S. & Meiler J. (2022). Sampling alternative conformational states of transporters and receptors with AlphaFold2. *eLife* 11:e75751. https://elifesciences.org/articles/75751
10. Monteiro da Silva G., Cui J. Y., Dalgarno D. C. & Rubenstein B. M. (2024). High-throughput prediction of protein conformational distributions with subsampled AlphaFold2. *Nature Communications* 15:2464. https://www.nature.com/articles/s41467-024-46715-9
11. Wayment-Steele H. K. et al. (2024). Predicting multiple conformations via sequence clustering and AlphaFold2 (AF-Cluster). *Nature* 625:832–839. https://www.nature.com/articles/s41586-023-06832-9
12. Wallner B. (2023). AFsample: improving multimer prediction with AlphaFold using massive sampling. *Bioinformatics* 39(9):btad573. https://academic.oup.com/bioinformatics/article/39/9/btad573/7274860
13. Kalakoti Y. & Wallner B. (2025). AFsample2 predicts multiple conformations and ensembles with AlphaFold2. *Communications Biology* 8. https://www.nature.com/articles/s42003-025-07791-9
14. Rajabi H. et al. (2026). NEFFy: a versatile tool for computing the number of effective sequences. *Bioinformatics* 42(6):btaf222. https://academic.oup.com/bioinformatics/article/42/6/btaf222/8155843

**재랭킹·모델 품질평가**

15. DeepRank-Ab: a scoring function for antibody–antigen complexes based on geometric deep learning (2026). *Communications Biology* (bioRxiv 2025.12.03.691974). https://www.biorxiv.org/content/10.64898/2025.12.03.691974v1
16. ABAG-Rank: Improving Model Selection of AlphaFold Antibody–Antigen Complexes by Learning to Rank (2026). *bioRxiv* 2026.03.17.712376. https://www.biorxiv.org/content/10.64898/2026.03.17.712376v1
17. ARID-sf: A physics-informed Deep Learning scoring function to improve Antibody–Antigen docking model ranking (2026). *bioRxiv* 2026.01.20.700530. https://www.biorxiv.org/content/10.64898/2026.01.20.700530v1
18. AntiConf: Confidence Scoring for AI-Predicted Antibody–Antigen Complexes (2026). *Briefings in Bioinformatics* 27(2):bbag137. https://academic.oup.com/bib/article/27/2/bbag137/8554115
19. Tam C. et al. (2021). NbX: Machine Learning-Guided Re-Ranking of Nanobody–Antigen Binding Poses. *Pharmaceuticals* 14(10):968. https://pmc.ncbi.nlm.nih.gov/articles/PMC8537642/
20. Zhu W., Shenoy A., Kundrotas P. & Elofsson A. (2023). Evaluation of AlphaFold-Multimer prediction on multi-chain protein complexes (pDockQ2). *Bioinformatics* 39(7):btad424. https://academic.oup.com/bioinformatics/article/39/7/btad424/7219714
21. Dunbrack R. L. Jr. (2025). Rēs ipSAE loquuntur: What's wrong with AlphaFold's ipTM score and how to fix it. *bioRxiv* 2025.02.10.637595. https://www.biorxiv.org/content/10.1101/2025.02.10.637595v1
22. Liu J., Guo Z., Wu T. & Cheng J. (2023). Combining pairwise structural similarity and deep learning interface contact prediction (MULTICOM_qa, CASP15). https://pubmed.ncbi.nlm.nih.gov/36945536/
23. Cheng J. et al. (2025). Highlights of Model Quality Assessment in CASP16. *Proteins*. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12750031/

**학습 편향·에피토프 과대표집**

24. Guan E. & Keating A. E. (2025). Training bias and sequence alignments shape protein–peptide docking by AlphaFold and related methods. *Protein Science* 34:e70331. https://pmc.ncbi.nlm.nih.gov/articles/PMC12518507/
25. Glukhov E., Vajda S. & Kozakov D. (2026). From memorization to generalization: Why physics will improve machine learning-based prediction of protein complexes. *Current Opinion in Structural Biology* 98:103288. https://www.sciencedirect.com/science/article/abs/pii/S0959440X26000709
26. SARS-CoV-2 antibodies recognize 23 distinct epitopic sites on the receptor binding domain (2023). *PMC10275037.* https://pmc.ncbi.nlm.nih.gov/articles/PMC10275037/
27. Liu C. et al. (2024). AsEP: Benchmarking Deep Learning Methods for Antibody-specific Epitope Prediction. *NeurIPS Datasets & Benchmarks* (arXiv:2407.18184). https://arxiv.org/abs/2407.18184
28. Resolving data bias improves generalization in binding affinity prediction (2025). *Nature Machine Intelligence.* https://www.nature.com/articles/s42256-025-01124-5

*주: 일부 인용(예: 26번 23-epitope 논문의 저자·서지, 그리고 bioRxiv 초안본의 세부 수치)은 제출 전 원문 PDF로 최종 확인 필요.*
