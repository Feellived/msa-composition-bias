# 실험 상세 로그 — MSA 조성 / 결합자리 (2026-07 아카이브)

세션 메모리가 너무 커져(563줄) 여기로 옮긴 원본이다. 현재 상태 요약은 메모리
`msa-bias-channel-experiment` 와 `msa-composition-phase2` 에 있고, 이 파일은 **날짜별 원문 기록**이다.

---

---
name: msa-bias-channel-experiment
description: "⭐**MSA 조성이 결합자리 선택을 좌우**; 중심 사례 2종 확정(8ulr_HL 1/30·8k3k_D 0/20, 둘 다 순열 p<5e-6 + 2×2 함께 유의) · 승격 후보 3종(대조 확대 12실행씩) · **전수 방향 5 대 5 → 선택 단계 필수** · 6.2 대상 선정 가설; ⛔p=3.9e-8 폐기; eval_dump_seedrep.py는 --data 필수"
metadata: 
  node_type: memory
  type: project
  originSessionId: be7059af-14d1-4717-a37d-c316fbae23f5
  modified: 2026-07-30T21:00:00.000Z
---

## ⭐⭐⭐ 2026-07-30 저녁 — 중심 사례(구 '앵커') 2종 확정 + 승격 후보 3종

**용어 확정:** 「앵커」→ **「중심 사례」**(통제를 모두 갖추어 집중 검증한 단일 복합체).
**「대조 확대」** = 원래 MSA(seedfull) 실행 횟수만 늘리는 것(조성 쪽은 불변). 0/8은 진짜 성공률이
20%여도 17% 확률로 나오지만 0/20이면 1.2% → **직접 맞대는 2×2 비교의 검정력을 올린다.**
순열(이질성) 검정은 조성 간 차이를 쓰므로 대조 크기에 덜 민감하다.

**확정 2종 — 9집단 이질성과 2×2 축약에서 함께 유의(검정방식 논쟁 없음).**
- **8ulr_HL**(Env, 자세·자리 양축): 조성 14/32(자세)·16/32(자리) 대 **원래 MSA 1/30**,
  통계량 41.85·49.52, 순열 **p < 5×10⁻⁶**, 2×2 Fisher 0.00016·0.000025. 성공 실행 recall 0.70(n=15)
  대 실패 0.24(n=47). 조성 내 sd 0.012 대 조성 간 0.259.
- **8k3k_D**(RBD 나노바디, 자리): 조성 18/24 대 **원래 MSA 0/20**, 통계량 53.99, 순열 p < 5×10⁻⁶,
  2×2 **1.3×10⁻⁷**. 자세 축은 0/24(무신호). **후보 생성(0.115→0.962)·데모 통과와 같은 복합체 = 척추.**

**승격 후보 3종 (이미 순열 유의 · 대조가 8회뿐이라 2×2만 미달).** 복합체당 12실행, 셋 36실행.
- 8tp5_HL 8/24 대 0/8(순열 0.0123) → 대조 20회면 2×2 ≈ 0.004
- 8k46_I 7/24 대 0/8(0.0047) → ≈ 0.009 · 8sis_HL 12/24 대 1/8(0.0308) → ≈ 0.005
- 8ume_HL은 확대해도 안 됨(≈0.26). ⚠️ **잘 나온 종만 늘리면 선택적 보강** → 조성 우세 5종 전부에
  동일 적용해야 방어 가능(2종 완료).

**⛔ 폐기 확정: `p = 3.9×10⁻⁸`** — 통제 실행 하나(자세 40개)를 독립 40표본으로 센 값. 통제 1회가
실패할 확률은 0.9이므로 0/40은 특이한 관측이 아니다. **다시 언급하지 말 것.**

**도구·사고** = `analyze_anchor_test.py`(레포 커밋 f7a85c5 · 4c6ef5d, 방향 열 포함).
⚠️ `eval_dump_seedrep.py`는 **`--data $DATA/compreps` 필수** — 빼면 옛 seedrep_cand(무효 boltz 3건)를
읽고 CSV를 덮어쓴다(2026-07-30 실제 사고, 예측 CIF는 무사해 재채점으로 복구). DockQ 없는 env에서는
파일을 쓰지 않고 멈춘다(정상 동작).

**전수 보강 분석(29종, 원래 MSA 포함 실행 단위 순열):** 자리 축 **유의 10/29인데 방향이 5 대 5**
(조성 우세 8k3k_D·8ulr_HL·8k46_I·8tp5_HL·8sis_HL / 원래 MSA 우세 9mqr_DE·8k5h_HL·9zdu_HL·8siq_HL·
9azt_HL). 자세 축은 8ulr 단독. → **"조성이 자리를 정한다"는 성립, "조성을 바꾸면 좋아진다"는 불성립**
→ 선택 단계 필수. **6.2 대상 선정 가설**: 데모 8종 중 조성 우세는 8k3k_D 하나였고 그 하나가 통과
(1/1 대 0/7, Fisher 0.125 = 가설이며 결과 아님). 검증 대상 = 조성 우세 4종 데모 미실시분.


## ⭐ 2026-07-29 본 검정 중간 결과 (18/29) — 사전 기준 세 개 다 통과

**결합자리(위치) 축에서 신호가 나왔다. 자세(DockQ) 축은 전멸.** 인수인계서 Ⅱ **6.6**에 기록,
빈 칸은 **6.7**(29종 완료 후 기입). 원자료 = 레포 `pipeline/results/`(커밋됨).

6.4 사전 기준 대조: **①효과** 이질성 유의 **5/18**(Fisher 결합 p=**0.00072**, 우연확률 0.0015) ·
**②재현성** 순열 p<0.05 **18/18**(실패 0) · **③후보생성** 정답자리 절반 이상 덮음 **16/18**.

⚠️ **2026-07-29 오후 정정** — 처음 집계의 `Fisher 0.00093 · 순열 17/18(8k5g_HL 실패)`는 **오염된 값**.
`8k5g_HL`(d58+d90)·`8q7s_C`(d35+d86) 두 종은 한 폴더에 **옛 arm과 본 검정 깊이가 같이** 있었고, 실행 폴더
이름이 두 깊이에서 같아(`c0_r0` 등) **한 실행의 자세가 5개→10개로 합쳐져** 최고값이 양쪽 다 부풀려짐;
반대로 결합자리 계산은 첫 깊이 폴더에서만 자세를 찾아 **조성 5~6가지만** 씀. 본 검정 깊이로 재채점 =
8k5g_HL 순열 **p 0.970→0.0125**(조성내 0.464/조성간 0.411), 8q7s_C **p 0.0165→0.0005**(0.745/0.661).
둘 다 이질성은 여전히 미유의 → **유의 5종 명단·최강 사례·기준③ 16/18은 불변.** 나머지 16종은 깊이 하나라
무관. 기록 = 인수인계서 Ⅱ **7.3.5**(사고) + 6.6.1(정정된 표·문단). 고침 커밋 `20a3695`·`1b1c43f`.
⭐**재사용 교훈(7.3.1·7.3.3과 같은 계열)**: 이 파이프라인의 실패는 멈추지 않고 **조용히 비거나 섞인다.**
채점 결과를 믿기 전에 **실행 수가 설계값과 맞는지(조성 수 × 반복 수 + 원래 MSA 횟수 = 32)** 먼저 볼 것.
같은 세션에서 발견한 것 하나 더: **DockQ conda 환경 밖에서 돌리면 dockq 열만 통째로 비고**
recall·overrep은 정상이라 겉보기엔 결과가 있다 → 이제 종료코드 5로 멈추고 파일을 안 쓴다.
→ 사전에 갈라둔 **갈래 1**(효과+재현성 확인) 조건 충족. 유의 = 8k5h_HL(.0003)·8k3k_D(.0013)·
9mqr_DE(.0028)·8k46_I(.031)·8siq_HL(.046); **BH 보정 시 앞 3종만 생존**.

⭐**가장 강한 사례** = 원래 MSA가 놓친 자리를 조성이 찾은 6종. **8k3k_D**(원래 0.12 → 조성 5개가
모인 후보 0.96, 조성내/조성간 0.82/0.22=3.7배, 순열 p=.0005) · **8sit_HL**(0.10 → 0.90, 조성 6개).
나머지 4종은 조성 1개만 모여 "운 좋은 추첨"에 가까움.

⚠️**반드시 함께 말할 것**: (a) **DockQ 576회 중 ≥0.49가 4회**, 절반은 원래 MSA 쪽. DockQ 이질성은
18종 전부 p=1.0(성공 0이라 무정보). → 정확한 서술 = **"자리는 찾되 자세는 못 맞힌다"** = 생성(조성
다양화)+선택(재랭커) 한 쌍 설계와 정합. (b) 기준③ 16/18은 관대 — 후보 1개뿐인 8종은 "전부 같은
자리로 감". (c) **8xsi_HL은 조성 7개가 전부 틀린 자리에 합의**(덮음 0.00) → 합의≠정답. (d) 유의 5종
중 4종은 조성 쪽 중앙값이 **더 낮음** — 명제와 모순 아님("조성이 성패를 가름"이지 "얕을수록 좋다"가
아님, 후자는 기각됨). (e) **원래 MSA vs 조성 직접비교는 깊이가 달라 교란** — 깨끗한 건 같은 깊이 안의
이질성·순열 검정뿐. (f) 분모가 29가 되면 우연확률 0.0136으로 올라가고 Fisher는 희석됨.

빠진 것: **9y0a_AB 채점 누락**(예측은 완료) · 8ulr_HL은 설계가 달라 별도 행.

## ⭐⭐⭐ 2026-07-30 본 검정 **최종 확정** (n=29, 예측·채점 30/30 완료) — 갈래 A

**세 기준 다 통과. 어젯밤 걱정한 Fisher 희석은 일어나지 않았다.**

| 기준 | 값 | 판정 |
|---|---|---|
| ①효과 이질성 유의 | **6/29** · 우연확률 **0.0027** | 통과 |
| ①효과 Fisher 결합 | **p = 0.002621** | 통과 |
| ①효과 BH(q=.05) 생존 | 3/29 | (보고용) |
| ②재현성 순열 p<0.05 | **25/29** (실패 8u44_ST·9azr_HL·9bdg_FI·9y0a_AB) | 통과 |
| ③후보생성 | **27/29** (관대함 — ★는 4종) | 통과 |
| DockQ ≥0.49 | **4회 / 928회** (8k5g_HL·8k5h_HL만) | "자리는 찾되 자세는 못 맞힌다" |

유의 6종 = **8k3k_D · 8k46_I · 8k5h_HL · 8siq_HL · 9mqr_DE · 9zdu_HL**.
⭐**★ 4종**(원래 MSA가 절반 못 덮었는데 후보가 넘김) = 8k3k_D(0.12→0.96) · 8sit_HL(0.10→0.90) ·
8ume_HL(0.23→0.95) · 8k46_I(0.24→0.72). **주장에 쓸 것은 이 4종.**

⚠️⚠️ **층별 결과가 사전 예상과 반대다 — 반드시 이대로 적을 것.**
`다양성없음(poor) 5/14 · 다양성있음(rich) 1/15` / 군별 `RBD 5/14 · C 1/5 · Env 0/2 · HA 0/8`.
§6.2.2의 사전 세 갈래는 ①다양성있음에서만 ②양쪽 다 ③양쪽 다 안 나옴이었는데, 실제는
**다양성없음에서만** = 사전에 없던 네 번째 경우. 그리고 **앵커 8ulr_HL은 Env·Neff80 3051 =
가장 다양한 층**에 있다 → **결정적 단일 사례가 전수 검정에서 신호가 거의 없는 층에 있다**는 긴장.
층 ≡ 항원 계열(다양성없음 14종 전부 RBD)이라 "조성 다양성이 낮아서"와 "RBD라서"를 분리 불가.
후자가 더 자연스럽다(RBD는 과대표집이 심한 항원). **이 교란을 빼고 층 해석을 쓰지 말 것.**

빈도 = 세트 3(분모 10) · 세트 4(분모 49, `sweep_targets.csv` 59행 = 49+10이므로 세트 3을 뺀다,
커밋 `5c7ebce`) · 층별. 도구 = `eval_criteria.py`(18종·24종·29종에서 검증).

## ⭐ 2026-07-30 데모 최종 — 8종 중 1종 통과, **3종은 우리 제약이 해로웠다**

| 복합체 | noconstraint | ours | 판정 |
|---|---|---|---|
| **8k3k_D** | 0.008/0.037 | **0.285/0.760** | ✅ 유일한 통과(네 팔 단조증가) |
| 8ume_HL | 0.021/0.231 | 0.068/0.692 | 개선되나 fullmsa·sizematch도 동일 = 자리 특이성 없음 |
| 8siq_HL | 0.146/0.579 | 0.159/0.579 | 무차이 |
| 8sit_HL | 0.818/0.893 | 0.787/0.923 | 무차이(모델이 원래 품) |
| 8t4a_PR | 0.093/0.433 | 0.058/0.300 | ❌ 해로움 |
| 8t4d_OQ | 0.092/0.357 | 0.013/0.143 | ❌ 해로움 |
| 9azt_HL | 0.199/0.706 | 0.015/**0.000** | ❌ 심하게 해로움 |
| 9mqr_DE | 0.007/0.000 | 0.008/0.000 | 둘 다 실패 |

(DockQ최고 / 자리겹침최고 · 네 팔 모두 자세 5개)

⭐**정확한 서술 = "현재 선택기로는 순 효과가 없다."** 1종 크게 개선 · 3종 악화. 8k3k_D는 다섯 팔
통제를 다 통과했으므로 **파이프라인이 작동할 수 있음의 존재증명**이지만, 적용 가능 여부를 미리
판별하지 못하면 실전에 못 쓴다. **9azt_HL이 선택기 실패의 명확한 사례**: 후보 7개 중 최고는 정답을
0.88 덮는데 `ncomp`는 0.00 덮는 것을 골랐다(서열 4개짜리라 '조성'이랄 게 없음).
→ 다음 과제(AbEpiScore-1.0로 선택기 교체)의 필요성이 추측이 아니라 **측정값**이 됐다.

## ⭐ 2026-07-29 밤 — 데모 실행 결과: 3종 중 1종만 통과 (기록 = 인수인계서 **Ⅰ §8**)

**⚠️ 8sit_HL 해석이 뒤집혔다. 낮에 "DockQ 0.787 대박"이라 한 건 틀렸다.**

| 복합체 | noconstraint | fullmsa | sizematch | ours | 판정 |
|---|---|---|---|---|---|
| 8k3k_D DockQ | 0.008 | 0.015 | 0.017 | **0.285** | ✅ 단조증가 |
| 8k3k_D 자리겹침 | 0.037 | 0.161 | 0.289 | **0.760** | ✅ |
| 8sit_HL DockQ | **0.818** | 0.023 | 0.019 | 0.787 | ❌ 제약없음이 더 높음 |
| 8ume_HL 자리겹침중앙 | 0.091 | 0.609 | 0.667 | 0.621 | ❌ 세 팔 동일 |

⭐**8k3k_D만 다섯 팔 설계 통과.** 8sit_HL은 **boltz가 원래 혼자 푸는 복합체**(제약 없이 0.818) →
우리 기여 0. 8ume_HL은 제약은 먹지만 **어느 자리를 주든 같음** → 우리 선택의 기여 0.
훑기 탈락 4종: 8siq_HL(0.579/0.159)·8t4a_PR(0.300/0.058)·8t4d_OQ(0.143/0.013)·9mqr_DE(0.000/0.008).

⚠️⚠️ **재사용 교훈: `ours` 팔만 보고 판단하지 말 것.** 통제 팔(noconstraint) 없이 본 숫자는
"모델이 원래 잘하는 것"과 구별이 안 된다. 훑기에 noconstraint를 넣도록 고침(`9a68a5e`).
**합격선 = ours > noconstraint 가 1차 관문**, 그다음 ours > sizematch > fullmsa 단조증가.

**protenix는 제약을 정상 적재하고도 안 움직인다**(로그 `#pocket:62` 확인, 네 팔 전부 동일).
→ 역할 분리 = **생성 protenix · 붙이기 boltz**. 근거는 실측 성능이지 아키텍처 추측이 아니다
(⚠️"boltz가 MSA에 덜 민감"은 a3m 사고로 **근거 무효** — 발표에 쓰지 말 것).

**AbEpiTope 트랙 = 9mqr_DE·8siq_HL**(정답 후보가 목록에 있는데 `ncomp`가 못 고름).
8t4d_OQ는 **선택이 옳았는데도 실패**라 성격이 다름(한계 슬라이드). **DiscoTope은 쓰지 않는다** —
항원 표면만 보고 인기 자리를 선호해 우리가 벗어나려는 편향을 선택 단계에 다시 넣는 순환 논리.

**발표 구성**: 근거(8ulr_HL 통제 + 본 검정 29종) · 데모 성공(8k3k_D) · 한계 3종.
**앵커를 더 만들지 말 것** — 타깃당 ~2.5h인데 n=29 체계적 증거가 이미 있어 한계효용 낮음.
**진행 중**: 본 검정 잔여 5종(9azr·9azt·9azv·9b7g·9bdg) ~04:00 완료 → `run_after_maintest.sh --apply SCREEN=1`이
채점(낡은 것 강제 재채점)·후보자리·noconstraint+ours 훑기까지 자동. **아침 할 일 = §6.7 통계 + Figure + 발표자료.**

### 사례를 늘릴 수 있는 경로 4개 (2026-07-29 밤 정리 — GPU가 비면 이 순서로)

본 검정이 ~04:00 끝나면 발표(7/31)까지 GPU 약 14시간. **아무것도 버리지 않는다** — 본 검정 실패
24종은 분모(우연확률·Fisher·빈도가 전부 분모에 의존), 데모 실패 6종은 **원인 4갈래 분류표**
(모델이 이미 잘함 / 제약이 자리를 안 가림 / 선택기가 틀림 / 유도가 안 먹힘).

1. ⭐**AbEpiTope 트랙 (9mqr_DE·8siq_HL, 각 ~40분) — 1순위.** 둘 다 **정답 후보가 목록에 이미 있는데**
   (천장 1.00·0.81) `ncomp`가 못 골랐을 뿐 → 선택기만 바꾸면 살아날 수 있음이 확인된 상태.
   성공 시 데모 +1~2 & "선택기 교체로 살아난다" 슬라이드. ⚠️관문 둘: (a)AbEpiTope이 **제약 걸린**
   구조를 변별하나(그 도구가 보던 입력과 분포가 다름) (b)골라도 `ours > noconstraint` 통과해야 함.
2. **8k3k_D를 두 번째 앵커로 (~2.5h).** 8ulr_HL과 같은 설계(조성 재추첨 + 예산 맞춘 통제).
   원래 0.12→후보 0.96·순열 p=0.0005라 성공 가능성 높음. 앵커 1건 → 2건이면 연구 주장이 단단해짐.
3. **boltz 본검정 (~1.5h/타깃) — 투기적.** 26종 중 **19종이 "고른 후보 = 원래 MSA 자리"라 데모에서 탈락**
   했는데, boltz로 조성 실험을 다시 하면 후보 구성이 갈려 일부가 ✅로 바뀔 수 있음. 될지는 돌려봐야 앎.
   ⚠️연구 주장으로 쓰려면 명단을 결과와 무관하게 정해야 함(데모용이면 "골랐다" 고지).
4. 8sit_HL·8k46_I 앵커(각 ~2.5h) — 2번이 되면 그다음.

도구: `eval_criteria.py`(6.7 판정, 커밋 `9ee2286` — 18종으로 검증: 유의 5·우연 0.0015·Fisher 0.000933·
BH 3·순열 17/18·후보생성 16/18·DockQ 4/576 전부 재현) · `run_maintest_boltz.sh`(3번용) ·
`prescreen_abepitope.sh`+`score_abepitope.py`(1번용).

## 2026-07-29 저녁 — 발표 데모(생성→선택→유도 재도킹) 코드 준비 완료

**사용자 결정: 연구 주장 아님. 잘 된 복합체만 골라 보여주고 "내가 골랐다"를 발표에서 구두 고지.**
목적 = 우리 설계가 이렇게 굴러간다는 **예상 흐름 시연**. 흐름 = 조성 다양화로 **생성** →
**선택**(재랭커) → 고른 에피토프를 제약으로 **도킹만 재실시** → 안 거쳤을 때 대비 **DockQ 상승**.

⭐**대조군 설계(중요)**: 대조를 "제약 없음"으로만 두면 보이는 건 "제약 주면 좋아진다"뿐 = 우리 기여
증명 안 됨. **세 팔** — (a)제약없음 (b)**원래 MSA가 간 자리로 유도** (c)우리가 고른 자리로 유도.
**판정은 (c)>(b).** 발표에서 반드시 찔리는 지점이라 (b)가 필수.

**선택 단계 근거**: 눈먼 선택(정답 안 보고 "조성이 가장 많이 모인 후보") = **13/18로 기준선(13/18)과 동률**,
천장은 16/18 → **생성 ○ · 선택 ✗**(Exp4 HADDOCK "생성 O 선택 X"와 같은 패턴). 그래서 재랭커 자리를
비워두고 붙일 수 있게 설계. 1순위 후보 = **DiscoTope-3.0**(항원 구조만·MSA 안 씀) · **AbEpiTope-1.0**
(계면 타당성·inverse folding) — **둘 다 미설치, 사용자 허가 대기**. 안 되면 manual(구두 고지).

**커밋**: msa-depth `e7a6eaa`(`analyze_site_reproducibility.py --dump-sites` → `results/sites_<타깃>.json`에 후보별
잔기 목록·모인 조성·원래MSA 포함여부; `run_analyze_target.sh`가 자동 생성) / antibody-ml `58fb379`
(`rank_sites.py`·`sites_to_pocket.py`·`run_demo_guided.sh`·`dockq_demo.py`, emit은 `zdock_to_pocket` 재사용).
잔기 키 = (항원 사슬 순번, 0-based) → 제약은 (사슬 id, 1-based). `rank_sites.py`는 `true_covered`·
`precision`을 **불러오며 버려** 눈먼 선택을 구조적으로 강제. 데모 후보 = 8k3k_D(0.12→0.96)·8sit_HL(0.10→0.90)·
8t4d_OQ(0.64→1.00) + 앵커 8ulr_HL(DockQ 성공·예산맞춘 통제 p=3.9e-8). **다음 = 본 검정 완료 후 1개로 3팔 스모크.**

## 2026-07-28 (최종) 본 검정 명단 확정 + 스모크 — 여기서 재개할 것

**상태: 판독·스모크 전부 끝났고(2026-07-28 16:15, 통과 59·실패 0), 본 검정(약 40시간)만 남았다.**

```bash
cd ~/projects/msa-composition-bias/pipeline && git pull
python prep_pick_depth.py --out maintest.csv
conda activate DockQ ; bash run_smoke_maintest.sh --gpu   # 통과해야 시작 (15~20분)
tmux new -s maintest ; HOURS=12 bash run_maintest.sh --apply
```

**⚠️ 스모크가 실제로 잡아낸 조용한 실패(고쳐서 커밋함 `04f0e41f`).** `eval_dump_seedrep.py`는
`seedrep_cand.csv`(예전 후보 5개)만 순회해서 **본 검정 29개 타깃을 한 번도 채점하지 않았고,
결과가 비면 파일을 안 쓰는 구조라 경고도 없었다.** `run_analyze_target.sh`도 같은 경로 →
그대로 뒀으면 40시간을 돌린 뒤 29개 전부 "자료 없음"을 받았다. 고침 = `--cand`에 없는 타깃은
`maintest.csv`(status=run)로 채점 · 대상 없으면 종료코드 2 · 채점 0개면 종료코드 3.
**교훈(반복됨): 이 파이프라인의 실패는 멈추는 게 아니라 조용히 비는 형태로 온다**
(7/27 a3m 질의행 사건과 같은 계열). 새 명단으로 도구를 처음 쓸 때는 반드시 스모크로 확인할 것.

**로그는 자동 저장된다**(`$DATA/logs/smoke_*.log`·`maintest_*.log`, 커밋 `f109ddf2`·`73b5d4a6`).
tmux가 죽어도 계산 결과는 디스크에 남고 `run_maintest.sh`가 완료분을 건너뛰므로 같은 명령으로
이어서 가면 된다. 판정 출력만 다시 보려면 `bash run_smoke_maintest.sh --recheck`(GPU 미사용).

**명단 = 54종 판독 → 본 검정 30개(8ulr 완료 → 신규 29개 × 32회 = 928회) · 무반응 24개.**
아무것도 고르지 않는다. 세트 3(10)·세트 4(49) 둘 다 결과와 무관하게 만든 명단이라 전수여야
빈도가 나온다. **"세트 4에서 7개를 고른다"는 이전 계획은 폐기**(예산이 없다고 본 판단이었음).
빈도 네 개를 각각 낸다 — 세트 3(분모 10) · 세트 4(분모 49) · 다양성있음 층(16) · 다양성없음 층(14).

**⭐ Neff80이 핵심 축이다(원 줄 수 아님).** 실측: 9zdu 27,944줄인데 **Neff80 28**(코로나 유행기
서열이라 거의 전부 중복) / 8k3k 141줄 Neff80 29(줄 수는 200배 적은데 유효 깊이는 같다) /
**8ulr 9,171줄 Neff80 3,051**(유일하게 진짜로 다양). 8ulr이 44종 중 1개였던 게 우연이 아니라
**입력의 성질** 때문일 수 있다.
→ Neff80(rung0) ≥ 50 = **다양성있음 층 16개**(Env 2·HA 8·C 5·8ulr), < 50 = **다양성없음 14개**(전부 RBD).
층은 성적이 아니라 입력으로 나눈 것이라 사후 기준 이동이 아니다. **제외하지 말고 층별 빈도를 따로 낼 것.**
⚠️ **교란: 층 ≡ 항원 계열**(다양성없음 14개가 전부 SARS-CoV-2 RBD). 층 안에 Neff 변이가 없어
구조적으로 분리 불가 → 층 차이를 기제로 쓸 때 반드시 교란을 함께 적을 것.

**세 갈래 결론(사전 확정):** ①다양성있음에서만 나옴 → "조성이 결합자리를 정한다, 단 조성에 실제
다양성이 있을 때"(+교란 명시) ②양쪽 다 → 조성 다양성이 기제가 아님 ③양쪽 다 안 나옴 → 8ulr 단일
사례로 축소, 헤드라인 음성.

**판독 규칙 보정 2건(둘 다 실행 불가능한 선택지를 막는 것 — 판정 기준 recall≥0.40은 불변):**
① **rung0 제외** — rung0은 make_composition_reps.sh가 그대로 대조군(seedfull.a3m)으로 쓰는 원래 MSA다.
거기서 조성을 재추첨하면 전체에서 전체를 뽑아 여섯 조성이 같아지고 조성군=대조군이 된다(9zdu에서 실제 발생).
② **Neff80 축·층 표시 추가** — 표시만. 규칙 불변.
✅ 규칙이 8ulr에 rung2·1,746서열·Neff80 702를 골랐다 = 확정 실험이 쓴 깊이와 정확히 일치(방증).

**실행 순서 Env → C → RBD → HA** (`run_maintest.sh` 기본, ORDER로 변경 가능).
①Env 2개(~6h)=8ulr과 같은 계열·같은 regime, 재현을 가장 직접 물음 ②C 5개(~8h)=적용 범위
(명제가 "조성이 자리를 정한다"로 바뀐 뒤 C군은 대조군이 아니라 범위 자료) ③RBD 14개(~11h)=
다양성없음 층 완결 + 세트 3 빈도 확보, 빠름 ④HA 8개(~17h)=가장 느림.

**⚠️ 스모크 없이 시작하지 말 것 (`run_smoke_maintest.sh`).** 군마다 코드 경로가 갈린다 —
항원 사슬 수(RBD·C 1개 / HA·Env 일부 2개), MSA 규모, 고른 칸 깊이(4줄~11,017줄), 과대표집 정의(C군 없음).
**8u44_ST가 본 검정 30개 중 유일한 다중사슬 항원(A|B)** → 반드시 스모크에 포함(자동 선정됨).
Phase A(GPU 없음) = 조성 개수·**md5 구별성**(같으면 이질성 검정 무의미)·질의행 일치·사슬별 폴더 /
Phase B(GPU 1회) = 산출물 5개·**입력 JSON의 항원 사슬 전부에 MSA 실림**·정제 후 질의행·run.log
MSA 경고·채점 값. 실패 시 exit 1. `make_composition_reps.sh`에 `GEN_ONLY=1`(조성만 만들고 GPU 전 종료) 추가됨.

**결과 해석 시 따로 표시할 복합체:**
- **9azt_HL** rung9·4줄·Neff80 3 → 사실상 단일서열. "MSA 조성"이 아니라 "어느 4개가 뽑혔나"
- **8siq_HL** 인기자리 겹침 0.91→0.05(최대 이탈). 다양성없음 층인데 크게 움직임 = 층 해석의 반례 후보
- **9azr_HL** 사다리 rung4에서 5/5인데 조성 재추첨 0/40 = 사다리 최고 칸 착시의 실증. 단 이번엔 rung1이라 같은 조건 재검증 아님
- **8p5m_GL** recall최고 0.34·DockQ최고 0.07(전 칸) → **역대조로 못 씀**(문턱 문제 아님을 두 지표로 확인). 실패로 계상
- **8k3k_D·8k46_I는 nanobody** → 세트 4를 "전부 Fab"으로 적은 이전 문서는 부정확
- 무반응 24/54 자체가 보고 대상(Protenix가 절반 가까이에서 결합자리를 전혀 못 찾음)

**판정:** 복합체별 이질성 정확검정 → 유의 개수. 실제 n으로 계산 — 세트 3(검정 6개)에서 2개 유의 = 0.033,
3개 = 0.002. 주 결과는 Fisher 결합. 빈도는 명단별 분모로 따로.

**신규 도구(전부 커밋됨):** `run_smoke_maintest.sh` · `run_maintest.sh` · `prep_pick_depth.py`
(Neff80·층·다지표 참고표·--only 누락 경고·--metric 민감도) · `make_composition_reps.sh GEN_ONLY`.
기록 = 노션 **인수인계서 Ⅱ** 6.2.1~6.2.4 · 6.5.1~6.5.2.

⚠️ **인수인계서가 2026-07-28에 둘로 분리됨** — Ⅰ Consensus Docking(모델 성능·채점·pose 선택·물리 도킹)
/ Ⅱ MSA Bias(이 실험). 세트 번호 1~5를 공유한다. 구판은 `[구판·보존용]`으로 제목이 바뀌었다.

---

Consensus docking(메인, [[bk-summer-2026-project]])의 "왜 특정 모델만 성공하나" 분석에서 파생된 서브 실험 방향. 2026-07-13 세션에서 구체화. 과대표집 prior 제안서(`Idea/연구제안_과대표집prior_novel항체도킹_BK2026.md`)의 MSA-채널 arm에 해당하나, 아래 검증으로 **방향이 바뀜**.

**가설 진화 (중요 — 반드시 이 상태로 이어갈 것):**
- 원래(제안서): "항원 PDB 과대표집 → 모델이 modal 에피토프로 끌려가 atypical 항체를 mis-dock. MSA-free(tFold)만 탈출."
- **2026-07-13 반증:** 9K6J(SARS-CoV-2 RBD + P5-1C8 Fab) 에피토프를 구조(9K6J.pdb, 4.5Å 접촉)에서 계산 → **modal RBM/ACE2 자리(Barnes class 1)로 판명**(24 접촉잔기 중 RBM 18·ACE2 접촉 13; K417·F486·Y489·Q493·Y505; cryptic core class4 = 0). → **9K6J로는 에피토프-prior 가설 검증 불가**(진짜 자리가 이미 modal이라 "modal로 끌림"이 실패 원인일 수 없음). 제안서 Phase-1 게이트의 "RBM인데도 실패=재설계" 가지에 걸림.
- **전환된(더 방어 가능) 가설 = 항체-쪽:** 과대표집 항원이라도 어려운 건 항원 에피토프가 아니라 **novel 항체(CDR-H3) 모델링**. MSA는 CDR에 장님(동족 없음), ESM/LM은 봄 → tFold(항체 MSA-free by design, ESM-PPI; PMC12815924)가 이김. 함의: consensus는 "**항체가 novel한 케이스에 MSA-free/LM 모델을 더 실어라**".

**뒷받침 데이터(2026-07-13):**
- 과대표집 실측(RCSB): 9K6J 항원 Spike만 압도적(항체복합체 776) — 나머지 9타깃 항원은 ≤22(대개 ≤7). 즉 과대표집 class = **n=1**(일반화 불가).
- tFold 성공 = 9K6J **단 하나**(Fab 0.491·Fv 0.525, 5/5 seed robust). MSA 모델 성공은 불안정(Boltz Fab 1/5 운, Protenix Fv만). Protenix는 희귀항원(B7-H3·Nectin-4·PlexinA1·uPAR)을 raw 정확도+inference-time scaling으로 이김 = 성공 이유가 tFold와 **직교**.

**결정적 다음 진단(먼저):** 서버 `runs_diverse/9k6j*` 예측 pose로 **모델별 예측 에피토프** 계산 → 실패한 MSA 모델이 RBM 근처(맞는 patch·틀린 pose = 항체/pose 문제 → Fork B)인지 엉뚱한 patch(에피토프 선택 문제 → Fork A)인지가 갈림.

**서브 실험 3군 설계:**
- G1 tFold-Ag (MSA-free + ESM 보상)
- G2 Chai-1 noMSA (MSA만 제거, 보상 없음 = ablation) — ⚠️ tFold와 같은 종류 아님(folding까지 약화; 실제 9K6J 더 못 맞힌 기록 있음). **G2가 더 나빠지면 편향이 추론-MSA가 아니라 weight prior에 있다는 신호 → MSA-depth로 못 벗기고 template injection이 맞는 개입.**
- G3 Chai(+MSA)·Protenix(+MSA) 대조
- 4복합체: Fork B면 9K6J 적합(축=항체 novelty). Fork A(에피토프-prior)면 9K6J 부적합 → over-rep 항원+atypical 항체(RBD class 3/4=S309·CR3022형) 검색 필요.

**순서:** ①9K6J 예측-에피토프 진단 → Fork 확정 → ②4복합체 재선정 → ③template injection·MSA-depth ablation.

**2026-07-14 재프레임 (Consensus Docking 재시작 — go/no-go 결정 국면, 반드시 이 상태로 이어갈 것):**
- **프로젝트 생사를 MSA-bias에 걸지 말 것.** MSA-free 이점은 9K6J **단 한 타깃**에서만 빛나고 나머지 9개서 tFold 전멸 → 아키텍처(MSA 유무)의 타깃-레벨 일반화는 데이터가 지지 안 함(n=1). MSA-bias는 "왜 실패하나"의 **논문 서사**로만 쓴다. **프로젝트 GO 근거 = ① 단일 모델 승자 없음(타깃마다 다른 모델 승) = consensus 여지 + ② ipTM 신뢰불가(9W43 ipTM0.854/DockQ0.067; 9K6J 정답pose ipTM<오답) = re-ranker 자리.**
- **grounded 핵심(레포 결과파일):** `diagnose_9K6J.txt` — 항체 내부 CA-RMSD가 Fv 전부 <1Å, Fab 최대 3.7Å = **어떤 모델도 항체를 '못 지어서' 실패하는 게 아님. 실패=배치(placement).** → 정답 pose가 앙상블 안에 도달가능 → consensus가 손댈 수 있음(단 n=1, 일반화는 Phase 0에서 확인). `where_docked_9K6J.txt` — 실패한 MSA 모델들은 면역우세 483-487 패치로 수렴, tFold만 진짜 RBM recall 0.95-1.0. `merged_4.txt` — 9K6J Boltz-Fab 0.733·Protenix-Fv 0.841·tFold 0.50-0.75 / 9Y0A Chai만 0.649 / 9LY5 Protenix만 0.874 / 9VXL 전멸.
- **지표 변경(2026-07-14 결정):** 주 성공기준 **DockQ ≥ 0.49 (CAPRI 'Medium', DockQ v2 legend 근거)** — 0.23(Acceptable)은 너무 낮음(동네만 맞음). **3-tier 리포트(0.23/0.49/0.80)**. **RMSD를 DockQ와 동등 지표로 이원화**: ①항체-내부 CA-RMSD(fold 품질) ②Ligand-RMSD(배치 품질) + ③epitope recall(region). 실패가 placement면 consensus 가능(GO), fold면 불가(NO-GO) → RMSD가 그 축을 가름. **⚠️0.49로 올리면 Protenix가 generalist화 → oracle−best_single 격차가 프로젝트 go/no-go 핵심 수치**(격차 미확정 = merged 6타깃 미계산 탓).
- **Phase 0 게이트(진행 중, 새 co-folder 예측 0개·서버안전):** `scripts/run_phase0.sh`(커밋 9563f76) = diagnose+where_docked+merge_dockq를 **10타깃×Fab/Fv 전부**에 실행. 산출 = diagnose_*·where_docked_*·merged_all.txt. **go/no-go 판정 기준:** GO = 0.49기준 oracle−best_single ≥2~3타깃 AND 실패다수가 placement(fold양호) AND ipTM-top이 oracle못고름 / NO-GO = 격차≤1(Protenix가 거의 다 먹음) OR 실패다수가 fold OR union-failure(9VXL·9MQR·9W43류)가 다수. **물리모델(HADDOCK)·MSA-free 신규모델·contact제약·Foldseek 정량은 전부 go 이후(Phase 1-3).**

**2026-07-14 Phase 0 결과 → ✅ GO 확정 (명제 = "재랭커가 망가진 ipTM을 이긴다"):**
- Phase 0 완주(10타깃×Fab/Fv, 5모델 각 5 pose). **merged oracle: 0.23→7/10, 0.49→3/10(9K6J·9Y0A·9LY5), 0.80→2/10.** best single(Protenix) 0.49→2/10 → consensus 격차 얇음(+1, 9Y0A=Chai단독). **"consensus>단일모델"은 약하므로 명제로 쓰지 말 것.**
- **결정타 = ipTM이 pose 선택기로 고장남(summary.csv 496 pose, ipTM 완전 채워짐).** cross-model ipTM-pick(global·pair **둘 다**): Fab 0.49 **1/10** vs oracle 3/10, 0.23 **2/10** vs 5/10, **regret 0.19–0.22**. per-model Spearman(ipTM,DockQ)=tFold 0.035 ~ Protenix 0.78(모델간 비교 불가). within-model best-pose 선택 1~5/10. 사례: 9K6J Boltz가 0.652·0.023 pose 둘 다 뽑는데 ipTM은 0.023 고름; 9Y0A 정답 Chai 0.662인데 ipTM은 Boltz 0.042 고름.
- **GO 근거 = ①ipTM 1–2/10 → oracle 3–7/10 격차 = 재랭커가 먹을 땅(크고 명확) ②실패 대부분 placement(fold 양호) → 정답 pose 앙상블에 존재 ③실패가 구조적·고칠 수 있음(글로벌 ipTM이 국소 상보성 못 봄 = 우리 설계가 겨냥한 것).** consensus=후보 풀 제공 / 재랭커=선택, 역할 분담.
- **한계(정직):** min-DockQ 기준(merged면 절대수↑, 랭킹결론 불변) / 천장 = 9VXL·9MQR·9W43 아무도 못 풂 / 9YC6 fold문제(abRMSD 8.5Å) / 재랭커가 실제 oracle 근접하는지는 Phase 1에서 증명. **9KKJ Fv merge 버그 fix 완료**(n_heavy를 chains[0](9KKJ선 항원)→항체 heavy 길이로, commit 102df6e) → Protenix Fv 0.642=Medium → **merged oracle@0.49 = 4/10**(9K6J·9Y0A·9LY5·9KKJ), 격차 여전 +1.
- **모델별 native 신뢰도 가용성(Phase 1 설계용, 정정됨):** Boltz·Chai·Protenix = global+pair ipTM+PAE. AF2-M = global ipTM + PAE(pair 없음). **tFold = global ipTM + pLDDT를 PDB REMARK 250 헤더에 냄(⚠️정정: "신뢰도 없음"은 틀림), 단 pair-ipTM·PAE 없음 → min-edge PAE만 불가.** 논문 tFold ipTM-DockQ Pearson r=0.77(그들 데이터)이나 우리 hard 세트선 Spearman 0.035(tFold가 9/10 실패→변량 없음). → **재랭커는 native 신뢰도에 안 기대고 구조-유래 피처(접촉그래프+PLM, 좌표만 있으면 5모델 균일) 중심**, min-edge PAE는 있으면 쓰는 피처.
- **역할배정 원칙(타입 무관 일관):** 사슬 heavy/light/항원을 길이·위치가 아니라 **chains.json 서열에 정렬해 배정**(diagnose·where_docked가 이미 이렇게 함 → 9KKJ서 안 깨졌음; merge_dockq만 길이기반이라 튀었음). Phase 1 빌더에서 merge·diagnose·재랭커 전부 이 서열기반으로 통일.
- **다음(Phase 1) = 재랭커 구축 + 데이터셋 빌더**: pose별 merged DockQ(라벨)+피처+ipTM 한 번에(서열기반 역할배정). 계면 접촉그래프 + PLM(ESM-2·ProstT5·SaProt) + min-edge PAE, **모델 간 정규화**로 ipTM 대체. 벤치 = diverse 10 + hackathon 5. 리포트 = 0.23/0.49/0.80 + epitope recall + RMSD 이원화.

**⭐ 2026-07-15 대전환 — RBD off-hotspot 세트 + MSA-깊이 인과실험 (반드시 이 상태로 이어갈 것; 위 재랭커 방향과 병행하는 새 서브축):**

세트 = RBD hotspot(RBM E484/F486/Q493, 잔기 483-496) **밖**에 붙는 항체 10개(9K6J류 off-hotspot; `scripts/select_rbd_subdominant.py`가 hotspot겹침 낮은 순으로 선별; 8SDH=다중카피[RBD2+Fab2, 사슬짝 어긋나 native 에피토프 n=0]로 드롭). 5 co-folder(Boltz·Chai·Protenix·AF2-M·tFold) × 5 pose 완주(`runs_rbd/`).

- **tFold KeyError('\x00') = 우리 스크립트 문제(tFold 무죄).** run_tfold의 a3m 선택 `find | head -1`이 colabfold 하위 `A_env/` 중간파일(117서열 불량행)을 집음 — 최종 `A.a3m`(143 clean) 대신. `-maxdepth 1` + A.a3m 우선으로 fix(commit 22bd8aa). 교훈: co-folder가 만든 a3m은 **최종 병합본만** 써야 함. ('같은 스크립트 다른 결과'의 답 = colabfold가 A_env/ 새로 만들며 find가 엉뚱한 파일 집던 것.)

- **결과① 5모델 전부 RBM hotspot 쏠림(where_docked+DockQ):** AF2-Multimer 최악(off-hotspot 10개 전부 true-recall 0·DockQ<0.06). Boltz·Protenix는 **core(369-385, CR3022형)** 회복(8SIT·9ML9·9ZDU DockQ 0.62~0.91, Protenix best). **cryptic back-face(462-520; 8SDF·8XSI·9ML8·9SBB)=5모델 전부 실패.** fold 정상(abRMSD<2-3Å)·placement만 틀림.

- **결과② tFold도 쏠린다(MSA-free-항체 안 통함):** tFold 예측 에피토프 n=44~76 diffuse(고분산, 5 seed 제각각), DockQ 0.02~0.22, true-recall≈0. **예측의 RBM% = 53~94%**(= AF3들과 같은 hotspot) + scatter, 진짜자리 부재. → "핫스팟+진짜 혼합"이 아니라 "핫스팟 쏠림(AF3와 동일)+노이즈". tFold 역할 = "항체 MSA-free만으론 편향 못 벗어남"을 보인 **증거(임무 완수)**, rescue 차량 부적합(약함+편향).

- **⭐결과③ MSA-깊이 인과실험 = 핵심 발견:** Boltz에서 다른 조건 고정하고 **항원 MSA 깊이만** 143→8→1(single-seq)로 낮춤(`subsample_a3m.py`·`make_boltz_msa_depth.py`·`run_msa_depth.sh`, 항체는 single-seq 상수, run-dir `runs_msad_<D>`). off-hotspot **rescue**: 8SIT(core) 0.00→**1.00**(d8/d1), 9ML8(cryptic) 0.00→**0.41**(d1, 462-469 찾음), 8XSI 근처이동, 8SDF 실패. **항원 MSA를 얕게/없애면 RBM 쏠림이 풀리고 진짜 에피토프를 찾음(4개 중 2-3; 전10 재현 실행 완료, 표 분석 중).** 비-단조(중간깊이 32·8 노이지 → "깊은 MSA=쏠림 / single-seq=진짜or노이즈" 이분법). **이게 tFold를 설명 = tFold는 '항체' MSA-free지 '항원'은 깊은 MSA 유지 → 그래서 tFold도 쏠림. 뺐어야 할 건 항원 MSA.**

- **재정립 결론(이전 결론 절반 교정):** 편향은 (상당부분) **항원 MSA**에서 온다 — 2026-07-14의 "편향=아키텍처, MSA 무관"은 **절반만 맞음**. single-seq 항원 MSA가 off-hotspot 배치를 rescue(**target-dependent** — 8SDF처럼 flat/polar/cryptic 계면은 실패=기하 문제). tFold의 "MSA-free"는 **엉뚱한 사슬(항체)**을 뺀 것.

- **문헌 대조(2026-07-15 웹서치; 반드시 인용 → [[guan-keating-msa-docking-bias]]):** ⭐**Guan & Keating 2025 (Protein Science 34(11):e70331; PMC12518507)** = 가장 가까운 선행. AF2-M/AF3/Boltz/Chai가 학습셋 과대표집 결합부위로 편향(AF3 실패의 **42%**가 native보다 학습셋-빈도 높은 계면에 도킹); **paired-MSA 셔플해도 도킹 무변**(공진화≠도킹), **unpaired MSA가 결합부위 정보 나름**(MSA 빼면 결합부위 바뀜). 단 **펩타이드**. MSA subsampling=**conformation**용(Del Alamo 2022 eLife/Wayment-Steele 2024 Nature). 항체-항원 MSA 무용(Gray lab AF3 항체 mAbs 2025/bioRxiv 2024.09.21.614257). MSA-free 항체모델(IgFold·tFold-Ab·RaptorX-Single PNAS 2024·xTrimoABFold)=**전부 항체 MSA만 뺌**(항원 유지). 물리하이브리드=AlphaRED·btaf129.
  - **우리 슬리버(novelty):** "**항원 MSA 깊이를 debiasing 노브로 써서 항체 off-hotspot 에피토프를 rescue**"는 미발표. Guan(펩타이드·MSA=결합부위 편향)+subsampling(conformation)의 **새 응용**. 방법 novelty 아님 → "Guan을 **항체-항원으로 확장** + 편향 나르는 게 **항원 MSA**임을 인과 분리"로 정직 포지셔닝.

- **✅ 확정 방향(2026-07-15, 사용자 동의) = "항원-MSA-깊이 debiasing" 메인 축.** 5단계: ①Guan식 편향 정량(over-rep 지도로 '예측계면 과대표집도 > native' 측정) ②MSA 분리 **2×2(항체×항원 full/empty)** + **unpaired MSA 셔플(Guan 방식)** = Phase B ③MSA-깊이=debiasing 노브(143→1, target-dependent 특성화) = Phase A ④기전 probe(항원 MSA 보존/공진화가 over-rep 패치에 몰리나) ⑤실용 레시피 + 물리 rescue 대비. 확장 모델 = Boltz(완료)+Protenix(best single, custom-MSA 포맷 확인 필요)+Chai.
- **HADDOCK = 사수 무조건 지시라 먼저/보조 arm.** gate-first(HADDOCK 짓기 전 오라클-ICF steering이 강한 모델[Protenix]을 진짜자리로 옮기나부터; ICF 텐서=zeros(L_H+L_L+L_A), icf[L_ab+r]=1, 0-based, tfold_ag_ppi). **critique 핵심(워크플로 2회):** 조건부 GO·전체구축 NO-GO; novelty=0(**BepiPocket/DiscoPocket bioRxiv 2025.09.17.676770**이 '에피토프→pocket제약→co-folder rescue' 직접 선점); 실패의 곱셈(P(후보에 진짜)×P(ICF override)×P(재랭커 선택)); HADDOCK ab-initio top-10 acc **0.31**; ⚠️**RBM=ACE2모티프=물리적으로도 stickiest면 HADDOCK도 hotspot 붕괴**(물리 de-confounder 무력화, modal≠sticky 계산 확인 필요); 제약은 강한 모델(Protenix/Boltz pocket·contact)에, tFold 단독 아님; ablation A0-A4(A2=DiscoTope=편향 양성대조, A3=native=천장); 지표 3종(DockQ 3-tier+epitope recall+**hotspot-drift 감소량**). ⚠️workflow synthesize/finalize는 정책 오탐(병원체+구축 조합)으로 재실행해도 실패 → 설계는 메인루프가 합침. Notion Study B에 52블록 종합 기록 완료.

- **⭐⭐ 2026-07-15(오후) 상위 프레이밍 확정 = "조건부 게이팅 에이전트"(사용자 명시, 반드시 이 상태로 이어갈 것):** 상위 서사 = "**모든 항원-항체 복합체에 균일 적용이 아니라, MSA 편향 받을 것 같은 케이스를 감지→추론 때 편향 덜 받는 방식(MSA 감소)으로 돌려 rescue하는 게이팅/에이전트**", 여러 docking 툴(consensus) 위에서. **MSA 감소 = 편향-의심 케이스를 어떻게 잘 rescue할지 정교화하는 생성-단계 '한 수'**(균일 X, 조건부). 방법 novelty는 사용자가 명시적으로 "상관없다"(CLAUDE.md 가중치와 일치).
  - **파이프라인 3단:** ①생성(여러 co-folder × {정상 MSA, 감소 MSA}) ②게이팅(이 케이스 편향 받나? → 감소버전 배치·플래그) ③선택(재랭커 = **후반부 숙제, 지금 핵심 아님**). 사용자 논리: "재랭커만 핵심이면 케이스를 나눌 이유가 없다" → **케이스 층화 실험(rescue 8SIT·9ML9 / 역방향 8P5M / 저항 8SDF·8XSI)이 곧 게이트 만드는 연구**.
  - **게이트의 정직한 한계 = triage/의심이지 최종판정 아님:** 추론시점 정답(진짜 에피토프) 모름 → 신호로만 의심(Foldseek 인기자리 겹침 / 정상 vs 감소-MSA 불일치 / 툴 간 불일치). **8P5M 문제**(진짜 RBM binder) = 게이트 신호가 8SIT(잘못 쏠림)와 구분 못 함 → 오탐 시 멀쩡한 것 망침(DockQ 0.74→0.06). → **의심 케이스엔 정상+감소 둘 다 생성, 최종선택은 재랭커로 연기.** 게이트=계산 어디 쓸지 triage + 불확실 플래그. 과대표집prior 제안서(`Idea/연구제안_과대표집prior_...md`)의 Foldseek 라우터와 수렴.
  - **연구 무게중심 = "칼을 날카롭게"가 아니라 "언제 칼 쓸지(게이트) + consensus".** MSA-수술 정교화(공진화 보존+편향만 제거)는 부서지기 쉬움 → 무딘 칼(depth 감소)이 이미 4/10 되니 그걸 생성 수로 쓰고, 지능은 게이트·오케스트레이션에.
  - **기법 탐색(2026-07-15 워크플로 9축) 결론 → [[msa-debias-technique-landscape]]:** 방법 novelty=0(전부 선점: IPW·LEACE·APC·robust-PCA·MACR·DeepRank-Ab). frozen+8주 개입지점 = 입력단 MSA / 출력단 점수 / 우리 재랭커뿐(co-folder 재학습 out). 최강 전이=IPW(recsys, 사후 pose 재가중, Guan 무관). 최저위험=APC식 인기배경 차감(출력단). 함정=LEACE(collinearity로 depth-down과 같은 병+원리적 위장으로 오히려 더 위험), BSS/ICA(co-folder가 covariance 안 씀+독립가정 깨짐 → 진단피처로만). **모든 MSA-수술 위협 = collinearity(인기⊥공진화 거짓)+Guan(편향이 MSA 아니라 가중치일 수도)+순환+"debiasing 아니라 confidence↓ 샘플링" 대안(미통제).**
  - **다음 실험(캐시 예측물, GPU 0, 밤샘 GPU와 병렬) = 게이트 연구 2갈래(Fork A/B를 재랭커 판정이 아니라 게이트 연구로 재배치):** **(A) 생성기 값어치** = 감소-MSA가 정상 앙상블에 **없던** 정답 pose를 만드나?(8SIT yes: 정상 MSA true-recall 0%) → 5모델 일반화 확인 = "감소-MSA를 파이프라인에 넣을 근거". **(B) 게이트 신호** = 정상 vs 감소-MSA 불일치 / Foldseek 겹침이 편향케이스(8SIT)를 튼튼케이스(9ZDU)와 구분하나? 8P5M 오탐 어떻게 뜨나? = 게이트 실현성 판정.
  - **5모델 MSA-depth 확장(진행/커밋):** Boltz(완료 d143/8/1) + Protenix(unpairedMsaPath로 항원 depth 주입, 스모크 OK N_msa=8, 전체 실행 중) + tFold(항원 a3m depth=제일 clean, d8/d1) + Chai(single-seq=d1 극단) + AF2(--msa-mode single_sequence=d1 극단). 러너 전부 커밋(`make_protenix_msa_depth.py`·`run_protenix_msa_depth.sh`·`run_tfold_msa_depth.sh`·`run_chai_single.sh`·`run_colabfold_single.sh`), 통합표=`build_msafree_summary.py`. HADDOCK 10타깃 = CPU(ncores 14)로 병렬 진행 중.

- **⭐ 2026-07-16 5모델 완주 결과 = 정직한 하향조정(Boltz-only 낙관 → 5모델 현실; 반드시 이 상태로).** `report_msa_depth.py`(RBM%/recall, results/msafree_summary.csv 150행). **어제 "Boltz 3-타깃 airtight rescue"는 유효하나 5모델로 보면 상당부분 착시.**
  - **발견1(제일 튼튼) = MSA 감소가 RBM 인력 약화, 전 모델 일관:** depth↓ → RBM% 하락(AF2 server 100%→d1 9~65%). **8P5M(진짜 RBM binder) 역방향 = 5모델 전부 recall 하락**(Boltz .67→.33·Prot .67→.17·Chai .42→.33·AF2 .58→.33·tFold .58→.08) → "항원 MSA가 RBM 인력" 메커니즘 모델 불문 지지.
  - **발견2 = recall rescue는 좁고 모델특이:** Boltz만 깨끗(8SIT d143→d8 0→1.00). **Protenix 노이지**(d8=sweet spot, d1 fold붕괴). **tFold rescue 없음**(d8 RBM 63~100%·recall≈0, d1 붕괴=사용자 "항원까지 빼면 못 맞춤" 확인). **Chai single-seq 오히려 손해**(8SIT server .95→d1 .05). **AF2 RBM%만 낮추고 recall 여전 ~0**.
  - **발견3(중요 교정) = "core rescue"는 착시:** 8SIT·9ML9·9ZDU(core)는 **server(full-both MSA)서 이미 풀림**(Boltz·Prot·Chai .82~1.00). 어제 "0.05→0.82 극적 rescue"는 **항체-single 베이스라인(d143)의 인위붕괴 대비**였음. **진짜 방어가능 rescue = server 실패한 cryptic을 reduced가 건짐: Boltz 9ML8(0→0.41)·Protenix 8XSI(0→0.46)·9ML8(0→0.29).** 8SDF back-face 전 조건 저항.
  - **함의:** (a) 효과 케이스·모델 딴판 = **게이팅 필요성 강화**(8P5M엔 독·tFold 무효·cryptic 값짐). (b) **server vs d143 = 항체MSA+파이프라인 교란** → **2×2(Phase B) 필수**(d143이 core 깬 게 항체MSA 탓인지 파이프라인 탓인지). (c) 깨끗한 비교는 msa_depth 내부(d143→d8→d1)뿐, server 대비 금지.
  - **HADDOCK**: 10/10 완료, pose=`haddock/<t>/run/6_seletopclusts/cluster_N_model_M.pdb.gz`(gz압축). `report_haddock.py`(per-cluster 점수·RBM%·recall + co-folder divergence, clustfcc tsv의 score/cluster_id 파싱, gunzip). Notion Study B에 22블록 종합 기록 완료.
- **⭐ 2026-07-16 HADDOCK(물리 ab-initio) 결과 = 물리도 de-confounder 아니나 좁은 밴드서 값짐(반드시 이 상태로).** ab-initio(CDR active+항원표면 passive) 클러스터 = FCC 접촉유사도로 물리가 찾은 후보 패치, 점수=HADDOCK score(음수 낮을수록 좋음). **① 물리도 RBM 쏠림**(대부분 타깃 상위 클러스터 RBM100%, 점수가 RBM 끈적임에 속음) = "RBM=ACE2가 물리적으로도 stickiest" 확인 → **물리는 bias-free 아님.** **② core 다 놓침**(8SIT·9ML9·9ZDU·8SIQ = co-folder가 푸는데 물리는 RBM으로 감). **③ 좁은 밴드서 진짜 값짐: 9SBB(#1클러스터 rec0.45)·9ML8(#1 rec0.29)이 물리 top-score인데 co-folder 전부 실패, 8XSI best_single 0.55.** **⭐수렴: 9ML8·8XSI·9SBB = reduced-MSA(Protenix d1) AND 물리 둘 다(독립 경로)가 건지고 8SDF는 둘 다 저항** → "co-folder가 놓친 cryptic/other 일부는 두 독립 경로로 발견"=consensus 근거. **Q2 패치-스캐닝 값 낮음**(ab-initio가 이미 패치스캔·HADDOCK점수가 RBM지배로 진짜패치 강제해도 선택X=생성O선택X·살릴 cryptic은 이미 top클러스터). → 기존 클러스터=생성기 + 계면 재랭커로 클러스터+co-folder 통합 재점수.
- **⭐ 2026-07-16 DockQ-확정 종합 + staged 파이프라인 정정(반드시 이 상태로).** `dockq_cofolder.py`(co-folder 전체 DockQ, H+L 병합)·`haddock_dockq.py`(HADDOCK)·`dockq_scorecard.py`(방법론별 분리 뷰). **recall이 여러 rescue를 부풀렸고 DockQ가 잡음.**
  - **⚠️ 파이프라인 프레이밍 정정(사용자 명시): 병렬 union 아니라 staged** — 전반부 **HADDOCK**(epitope·pose 후보 생성, MSA편향 없이 표면탐색) → 후반부 **Protenix/tFold**(후보 받아 정교 refine, 편향회피). 방법론별 분리 저장 후 조립. **"유니온 9/10"은 병렬 오해.**
  - **co-folder server DockQ(blind): Protenix가 core 챔피언**(8SIT 0.83·9ML9 0.87·9ZDU 0.91·8SIQ 0.84). **cryptic/other/8P5M/8SIS = 전 co-folder blind로 실패(<0.25).**
  - **MSA-depth DockQ: 좁고 Boltz-특이.** Boltz 8SIT(d8 0.82)·9ML9(d1 0.83>server 0.62 개선)·**8SIS(d1 0.44, 아무도 못 푸는데 신규)**; 역방향 8P5M(d143 0.74→d8/d1 0.06). **⚠️ Protenix msa_depth는 server보다 다 나쁨**(항체-single이 core해답 깸) → **rescue는 Protenix엔 안 나타남.** cryptic(9ML8·8SDF·8XSI) 감소-MSA로 안 풀림(recall 착시).
  - **DockQ로 살아남은 off-hotspot 감소-MSA rescue = Boltz 8SIS d1(0.44) + AF2 9ML8 d1(0.247 경계)뿐.**
  - **HADDOCK(전반부 후보) — rigidbody DockQ 10타깃 완결(results/haddock_dockq_rigidbody.csv)**: raw 풀 near-native = **8P5M 0.98·9ML9 0.99·8XSI 0.95·9SBB 0.95·8SIS 0.86**(cryptic/other는 co-folder 전부 실패인데!), 8SIQ 0.33·9ML8 0.22(경계). **seletop 점수-컷 손실 큼: 9ML9 0.99→emref 0.41·8SIQ 0.33→0.06·9ML8 0.22→0.06 = near-native 대량 학살(생성O선택X 확정)**. core는 co-folder(Protenix)가 blind로 품. **8SDF만 raw도 0.16=진짜 hard(아무도).** ⚠️9ML8 정정: emref 0.06라 착시라 했으나 raw 0.22=seletop이 죽인 것(완전 착시 아님).
  - **⚠️ 진짜 파이프라인(HADDOCK 후보 → co-folder guided refine)은 아직 미실험 = 다음 단계.** 지금까진 co-folder blind + HADDOCK blind만. Notion Study B 17블록 기록. patch-scan 도구(partition_surface·make_haddock_patch_restraints·run_haddock_patches) 커밋됨(POC/큰항원용, 이 세트엔 불필요=샘플링 아닌 점수 문제).

- **⭐ 2026-07-20 데이터셋-우선 전환 + A-vs-B 척추 확정(반드시 이 상태로 이어갈 것):** 사용자 지시 = "복합체 데이터셋 제대로 다 모으고 검증까지 마친 뒤 전 모델 한 번에 학습." → 생성 전 데이터셋 lock. **143은 앵커 아님**(한 항원의 full 행수일 뿐) → 사다리는 **항원별 full→single 상대정의, 축=Neff80(θ=0.8 per-residue Neff median, AF 방식)**.
  - **적대검증(워크플로)이 3그룹 설계 4구멍 지목:** ①MSA-깊이=공진화(fold)채널이지 prior채널(template+training membership) 아님 → depth-response가 fold붕괴 측정 위험(→ fold-안정성 control로 반박: co-folder 템플릿-프리 + Phase0 abRMSD) ②"dominant를 PDB빈도로 라벨"=순환(→PDB-독립 라벨) ③A(쉬움)vs B(어려움) 난이도 confound(→baseline 매칭) ④C그룹=novel-fold+얕은MSA+저해상도 교란덩어리+sweep range 없음.
  - **재설계 척추 = 같은 과대표집 항원 내 A(지배)vs B(비지배)** — 항원·항원MSA·fold·Neff range를 설계상 상수화 → confound 대부분 제거. **항원 3종 = RBD·HA·HIV Env**(사용자 확정). PDB-독립 dominance: RBD=RBM/ACE2·Barnes class, HA=head/stem(교과서 immunodominance=가장 깨끗), Env=bnAb supersite(⚠️PDB-빈발과 감염-면역우세가 갈림 → ENV_DOMINANT_MODE=structural 기본+한계 명시). C는 약한 보조로 강등.
  - **SAbDab 카운트(Mac scratchpad, `date`=진짜 구조일; SABDABdepo_date=DB재구축일2026 폐기):** post-2023-06(=Boltz-2 최신컷, 3 co-folder 공통 clean)·res≤3.5 distinct PDB. 과대표집 풀 천장: RBD 214·HIV Env 140·HA 115·Ebola 48·RSV 30 / lysozyme·PD-1·TNF는 post-cut 풀 2~10=탈락. **⚠️GPCR용 G-protein/Fab-fiducial 대량건 제외.**
  - **✅ 데이터셋 파이프라인 구축·검증·커밋(repo `consensus_docking/dataset/`):** build_manifest→**후보 375 PDB(RBD169[141paired+33nb]·HA115·Env91)** → fetch_structures(RCSB)→classify_epitope(접촉4.5Å→에피토프클래스→A/B). RBD=Wuhan auth직접(RBM 437-508 overlap≥0.5=A), HA=P03437 참조정렬(head range 56-306 비율≥0.6=A), Env=HXB2(P04578) bnAb class footprint. **알려진-에피토프 참조 6/6 정답**(6W41 CR3022→class4/B, 3GBN/4FQI stem→B, 2VIR/4FP8 head→A, 3NGB VRC01→CD4bs/A; Mac biopython 1.85). "no-RBD-contact"로 NTD/S2 결합 자동배제, HA head/stem 경계=mixed flag.
  - **다음:** 서버서 fetch(375, ~1GB)+classify → 실제 A/B 카운트 → 검증게이트(fold-stability·난이도매칭·homolog leakage·Neff80 사다리 빌더[143폐기]) → lock → 전모델 depth-sweep 배치. Chai depth 러너 신규 필요(make_chai=FASTA-only), Protenix=base(2021-09, 2025 leaky 금지). must-cite=[[guan-keating-msa-docking-bias]] NEFFy AsEP Barnes2020.

- **⭐⭐ 2026-07-22 Boltz 12단 depth-sweep 완주 + 정직한 재프레임(반드시 이 상태로 이어갈 것). 별도 자기완결 레포 [[bk-summer-2026-project]]의 `Feellived/msa-composition-bias`(guided는 consensus_docking 소관).** 47복합체×12rung(geomspace full→single-seq, Neff80축)×best-of-5, Boltz만 완주. 채점=`eval_dockq_sweep.py`(다중copy 오짝 버그 fix: native를 src_chains로 추출). 리포트·CSV 커밋 = `report/boltz_depth_sweep_analysis.md` + `report/dockq_sweep_boltz.csv`(재-투입용).
  - **블랭킷 A/B/C 가설 = 깨끗한 null(2 워크플로 13+7에이전트 적대검증).** Δmean(rung1-9−rung0) A−0.004·B−0.022·C+0.018 전부 Wilcoxon NS; B의 −0.022는 단일 8wpy 붕괴 산물(빼면 +0.014); **argmax-deep A19/20·B15/18·C8/9 동일 = best-of-N 표집 아티팩트(음성대조 C가 B와 동일 재현, 특히 9y0a rung1/6/10). B vs C 구별불가(P=0.92).** MSA축소는 정상작동(Neff 119→2). → "off-site rescue 블랭킷"은 지지 안 됨.
  - **✅ 재프레임(사용자 지시 = ipTM/pose/consensus 앵커 버리고 MSA-깊이 현상 자체를 연구 대상으로): "특성화 연구 + 방법론 발견".** 건질 것: ①재현성 있는 깊이-반응은 드묾(≥0.49=14/47, 진짜 지속전이 **1개**), ②방향이 복합체별로 갈림(깊은MSA 필수 vs 방해), ③**핵심발견=원인은 양(깊이) 아니라 조성(무엇이 빠지나)** — 8wpy 절벽이 Neff 19.4→18.0(거의 불변)에서 일어남 + 전이Neff 1~688 산발(특정값 비수렴).
  - **⭐ 기억할 case study 3(+대조):** **8wpy_AB(RBD,B)=간판·유일확정 지속전이**(rung0-2 0.85→rung3~ 0.04 절벽 후 7rung 고착=깊은MSA 필수). **9y0a_AB(C)=반전·유력**(full 0.06, neff688·28·2서 0.57/0.72/0.71 반복회복=과대표집 prior가 오답 유도, 인접지지0라 seed재현 필요). **8t4d_OQ(Env,A)=중간깊이최적**(neff319·154 인접 반등 0.39·0.19). 대조=8vye_AD(깊이완전무관 0.91~0.93). ⚠️경고사례=8y6a·9b7g(진폭 크나 neff1~2 단발스파이크=운, "진폭에 속지마라").
  - **8wpy 후속(사용자 지정 핵심 분석) = `ladders/8wpy_AB/A/rung2.a3m` vs `rung3.a3m` 서열 diff → 절벽서 빠진 동족서열군 + 에피토프 열별 보존도.** "조성 가설" 직접 검증. (a3m 서버에 있음, Protenix·Chai 다 돌린 뒤 feature와 함께.)
  - **Boltz MSA-의존 정확화:** Boltz는 아키텍처상 MSA기반이나 **이 데이터서 "항원 MSA 깊이"에 둔감**(8ulr neff3050→1에도 0.73~0.80; A그룹 Δ−0.004). Protenix는 강의존(consensus data: core full 0.83→d143 0.037 붕괴). → **Protenix 12단 스윕이 진짜 시험대**(A를 positive control로 "이 파이프라인이 깊이효과 잡나"부터; B가 flat이면 강기각, 중간peak면 첫 신호). Chai도.
  - **⚠️ leaky Protenix 정정(중요):** consensus_docking의 Protenix(Exp1~4 blind·depth·guided 전부)가 **`protenix_base_20250630_v1.0.0`(학습컷오프 2025-06-30=leaky)** 사용. 타깃 10개 중 8개가 2025-06-30 이전 공개=오염(9SBB·9ZDU만 clean). leakage-free=`protenix_base_default_v1.0.0`(2021-09-30) 또는 미공개 `protenix-v2`(다운로드 403). **run_protenix_*.sh 기본값 v1.0.0으로 수정·커밋(Feellived/epitope-guided-docking).** 사용자 결정: 학부프로그램·비출판이라 옛것 재실행 안 함(단 기록엔 실제모델 정직표기, 안 쓴 걸 v1.0.0이라 바꿔쓰기 금지). 새 3타깃 guided(8P5M·8SIQ·9ZDU)만 v1.0.0로.
  - **다음 우선순위:** ①Protenix 12단(A=pos control) ②Chai 12단 ③seed 5→25~50(best-of-N 분리, 분포기반) ④8wpy 조성 diff ⑤pose-level 에피토프 recall. **feature 실현성은 Protenix·Chai 후로 연기**(DockQ-반응=정답필요 라벨용, GT-free 대리=깊이 구조발산/신뢰도 스윙을 pose원본서 별도계산, pose선택기 아니라 불확실성 triage 플래그로만).
  - **Notion:** Results 페이지(3a5ac1a3-a28a-80c3-9a69-df5a2f0e4dab)=결과·해석, Progress 페이지(3a3ac1a3-...863836)=날짜별 하위페이지 로그. (이 depth-sweep 결과·Progress 기록은 아직 미반영, 다음 작업.)
  - **⭐ 2026-07-22(오후) 사용자 방향 재정립 = "깊이 자체 아니라 특정 MSA 서열(조성)" + sweet-spot 렌즈.** 사용자가 "8wpy만" 강조에 반발(맞음) → EDA 완화: **깊이-반응 복합체 16개**(sweet-spot=reduced>full 14 + deep-required 8wpy). **EDA 그림 커밋**(`report/figures/`): heatmap_dockq_depth·**heatmap_normalized_shape**(행정규화=sweet-spot 위치 드러남, peak가 r1~r10 산발=단일최적깊이 없음=조성 방향 재확인)·lines_rung_unified(x=rung 통일)·small_multiples·sweetspot_scatter·peakrung_hist. 추천 미완 그림=**epitope-recall vs depth**(pose 필요=2a 위치편향 증명) + seed-복제 에러바.
  - **⭐ 통제 실험 계획 확정·커밋(`report/msa_composition_control_plan.md`) = "깊이(개수) vs 조성(특정서열)" 인과 분리:** ⚠️현 사다리=rung별 독립 랜덤 subsample이라 개수/조성 못 가름. **Exp1 seed복제**(같은 개수 N seed→분산 작음=깊이/큼=조성, 첫 필터) → **Exp2 nested 사다리**(rung_{k+1}⊂rung_k, 전이=빠진 특정서열) → **Exp3 LOCO/AOCI**(클러스터 leave-out/add-in=인과 클러스터 ID) → **Exp4 특성화**(taxa·%id·에피토프 열·공변→mutation/DMS 표적). Readout=DockQ+epitope recall. **★사용자 지시: 통제 실험은 Boltz/Chai/Protenix 중 반응성 최고 1모델만**(Exp0=3모델 스윕 후 sustained전이수·진폭·스파이크비율로 선택; Protenix 유력). 구현대상=build_ladder `--nested`·LOCO/AOCI 서브샘플러·epitope recall 파서.

- **⭐ 2026-07-23 iDist 과대표집 채점 + 검정2 + 적대적 재검증 = "과대표집→깊이-rescue" 명제 하향(반드시 이 상태로).** iDist(PPIRef; 계면 6Å heavy-atom, 20-dim AA조성 임베딩, near-dup 0.04) 채점 완료(NaN 버그=embed_parallel 실패계면 NaN행 → 필터+반경 실측화(0.02~0.06)+min_dist·mean_knn·frac_ndup·nearest+npz캐시; `overrep_idist.csv` 40행, n_ref RBD442·HA132·Env209).
  - **서술적 발견(과대표집≠면역우세):** 신규(과소표집) 계면 = 거의 전부 Env 측면(MPER·FP·gp120gp41: 8vfv 0.09·8tq1·8tkc·8tjr 0.27). RBD 포화(A≈B≈0.92=최다연구+20-dim 조성 해상도한계). HA 역전(stem B가 head A만큼 과대=범용백신 표적이라 많이 침착). → PDB 과대표집=연구량 반영.
  - **검정2(과대표집 vs 깊이지표) = null.** pooled(n=38) 전부 |rho|<0.24. family별 HA만 겉보기 유의(frac vs full-recall −0.56 p0.02).
  - **⚠️ 적대검증(9에이전트 워크플로 stats/confound/generalize, 잔여신뢰 0.30 + 결정테스트)으로 HA 신호 대부분 인공물:** ①변화점수(gain=max−full_recall이 −0.56을 +0.61로 기계증폭; full_recall 통제 시 +0.44) ②best-of-N(11rung×5pose max) ③**고정-rung LEVEL 테스트=과대표집 모든 깊이서 recall 낮음(−0.41~−0.53) = "깊이-rescue"가 아니라 "어느 깊이든 어렵다"** ④다중비교(~150 test, Bonferroni·BH 미통과; n14; 3가족 중 1개) ⑤stem/head 부분교란(HA A/B=head/stem; stem 통제 후 −0.56→−0.43 약화·생존; n_ag 통제 강화 +0.64=suppressor 신호). LOO-robust·n_ag통제는 인공물서도 나오는 성질=근거 아님.
  - **판정 = REFINE: iDist=탐색적 공변량, "과대표집→깊이-rescue" 미지지.** 대안평가: Foldseek=**skip**(가족내 더 거침), iAlign=RBD포화 해상용 **do-scoped(단 rescue 실재 확인 후)**, positional epitope Neff=cheap(neff_col 슬라이스), direct-count(Guan식 epitope-freq)=defer. novelty=0(iDist=Bushuiev2024·iAlign=Gao&Skolnick2010·Foldseek=vanKempen2024).
  - **⭐ 더 큰 함의 = rescue 실재 자체 미확인:** 실패복합체가 11rung 중 노이즈로 0.3 넘을 P≈0.73, 관측 rescue 8/11=0.73 = 노이즈밴드 안. → 과대표집이든 Neff든으로 rescue 설명 전에 **rescue가 진짜인지부터**(seed-복제). 통제실험(특히 Exp1)이 **linchpin으로 승격.**
  - **★2026-07-23(후속) 순열검정(Ojala/Porter식, `analyze_perm_null.py` 커밋)으로 블랭킷 rescue = best-of-N 잡음 확정**(rescue 11≈우연 9.5 p0.30 · mean-gain null 밴드 안 · 평균 rescue 소멸 1<2.5). 문헌 정합: Porter그룹이 AF-cluster를 noise+학습셋암기로 반박(RfaH 무작위 24%>클러스터 4%) + McCoy2024가 Ab-Ag서 Neff-DockQ 무상관 = 같은 결. 왜 재현 안 되나=논문은 '다양성/best-of-N'을 봤지 '정답 에피토프 정조준'이 아님(AFsample2: N=1이면 오히려 나쁨).
  - **개별 case = 9y0a(bistable) 단일사례를 seed-복제 재현게이트(무작위-subsample 밴드 넘나 + 진짜 에피토프로 가나)로만** 진행 → 통과 시 nested/LOCO 기전 + 3모델 전이 + 암기배제, 실패 시 정직한 negative. arc = 선행 null 재현 > 단일복합체 존재증명 > 기전 > 조건부 게이팅(일반화·보편 MSA감소 주장 금지). de-noising 프로토콜 = `report/denoising_protocol.md`. must-cite: Chen2021(pass@k)·Porter2024/2025·McCoy2024·delAlamo2022·AFsample2(Kalakoti&Wallner2025)·Guan&Keating2025.
  - **⚠️ Proposal(Notion, v1 2026-07-20) 프레이밍 전환 필요:** v1의 H1-H4(항원 MSA↓→rescue·sweet-spot·depth-response=선택신호)가 **우리 분석으로 사실상 falsify됨**(블랭킷 rescue=잡음). 새 프레이밍 = **"MSA-의존이 가변성 큰 Ab-Ag 필드에선 편향/취약성의 통로 = 해소 대상, MSA-비의존 구조기반 새 피처 제안"** + **헤드라인=null(블랭킷)** + case-study(9y0a)=앵커된 존재증명. v1이 이미 인용한 McCoy2024(Neff-DockQ 무상관)·Guan&Keating(암기)가 새 프레이밍의 발판.
  - **통제 입력생성기 커밋(msa-composition-bias/pipeline, GPU불필요):** `run_seed_replicate.py`(Exp1)·`prep_ladder_nested.py`(Exp2 dropped.tsv)·`analyze_loco_aoci.py`(Exp3 %동일성클러스터 LOCO/AOCI). 순서=Protenix완주→Exp0(모델)→Exp1(rescue 실재+깊이vs조성). 대상=지지 8y6a·8kdm / 비회수대조 9evz·8tl5 / bistable 9y0a·8t49_PR. Results §7·§7-1 + Progress 2026-07-23 기록 완료.

- **⭐ 2026-07-23(저녁) eval_epitope_shift.py 확장(순위점수·중심거리) + '진짜 편향' 필터링 = 3중 확인으로 null 재확인·강화(반드시 이 상태로).** `eval_epitope_shift.py`에 DiscoTope식 순위점수(threshold-free)·DCC식 중심거리(Å, pose 내부 비교) 추가, oracle(recall최고 pose)/mean(5pose평균) 명시분리 후 Boltz 47타깃 재채점(`results/epitope_shift.csv` 504행).
  - **recall Δ뿐 아니라 관대한 연속지표(true_rank·dcc_true)로도 신호 없음**(RBD/B true_rank Δ=−0.063·dcc_true Δ=+1.45Å로 오히려 정답에서 멀어짐; family/B별 corr(이탈량,정답도달량)=RBD/B−0.43·HA/B−0.18·Env/B+0.32=뚜렷한 연동 없음).
  - **핵심 방법론 정정(사용자 지적) = A/B는 '모델 편향 여부'가 아니라 '정답이 인기자리인지'의 구조적 분류.** B라벨(정답=인기자리 아님) 18개를 전부 편향된 것으로 취급하면 안 됨 → full-depth에서 실제로 native_overrep보다 뚜렷이(excess≥0.3) 더 인기자리에 쏠린 **'진짜 편향' 케이스만 재필터 = 5개뿐**(8txu_HL·8q7s_H·8tx3_FK·9evz_HL·8tl5_IJ; 나머지 13개는 정답이 인기자리 밖이어도 모델이 처음부터 그렇게까지 안 쏠려있었음=벗어날 편향 자체가 없었음).
  - **⭐가장 강한 결과: 그 5개 중 4개는 MSA 깊이를 60배 이상 줄여도(neff 50~200대→1~3) 인기자리 겹침%가 사실상 무변화**(8q7s_H 0.990→0.989, 8tx3_FK 1.0→0.989(전구간 거의 항상 1.0), 8tl5_IJ 0.756→0.761, 8txu_HL 0.657→0.656). 유일하게 움직인 9evz_HL(0.793→0.721)도 recall은 끝까지 0. **무작위 대조 패치 불필요**(애매하게 흔들린 게 아니라 아예 안 움직여서 비교할 대상이 없음). recall(§7-2)·관대한지표·진짜편향필터 3가지 독립 방식이 전부 같은 결론 = 견고.
  - **Notion 반영 완료**(Results §7-3 + Progress 2026-07-23, 그래프 2종 fig_family_ab_delta.png·fig_biased_flat.png 첨부). 로컬 재사용 스크립트: `analyze_shift.py`·`analyze_shift_full.py`·`find_biased_cases.py`(scratchpad/analysis, Protenix 도착 시 그대로 재사용).
  - **⭐⭐ 사용자 확정 다음 계획(3단계, 반드시 이 순서로):**
    1. **Protenix 도착 시** — DockQ 다각도 분석(단순 추이→family요약→Spearman→oracle-vs-ipTM regret→순열검정→pass@k→짝Wilcoxon→혼합모델→TOST, 쉬운 것부터 어려운 것까지 전부) **+** 에피토프 위치 분석(recall·over-rep·순위점수·dcc·진짜편향필터링 전부, 위 Boltz와 동일 스크립트 재사용) → **이 결과로 Chai를 돌릴 가치가 있는지 판단**(Protenix가 Boltz처럼 flat이면 Chai도 별 의미 없을 가능성↑; Protenix가 진짜 깊이반응 보이면 Chai로 삼각확인 가치↑).
    2. **case-study 전환** — Protenix·Boltz 양쪽에서 가설과 비슷하게 움직이는 대표 항원-항체 복합체를 찾아, 이미 커밋된 통제실험 생성기(`run_seed_replicate.py`=Exp1 seed복제·`prep_ladder_nested.py`=Exp2 nested사다리·`analyze_loco_aoci.py`=Exp3 LOCO/AOCI) + MSA 돌연변이 실험으로 실재 가능성을 재확인 시도.
    3. **그래도 안 되면(=지금 Boltz 5/18 필터 결과가 시사하듯 낮은 히트율 예상)** — 전체 결과 정리 후 "**재랭커가 답이다**"로 결론 확정 → **데모**: 재랭커를 붙였을 때 특정 항원 세트의 예측이 어떻게 달라지는지 시연 + 기존 재랭커가 못 잡는 근본 문제(모델 아키텍처의 편향·leakage 등)를 명시적으로 개선하는 방안까지 시연.
  - Why 이 순서: 8/5 케이스 결과가 이미 히트율이 낮음을 시사하지만(0/5 완전이탈·1/5 부분이탈-정답미도달), 사례연구 표준 관행(KaiB/AF-cluster식)은 소수 사례+통제실험으로 존재증명하는 것이므로 시도할 가치는 있음 — 단 기대치를 낮게 잡고 시간을 크게 태우지 않을 것.

- **⭐ 2026-07-23(밤) "MSA 깊이 변화로 편향을 확증할 수 있나" 4각도 설계+적대검증 워크플로 = 전부 같은 벽(반드시 이 상태로).** 사용자가 "안 된다고만 하지 말고 설계해보라" 요청 → 4개 독립 설계(①germline/CDR 서열 축 ②대조쌍 데이터셋 확장 ③깊이-반응 모양/분산 분석 ④통계 실현가능성 감사)를 각각 만들고 서로 적대검증(웹검색·레포 실데이터 대조 포함, 8에이전트).
  - **공통 결론(4개 독립 도달) = A그룹(진짜 정답=인기자리)도 B그룹(편향)도 똑같이 안 움직이는 것은 "표본부족(검정력)" 문제가 아니라 "식별불가능성(identifiability)" 문제 — 표본을 아무리 늘려도 안 풀림.** ①은 "이미 아는 정답이 서열로도 그럴듯한지"만 확인(모델 행동 검증 아님, 순환·클론중복·재현안됨 문제도 발견) → 정성적 대조표로만 격하. ②는 관측 효과크기(≤0.02)로 8주 내 가능 표본이 필요량의 10~15%뿐, 게다가 "공짜 사전점검"이 이미 순열검정으로 부정적 결론 나 있음(우선순위가 seed-복제 통제보다 뒤처짐). ③은 "암기=뚝끊김/증거=매끄러움" 전제 자체가 미검증 가정이고 분산추정이 평균추정보다 원래 더 잡음 큼. ④는 원안(이분법) 기각은 타당하나 대안(연속상관)도 깊이조작을 아예 빼버린 별개 질문(동어반복 위험, r=1.0이어도 인과증거 안 됨)이라 기각.
  - **⚠️✅ 인용 오류 확정 + 해결(원문 PMC10349958 직접 확인 완료, 2026-07-23 밤): "McCoy2024=Ab-Ag Neff-DockQ 무상관"은 틀린 인용.** 실제 논문 = Yin & Pierce 2024(Protein Science, PMID 38073135; McCoy et al. 2024=PMC11337930=TERM/CIM 논문과는 별개, 혼동 금지). 원문(Methods: Neff=CD-Hit 80% identity, 항체+항원 합친 복합체 MSA 전체 기준; Results): **"Neff 분포가 Incorrect vs Medium(p≤0.01) · Incorrect vs High(p≤0.01) 유의하게 다름 — MSA 깊을수록 DockQ 높은 경향"** — 무상관 아니라 **유의한 정상관**.
    - **⭐ 근데 이게 우리 결론과 모순 아님 — 다른 질문(반드시 이 구분 유지):** Yin&Pierce = **복합체 간(between-complex) 관찰적 상관**("자연적으로 MSA 깊은 복합체가 얕은 복합체보다 대체로 잘 맞는다" — confound 가능, 그냥 쉬운 항원일 수도). 우리 프로젝트 = **복합체 내(within-complex) 인과적 개입**("같은 복합체 하나를 잡고 MSA 깊이를 인위조작하면 그 편향이 풀리나" — 안 풀림). 서로 다른 질문이라 둘 다 참 가능(유비: "키 큰 사람이 농구 잘한다"[관찰] vs "이미 큰 선수 키 깎아도 실력 안 변한다"[개입]).
    - **새 프레이밍(보고서·Proposal에 반영할 것): "선행연구는 복합체 간 자연적 깊이-정확도 상관을 보였으나, 특정 복합체 안에서 깊이 조작이 편향을 없애는지는 테스트한 적 없음 — 우리가 그 구체적 개입 질문을 처음 테스트했고 안 됨을 보임."** "문헌과 무상관 정합" 식 서술은 폐기.
  - **유일하게 남는 조각 = germline 정성적 교차확인표**(5개 확정편향 케이스, IGHV3-53/3-66=RBD RBM·IGHV1-69=HA stem 등 문헌 대조) — p-value 없는 랩노트 부록용으로만.
  - **다음 우선순위(3/4 설계안이 독립적으로 지목) = 새 타깃 확장·새 통계 전에 `run_seed_replicate.py`부터 실행.** 지금 새 설계들이 암묵적으로 전제하는 "rescue 현상 자체가 진짜"가 이미 순열검정으로 잡음(p=0.30) 판명 났으므로, 그 전제부터 seed-복제로 재확인 안 하면 위에 뭘 쌓아도 모래성.
  - **⚠️ 미래 세션: "MSA 깊이 변화로 편향 확증" 시도를 다시 밀지 말 것 — 4각도에서 이미 막힘(식별불가능성, 표본으로 안 풀림). 다시 하려면 깊이-반응이 아닌 완전히 다른 축(예: wet-lab 비닝 데이터, Foldseek/TERM 계면 genericness)을 독립 진실축으로 결합해야 함(이미 앞 대화에서 논의됨).**

- **⭐⭐ 2026-07-23(밤) 사용자 명시 지시 = Protenix 결과 나온 뒤 "가설 살리기(save the hypothesis)" 단계 — 반드시 이 지시대로 진행할 것.** 트리거: Protenix DockQ 다각도 채점(eval_dockq_sweep.py, `--watch 300`으로 생성과 병행 채점 중, 528/588≈90% 완료) + 에피토프 분석(eval_epitope_shift.py) 결과가 나오고 **사용자가 "가설 살리기 하자"고 말하면** 착수.
  - **핵심 지시(사용자 원문 취지): "재채점기(재랭커)로 실패 결론 내리지 말 것. 이렇게 끝내면 일주일이 날아가고 그래서는 안 됨. 최대한 우리 가설을 살릴 방향을 어떻게든 찾아내라."** → 블랭킷/모집단 수준이 아니라 **단일 복합체 하나에 가능한 모든 실험을 다 퍼부어서라도** "이 복합체의 이런 특성 때문에 MSA 축소(또는 특정 서열 제거+추가)로 정확도를 높일 수 있다" 또는 "정확히 어떤 경우에 못 찾는지"를 존재증명으로 뽑아내는 방향(KaiB/AF-cluster식 단일사례 심층).
  - **동원할 도구(이미 커밋됨, msa-composition-bias/pipeline): `run_seed_replicate.py`(같은 조건 seed만 바꿔 반복=깊이반응 vs seed잡음 분리) · `prep_ladder_nested.py`(rung_{k+1}⊂rung_k, 전이=빠진 특정서열) · `analyze_loco_aoci.py`(클러스터 leave-out/add-in=인과 클러스터 ID) + MSA 돌연변이/특정서열 제거·추가 실험.** 후보 복합체 = Boltz·Protenix 양쪽에서 가설과 비슷하게 움직이는 것(8wpy_AB=유일 지속전이·9y0a_AB=bistable 반전·8t49_PR=중간깊이최적 등, 단 Protenix 결과 보고 재선정).
  - **⚠️ 균형 주의(과거 교훈과 충돌 관리): 이 "살리기"는 지금까지 확인된 것(블랭킷 null 확정·깊이확증 4각도 식별불가능성·진짜편향5개중4개 무변화)을 부정하는 게 아니라, "모집단은 안 되지만 특수 단일사례는 되는가"를 끝까지 확인하는 것. 억지로 유의성 만들지 말고(best-of-N·garden-of-forking-paths 재발 금지, seed-복제 게이트 필수), 진짜 안 되면 "정확히 이런 특성이면 안 된다"는 negative characterization 자체를 산출물로.** 사용자는 "진짜 실패한다면"의 출구도 열어둠 — 단 재랭커로 곧장 도망가는 건 금지.
  - **지금은 대기: Protenix 채점·에피토프 분석 완료 + 사용자 "가설 살리기" 신호를 기다린다. 신호 오면 이 항목대로 단일복합체 심층 워크플로로 착수.**

- **⭐⭐ 2026-07-24 Protenix DockQ 채점 거의 완료 + 결정적 재프레이밍(반드시 이 상태로): "Boltz는 MSA-둔감 아키텍처라 애초에 이 실험의 적절한 테스트베드가 아니었다. Protenix가 진짜 테스트베드."** eval_dockq_sweep.py로 boltz+protenix 동시 채점(results/dockq_sweep.csv, 타깃마다 flush, ~95% 완료). 파이프라인 실패 0건(생성 boltz 504·protenix 528 조합, best_dockq 빈행 0). Boltz single-seq rung은 MIN_MSA=2로 정상 skip(504<588 설명).
  - **⭐ 핵심 관찰 = Boltz depth-invariant(평평) vs Protenix depth-sensitive(요동).** 예: **8ulr_HL(Env/A CD4bs) Boltz 전rung 0.75~0.80 평평 / Protenix rung0~1=0.01(실패)→rung2~4(neff700~160)=0.66/0.61/0.59 살아남→저깊이 붕괴 = 중간깊이 sweet-spot.** 8k5g_HL Boltz 0.81평평/Protenix rung2만 0.34. 9b7g_QP Protenix full 0.28→감소. 8t4d_OQ Boltz rung3(neff319)만 0.39. 채점 아티팩트 아님(동일 native·병합·DockQ CLI, pose 자체가 다름).
  - **⭐ 문헌 뒷받침(왜 Boltz가 둔감한가): Boltz-2는 Pairformer에서 MSA column-wise attention을 제거(효율화)** = MSA 서열간 공진화 처리 경로를 아예 뺌(AF3/Protenix는 유지) + "공진화 신호 의존 않고 내부 학습 구조 prior를 MSA 없이도 꺼내씀"(biorxiv 2026.01.23.701250). → Boltz는 설계상 MSA 깊이에 덜 반응. **함의: "Boltz에서 안 움직였다"는 우리 가설 반증이 아니라 "MSA 채널 약한 모델로는 이 실험 자체가 불가"였다는 뜻(잘못된 도구로 측정).**
  - **⭐⭐ 재프레이밍(가설 살리기의 유력 출발점, 방어가능+문헌뒷받침): "MSA-깊이 편향 실험은 MSA-의존 모델에서만 유효. Boltz는 아키텍처상 depth-invariant라 부적절한 테스트베드. Protenix가 진짜 테스트베드이고 거기서 depth-response(특히 중간깊이 sweet-spot)가 실제 관찰됨."** → 지금까지 "MSA로 편향 못푼다" 결론이 거의 Boltz기반이었다는 점 재고. **단 seed-복제로 sweet-spot이 진짜인지(best-of-N/seed잡음 아닌지) 먼저 게이트 — 특히 8ulr_HL·8t4d_OQ·9y0a_AB(bistable, Boltz rung1/6/10 반복회복)·8t49_PR 후보.** Chai(MSA-의존)도 삼각확인 가치↑.
  - **다음(사용자 "가설 살리기" 신호 시 착수): Protenix 다각도 분석(추이·family·Spearman·oracle-vs-ipTM·순열·pass@k) + 에피토프 분석 → depth-response 강한 복합체 골라 단일사례 심층(seed_replicate/nested/loco + MSA 특정서열 제거·추가)으로 "이 복합체 이 특성이면 MSA조작으로 정확도 오른다/못찾는다" 존재증명.**

- **⭐⭐ 2026-07-24 (오후) Protenix full-data(44타깃) 완주 + 정직한 이원 결론(반드시 이 상태로 이어갈 것; 위 07-24 오전 항목을 이걸로 갱신).** dockq_sweep 채점 완료(boltz 47·protenix 44 타깃) + epitope_shift 양모델 완주 + 통합 분석기 2종 신규(`analyze_dockq_sweep.py` 8분석·`analyze_epitope_shift.py` 6분석, 순수 stdlib) + `eval_epitope_shift.py` model컬럼 버그 fix. 결과가 **두 갈래**로 갈림:
  - **(A) DockQ = 블랭킷 rescue 여전히 잡음(양모델).** 순열검정 full-data: boltz gain 관측0.110/null0.116 p=0.670 · protenix 0.064/0.070 p=0.773(rescue p 0.77/0.69). **⚠️ 07-24 오전에 낙관했던 "Protenix sweet-spot이 진짜 테스트베드" 중 DockQ 부분은 하향** — 7타깃 파일럿 protenix p=0.030은 소표본 착시, 44타깃서 소멸. Protenix 절대성능도 약함(full-depth acceptable **1/44** vs Boltz 13/44). 단 개별 rescue 강후보 존재(8y6a_CD b 0.04→0.80·8ulr_HL p 0.01→0.66·9y0a_AB 0.06→0.72·9azr_HL).
  - **(B) 에피토프 위치(over-rep) = 가설이 처음으로 방향 맞는 신호(핵심, 지킬 것).** Protenix RBD/A over-rep Δ=−0.457 **Spearman +0.75(7/7 만장일치)** = 깊이↓→인기자리 이탈; **거울상 A/B**: RBD/A recall −0.201(정답=인기라 이탈시 나빠짐) vs RBD/B recall +0.024(정답≠인기라 이탈시 좋아짐). 진짜편향(excess≥0.3) 7타깃 중 **5개 이탈(4개 |Δ|>0.3: 8k3k_D −0.598·8k46_I −0.579·8tx3_FK −0.486·8wpy_AB −0.337)**, Boltz는 5중 1(4개 flat). **Boltz 무반응=아키텍처(MSA column-attn 제거) 탓=07-24 오전 "잘못된 테스트베드" 재확인, 단 방향은 DockQ 아니라 에피토프-위치에서 잡힘.**
  - **모델간 교차확인: 편향 10타깃 중 8개가 한쪽 모델만 = 모델-특이적(재랭커/consensus 여지), 2개(8tx3_FK HA-stem·8txu_HL)만 양쪽 = 근본난제.**
  - **정직 경계: (B)의 위치이탈이 (A)의 DockQ정답으로 안정전환 안 됨(recall Δ 작음, 개별 corr 엇갈림 protenix RBD/B +0.39=반대), family당 n=7.** = 관대하게 "첫 신호"지만 과장 금지·아직 "가설 살리기" 착수 아님.
  - **Notion 기록 완료:** Results **§7-4**(그래프3: RBD거울상·rescue6후보·인기자리이탈 + 표5 + 원자료 CSV 9개) + Results 맨 위 **현재 상태 보드**(H1 뒤 삽입, 정착/살아있는신호/다음한수) + **Progress 2026-07-24** 신규(단순 로그). 그래프 재생성 = `scratchpad/analysis/fig_protenix_final.py`(폴더 CSV 직접 읽음), 데이터 = `~/Downloads/'Notion '`(끝 공백 주의). 서버 figure 2종은 한글깨짐이라 미사용.
  - **⭐ 다음 한 수(우선): RBD 편향 강후보 8k3k_D·8k46_I·8wpy_AB를 `run_seed_replicate.py`(같은 깊이×다른 seed)로 → 인기자리 이탈이 seed잡음 넘어 재현되면 case-study, 아니면 best-of-N. 그다음 Chai smoke-test(3번째 아키텍처; PLM 보정이라 얕은 MSA로도 Ab-Ag 성공 2배[bioRxiv 2025.09.17.676770]=깊이신호 관측 가능). ⚠️ 재랭커로 도망 금지(관대 해석 지시 유지).**

- **⭐⭐ 2026-07-24(밤) 다음 단계 계획 확정 — 재랭커를 "학습 모델"이 아니라 "ipTM 대체 선택기"로 재정의(반드시 이 상태로 이어갈 것).** 사용자가 서사 초안(진단→블랭킷 기각→MSA섭동은 치료제 아닌 샘플러→계면 특화 plug-in 모듈; DeepSCFold식 구조이되 **채점기가 슬리버**)을 가져와 계획 수립. 사용자 명시 지시 = **이전 계획(PLM+접촉그래프+min-edge PAE 재랭커 학습)에 매몰되지 말 것.** 탈출한 함정 2개:
  - **함정① 서열전용 PLM은 pose를 못 가름**: 같은 복합체의 pose들은 **서열이 전부 동일** → ESM-2 임베딩이 pose마다 같은 값 → within-complex 순위에 기여 0. pose를 가르는 건 **기하**(PAE·접촉수·매몰면적·dcc·rank). SaProt/ProstT5(3Di 구조토큰)는 pose마다 바뀌므로 예외로 유지.
  - **함정② n=47 과적합**: 유효표본 = 복합체 47개(그룹)라 처음부터 학습한 GNN/LTR은 무리(LTR은 통상 수백~수천 쿼리). → **재랭커 = "ipTM을 대체하는 선택기"로 재정의**하고 **Phase 1~2를 무학습으로** 설계해 결과가 먼저 나오게. 학습 모델은 맨 위층(선택).
  - **⚠️ 데이터 정정(중요)**: "학습 라벨 이미 있다"는 **틀림**. `dockq_sweep.csv`·`epitope_shift.csv`는 (target,model,rung) 단위 **best-of-5 / 평균 집계**라 **pose 단위 라벨이 없음** → pose 재채점이 선행 필수. 또 "min-edge PAE가 ipTM 압도(Oracle 1.00 vs 0.20)"는 **사수 초안(Epitope_binning, AF3 69 PDB×75모델)** 결과지 이 depth-sweep 결과가 아님 → sweep 출력에 **confidence(ipTM/PAE) 보존됐는지 게이트 확인** 필요.
  - **5단계 계획**: **Phase0** pose 단위 라벨·피처(`lib_pose_features.py` 커밋 `800d78b4`, CPU~3h, 이어달리기; 출력 마지막 줄에 confidence 회수율=게이트) → **Phase1** 무학습 "ipTM vs 기하 피처" 복합체-내부 순위 비교(1~2일, **make-or-break**; best-of-N null 대비) → **Phase2** 사전학습 계면 채점기(DeepRank-GNN-esm·VoroIF-GNN) 적용(2~3일, 이기면 **그게 곧 모듈**) → **Phase3** 깊이-안정성 피처 + matched-N(3~5일) → **Phase4** seed-복제 게이트(8y6a·9y0a) → **Phase5**(선택) 작은 학습 랭커.
  - **사용자 결정 2건**: 라벨 = **DockQ + 에피토프 recall 둘 다** / matched-N 범위 = **반응 서브셋 ~10개**(전체 47 아님).
  - **⭐ 전략 통합(가장 중요)**: 이 재랭커 = **consensus docking Phase 1과 같은 모듈**. 두 개 만들지 말고 하나로, 평가 벤치 2개 — **diverse 10**(co-folder 5종+물리 ZDOCK·HADDOCK = 모델·엔진 다양성 축) / **msa-depth 47**(같은 모델×깊이 섭동 = 입력 섭동 축). 같은 채점기가 양쪽서 ipTM·물리점수를 다 이기면 훨씬 강한 주장이고 두 프로젝트가 한 서사로 묶임("생성은 다양화, 선택은 계면 재랭커").
  - **깊이-안정성 피처 = 이 실험만의 자산**: 깊이 섭동 풀이 있으니 pose마다 "이 에피토프가 섭동 전반에서 얼마나 안정적으로 나왔나"를 공짜로 계산 가능 — 범용 MQA는 못 쓰는 신호. **matched-N에서 MSA≈seed로 나와도 이 피처 때문에 MSA 축의 가치는 남음**(oracle이 아니라 선택기 피처로 재프레임). ⚠️ 단 인기자리가 오히려 섭동에 안정적이므로(Boltz 4/5 flat) "합의=정답"이 아니라 **합의도와 이탈 둘 다** 피처로 줄 것.
  - **상태 = 기록만. 당분간 Consensus Docking(물리 arm·guided)에 집중**(사용자 지시 2026-07-24). 재개 시 Phase 0부터.

- **⭐ 2026-07-26 "가설 살리기(존재증명)에 깊이 감소가 최선 개입인가" 적대검증(11에이전트 5렌즈+5반박+종합) = 아니오/conditional(반드시 이 상태로).** 5렌즈·5반박 만장일치 두 결론: ①깊이 감소는 **시연자 아님, 값싼 스캔 역할로만 정당**(단 유일한 살아있는 신호 +0.75를 낸 도구라 복권=꼴찌 아님) ②`seed_replicate`가 **협상불가 게이트**(2026-07-23 linchpin, nested/LOCO보다 논리적 위). **도구 서열 = seed_replicate(0게이트) → 깊이감소(0.5 스캔, 사전등록 고정rung·적응형금지) → nested(1 조성검증, ⚠️커밋된 build_nested_ladder는 geomspace라 깊이도 같이 줆=순수격리 아님) → 크기매칭 LOCO/AOCI(2 국소화, 신호 게이트통과 후만·⚠️collinearity로 깨끗한 인과귀속은 거짓약속) → template injection(3, ⚠️시연자 아니라 weight-prior 반증자 — 양성이면 'MSA가 편향 나른다' 청구가 죽음).** **⭐결정적 전환 3: (a)채점=DockQ 아니라 epitope-recall(블랭킷 DockQ=잡음 p0.67/0.77, 정직한 종점=예측 에피토프가 진짜 쪽 이동) (b)앵커=에피토프-이탈 실재 사례 8wpy_AB·8k3k_D·8k46_I, 9y0a_AB는 DockQ 교차확인만 (c)Chai 삼각확인=옵션 아니라 필수(Boltz 제외=둔감·잘못된 시험대, Protenix 단독이면 'Protenix 특이 인공물' 반박에 즉사).** 정직한 천장 낮음: novelty=0(전부 선점), 존재증명 장르(KaiB/AF-cluster)=Porter가 noise+암기로 debunk, 조성>깊이는 단일 8wpy 절벽 의존(n=1 인공물후보), n1~3×seed5로 best-of-N 밴드 특성화 못함. unpaired셔플=Ab-Ag서 near-noop 제외·IPW/APC=개입 아니라 선택→Track B(재랭커). **Stage 0 드라이버 `run_track_a_seedrep.sh` 커밋(6b22f665, msa-composition-bias/pipeline; 앵커 항원 a3m을 neff.tsv 기반 얕은~중간 3깊이×5seed 재추첨, CPU). 다음=make_input→Protenix/Chai 예측(GPU, tFold 후)→epitope-recall.** 워크플로 원본 `tasks/wn8tlq0xs.output`. 재랭커=별개 Track B(Phase0 lib_pose_features.py `--models boltz protenix --data /mnt/data/msadepth`, CPU~3h, 마지막줄 confidence 게이트).

- **⚠️ [이 항목은 2026-07-27(밤) 항목이 정정함 — p=3.9e-8과 "best-of-N 배제" 논리는 철회됨. 아래 밤 항목을 먼저 읽을 것]** 2026-07-27 seed-복제 + 예산맞춤 통제.
  - **확정 사례 = 8ulr_HL / Protenix (HIV Env + Fab).** peak rung2 깊이(원서열 1746개, neff 702)에서 조성만 8번 재추첨(개수 고정) → **조성 4/8 성공, 자세 20/40 (DockQ≥0.49), 최고 0.64**. 대조로 **full MSA에 자세 예산을 8배(5→40)** 준 통제 = **0/40, 최고 0.05**. **Fisher 단측 p = 3.9×10⁻⁸.**
  - **best-of-N 구조적 배제(가장 강한 근거):** 분산 분해에서 **조성 안 sd 0.012 vs 조성 사이 sd 0.259 (21.7배)**. 한 조성 안에선 자세 5개가 거의 동일 → full MSA에서 자세를 아무리 뽑아도 같은 답. 실제로 실패 조성(seed1·seed5) 값 0.009~0.012 = **full MSA 값 0.013과 동일**(일부 조성은 full의 실패를 그대로 재현).
  - **정보부족이 아님(핵심 논리):** 성공한 조성은 전부 full a3m에서 뽑았으므로 **full MSA ⊃ 모든 성공 부분집합**. full은 필요한 서열을 다 갖고도 진다 → **일부 서열이 모델을 엉뚱한 자리로 적극적으로 끈다.** "MSA는 깊을수록 좋다"의 깨끗한 반례.
  - **성패는 결합자리에서 갈림:** 성공 자세 recall **0.71** vs 실패 자세 **0.27**. (Boltz 사례들은 더 극단적 — 실패 자세 recall이 정확히 0: 9y0a 33/35, 8y6a 36/37 = 맞히거나 아예 딴 데.)
  - **9azr_HL/protenix 0/40(최고 0.158) = 유효한 음성.** 조성을 다시 뽑아도 안 살아남. 원 관측 0.492는 사다리 11~12칸 중 **최고칸을 고른 winner's curse**(각 칸이 추첨 1회). → **교훈: 단일 추첨의 최고점으로 후보를 고르지 말 것.** (단 "드문 rescue"까지 배제는 못 함 — 1/8 비율이면 8번에 0번 나올 확률 34%.)
  - ⚠️ **8txu·9y0a·8y6a(boltz) 결과는 전부 무효** — 아래 a3m 사건으로 MSA가 아예 안 들어감. **재현 실패가 아니라 미시험.** 2026-07-27 밤 재실행이 첫 시험.

- **⚠️⚠️ 2026-07-27(오후) a3m 질의행 오염 사건 — 범위 확정과 주장 철회(반드시 이 상태로 이어갈 것).**
  - **원인:** `prep_ladder_neff.py`의 `read_raw`가 ColabFold a3m 첫 줄 메타주석(`#<len>\t<card>`, 6자)을 **질의 서열 앞에 붙여** 저장. 첫 `>`를 만날 때 `h`가 None이라 `cur=[]`가 실행되지 않아 주석이 첫 서열에 합쳐짐. `run_seed_replicate.py`도 같은 함수를 import → 사다리·seed a3m **780개 전부** 오염. 49타깃 61사슬 전수 확인(`prep_a3m_check_match.py`), **진짜 서열 불일치는 0건**(전부 머리말만).
  - **⭐ 범위(가장 중요 — 과대평가하지 말 것):** `make_input.py`가 **2026-07-22부터 protenix·chai 분기에 `clean_a3m` 적용**(정규식 `^#\d+\s+\d+\s*`이 붙어버린 경우까지 제거). → **Protenix·Chai는 처음부터 정상 MSA를 받았고 결과 전부 유효**: **8ulr 확정 결론(p=3.9e-8) 그대로 성립**, 9azr 음성 유효, **Protenix 44타깃 sweep·2026-07-24 에피토프 이동(RBD Spearman+0.75) 전부 유효**. **boltz 분기만** 원본 경로를 그대로 넘김(코드 주석 "Boltz는 a3m 쿼리줄 무시"가 **틀린 가정**) → boltz는 불일치를 감지해 **MSA를 통째로 버리고 단일서열로 예측**하며 **경고만 남기고 종료코드는 정상** → 아무도 눈치 못 챔. **boltz sweep 496건 + boltz 후보 3개 무효.**
  - **철회하는 주장 3개:** ①"8txu 재현실패 = winner's curse" → **미시험** ②"조성 안 흔들림 Protenix 0.012 vs Boltz 0.089~0.119"라는 **모델 간 비교** → boltz 수치는 MSA 없는 상태라 비교 불성립(단 Protenix가 MSA 주면 거의 결정적이라는 것 자체는 유효) ③**"Boltz는 깊이에 무반응 = MSA column-attention 없는 아키텍처라 잘못된 시험대"**(2026-07-15·07-24 기록) → **근거 없음. boltz는 깊이 실험을 한 번도 받은 적이 없다.** 2026-07-27 밤 boltz sweep이 첫 검증.
  - **수정 커밋(msa-composition-bias/pipeline):** `prep_ladder_neff.py`(read_raw 근본수정: 주석 줄 처리 + 버퍼 항상 초기화, 손상 파일도 서열 보존) · `make_input.py`(boltz 분기도 `clean_a3m`, 절대경로) · `prep_a3m_fix_query.py`(기존 파일 복구, **기본 dry-run**·멱등) · `prep_a3m_check_match.py`(**정상/머리말오염/서열자체다름 3분류**).
  - **⭐ 교훈(앞으로 규율로):** 모델이 **경고만 내고 조용히 성능을 떨어뜨리는 경우가 있다**(종료코드 정상, 결과도 그럴듯). → **실행 전 `prep_a3m_check_match.py`를 게이트로**, **실행 후 `grep -rl "does not match input sequence"`로 로그 확인**. "오류 없이 돌았다"는 정상 작동의 증거가 아니다.
  - **⚠️ 단어 규율(엄수):** ❌"깊이를 줄이면 rescue"(블랭킷 기각 유지 + 8ulr은 깊이를 **상수로 고정**하고 조성만 바꾼 실험이라 "적을수록 좋다"는 미검증) ❌"에피토프/과대표집 편향"(2026-07-23 인공물 판정) → ✅ **"MSA 조성이 결합자리 선택을 좌우한다"**, 깊이 감소는 조성을 바꿀 수 있게 하는 **전제조건**일 뿐. "인기자리 이탈"은 2026-07-24 RBD over-rep Spearman +0.75(7/7, **탐색적**)를 **별도 문장·별도 근거**로만 — p=3.9e-8을 여기 붙이면 7/23 지적이 그대로 돌아옴.
  - **처방 = 생성 + 선택 한 쌍.** 조성 8개 중 4개만 성공 → 조성 재추첨은 **정답이 든 후보군을 열 뿐** 답을 주지 않음. 재랭커·consensus(Track B)가 반드시 붙어야 함 = 두 갈래가 여기서 만남.
  - **선점 확인(2026-07-27 문헌조사, 반드시 명시):** 기법 novelty = **0**. MSA 부분추출은 표준이고 아키텍처 수술 불필요 — **AF2가 이미 내부에서 MSA 클러스터를 무작위 추출**(seed마다 조성이 바뀜), 손잡이 = ColabFold `--max-msa`·AF2 `max_msa_clusters`. 선행: **Subsampled AF2**(Nat Commun 2024, 구조분포) · **AF-Cluster**(Nature 2024, 접힘전환) · **SPEACH_AF**(PLoS CB 2022, MSA 열 변이) · **AFsample2**(Commun Biol 2025) · **SF-Cluster**(arXiv 2026, MSA 부분추출을 "표현기반 재가중 문제"로 정식화). **전부 한 사슬의 형태(conformation)용.** ⚠️ **Porter PNAS 2024 "AF2가 너무 많이 외운다"** = MSA 부분추출은 **일부 단백질에서만 작동하고 이유 미규명** → **우리 5개 중 1개 확정이 이 관찰과 정합**(실패 사례는 변명이 아니라 선행문헌 부합). **AF3 항체논문(mAbs 2025)**은 "MSA 모듈이 여전히 관건, MSA 입력 질 개선이 매우 가치 있다"고 명시 = 우리가 겨눈 빈틈. **우리 슬리버 = (a)형태가 아니라 에피토프/계면 선택에 적용 (b)예산 맞춘 통제(full이 자기 부분집합에 짐) — 선행은 대개 "부분추출하면 다양해진다"까지.**
  - **📝 확정 결론 문단(보고서·발표용 앵커):** "항체-항원 복합체 예측에서 co-folder가 어느 결합자리를 고를지는 MSA에 어느 서열이 들어 있느냐에 좌우된다. Protenix·8ulr 사례에서 full MSA는 성공한 부분집합들을 모두 포함하고도 틀린 자리를 고르며, 자세 예산을 8배로 늘려도 복구되지 않는 반면(0/40), 같은 서열 수로 조성만 다시 뽑으면 조성 8개 중 4개가 정답에 가까운 자세를 만든다(자세 20/40, DockQ≥0.49, p=3.9×10⁻⁸). 이때 성공한 자세의 결합자리 회복률은 0.71, 실패한 자세는 0.27로, 성패가 결합자리 선택에서 갈린다. 별도로, 여러 복합체를 대상으로 한 위치 분석에서는 MSA를 사용하는 Protenix가 깊이를 줄일수록 자주 관측되는 결합자리에서 벗어나는 경향이 관찰되었다(RBD 계열 7개 일치, 탐색적). 종합하면 MSA 조성 재추첨은 항체-항원 결합자리 예측의 실패를 교정할 수 있는 통로를 열어주지만, 정답을 보장하는 것이 아니라 후보군이 열릴 뿐이므로 선택 단계가 함께 필요하다."
  - **다음 실험 = 빈도 측정(사전등록형, winner's curse 회피).** 규칙을 **먼저** 확정: `Protenix full-MSA best-of-5 < 0.23(실패) AND 정답 도달가능(다른 모델이나 oracle ≥0.49)` → 해당 타깃 **6~8개** × **한 중간 깊이** × **조성 8회**(처음부터 복제, 최고칸 선별 금지) → 편향 없는 **M/N**("N개 중 M개에서 조성 재추첨이 rescue"). 지금 제일 없는 숫자가 이것. 층위 3(깊이별 성공확률 곡선)·층위 4(어느 서열이 해로운가=여기부터 "방법")은 그 다음.
  - **커밋(msa-composition-bias/pipeline):** `run_seedrep_cand.sh`(⚠️**다중사슬 버그 fix `456e78e0`** — 사슬별 서열수가 다르면(8txu A=413/B=579) 첫 사슬 깊이 폴더명으로 전부 skip되던 문제; 이제 사슬별 깊이 기억+seed축 순회) · `eval_dump_seedrep.py`(자세 단위 전체 덤프 = full/사다리/seed복제 3층 나란히; `--only` 쓰면 원자료 별도 파일로 분리) · `run_fullmsa_control.sh`(full MSA × N자세 예산맞춤 통제) · `eval_fullmsa_control.py`(자세단위 맞대결 + Fisher 단측검정, stdlib 구현·scipy 대조 검증) · `plot_export_data.py`(그림용 데이터 추출).

- **⭐⭐⭐⭐ 2026-07-27(밤) 올바른 단위로 재측정 → 가설 확인. 이 항목이 최신이고 앞의 두 항목을 정정한다.**
  - **⚠️ 철회 1 — p=3.9×10⁻⁸.** 자세 40개를 독립 표본으로 셌으나, **한 실행 안 자세 5개는 서로 상관**되어 유효표본은 **실행 수**다. 내가 "조성 안 sd 0.012"를 *결정성의 증거*로 읽은 게 반대였다(=자세들이 독립이 아니라는 뜻). 
  - **⚠️ 철회 2 — 재현성.** 바이트 동일 입력(input.json·ag_A_clean.a3m·항체 a3m 4종 md5 일치)으로 재실행 시 8ulr seed0 **0.588 → 0.011**. 확산 모델이라 실행마다 답이 갈릴 수 있다. 대응은 포기가 아니라 **여러 번 돌려 성공률을 재는 것**.
  - **⭐ 올바른 설계 = 조성 × 반복** (`make_composition_reps.sh`): 조성 8개 × 4회 + full MSA 10회 = 42회. **성공 = DockQ ≥ 0.49, 단위 = 실행 1회(자세 5개 중 최고).**
  - **⭐⭐ 결과(8ulr_HL / Protenix, 깊이 rung2=1746서열):** 조성별 성공률 **seed3 4/4 · seed4 4/4 · seed7 3/4 · seed1·5·6 각 1/4 · seed0·seed2 각 0/4 · full MSA 1/10.**
    **→ 조성 간 이질성 정확검정 p = 0.0025 (사후선택 아님, 옴니버스).** 이게 **핵심 결과**이자 원래 가설의 확인이다.
    - 얕은 전체 14/32(43.8%) vs full 1/10(10%): **p = 0.054(경계)** — "깊이를 줄이면 좋아진다"는 **약함**. 게다가 줄인 조성 중에도 0/4가 둘(seed0·seed2)이라 **full보다 나쁜 조성도 있다.**
    - ⚠️ seed3+4만 뽑으면 8/8 vs 1/10 → p=0.0002지만 **사후선택이라 주장 금지.**
  - **⭐ 기제(`analyze_epitope_cluster.py`, 외부 PDB 인기도 안 씀):** 성공 실행끼리 예측 접촉면 겹침 **0.711**(우연 0.045의 16배), 실패끼리 **0.386**(우연 0.111의 3.5배 = **공통 선호 영역 존재**), 성공 vs 실패 **0.231**(서로 다른 자리). 합의 접촉면 = 성공 **34잔기·진짜 30개 중 27개(90%) 포함·정밀도 0.79** / 실패 **72잔기·정밀도 0.28**. 접촉면 평균 크기 **성공 38 vs 실패 88잔기**(항원 440잔기의 20% = 물리적으로 말 안 되는 넓이) → **접촉면 크기 자체가 정답 없이 쓸 수 있는 재랭커 피처**(Arm A와 연결).
  - **⚠️ 용어 정정(중요):** `epitope_shift.frac` = `|pred ∩ region| / |pred|` = **정밀도(precision)**이지 recall이 아니다. 지금까지 "결합자리 회복률"로 부른 값은 전부 **"예측 접촉 중 진짜인 비율"**이다. 보고서에 "회복률"로 쓰면 틀림.
  - **⭐ 후보 재선별(`analyze_screen_candidates.py`) — 기준 변경:** full 실패 + 얕은 깊이 **연속 2칸 이상** 성공(한 칸만 튄 건 실행 운). **Protenix 44타깃 중 ★강후보 = 8ulr 하나(1/44)**, 9azr는 1칸뿐이고 실제로 재현 실패(0/40) → 기준이 작동함을 입증.
  - **⭐ boltz는 반대 — 가설 기각(제대로 된 첫 측정):** full MSA가 더 좋거나 같음. **8y6a 15/40 vs 얕은 1/40 · 9y0a 12/40 vs 10/40 · 8txu 0/40 vs 1/40.** 즉 MSA는 도움이 되고 줄이면 손해. (8y6a의 "얕은 깊이"는 서열 3개라 사실상 MSA 제거.)
  - **✅ 쓸 수 있는 문장:** *"8ulr에서 Protenix의 성공 여부는 MSA에 어느 서열이 들어 있느냐에 좌우된다. 서열 수를 1746개로 고정한 채 조성만 8가지로 바꿔 각 4회 예측했더니 성공률이 0/4~4/4로 갈렸다(이질성 p=0.0025). 실패 시 항체가 넓게 퍼져 붙고(88잔기, 정밀도 0.28) 그 잘못된 영역은 실행 간 공유되며(겹침 0.386=우연의 3.5배), 성공 시 좁고 정확한 접촉면(38잔기, 정밀도 0.79, 진짜의 90%)으로 수렴한다. Protenix 44개 중 이런 사례는 1개였다."*
  - **❌ 쓰면 안 되는 것:** ①"깊이를 줄이면 정확도가 오른다"(p=0.054·역전 조성 존재) ②"일반적으로"(1/44) ③**"PDB에서 인기 있는 자리에서 벗어난다"** — 8ulr의 진짜 에피토프 CD4bs가 **곧 인기 자리**라 데이터가 반대다(성공할수록 인기자리와 더 겹침: over-rep 0.352→0.703, `chains.json` `"AB":"A"`=on-site). 2026-07-23 과대표집 인공물 판정과도 겹침.
  - **표준 분석 절차(모든 타깃 동일):** `bash run_analyze_target.sh <타깃>` = dump_seedrep_full → score_compreps(성공률+이질성) → epitope_cluster(기제). ⚠️ **이질성 검정은 조성당 반복 2회 이상일 때만** 나온다(예비검정은 1회라 불가 → 유망하면 반복 늘려 재실행).
  - **진행 중(2026-07-27 밤):** 예비검정 4개 = 9azr_HL(RUNG4) · 8k5g_HL(RUNG2) · 8q7s_C(RUNG3) · 8ume_HL(RUNG1), 각 얕은 5조성×1회 + full 5회. 그다음 boltz 47타깃 sweep(`run_sweep.sh boltz`, 11h 규모, `$DATA/boltz` mv 선행 필수) → 내일 `analyze_screen_candidates.py --model boltz`로 boltz 빈도.

- **⭐⭐ 2026-07-27(밤, 후속) 실험계획 [⚠️ 아래 '밤, 최종' 항목이 이것을 갱신함 — 본 검정 27→20, 거울상 기각]**
  - **상위 서사(사용자 확정):** *"MSA 조성(깊이 포함)을 바꾸면 통계적으로 유의하게 예측 결합자리가 이동한다. DockQ 상승까지 바로 이어지는 경우는 드물지만, 기존 연구가 MSA 부분추출을 **한 사슬의 접힘/형태**에만 쓴 것과 달리, **에피토프 비닝을 위한 유의미한 후보 생성**이 가능함을 확인했다."* → **주 지표 = 결합자리(위치), DockQ는 보고만.**
  - **⭐ 새 데이터 원천 발견 = `runs_rbd`** (`~/projects/epitope-guided-docking/consensus_docking/runs_rbd`, **다른 레포**=Consensus Docking 쪽). 노션 인수인계서 1.3절 "결합부위 밖 항체 10종". **10종 전부가 인기 부위를 피해 붙는 항체**(44개에서 겨우 8개 찾던 조건 B군을 10/10이 만족). 8sdh는 노션 표에 없고 분류 미기록이라 **사용자 지시로 제외**.
    - 구성: **B군 9** = 8sdf_HL·8siq_HL·8sis_HL·8sit_HL·8xsi_HL·9ml8_HL·9ml9_HL·9sbb_HL·9zdu_HL (숨은면 3·코어 4·코어와그밖 1·그밖 1) / **A군 1** = **8p5m_GL**(인기 부위 가장자리 = **반대방향 대조군, 필수 포함**).
    - **깊이 사다리는 없음**(out_protenix 등은 원래 MSA 1조건). 하지만 **항원 a3m이 이미 있음**(`msa_<pdb>/A.a3m`) → 가장 비싼 단계 생략. 없는 것은 native.cif뿐(RCSB에서 받음).
    - ⭐ **10종 전부 같은 a3m = 143개 서열, 항원 195잔기 내외(같은 RBD).** 8ulr(1746서열·440잔기)보다 훨씬 빠름.
    - ⭐⭐ **설계상 강점:** 항원·MSA가 완전히 동일하고 **항체만 다름** → 항원 차이·MSA 크기 같은 교란이 원천 차단. 같은 MSA에서 조성만 바꿔 **각 항체의 서로 다른 진짜 자리로** 예측이 옮겨가면 항원 요인으로 설명 불가. ⚠️ **반대급부: 항원은 1개다. "항원 10개에서 확인"이라고 쓰면 안 됨** — 항원 다양성(Env·HA·RBD·C)은 기존 44가 담당.
  - **⭐ 최종 데이터셋(확정):**
    | 검정 | 대상 | 개수 | GPU |
    |---|---|---|---|
    | ① 위치이동 전수 | 기존 44 + RBD 10 | **54** | 0(44는 예측물 있음) |
    | ② 본 검정(조성8×반복4 + full 10 = 42회/타깃) | RBD 10 + 44에서 16 + 8ulr | **27** | 실측 후 조정 |
    - 본 검정 27 내역: **B군 17**(RBD 9 + 44의 8) · **A군 10**(8p5m + 44의 8 + 8ulr).
  - **⭐ 판정 기준(사전 확정 — 결과 보고 고르지 않기 위해 미리 적음):**
    - **성공 = 결합자리 회복률 ≥ 0.4**(DockQ 아님). `run_analyze_target.sh`가 이미 세 지표(DockQ·회복률 th0.4·인기자리겹침 `--lower-better` th0.3)를 찍음.
    - **검정① 조성 효과:** 타깃별 이질성 정확검정 → **27개 중 유의한 타깃 수**. 귀무 기대 1.35개. **P(X≥4)=0.044(경계) · P(X≥5)=0.010 · P(X≥6)=0.0019 · P(X≥7)=0.0003.** → **기준선 = 5개.**
    - **검정② 거울상:** 타깃별로 "조성 축소 시 인기자리 겹침이 줄었나" 부호 → A군 vs B군 부호검정(`analyze_shift_direction.py`). 최소 p: B군 n=17 → 1e-5, A군 n=10 → 0.00098.
    - **결론 두 갈래(미리 정함):** (a)유의≥5 **AND** 거울상 성립 → **"MSA 조성이 인기 결합자리 편향을 만들고 조성 변경으로 풀린다"**(강한 주장) (b)유의≥5, 거울상 실패 → **"MSA 조성이 결합자리 선택을 좌우하며, 조성 다양화가 재현성 있는 에피토프 후보 생성기가 된다"**(약하지만 안전, 선행연구와 여전히 비중복) (c)유의<4 → 8ulr 사례연구로 축소, 헤드라인은 음성.
  - **⚠️⚠️ 거울상은 이미 반례 1개를 안고 출발한다 — 잊지 말 것.** **8ulr은 A군**(진짜 자리 CD4bs = 인기 자리)인데 조성 변경으로 **좋아졌고** 인기자리 겹침도 **0.35→0.70으로 늘었다**. 즉 8ulr의 기제는 "편향 해소"가 아니라 **"넓게 퍼진 접촉이 좁게 수렴"**. 거울상 서사를 데이터 확인 없이 밀면 자체 데이터에 반박당함.
  - **⭐ 다른 갈래 처리(확정):** ①**예비검정 4개**(9azr_HL·8k5g_HL·8q7s_C·8ume_HL) = 끝까지 돌리고 채점하되 **용도 전환 → DockQ 기준으로 뽑은 것들이라 "DockQ로 뽑으면 재현 안 된다"는 음성 대조군**. ②**`candidates.csv` 6개**(8k3k_D·8tx3_FK·8k46_I·8xsj_HL·8y6a_CD·8wpy_AB) = 지금 안 돌리고 **본 검정 후보 pool로 대기**(검정① 결과로 44에서 16개 확정할 때 함께 고름). ③**boltz 47종 sweep = 3순위로 미룸** — boltz는 제대로 돌린 3개에서 full MSA가 더 좋다고 **반대로** 나왔으므로 테스트베드로 약함. 논문에는 "모델에 따라 다르다" 한 문단으로 충분하고 그건 지금 3개로 됨. GPU 순서 = **예비검정(진행중) → RBD 사다리 → 본 검정 → (여유시) boltz**.
  - **⭐ 새 커밋 `prep_import_rbd.py`(`3d1357ae`, msa-composition-bias/pipeline):** runs_rbd → 깊이실험 형식 어댑터. chains.json에서 항원/항체 판별(표기 없으면 길이로 추정) · 기존 a3m을 `$DATA/ladders/<타깃>/<사슬>/rung0.a3m`으로 복사(**머리말 오염 자동 제거**) · native.cif를 RCSB에서 받음 · sweep_targets.csv에 행 추가. **기본 dry-run, 기존 파일 절대 안 덮어씀.** 자가검증 4항목 통과.
    - ✅ **열 매핑·8sdh 제외 수정 완료(`4706e28625`).** 실제 열 = `target,pdb,group,ab,dirtype,ag_chains,label` — **`ab`가 A/B 표식이고 항체 사슬이 아님**(8q7s_O=A, 8y6a_CD=B), `label`=붙는자리 분류(기존은 `class1/2(RBM)`·`class3` 같은 Barnes 분류; 우리 것은 `offhot:숨은면` 형태로 구분), `dirtype`=`targets` 자동 승계. 기존 `targets/*/chains.json`을 본보기로 찍고 최상위 키 누락을 경고(사슬 ID는 비교 제외). 가짜 환경 종단 자가검증 통과(dry-run/apply, 143서열 머리말 제거, 8sdh 건너뜀, A/B 표식·항원사슬 정확).
    - ⚠️⚠️ **2026-07-27 밤 형식 사고와 재작성(`634717aea8`) — 반드시 이 상태로.** 첫 판이 runs_rbd의 chains.json을 **그대로 복사**했으나 이 저장소 형식은 완전히 다름: `prep_targets.py`가 **사슬을 항원=A·중쇄=B·경쇄=C로 개명**하고 `{pdb_id,target,antigen_grp,AB,label,antigen,antibody,chains:[{id,role,seq,src,crop}],src_chains}` 구조를 만들며 `native.cif`는 `structures/<pdb>.cif` **심링크**. 게다가 **사다리 폴더는 개명 ID로 만들어야 하는데**(sweep_targets.csv `ag_chains`가 전부 `A`) 원본 ID(I·Z·B·R)로 만들어 4개 타깃이 어긋남. → **chains.json을 손으로 흉내내지 말고 `prep_targets.py` 정규 경로에 태울 것.**
    - ✅ **RBD 크롭 걱정 없음:** `prep_targets.py:66`은 RBD 그룹이라도 **400잔기 초과 사슬만** 319-541로 자름. 우리 항원은 194~197잔기라 `crop=null` → runs_rbd a3m(195잔기 기준) 질의서열이 그대로 일치.
    - **재작성 3단계(각각 기본 dry-run):** `--stage undo`(첫 시도 산출물 정리, 우리가 만든 이름만) → `--stage csv`(prep용 CSV `pdb,Hchain,Lchain,antigen_chain,antigen,AB,label` + `structures/<pdb>.cif`; **중/경쇄는 J영역 모티프 중쇄 WGxG·경쇄 FGxG로 판정**, 실패 시 길이) → **사용자가 `prep_targets.py` 실행** → `--stage msa`(개명된 항원 사슬에 a3m 연결 + sweep_targets.csv 등록). 가짜 환경 3단계 종단 자가검증 통과.
    - **▶ 바로 실행할 명령:** `cd ~/projects/msa-composition-bias && git pull && cd pipeline` → `python prep_import_rbd.py --stage undo`(확인 후 `--apply`) → `python prep_import_rbd.py --stage csv --apply` → `python prep_targets.py --csv rbd_offhot.csv --struct structures --outdir targets` → `python prep_import_rbd.py --stage msa --apply` → **`python prep_a3m_check_match.py`(게이트, 필수)** → `python prep_ladder_neff.py` → `bash run_sweep.sh`. group 값 = **RBD**(기존 목록 C·Env·HA·RBD).
  - **⭐ 실행 순서:** ①import_rbd --apply → neff_ladder(143서열 → 8칸) ②사다리 예측 10×8칸×1회=80회(**~2시간**, 깊이 고르기 + 위치이동 관찰) ③`lib_pose_features.py` → `analyze_shift_direction.py`로 **54 전수 거울상**(GPU 0) ④본 검정 RBD 10 전수 ⑤본 검정 44 중 16. **②의 첫 타깃 실소요를 재서 ④⑤ 개수 재조정**(빠르면 44에서 16→20으로 확대).
  - **선점 위치(변함없음):** 주장은 "방법이 새롭다"가 아니라 **"적용 대상이 새롭다(형태가 아니라 결합자리 선택) + 예산 맞춘 통제"**.

- **⭐⭐⭐ 2026-07-27(밤, 2차) 거울상 검정 결과 = 기각 → 서사 확정(결론 갈래 b).** `analyze_shift_direction.py --model protenix`, 44타깃(A 20·B 15, pose_features 2640행 완비).
  - **거울상 없음(확정):** 진짜자리 겹침 변화 평균 A **−0.058** vs B **−0.018** — 둘 다 음수, 부호검정 A p=0.94·B p=0.30, 흔한자리 변화도 A p=0.41·B p=0.27. **"MSA 깊이를 줄이면 B군이 편향에서 풀린다"는 지지되지 않음.** (깊이 축 블랭킷 기각과 정합 — shift_direction은 rung>0 전체 평균이라 조성-특이 효과는 원래 희석된다.)
  - ⚠️ **결합(coupling)은 거울상처럼 보이나 증거로 쓰면 안 됨:** Spearman(진짜자리변화, 흔한자리변화) = A **+0.52** / B **−0.60**. 그러나 A군은 정의상 진짜자리≈흔한자리라 **부호가 기하학적으로 거의 자동**이다. 같은 이유로 "진짜↑ & 흔한↓ 동시" B 6/15 vs A 2/20(Fisher 단측 p=0.046)도 **A군에서는 구조적으로 일어나기 어려운 조합**이라 독립 증거가 아니다. 발표·보고서에 인과 증거로 쓰지 말 것.
  - ⭐ **독립적으로 진짜인 것 2가지(= 새 서사의 근거):** ①**실패가 재현된다** — 실패 실행끼리 예측 접촉면 자카드가 우연의 **3.5~7.9배**(8ulr 3.5 · 8k5g 3.9 · 9azr 6.3 · 8q7s 7.9 · 8ume 5.9), 즉 모델이 반복해서 같은 엉뚱한 자리로 간다. ②**그 자리는 PDB 인기 자리가 아니다** — 예비검정 4개 전부 "실패 예측 ↔ 흔한 자리" 겹침이 진짜 자리보다 **낮다**. → **"PDB 과대표집/인기자리 편향" 서사는 우리 데이터가 지지하지 않음**(2026-07-23 인공물 판정과 정합). 대신 **"모델 내부의 재현되는 잘못된 선호 자리"**.
  - **✅ 사용자 승인 서사(2026-07-27 밤, 표현 확정):** *"MSA 조성이 항체-항원 결합자리 선택을 좌우한다. 조성을 다시 뽑으면 **믿을 만한(재현되는) 다른 결합자리 후보**를 만들 수 있으므로, 에피토프 비닝의 후보 생성 단계에 쓸 수 있다."* — 헤드라인은 **위치(결합자리)**, DockQ는 보고만. ⚠️**"잘못된 자리로 수렴한다"를 강하게 단정하지 말 것**(사용자 명시). 실패의 재현성(우연의 3.5~7.9배)·인기자리 비해당은 **보조 관찰로만** 서술하고 주장의 축으로 올리지 않는다.
  - ⚠️ **A/B 표식이 실제 겹침과 안 맞음:** 예비검정 4개 전부 표식 `A`인데 계산된 진짜↔흔한 겹침은 0.068·0.111·0.189·0.324로 낮다. 표식은 다른 기준(항원 과대표집 등)으로 붙은 듯 → **거울상류 분석을 표식으로 나누지 말 것.**
  - **예비검정 4개 판정:** 9azr_HL·8k5g_HL **음성 확정**(음성 대조군으로 사용) · 8q7s_C **탈락**(DockQ 점수순위 p=0.016이나 절대값 0.089→0.132로 둘 다 실패권이고 **진짜자리 겹침은 0.84→0.68로 악화**) · **8ume_HL 승격**(세 지표 동시 방향 일치: DockQ 0.023→0.145 p=0.075 · 진짜자리 0.33→0.63 p=0.052 · 흔한자리 0.27→0.00). 조성당 1회라 이질성 검정 불가 = 검정력 없음.
  - ⭐ **본 검정 B군 후보 6개가 두 경로에서 일치:** `analyze_screen_epitope_shift.py`가 뽑은 목록과 `shift_direction`의 "진짜↑ & 흔한↓ 동시" 목록이 **정확히 같음** = **8xsj_HL · 8k3k_D · 8y6a_CD · 8tx3_FK · 8wpy_AB · 8k46_I**(= `candidates.csv` 그대로). 여기에 8ume_HL 추가.
  - **RBD 10종 실측 정정:** a3m이 143개가 아니라 **21,737~27,946개**(출처 = `out_colabfold/<pdb>_env/uniref.a3m` 또는 boltz `..._unpaired_tmp_env/uniref.a3m`, 질의행 = RBD 서열 확인). 항원 길이도 194~197로 제각각이라 **"10종이 항원·MSA가 완전히 동일해 교란이 원천 차단"이라는 앞선 서술은 과장** — 같은 단백질(RBD) 계열이라는 정도로만. 깊이는 오히려 44종과 비슷해져 사다리가 제대로 퍼짐. 예측 소요도 재추정 필요(하룻밤).
  - **커밋:** `prep_import_rbd.py`(3단계 재작성 `634717aea8`) · `prep_ladders.sh`(`b6105be5`, rung1 있으면 건너뛰어 기존 61사슬 조성 보호, 멱등).

- **⭐⭐⭐⭐⭐ 2026-07-27(밤, 최종) 연구 전체 정리 — 세션이 바뀌면 이 항목 하나만 읽고 이어가면 된다. 노션 인수인계서 6장과 같은 내용.**
  - **📖 핵심 이야기(사용자 확정 문장):** *"MSA에 어느 서열이 들어 있느냐가 모델이 항체를 항원의 **어느 자리**에 붙일지를 좌우한다. 그래서 서열 목록을 다르게 뽑으면 **믿을 만한(재현되는) 서로 다른 결합자리 후보**를 만들 수 있고, 이것을 **에피토프 비닝의 후보 생성 단계**에 쓸 수 있다."* — 두 단어가 핵심: **재현성**(같은 조성이면 같은 자리 = 무작위 흔들림과 구별) · **후보**(정답 보장 아님, 고를 거리를 만듦 → 선택 단계는 재랭커·consensus 몫). **헤드라인 = 결합자리 위치. DockQ는 보고만**(같이 오르는 경우는 드묾).
  - **✅ 확인된 것:** ①조성이 자리를 바꾼다(8ulr, 서열수 1746 고정·조성 8가지×4회 → 성공률 0/4~4/4, **이질성 정확검정 p=0.0025**) ②개수 축소 자체는 답이 아님(줄인 쪽 전체 vs full은 p=0.054 경계 + 줄인 조성 중에도 full보다 나쁜 것 존재). **보조 관찰(강하게 단정 금지)**: ③실패 실행끼리 예측 자리가 우연의 **3.5~7.9배**로 겹침(8ulr 3.5·8k5g 3.9·9azr 6.3·8q7s 7.9·8ume 5.9) ④그 자리는 **PDB 인기 자리가 아님**(예비검정 4개 전부 실패 예측↔흔한자리 겹침이 진짜 자리보다 낮음).
  - **❌ 검정하고 기각한 것(보고서에 숨기지 말 것):** 설계 가설 = "B군(진짜 자리가 흔한 부위 밖)은 조성 축소로 편향이 풀려 좋아지고 A군은 나빠진다." **44개 전수 `analyze_shift_direction.py` 결과 = 기각.** 진짜자리 겹침 변화 A **−0.058**(부호검정 p=0.94) vs B **−0.018**(p=0.30), 흔한자리 변화도 A p=0.41·B p=0.27. **거울상 없음.**
    ⚠️ **솔깃하지만 근거로 쓰면 안 되는 것:** 두 변화량 Spearman이 A **+0.52** / B **−0.60**으로 정확히 반대이고 "진짜↑&흔한↓ 동시"가 B 6/15 vs A 2/20(Fisher p=0.046)이지만, **A군은 정의상 진짜자리≈흔한자리라 부호가 기하학적으로 자동**이다. 인과 증거 아님.
    ⚠️ **A/B 표식이 실제 겹침과 불일치**(예비검정 4개 모두 표식 A인데 계산된 겹침 0.068·0.111·0.189·0.324). **앞으로 A/B로 나누어 분석하지 말 것.** 본 검정도 A·B 균형을 맞출 이유 없음 → **효과 있는 것 위주로 선정.**
  - **⭐ 본 검정 = 20 복합체(확정, 27에서 축소).** 타깃당 조성 8×반복 4 + 원래 10 = 42회. **성공 기준 = 결합자리 회복 ≥0.4(DockQ 아님).**
    | 세트 | n | 역할 |
    |---|---|---|
    | RBD 인기부위 밖 항체 10종 | 10 | ⭐**거르지 않고 전수** → 편향 없는 빈도 M/10 = 헤드라인 숫자 |
    | 44에서 자리이동 큰 것 | 7 | 기제 설명(선별했으므로 빈도 불가) |
    | 음성 대조군 | 2 | 9azr_HL·8k5g_HL(완료) |
    | 완료 | 1 | 8ulr_HL |
    **44측 7개 = 8xsj_HL · 8k3k_D · 8y6a_CD · 8tx3_FK · 8wpy_AB · 8k46_I · 8ume_HL** (앞 6개는 `screen_epitope_shift`와 `shift_direction`의 "진짜↑&흔한↓" 목록이 **정확히 일치** = `candidates.csv` 그대로; 8ume_HL은 예비검정에서 세 지표 동시 방향 일치로 승격).
    ⚠️ **RBD 10종을 걸러 돌리지 말 것** — 걸러내면 빈도를 못 잰다. 지금 제일 없는 숫자가 편향 없는 M/N.
  - **⭐ 판정 기준(사전 확정, 결과 보고 바꾸지 말 것):** ①**복합체별 효과** = 이질성 정확검정 유의 타깃 수, **20개 중 4개 이상이면 우연확률 0.016 · 5개 0.003 · 3개는 0.076(불충분)** → **기준선 4개** ②**재현성** = (조성 내 결합자리 겹침)/(조성 간), 조성 딱지 뒤섞기 검정 **p<0.05** ③**후보 수** = 조성별 합의 자리를 겹침 0.5로 묶은 개수, **하나라도 진짜 자리 절반 이상 덮으면 성공**. **결론 갈래**: (a)효과+재현성 → 위 핵심 이야기 그대로 (b)효과만, 재현성 낮음 → "자리가 흔들린다"까지만, 후보 생성 주장 뺌 (c)유의 3개 이하 → 8ulr 사례연구로 축소, 헤드라인 음성.
  - **⭐ 새 커밋:** **`analyze_site_reproducibility.py`(`490c74db`) = 주장의 핵심 통계** — 조성 내/간 결합자리 겹침 비율 + 뒤섞기 검정 + 서로 구별되는 자리 후보 개수와 각 후보의 진짜자리 덮음/정밀도. 자가검증(대립 p=0.0003 · 귀무 p=0.53) 통과. **`run_analyze_target.sh`에 ④단계로 편입**(`12683daa`) → 타깃마다 자동. `comp_x_reps`를 거친 타깃에만 성립(44 사다리는 칸당 1회라 조성 내 편차가 없음).
  - **선점(반드시 명시):** 기법 novelty 0. AF-Cluster(Nature 2024)·Subsampled AF2(Nat Commun 2024)·SPEACH_AF(PLoS CB 2022)·AFsample2(Commun Biol 2025)·SF-Cluster(arXiv 2026) **전부 한 사슬의 형태(conformation)용**. **결합자리 선택에 적용한 사례 없음 = 우리 슬리버.** 주장은 "방법이 새롭다"가 아니라 **"적용 대상이 새롭다 + 예산 맞춘 통제"**. Porter PNAS 2024(MSA 부분추출은 일부에서만 작동)와 우리 낮은 빈도가 정합.
  - **📍 실행 위치:** 노션 인수인계서(page 3a7ac1a3-a28a-8151-a962-e10474fa91ea) **5장=8ulr 사례 · 6장=현재 결론과 계획**. 코드 = `Feellived/msa-composition-bias/pipeline`.

- **⭐⭐⭐⭐⭐ 2026-07-28 설계·통계 확정 + protenix 폴더 사고 수습 — 여기가 최신. 위 '밤, 최종'과 함께 읽으면 현재 상태 전부.**
  - **⭐ 통계 구조 확정(사용자 지적으로 교정 — 내가 처음에 틀림):** 모든 복합체가 **똑같은 실험**(조성 6 × 반복 4 + 원래 8 = 32회)을 받고 결과 형태도 같다. 세트별로 다른 것은 **어떻게 명단에 들어왔나 하나뿐**이고, 그건 **비율에만** 영향을 준다.
    - **주 결과 = 검정한 전부를 모아 합친 p(Fisher).** ⭐**선별해도 타당하다**: 이질성 정확검정은 그 실험의 **총 성공 횟수를 조건으로 삼아** 그 깊이의 성공률을 상쇄하고, 본 검정은 조성·반복을 **새로 뽑으므로** 선별에 쓴 사다리 자료가 새 자료에 들어오지 않는다 → 귀무에서 각 타깃의 p가 균등분포 → Fisher 결합이 유효. **N을 늘리는 게 실제로 이득**(각 p=0.20이어도 8개면 0.058, **16개면 0.016**). ⚠️단 **신호 없는 타깃을 넣으면 희석**(강한 4개+무반응 6개=0.006 vs +무반응 12개=0.035) → 아무거나 많이가 아니라 **신호 있을 만한 것을 골라 넣는 게 맞다**.
    - **하한** = "44개(또는 54개) 중 **최소 M개**" — 명단이 먼저 정해진 세트를 분모로 쓰고 **안 돌린 것을 실패로 계상**. 밑으로만 틀리는 안전한 숫자.
    - **비율** = 코로나 10종 안에서 M/10(선별 없이 명단 확정된 유일한 세트). **곁가지 — 안 써도 됨.** 사용자가 이 표현을 반복해서 헷갈려했으니 앞으로 간단히만 언급할 것.
    - ⚠️ **"마음대로 뽑았다"고 쓰지 말 것**(사용자 제안이었으나 사실과 다르고 스크립트가 저장소에 남아 있음). **"골랐다"고 그대로 쓰고 왜 검정이 무해한지 설명**하는 쪽이 더 강하다.
    - ⚠️ 내가 붙였던 "사례 증명/본 검정 후보/편향 없는 빈도" 3분류는 **과한 이름 붙이기**였음 — 실제로는 **하나의 실험 + 두 종류의 명단**.
  - **⭐ 설계 축소(검정력 계산 근거):** 이상적 형태에서 8×4=p<1e-5 · **6×4=1e-5** · 6×3=4e-4 · 5×3=2e-3. 8ulr 실제 형태를 축소하면 8×4 p=0.0025 · **6×4 p=0.0011** · 6×3 p=0.034(경계). → **반복은 4회 유지, 조성만 8→6.** **6조성 × 4반복 + 원래 8회 = 32회**(24% 절감, 검정력 유지). 8ulr은 이미 8×4+10으로 돌아가 있어 그대로 사용.
  - **⭐ 깊이는 결과를 보고 골라도 된다(내가 처음에 과하게 보수적이었음).** 근거는 위 조건부 검정 논리. **금지되는 것은 ①타깃을 결과로 고르는 것(비율이 부풀려짐) ②"그 깊이에서 성공률 X%"·"얕은 쪽이 원래보다 낫다"는 절대·비교 주장.** 오히려 성공률이 0도 1도 아닌 **중간 지대**를 고르는 게 검정력이 최대.
  - **⭐ 판독 규칙(2026-07-28 결과 보기 전 확정, `prep_pick_depth.py` `a03c7543`):** 칸마다 구조 5개 중 **결합자리 덮음 ≥0.4** 개수를 세어 — ①60개 구조에서 **하나도 없음 → 본 검정 안 함, 분모에 실패로 계상**(GPU 절약 + 안전한 방향) ②**모든 칸 5/5 → 안 함, '구제 불필요'로 분리** ③그 밖 → **본 검정 대상, 깊이 = 1~4개인 칸 중 서열이 가장 많은 칸**(없으면 서열 **1746개**에 가장 가까운 칸). 출력 `maintest.csv`. 가짜 사다리 5경우 검증 완료.
  - **⭐ 본 검정 규모(유동):** **코로나 10종 = 전수 유지 권장**(항원 195잔기라 타깃당 48분, 10개 8시간 — 싸고 유일하게 선별 없는 세트). **비싼 쪽은 44종**(타깃당 2.7시간). 44종은 **개수 고정하지 않고 우선순위 큐**: 8k3k_D → 8wpy_AB → 8tx3_FK → 8k46_I → 8xsj_HL → 8y6a_CD → 8ume_HL. 음성 대조 2(9azr_HL·8k5g_HL)는 유지. **4개까지면 하루, 7개까지면 이틀.**
  - **⚠️⚠️ protenix 예측 폴더 사고(2026-07-27~28) — 반복 금지:** 07-27 15:05 a3m 수습 중 `$DATA/protenix`를 stale로 옮겼으나 **Protenix는 그 사고와 무관**(clean_a3m이 07-22부터 적용, MSA를 잃은 건 boltz뿐). 폴더가 비니 `run_sweep.sh`의 self-heal skip이 무력화되어 07-27 밤 sweep이 **44종 중 26종을 재계산하고 12h 예산 소진**, RBD 10종(목록 맨 뒤)은 미실행. → `run_fix_protenix_dirs.sh`(`612c7802`)로 중복 26 삭제 + 원본 49 복귀 + `$DATA/README_protenix.txt` 기록. **교훈: 출력 폴더를 옮기면 skip이 무력화되어 목록 앞쪽부터 전부 재실행된다. 새 타깃만 돌리려면 `LIST=`로 목록을 제한할 것.** 정본 = `results/pose_features.csv`, **`--rescore` 금지**.
  - **RBD 10종 실측:** a3m은 **21,737~27,946개**(`out_colabfold/<pdb>_env/uniref.a3m` 또는 boltz `..._unpaired_tmp_env/uniref.a3m`), 항원 194~197잔기. **칸당 1~1.5분**(44종은 5~6분)이라 사다리 10종 12칸이 **약 3시간**. 깊이가 깊을수록 느린 것 = MSA를 실제로 읽고 있다는 증거.
  - **RBD 해석 주의:** 코로나 세트에서 더 잘 나오면 "인기 부위 밖 항체라서"로 읽고 싶겠지만, 그 세트는 **항원 계열·MSA 크기·항원 길이도 함께 다르다** → **가설로만** 서술하고 단정 금지.

**Why:** 사용자가 이 방향에 감을 잡았고, 9K6J-modal 반증 + 2026-07-14 go/no-go 재프레임 + Phase 0 결과(GO) + **2026-07-15 MSA-깊이 인과실험(항원 MSA가 편향 나름)** + **2026-07-20 데이터셋-우선/A-vs-B 척추**가 가설·프로젝트 방향을 확정한 핵심 분기점이라 정확히 남겨야 함(미래 세션이 9K6J에 에피토프-prior를 다시 밀거나, 프로젝트 생사를 n=1 MSA 스토리 또는 얇은 consensus-격차에 걸면 안 됨).
**실험 로그 위치:** ⭐**2026-07-20 사용자 지정 전용 페이지 = Notion "MSA Bias"**(page **3a3ac1a3-a28a-80ce-b19a-c6f145cd6762**, parent 388ac1a3-a28a-809b-80d7-eabf91ba17b5). **앞으로 MSA-bias 데이터셋·depth-response 작업은 전부 여기 기록**(이전 Study B page 39dac1a3는 초기 로그). 토큰=scratchpad/.notion_tok, append는 end-append만([[notion-api-gotchas]]). 실험 = E1(Neff+Foldseek/PDB 대조, 착수) · E2(과대표집 항원 벤치 tFold vs MSA모델 n>1=척추) · E3(2×2 epitope-only×Fab/Fv) · E4(HADDOCK crop, E3 게이트). 방법론: Neff=가설지표(tFold와 MSA모델 가르는 유일 채널), Foldseek·PDB=대조군(deconfound; Spike는 세 축 공변→divergent 타깃 필요). ablation·다른 MSA-free 복합체모델=폐기/없음.

**How to apply:** ⭐⭐⭐**재개 = 맨 위 두 항목만 읽으면 된다** — ①'2026-07-27(밤, 최종) 연구 전체 정리'(핵심 이야기·기각된 것·판정 기준, = 노션 6장) ②'2026-07-28 설계·통계 확정'(설계 6×4+8·깊이 선택 규칙·Fisher 결합·폴더 사고). 그 아래는 이력. **현재 상태(2026-07-28 오전)** = RBD 10종 사다리 예측 진행 중(3h) → `python lib_pose_features.py --models protenix`(CPU 1h, **--rescore 금지**) → `python prep_pick_depth.py`(본 검정 대상·깊이 자동 선택) → `run_maintest.sh`(아직 미작성) → `run_analyze_target.sh` → 판정. **실험 종료 예상 = 2026-07-29 저녁.** 첫 할 일 = `head -3 sweep_targets.csv`로 열 형식 확인 → `prep_import_rbd.py` 열 매핑 수정 + 8sdh 제거 → `--apply --group RBD`. 그 밖 지속 규율: (a)결론은 **"MSA 조성이 결합자리 선택을 좌우"**(존재증명 확보, 8ulr/Protenix p=3.9e-8) — 단어 규율 엄수(깊이↓rescue·과대표집 편향 금지) (b)다음 한 수 = **사전등록 규칙으로 타깃 6~8개 빈도 M/N 측정**, Protenix로 좁혀서 (c)처방은 항상 **생성+선택 한 쌍**으로 서술 (d)novelty 주장 시 선점 5종(Subsampled AF2·AF-Cluster·SPEACH_AF·AFsample2·SF-Cluster)+Porter PNAS2024 반드시 명시하고 슬리버는 "에피토프 선택 적용 + 예산맞춘 통제"로만.
이전 이력 참고: (a)9K6J=modal(에피토프-prior 부적합) (b)프로젝트 GO 근거=consensus 여지+ipTM 대체 (c)지표=DockQ 0.49+RMSD 이원화. 선점 추가: peptide training-bias(2025, PMC12518507), ab-ag systematic-biases/TERM(2024, bioRxiv 2024.03.15.585121), depth-over-pairing(2026). [[hackathon-consensus-docking]]와 모델 풀 공유, [[guan-keating-msa-docking-bias]]·[[msa-debias-technique-landscape]]·[[no-invented-jargon]] 참조.
