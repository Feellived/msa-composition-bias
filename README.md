# msa-composition-bias

서울대학교 의과대학 항체·면역학 실험실(이창한 교수 연구실) **BK21 여름 학부연구 2026** 코드 저장소.

---

## 무엇을 하는 저장소인가

**항원 MSA(다중서열정렬)의 조성이 항체–항원 co-folding 모델의 결합 자리 선택을 좌우하는지
검정하고, 그 성질을 결합 자리 후보 생성에 이용하는 실험 코드.**

MSA의 **깊이(서열 수)를 고정한 채 구성원만 다시 뽑으면**, 같은 복합체·같은 모델·같은 예산
아래에서도 모델이 지목하는 항원 표면이 바뀐다. 이 저장소는 그 현상을 30개 복합체 × 960회
실행으로 측정하고, 조성별 결합 자리를 모아 후보를 만드는 데까지를 담는다.

> **상세 인수인계 문서는 Notion 「인수인계서 Ⅱ · MSA Bias」에 있다.**
> 연구 경과·복합체 명세·통계 설계·결과 해석·알려진 함정이 그쪽에 정리돼 있다.
> 이 README는 저장소를 열었을 때 길을 잃지 않기 위한 최소 안내다.
>
> 후보 자리를 **고르고**(선택) 제약으로 주어 **다시 접는**(유도 재도킹) 단계는 형제 저장소
> **`epitope-guided-docking`** 과 「인수인계서 Ⅰ · Consensus Docking」에 있다.

---

## 확정된 결론

| 주장 | 상태 |
|---|---|
| 항원 MSA 조성이 결합 자리 선택을 좌우한다 | **성립.** 조성 내 재현성이 조성 간 재현성보다 높다(자카드 0.503 대 0.413, 부호검정 p = 1.9×10⁻⁹). 복합체별 순열검정에서 30종 중 26종 유의(이항 p = 3.4×10⁻³⁰) |
| 조성 재추첨이 시드 변경보다 자리를 크게 움직인다 | **성립.** 30종 중 24종에서 우세, 중앙값 차이 +0.47 (p = 1.8×10⁻⁴) |
| 조성을 바꾸면 예측이 좋아진다 | **불성립.** 방향이 갈린다. 그래서 **후보 중 하나를 고르는 선택 단계가 반드시 필요하다** |
| 조성 재추첨으로 자세(pose)까지 맞는다 | **불성립.** 960회 중 DockQ ≥ 0.49는 4회. 정확한 서술은 **"자리는 찾되 자세는 못 맞힌다"** |
| 정답 없이 후보 자리를 고를 수 있다 | **현재로선 불가.** 선택기 다섯 개 전부 무작위 선택과 구별되지 않는다(순열 p ≥ 0.074) |

### ⛔ 되살리지 말 것 — 검정으로 기각된 주장

- **"MSA 깊이를 낮추면 정답 자세가 되살아난다"**(전면 적용). 순열검정에서 잡음으로 확정됐다.
  깊이 감소는 조성을 바꿀 수 있게 하는 **전제조건**일 뿐, 그 자체가 효과가 아니다.
- **"과대표집 항원일수록 깊이 감소 효과가 크다"**. 적대적 재검증에서 인공물로 판정됐다.
- **`p = 3.9×10⁻⁸`**. 통제 실행 1회(자세 40개)를 독립 표본 40개로 센 값이다. 폐기.
- **"Boltz-2는 구조상 MSA에 둔감하다"**. 근거였던 실행이 a3m 사고로 무효였고, 예산을 맞춘
  정식 시험에서 3종 전부 기각됐다. 조성 효과는 현재까지 Protenix 계열에서 확인된 것이다.
- **'조성 다양성'(Neff80 층) 축**. 층이 항원 계열과 사실상 같아 교란을 분리할 수 없다.

---

## 실험 설계

### 조성 재추첨

항원 MSA에서 **깊이(서열 수)를 고정한 채 구성원만 다시 뽑는다.** 질의 서열은 항상 포함한다.
복합체 하나당 **조성 6가지 × 반복 4회 + 원래 MSA 8회 = 32회** 실행하고, 실행 1회마다
자세 5개를 얻는다. 30개 복합체 전체로 **960회 실행**이다.

실행 1회의 **결합 자리**는 그 실행의 자세 5개 중 **과반에서 관찰된 항원 잔기 집합**으로
정의한다. 자세 하나가 아니라 실행 단위로 정의해야 자세 순위의 흔들림에 휘둘리지 않는다.

### 세 가지 판정 기준

같은 데이터를 세 가지 기준으로 읽는다. 결론이 기준에 따라 달라지는지 보기 위해서다.

| 기준 | 결과 |
|---|---|
| 조성별 결합 자리가 서로 다른가 (엄격) | 30종 중 **7종** |
| 조성 내 재현성이 조성 간 재현성보다 높은가 | 30종 중 **26종** |
| 조성 효과가 시드 효과보다 큰가 | 30종 중 **23종** |

다중비교는 **Bonferroni** 보정을 쓴다(선택기 9개 비교 시 문턱 0.0056).

### 지표

- **자카드 지수** — 두 결합 자리 집합의 겹침. 재현성 비교의 주 지표
- **에피토프 F1** — 예측 자리와 정답 자리의 조화평균. 덮음(recall)만 보면 자리를 크게
  잡을수록 유리해지므로 F1으로 바꿨다
- **DockQ** — 자세 품질. CAPRI 기준 0.23(acceptable)·0.49(medium)·0.80(high)
- **Neff80** — 80% 서열동일성 기준 잔기별 유효 서열 수의 중앙값. 깊이 칸을 정하는 데 쓴다

---

## 디렉토리

| 경로 | 내용 |
|---|---|
| `pipeline/scripts/` | 코드 97개. 명명 규칙은 다음 절 |
| `pipeline/*.csv` | 타깃 명단·매니페스트 19개. 실행은 `pipeline/` 안에서 하므로 상대경로로 읽힌다 |
| `pipeline/results/` | 커밋되는 집계 산출물. `.log`·`.bak`은 추적하지 않는다 |
| `report/` | 최종 보고서 본문·그림 스크립트(`plot_*.py`)·그림 |
| `report/data/` | 재분석용 CSV 약 300개. **GPU 없이 결과를 다시 확인할 수 있는 유일한 자산** |
| `plan/` | 사전 등록 문서와 실행 계획. 시점 기준 문서이므로 상단 배너를 확인할 것 |

대용량 산출물(구조·MSA·예측 자세)은 저장소에 넣지 않는다. 데이터 루트는 환경변수 `DATA`로
지정하며, 하위 구조는 `pipeline/DATA.md`에 있다.

---

## 스크립트 명명 규칙

파일명 앞부분만 보고 파이프라인의 어느 단계인지 알 수 있도록 여덟 개 접두사로 통일했다.
형제 저장소 `epitope-guided-docking`과 같은 규칙이다.

| 접두사 | 의미 | 예 |
|---|---|---|
| `prep_` | 데이터셋 조립·구조/MSA 전처리 | `prep_manifest.py` `prep_ladder.py` `prep_a3m_check_match.py` |
| `make_` | 모델 입력·MSA 생성 | `make_msa.sh` `make_composition_reps.sh` |
| `run_` | 실행 오케스트레이션 | `run_maintest.sh` `run_seed_replicate.py` |
| `eval_` | 정답을 보고 채점 | `eval_dockq_sweep.py` `eval_compreps.py` |
| `select_` | 정답을 보지 않고 후보 선택 | `select_eval_selectors.py` `select_refine_sites.py` |
| `analyze_` | 집계·통계·진단 | `analyze_perm_null.py` `analyze_site_reproducibility.py` |
| `plot_` | 그림 생성 (`report/`) | `plot_agg.py` `plot_refine.py` |
| `lib_` | 공유 라이브러리(CLI 없음) | `lib_pose_features.py` `lib_epitope_recall.py` `lib_epitope_defs.py` |

`lib_pose_features.py`가 채점의 중심 모듈이다. DockQ·결합 자리 추출·에피토프 지표를 담고
있으며, 형제 저장소 `epitope-guided-docking`도 이 모듈을 불러 쓴다.

---

## 실행 순서

**모든 명령은 `pipeline/` 안에서 실행한다.** 타깃 명단 CSV가 상대경로로 읽히기 때문이다.

```bash
cd ~/projects/msa-composition-bias/pipeline
export DATA=/path/to/large-disk/msadepth

# 1) 데이터셋 조립 — SAbDab → 매니페스트 → 에피토프 라벨 → 파일럿 선별
python scripts/prep_fetch_sabdab.sh
python -u scripts/prep_manifest.py
python -u scripts/prep_classify_epitope.py
python -u scripts/prep_select_pilot.py

# 2) MSA 생성과 깊이 사다리
bash scripts/make_msa.sh
python -u scripts/prep_ladder_neff.py
python -u scripts/prep_pick_depth.py        # 복합체마다 깊이 칸 선택

# 3) ⭐ 실행 전 게이트 — a3m 질의행 오염 점검 (건너뛰지 말 것)
python -u scripts/prep_a3m_check_match.py

# 4) 조성 재추첨 실행 (조성 6 × 반복 4 + 원래 MSA 8)
bash scripts/make_composition_reps.sh
bash scripts/run_maintest.sh

# 5) 채점
python -u scripts/eval_compreps.py
python -u scripts/eval_fullmsa_control.py

# 6) 통계·집계
python -u scripts/analyze_site_reproducibility.py     # 조성 내 대 조성 간 자카드
python -u scripts/analyze_perm_null.py                # 복합체별 순열검정
python -u scripts/analyze_seed_vs_comp.py             # 조성 효과 대 시드 효과

# 7) 그림
cd ../report && python -u plot_agg.py
```

전체 인자는 `--help`로 확인한다.

---

## conda 환경

**환경은 `boltz` 하나로 충분하다.** DockQ·biopython·scipy가 모두 여기 있다.
**채점과 통계 분석도 반드시 이 환경에서 실행한다** — 환경 밖에서 돌리면 `dockq` 열만 조용히
비고 나머지 열은 정상이라 결과가 있는 것처럼 보인다.

---

## ⚠️ 이 파이프라인은 멈추지 않고 조용히 빈다

실패해도 종료 코드가 정상이고 결과 파일도 생기는 사고를 겪었다.
**채점 결과를 믿기 전에 실행 수가 설계값(복합체당 32회)과 맞는지 먼저 확인한다.**

| 유형 | 증상 |
|---|---|
| a3m 질의행 오염 | 모델이 MSA를 통째로 버리고 경고만 남긴 뒤 정상 종료한다. 조성 조작이 아무 효과도 못 낸다 |
| 깊이·조건 폴더 혼입 | 한 폴더에 두 조건이 섞여 자세 수가 배로 늘고 최고값이 부풀려진다 |
| DockQ 환경 밖 실행 | `dockq` 열만 통째로 비는데 나머지 지표는 정상이라 결과가 있어 보인다 |
| 비교군 예산 불일치 | 조건마다 자세 예산이 달라 비교가 성립하지 않는다 |

- **실행 전** — `scripts/prep_a3m_check_match.py` 게이트를 통과시킨다
- **실행 후** — `run.log`에서 경고를 검색한다. 다만 항체 사슬에 대한
  `Found explicit empty MSA for some proteins`는 **정상 메시지**다(항원에만 MSA를 준다)

---

## 형제 저장소와의 관계

이 저장소는 **① 생성**(조성 재추첨으로 결합 자리 후보를 넓힌다)과 **② 후보 조립**(조성 간
투표로 후보 자리를 세운다)을 담당한다. **③ 선택**과 **④ 유도 재도킹**은 형제 저장소
[`epitope-guided-docking`](https://github.com/Feellived/epitope-guided-docking)에 있다.

의존은 양방향이다.

- [`epitope-guided-docking`](https://github.com/Feellived/epitope-guided-docking)의 데모·honest·전체후보
  흐름이 이 저장소를 환경변수 `MSAD`(기본값 `~/projects/msa-composition-bias/pipeline`)로 참조한다.
  공유 모듈 `lib_pose_features`·`lib_epitope_recall`을 `$MSAD/scripts`에서 불러오고,
  결합 자리 후보 `$MSAD/results/sites_<타깃>.json`을 읽는다
- 이 저장소의 재선택·후보천장 흐름은 `epitope-guided-docking`의 스크립트를 호출한다

**한쪽 저장소의 파일명을 바꾸면 반대쪽도 함께 고쳐야 한다.**

---

## 비고

- 예측 구조·MSA·weight·대용량 CSV는 저장소에 올리지 않는다(`.gitignore`).
- 서버 접속 정보(IP·포트·계정·비밀번호)와 API 토큰은 문서·코드·커밋에 남기지 않는다.
- 코드는 붙여넣기가 아니라 커밋 후 서버에서 `git pull`로 전달한다.
