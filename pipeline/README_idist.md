# iDist 계면 과대표집 정량 — 파이프라인

**목적.** 각 테스트 항체-항원 **계면(에피토프+파라토프)**이 "모델이 학습한 컷 이전 PDB(SAbDab)"에 얼마나
흔한가 = **과대표집 점수**를 계면 단위로 잰다. A/B/C 손라벨을 **데이터 기반 연속 점수**로 바꾸고,
이 점수가 우리 47복합체의 **깊이-취약성**(epitope recall/DockQ의 MSA 깊이 반응)과 상관되는지로
**위치편향 인과 축**을 검증한다.

> ⚠️ 왜 Foldseek이 아니라 iDist인가: Foldseek은 "복합체/항원 fold 전체" 유사도라 **계면 국소를 놓친다**
> (Bushuiev et al. 2024가 명시). A/B는 **같은 항원의 다른 에피토프**(같은 fold) → Foldseek은 A vs B를
> 못 가름. 계면 과대표집의 정석 = **iAlign(IS-score)** 또는 그 스케일업 근사 **iDist(PPIRef)**. Foldseek은
> 항원 fold coarse dedup 보조로만.

## 방법 novelty 없음 (선점 도구 인용)
iDist/iAlign은 기성 도구. 여기선 **"항체-항원 계면 과대표집·누수 정량"의 재현·적용**으로만 씀.
- iDist / PPIRef — Bushuiev et al., ICLR 2024, "Revealing data leakage in PPI benchmarks" (github.com/anton-bushuiev/PPIRef)
- iAlign(IS-score) — Gao & Skolnick, Bioinformatics 2010 · PINDER(계면 non-redundant split) · AsEP(epitope-group split)

## 설치 (서버, 사수 허가 후 — 순수 파이썬, 외부 바이너리 불필요)
```bash
conda create -n ppiref python=3.11 -y && conda activate ppiref
pip install git+https://github.com/anton-bushuiev/PPIRef.git
pip install biopython scipy pandas matplotlib
```
계면 추출+iDist 임베딩은 biopandas 접촉 기반이라 DSSP/reduce/foldseek이 **안 걸린다**. (iAlign 교차검증을
쓸 때만 Perl iAlign 별도.)

## 실행 순서 (`pipeline/`에서, ppiref env)
```bash
# 0) 레퍼런스: 컷 이전 SAbDab 항체-항원 매니페스트 → 구조 다운
python idist_ref_manifest.py                                   # → ref_manifest.csv (RBD/HA/Env, date<2023-06-01)
python fetch_structures.py --manifest ref_manifest.csv --outdir ref_structures
# (테스트 구조 targets/ 는 기존 재생성: fetch_structures + prep_targets)

# 1) 계면 추출(테스트+레퍼런스) — 6Å 접촉, 캐시 재사용
python idist_overrep.py --stage extract

# 2) 임베딩 → 이웃수(과대표집) 채점
python idist_overrep.py --stage score                          # → results/overrep_idist.csv

# 3) 가설 검정: A>B? 과대표집 vs 깊이-취약성 상관 + 그림
python idist_analyze.py                                        # → results/idist_overrep.png
```
⚠️ **먼저 2개 복합체로 스모크** 권장(코드가 로컬 미검증): `sweep_targets.csv`를 2행으로 줄여 --stage all 후
계면 파일·이웃수가 나오는지 확인.

## 핵심 설계 (근거)
- **병합**: native → 2체인 PDB(항원=A, 항체 H+L=B). `dockq_sweep.merged_pdb` 재사용 → 내부 VH-VL 계면 제외
  (H·L을 한 파트너로 안 묶으면 모든 항체 공유 프레임워크 계면이 섞여 과대표집이 편향됨).
- **추출**: `PPIExtractor(kind='heavy', radius=6.0)`. **radius=6.0 ↔ near-dup 임계 0.04는 반드시 짝**
  (기본 radius=10.0이라 명시 지정 필수).
- **점수**: `IDist` 임베딩(비파라미터 1-step message passing) → 테스트 계면의 레퍼런스 내 유클리드 이웃수.
  `n_0.04`=보정된 near-duplicate(거의 동일), `n_0.08/0.15/0.30`=탐색용(**비보정** — 등급 참고만, 임계 주장 금지).
- **매칭**: 항원군(RBD/HA/Env) 단위. 테스트 A/B ↔ 같은 항원군 컷 이전 레퍼런스.

## 해석 (가설)
- **검정1**: 과대표집 **A(우세) > B(비우세)** — 성립하면 A/B 라벨 + "B=희귀 에피토프" 검증.
- **검정2**: 과대표집 클수록 **full-MSA에서 잘 붙고**(prior가 modal 자리 도움), 낮을수록 **깊이-취약**
  → 위치편향(과대표집 prior가 배치를 좌우)을 데이터로 확정.

## 한계 (정직)
- 과대표집=raw 이웃수라 **레퍼런스 크기에 스케일 의존** → A/B **대비(비율)**로 해석. `--max-per-family` 상한(초과 로그).
- leakage 기준 = SAbDab `date`(release). **최종 lock은 RCSB release date 재확인** 권장(현재 1차 필터).
- C군(신규 항원)은 RBD/HA/Env 분류 밖 → v1에서 항원군-매칭 레퍼런스 없음(점수 미측정). C는 항원명-매칭 레퍼런스로 v2.
- 컷오프 통일: 우리 접촉 코드가 5.0/4.5Å 혼재였으나 iDist는 **PPIExtractor 6.0Å로 일원화**(임계와 짝).
- '유사(관련)' 느슨한 임계는 공식 보정 없음 → 경계 사례는 iAlign IS-score/P-value로 교차확인.