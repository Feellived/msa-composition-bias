# 환경 · 의존성

무거운 모델은 각각 **별도 conda env**로 설치(가중치 자동 다운로드). 스크립트별로 필요한 env가 다르다.

## 1. 데이터셋·준비·채점 (가벼움, CPU)

`build_manifest`·`classify_epitope`·`prep_targets`·`build_ladder`·`dockq_sweep` 등.

```bash
conda create -n msadepth python=3.10 -y && conda activate msadepth
pip install biopython pandas numpy DockQ
# (build_ladder의 Neff80 계산엔 numpy만, classify_epitope엔 biopython)
```
- **DockQ** = pose 채점(`eval_dockq_sweep.py`)에 필수. `which DockQ`로 확인.
- 서버 실무에선 `boltz`/`DockQ` 깔린 env에 biopython만 추가해 써도 됨(전에 boltz env에서 채점함).

## 2. MSA 생성 (CPU/네트워크)

`make_msa.sh` — colabfold로 항원 MSA 생성.
```bash
# ColabFold 설치 (예: localcolabfold) → colabfold_batch 사용
#   make_msa.sh 의 MSA_CMD 기본값 = "colabfold_batch --msa-only"
```
서버 colabfold 명령이 다르면 `MSA_CMD` env로 조정.

## 3. co-folder 생성 (GPU, 각각 별도 env)

`run_sweep.sh` — 모델별 env를 `PROT_ENV` 등으로 지정. 가중치는 대개 자동 다운로드.

| 모델 | env | 설치 | 비고 |
|---|---|---|---|
| Boltz-2 | `boltz` | `pip install boltz` | `--use_msa_server` 내장 / 여기선 사다리 a3m 주입 |
| Protenix-v1 | `protenix` | `pip install protenix` | 기본 `protenix_base_default_v1.0.0`(컷오프 **2021-09-30**=leakage-free, 공식 권장·AF3 능가). ⚠️`protenix_base_20250630`=2025컷오프=**leaky 금지**; v2(동일 컷오프·더 강함)는 공개 다운로드 403이라 보류. `LAYERNORM_TYPE=torch`+`--trimul/triatt_kernel torch`로 CUDA 커널 JIT 회피 |
| Chai-1 | `chai` | `pip install chai_lab` | `CHAI_ENV`로 지정. MSA=항원 `aligned.pqt`(make_input이 chai_lab 변환기로 생성), 항체 single-seq. 첫 실행은 스모크 권장 |

- Protenix 실행 env: `PROTENIX_ROOT_DIR`(가중치 경로)·`LAYERNORM_TYPE=torch` 설정(`run_sweep.sh`가 처리).
- 물리 생성기(HADDOCK·SnugDock)는 이 파이프라인 범위 밖(계획서상 decorrelated 생성기, 보조).

## 스크립트 → 필요 env 요약

| 스크립트 | env |
|---|---|
| fetch_sabdab · fetch_structures · prep_targets · classify_epitope · build_ladder | msadepth (biopython/numpy) |
| gen_msa | colabfold |
| run_sweep boltz | boltz |
| run_sweep protenix | protenix (PROT_MODEL=base) |
| dockq_sweep | DockQ 설치된 env (msadepth 또는 boltz+DockQ) |
