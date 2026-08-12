# depth-sweep 생성 — 서버 실행 흐름 (lean 파일럿 49)

세트 = `pilot_lean.csv` (RBD 7A+7B, HA 7A+7B, Env 6A+6B, C 9 재활용). 순서 = `sweep_targets.csv`(round-robin).
출력 = **`/mnt/data/admuser/msadepth/`** (홈 디스크 97%→금지). co-folder 순서 = **Boltz → Protenix → Chai**.
depth = rung 6개(full→single, 사슬별 Neff80 라벨). pose = diffusion sample 5.

## 0. 착수 전 확인 (사수 GPU 허락 + 서버값)

- **GPU 며칠 점유** = 사수 허락 필수. 반나절(밤샘)/주말 청크로.
- 아래 3개를 서버에 맞게 확인:
  - `MSA_CMD` (기본 `colabfold_batch --msa-only`) — 서버 colabfold 명령.
  - `PROT_MODEL` = **Protenix base(2021-09) 체크포인트명** (⚠️ 2025 `protenix_base_20250630`는 leaky, 금지).
  - `DIVERSE` = C 재활용 위치(기존 `runs_diverse`). `PROTENIX_ROOT_DIR` = weights 경로.

## 1. 준비 (GPU 불필요)

```bash
cd ~/projects/bk21-antibody-ml && git pull
cd pipeline
conda activate boltz     # biopython 필요

# (a) A/B 구조 fetch + 입력 준비 (C는 runs_diverse 재활용이라 제외)
python prep_fetch_structures.py --manifest pilot_lean.csv --outdir structures   # A/B 40개(+C skip)
python prep_targets.py --struct structures --outdir targets                # targets/<id>/chains.json

# (b) 항원 MSA + depth 사다리 (CPU/네트워크, self-heal)
ONLY=8q7s_O bash make_msa.sh      # ⚠️ 먼저 한 타깃 smoke — colabfold 동작·a3m 생성 확인
bash make_msa.sh                  # 전체 → /mnt/data/admuser/msadepth/ladders/<target>/<chain>/rung{0..5}.a3m
```

## 2. smoke test (1건, 타이밍 캘리브레이션)

```bash
SMOKE=1 bash run_sweep.sh boltz 1
# → /mnt/data/admuser/msadepth/boltz/<target>/rung0/results/*.cif 확인 + 소요시간 측정
#   1건 소요시간 × 49 × 6 으로 청크 크기 재보정.
```

## 3. 청크 실행 (시간-박스, 재개형)

```bash
tmux new -s sweep
# 평일 퇴근(밤샘 ~11h):
bash run_sweep.sh boltz 11
# 금요일 퇴근(주말 ~54h):
bash run_sweep.sh boltz 54
# Boltz 다 끝나면 Protenix:
PROT_MODEL=<base모델명> bash run_sweep.sh protenix 54
```
- 예산 소진 시 자동 정지, **재실행하면 완료분 skip하고 이어감**(kill돼도 안전).
- 복합체 단위 완결 순서라, 중간에 멈춰도 "완전히 스윕된 복합체 K개"(그룹 균형) = 바로 분석 가능.
- Chai는 추후(MSA 주입 러너 신규 필요).

## 산출물 구조

```
/mnt/data/admuser/msadepth/
  ladders/<target>/<chain>/rung{0..5}.a3m + neff.tsv   # 사슬별 depth 사다리(Neff80 축)
  boltz/<target>/rung{0..5}/results/*.cif              # pose (5 sample)
  protenix/<target>/rung{0..5}/results/*.cif
```

## 다음 (생성 후)

- `eval_dockq_sweep.py`(추후): pose별 DockQ(다사슬 항원 native 대비, RBD는 크롭 native) → depth-response 곡선.
- 분석: rung(=Neff80)별 pose 이동 / off-site rescue 비율 A vs B vs C / 지배-centroid 방향성 shift.
- HADDOCK·SnugDock = CPU 병렬(별도).

## 파일

| 파일 | 역할 |
|---|---|
| pilot_lean.csv / sweep_targets.csv | 49 세트 / round-robin 순서 |
| prep_targets.py | 구조→chains.json(다사슬 항원·RBD크롭)+native |
| make_msa.sh + prep_ladder.py | 항원 MSA + rung 사다리(Neff80) |
| make_input.py | chains.json+rung a3m → Boltz YAML / Protenix JSON |
| run_sweep.sh | 시간-박스 재개형 dispatcher (boltz\|protenix) |
