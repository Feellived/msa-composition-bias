# 데이터 획득·저장 (provenance)

모든 데이터가 어디서 와서 어디에 쌓이는지. **다운로드 데이터**와 **생성 데이터**를 구분한다.

## 흐름: 원천 → 확정 세트 → 생성물

```
[다운] SAbDab2 summary  ──prep_fetch_sabdab.sh──▶  sabdab2_summary.csv
                                                   │
                          build_manifest → classify_epitope → select_pilot   (Mac에서 1회 수행, 결과 committed)
                                                   ▼
                                    ★ pilot_lean_full.csv (49, FROZEN 확정 세트)   ← 재현은 이걸 그대로 사용
                                                   │
[다운] RCSB 구조  ──prep_fetch_structures.py──▶  structures/         (49 PDB)
                          prep_targets.py ──▶  targets/          (chains.json·항원 fasta·native)
[생성] colabfold  ──make_msa.sh + build_ladder──▶  ladders/       (항원 MSA → Neff80 사다리)
[생성] co-folder  ──run_sweep.sh──▶  boltz/ protenix/ ...        (pose)
```

⚠️ **SAbDab은 live DB** — `prep_fetch_sabdab.sh`로 다시 받아 `build_manifest`를 재실행하면 세트가 **달라진다**.
확정 세트를 재현하려면 committed `pilot_lean_full.csv`를 그대로 쓰고 구조만 받으면 된다(아래 최소 재현).
`prep_fetch_sabdab.sh` + `build_manifest` 계열은 **"어떻게 만들었나"(provenance) 재현·갱신용**.

## 저장 위치 (다운로드 vs 생성)

| 데이터 | 종류 | 위치(권장) | 대략 크기 | 만든 것 |
|---|---|---|---|---|
| SAbDab summary | 다운 | `pipeline/` | ~11 MB | `prep_fetch_sabdab.sh` |
| 구조 (cif) | 다운 | `structures/` 또는 `/mnt/data/.../structures` | ~수백 MB | `prep_fetch_structures.py` (RCSB) |
| chains·fasta·native | 생성 | `targets/` | 작음 | `prep_targets.py` |
| MSA 사다리 (a3m) | 생성 | **`/mnt/data/msadepth/ladders`** | ~GB | `make_msa.sh` |
| pose (cif) | 생성 | **`/mnt/data/msadepth/{boltz,protenix}`** | 수십~수백 GB | `run_sweep.sh` |

## ⚠️ 홈 디스크 보호

서버 홈(`/`)은 여유가 적다(예: 97% 사용). **큰 것은 전부 `/mnt/data`로.** MSA·pose는 이미 `DATA` env로 그쪽에 간다.
**구조도 홈에 안 쌓으려면** `--outdir`을 `/mnt/data`로:

```bash
DROOT=/mnt/data/msadepth
python prep_fetch_structures.py --manifest pilot_lean_full.csv --outdir "$DROOT/structures"
python prep_targets.py     --csv pilot_lean_full.csv --struct "$DROOT/structures" --outdir targets
#   targets/ 는 작아서 레포 폴더에 둬도 됨(gen_msa·run_sweep이 $HERE/targets 를 봄).
```

## 최소 재현 (확정 세트 그대로)

```bash
cd pipeline && conda activate boltz
python prep_fetch_structures.py --manifest pilot_lean_full.csv --outdir structures   # 구조만 받으면 됨
python prep_targets.py     --csv pilot_lean_full.csv --struct structures --outdir targets
bash make_msa.sh                                                                  # MSA·사다리 → /mnt/data
bash run_sweep.sh boltz 11                                                       # pose → /mnt/data
```

## 원천 URL·참조

- SAbDab2 summary: `https://sabdab.opig.stats.ox.ac.uk/api/download/all-summary` (콤마 CSV, GET only)
- 구조: RCSB `https://files.rcsb.org/download/<PDB>.cif`
- 참조 서열(HA/Env 정렬): UniProt `P03437`(H3 HA)·`P04578`(HIV Env) — `prep_classify_epitope.py`가 런타임 fetch·캐시
