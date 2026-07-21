#!/usr/bin/env bash
# 시간-박스 재개형 depth-sweep dispatcher (co-folder 하나씩).
# 복합체 단위로 rung 전부 → 다음 복합체(round-robin 순서=중간에 끊겨도 그룹 균형).
# 예산 소진 시 깔끔 정지, 재실행 시 self-heal(출력 있으면 skip).
#   bash run_sweep.sh boltz 11     # Boltz, 11시간(반나절/밤샘)
#   bash run_sweep.sh protenix 54  # Protenix, 주말
#   bash run_sweep.sh chai 11      # Chai-1 (CHAI_ENV=chai)
#   SMOKE=1 bash run_sweep.sh boltz 1   # 1건만
set -uo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"          # .../consensus_docking/dataset
cd "$HERE" || exit 1
COF="${1:?boltz|protenix|chai}"; HOURS="${2:-11}"
DATA="${DATA:-/mnt/data/admuser/msadepth}"      # ⚠️ 대용량 출력 = /mnt/data(홈 금지)
LADDIR="$DATA/ladders"; RUNGS="${RUNGS:-6}"; SAMP="${SAMP:-5}"; SEED="${SEED:-0}"; SMOKE="${SMOKE:-0}"
DIVERSE="${DIVERSE:-$HERE/../runs_diverse}"      # C 재활용 위치(서버 확인)
LIST="${LIST:-$HERE/sweep_targets.csv}"
say(){ echo "[$(date '+%m-%d %H:%M:%S')] $*"; }
DEADLINE=$(( $(date +%s) + HOURS*3600 ))
source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null

case "$COF" in
  boltz)
    conda activate boltz 2>/dev/null
    command -v boltz >/dev/null || { say "!! boltz 없음"; exit 1; }
    EXT="yaml"
    RUN(){ ( cd "$2" && boltz predict "$1" --out_dir results --no_kernels --diffusion_samples "$SAMP" >"$3" 2>&1 ); }
    DONE(){ find "$1/results" -name '*_model_*.cif' 2>/dev/null | grep -q .; } ;;
  protenix)
    conda activate "${PROT_ENV:-protenix}" 2>/dev/null
    export PROTENIX_ROOT_DIR="${PROTENIX_ROOT_DIR:-/mnt/data/admuser/protenix_weights}" LAYERNORM_TYPE=torch
    command -v protenix >/dev/null || { say "!! protenix 없음"; exit 1; }
    # 모든 Protenix 공개 체크포인트 = 학습컷오프 2021-09-30(AF3정렬) → post-2023-06 세트에 leakage-free.
    # 표의 2025/2026은 출시일일 뿐(컷오프 아님). 기본=권장·최강 v1.0.0.
    PROT_MODEL="${PROT_MODEL:-protenix_base_default_v1.0.0}"
    EXT="json"
    RUN(){ ( cd "$2" && protenix pred -i "$1" -o results -n "$PROT_MODEL" -s "$SEED" -e "$SAMP" \
             --trimul_kernel torch --triatt_kernel torch >"$3" 2>&1 ); }
    DONE(){ find "$1/results" -name '*sample*.cif' 2>/dev/null | grep -q .; } ;;
  chai)
    conda activate "${CHAI_ENV:-chai}" 2>/dev/null
    command -v chai-lab >/dev/null || { say "!! chai-lab 없음 (conda activate ${CHAI_ENV:-chai})"; exit 1; }
    EXT="fasta"    # make_input(chai)가 FASTA + $out/msa/*.aligned.pqt(항원, 파일명=서열해시) 생성. 항체=pqt없음→single-seq
    # chai-lab fold는 출력 폴더가 이미 있으면 실패 → 미완성 results 제거 후 실행(DONE이면 여기 안 옴)
    RUN(){ rm -rf "$2/results"; ( cd "$2" && chai-lab fold --msa-directory "$2/msa" "$1" results >"$3" 2>&1 ); }
    DONE(){ find "$1/results" -name '*.cif' 2>/dev/null | grep -q .; } ;;
  *) say "지원: boltz|protenix|chai"; exit 1;;
esac

say "=== sweep $COF | 예산 ${HOURS}h | rungs $RUNGS | samples $SAMP | out $DATA/$COF ==="
n_run=0; n_skip=0
# process substitution(파이프 아님) → 시간-박스 exit·카운터가 메인 셸에서 동작
while IFS=, read -r target pdb group ab dirtype ag_chains label; do
  [ -z "$target" ] && continue
  if [ "$dirtype" = "diverse" ]; then tdir="$DIVERSE/$pdb"; else tdir="$HERE/targets/$target"; fi
  cj="$tdir/chains.json"
  [ -f "$cj" ] || { say "skip $target (chains.json 없음: $cj)"; continue; }
  IFS='|' read -ra AGC <<< "$ag_chains"
  for r in $(seq 0 $((RUNGS-1))); do
    out="$DATA/$COF/$target/rung$r"
    if DONE "$out"; then n_skip=$((n_skip+1)); continue; fi
    # 항원 사슬별 이 rung의 a3m 매핑
    map=""; ok=1
    for c in "${AGC[@]}"; do
      a3m="$LADDIR/$target/$c/rung$r.a3m"
      [ -f "$a3m" ] || { say "  !! $target $c rung$r a3m 없음($a3m) — gen_msa 먼저"; ok=0; break; }
      map="${map:+$map,}$c=$a3m"
    done
    [ "$ok" = 1 ] || continue
    mkdir -p "$out"
    inp="$out/input.$EXT"
    python "$HERE/make_input.py" --cofolder "$COF" --chains "$cj" --ag-a3m "$map" --dir "$out" --out "$inp" >/dev/null \
      || { say "  !! $target rung$r 입력생성 실패"; continue; }
    say "run $COF $target rung$r ($group-$ab, ag=$ag_chains) ..."
    RUN "$inp" "$out" "$out/run.log" && { say "  OK $target rung$r"; n_run=$((n_run+1)); } \
      || say "  !! $target rung$r 실패 → $out/run.log (tail: $(tail -1 "$out/run.log"))"
    [ "$SMOKE" = 1 ] && { say "SMOKE=1 → 1건 후 종료. 확인: $out/results"; exit 0; }
    if [ "$(date +%s)" -ge "$DEADLINE" ]; then say "예산 소진 → 정지 (run=$n_run skip=$n_skip). 재실행하면 이어감."; exit 0; fi
  done
done < <(tail -n +2 "$LIST")
say "완료 스캔 끝 (run=$n_run skip=$n_skip). 채점: python dockq_sweep.py --models $COF (DockQ env)."
