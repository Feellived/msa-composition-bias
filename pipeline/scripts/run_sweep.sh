#!/usr/bin/env bash
# 시간-박스 재개형 depth-sweep dispatcher (co-folder 하나씩).
# 복합체 단위로 rung 전부 → 다음 복합체(round-robin 순서=중간에 끊겨도 그룹 균형).
# 예산 소진 시 깔끔 정지, 재실행 시 self-heal(출력 있으면 skip).
#   bash run_sweep.sh boltz 11     # Boltz, 11시간(반나절/밤샘)
#   bash run_sweep.sh protenix 54  # Protenix, 주말
#   bash run_sweep.sh chai 11      # Chai-1 (CHAI_ENV=chai)
#   SMOKE=1 bash run_sweep.sh boltz 1   # 1건만
set -uo pipefail
# 스크립트는 pipeline/scripts/ 에 있고, 기준 디렉토리는 상위 pipeline/ 이다
HERE="$(cd "$(dirname "$0")/.." && pwd)"          
cd "$HERE" || exit 1
COF="${1:?boltz|protenix|chai}"; HOURS="${2:-11}"
DATA="${DATA:-/mnt/data/msadepth}"      # ⚠️ 대용량 출력 = /mnt/data(홈 금지)
LADDIR="$DATA/ladders"; RUNGS="${RUNGS:-12}"; SAMP="${SAMP:-5}"; SEED="${SEED:-0}"; SMOKE="${SMOKE:-0}"   # 12단(깊이 촘촘)
DIVERSE="${DIVERSE:-$HERE/../runs_diverse}"      # C 재활용 위치(서버 확인)
LIST="${LIST:-$HERE/sweep_targets.csv}"
say(){ echo "[$(date '+%m-%d %H:%M:%S')] $*"; }
DEADLINE=$(( $(date +%s) + HOURS*3600 ))
# ⚠️ 홈(/) 보호: 기본이 홈(~/.cache, ~/.boltz)인 캐시들을 전부 /mnt/data로 — 홈 디스크 full 재발 방지
export XDG_CACHE_HOME="${XDG_CACHE_HOME:-$DATA/.cache}"          # ~/.cache 전반(torch_extensions 포함 다수)
export TORCH_EXTENSIONS_DIR="${TORCH_EXTENSIONS_DIR:-$DATA/.cache/torch_ext}"
export HF_HOME="${HF_HOME:-$DATA/.cache/hf}" PIP_CACHE_DIR="${PIP_CACHE_DIR:-$DATA/.cache/pip}"
mkdir -p "$XDG_CACHE_HOME" "$TORCH_EXTENSIONS_DIR" "$HF_HOME"
source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null

case "$COF" in
  boltz)
    conda activate boltz 2>/dev/null
    command -v boltz >/dev/null || { say "!! boltz 없음"; exit 1; }
    export BOLTZ_CACHE="${BOLTZ_CACHE:-$DATA/boltz_cache}"; mkdir -p "$BOLTZ_CACHE"   # 가중치·CCD: 홈(~/.boltz) 대신 /mnt/data
    MIN_MSA="${MIN_MSA:-2}"   # ⚠️ Boltz만: 순수 single-seq(1행)서 데이터로더 폭주(수십GB·stall) → <2서열 rung skip
    EXT="yaml"
    RUN(){ ( cd "$2" && boltz predict "$1" --out_dir results --cache "$BOLTZ_CACHE" --no_kernels --diffusion_samples "$SAMP" >"$3" 2>&1 ); }
    DONE(){ find "$1/results" -name '*_model_*.cif' 2>/dev/null | grep -q .; } ;;
  protenix)
    conda activate "${PROT_ENV:-protenix}" 2>/dev/null
    export PROTENIX_ROOT_DIR="${PROTENIX_ROOT_DIR:-/mnt/data/protenix_weights}" LAYERNORM_TYPE=torch
    command -v protenix >/dev/null || { say "!! protenix 없음"; exit 1; }
    # ⚠️ leakage: protenix_base_20250630_v1.0.0 = 학습컷오프 2025-06-30(post-2023-06 테스트셋 암기) → 금지.
    # 기본=protenix_base_default_v1.0.0(컷오프 2021-09-30=leakage-free, 공식 권장·AF3 능가).
    # protenix-v2(더 강함, 동일 컷오프)는 공개 다운로드 403(미공개)이라 보류 — 받히면 PROT_MODEL=protenix-v2로 교체.
    PROT_MODEL="${PROT_MODEL:-protenix_base_default_v1.0.0}"
    MIN_MSA="${MIN_MSA:-1}"   # single-seq도 실행(공진화 0 극단점 = 편향 완전제거)
    EXT="json"
    # 커널 JIT 회피 3중: LAYERNORM_TYPE=torch(env) + 삼각커널 torch + fusion off → 무인 실행 시 컴파일 실패-루프 방지
    RUN(){ ( cd "$2" && protenix pred -i "$1" -o results -n "$PROT_MODEL" -s "$SEED" -e "$SAMP" \
             --trimul_kernel torch --triatt_kernel torch --enable_fusion False >"$3" 2>&1 ); }
    DONE(){ find "$1/results" -name '*sample*.cif' 2>/dev/null | grep -q .; } ;;
  chai)
    conda activate "${CHAI_ENV:-chai}" 2>/dev/null
    command -v chai-lab >/dev/null || { say "!! chai-lab 없음 (conda activate ${CHAI_ENV:-chai})"; exit 1; }
    export CHAI_DOWNLOADS_DIR="${CHAI_DOWNLOADS_DIR:-$DATA/chai_downloads}"; mkdir -p "$CHAI_DOWNLOADS_DIR"   # 가중치: 홈 대신 /mnt/data
    MIN_MSA="${MIN_MSA:-1}"   # ESM-2 PLM 내장 → single-seq(MSA-free)가 설계상 유효, 실행
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
    # 최소 서열 수(항원 사슬 중 최소) < MIN_MSA → 이 모델은 skip. Boltz(MIN_MSA=2)만 single-seq rung 건너뜀; Chai/Protenix(1)는 실행.
    minseq=1000000
    for c in "${AGC[@]}"; do n=$(grep -c '^>' "$LADDIR/$target/$c/rung$r.a3m" 2>/dev/null || echo 0); [ "$n" -lt "$minseq" ] && minseq=$n; done
    if [ "$minseq" -lt "$MIN_MSA" ]; then say "  skip $target rung$r ($COF: ${minseq}서열 < MIN_MSA=$MIN_MSA)"; n_skip=$((n_skip+1)); continue; fi
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
say "완료 스캔 끝 (run=$n_run skip=$n_skip). 채점: python eval_dockq_sweep.py --models $COF (DockQ env)."
