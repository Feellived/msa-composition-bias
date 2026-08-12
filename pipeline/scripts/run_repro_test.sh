#!/usr/bin/env bash
# ══════════════════════════════════════════════════════════════════════════════
# [재현성 측정] 완전히 같은 입력을 N번 돌려 '실행 간' 편차를 잰다.
#
# 왜: 우리가 "조성 안 흔들림 0.012"라고 부른 값은 **한 실행 안의 자세 5개** 편차였다.
#     같은 입력을 두 번 돌렸을 때의 편차는 한 번도 측정한 적이 없다.
#     그런데 2026-07-27 8ulr에서 md5가 같은 MSA로 다른 결과가 나왔다.
#     이 값이 크면 "조성이 답을 정한다"는 결론은 실행 잡음과 구분되지 않는다.
#
# 방법: 기준 실행의 input.json을 그대로 복사해 N번 재실행(같은 시드·같은 MSA·같은 모델).
#       출력 구조를 seedrep_cand와 똑같이 만들어 eval_dump_seedrep.py로 바로 채점.
#
# 사용:
#   bash run_repro_test.sh                       # 8ulr_HL protenix seed0 기준 3회
#   TARGET=9y0a_AB MODEL=boltz SRCSEED=seed0 REPS=3 bash run_repro_test.sh
#   # 채점: python eval_dump_seedrep.py --data $DATA/repro --only 8ulr_HL
# ══════════════════════════════════════════════════════════════════════════════
set -uo pipefail
cd "$(cd "$(dirname "$0")" && pwd)"
DATA="${DATA:-/mnt/data/msadepth}"
TARGET="${TARGET:-8ulr_HL}"; MODEL="${MODEL:-protenix}"; SRCSEED="${SRCSEED:-seed0}"
REPS="${REPS:-3}"; SAMP="${SAMP:-5}"; SEED="${SEED:-0}"
say(){ echo "[$(date '+%m-%d %H:%M:%S')] $*"; }
export XDG_CACHE_HOME="${XDG_CACHE_HOME:-$DATA/.cache}" TORCH_EXTENSIONS_DIR="${TORCH_EXTENSIONS_DIR:-$DATA/.cache/torch_ext}"
export HF_HOME="${HF_HOME:-$DATA/.cache/hf}" PIP_CACHE_DIR="${PIP_CACHE_DIR:-$DATA/.cache/pip}"
source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null

# 기준 실행 찾기(현재 → 없으면 stale)
SRC=""
for base in "$DATA/seedrep_cand" $DATA/seedrep_cand_stale_*; do
  c=$(ls -d "$base/$MODEL/$TARGET"/d*/"$SRCSEED" 2>/dev/null | head -1)
  [ -n "$c" ] && [ -f "$c/input.json" -o -f "$c/input.yaml" ] && { SRC="$c"; break; }
done
[ -n "$SRC" ] || { say "!! 기준 실행을 못 찾음: $MODEL/$TARGET/*/$SRCSEED"; exit 1; }
DEPTH="$(basename "$(dirname "$SRC")")"
say "기준 = $SRC  (깊이 $DEPTH)"

case "$MODEL" in
  protenix)
    conda activate "${PROT_ENV:-protenix}" 2>/dev/null
    export PROTENIX_ROOT_DIR="${PROTENIX_ROOT_DIR:-/mnt/data/protenix_weights}" LAYERNORM_TYPE=torch
    PROT_MODEL="${PROT_MODEL:-protenix_base_default_v1.0.0}"; INP="input.json"
    RUN(){ ( cd "$1" && protenix pred -i "$INP" -o results -n "$PROT_MODEL" -s "$SEED" -e "$SAMP" \
             --trimul_kernel torch --triatt_kernel torch --enable_fusion False >run.log 2>&1 ); }
    DONE(){ find "$1/results" -name '*sample*.cif' 2>/dev/null | grep -q .; } ;;
  boltz)
    conda activate boltz 2>/dev/null
    export BOLTZ_CACHE="${BOLTZ_CACHE:-$DATA/boltz_cache}"; INP="input.yaml"
    RUN(){ ( cd "$1" && boltz predict "$INP" --out_dir results --cache "$BOLTZ_CACHE" --no_kernels --diffusion_samples "$SAMP" >run.log 2>&1 ); }
    DONE(){ find "$1/results" -name '*_model_*.cif' 2>/dev/null | grep -q .; } ;;
  *) say "!! MODEL은 protenix|boltz"; exit 1;;
esac
[ -f "$SRC/$INP" ] || { say "!! $SRC/$INP 없음"; exit 1; }

# seedrep_cand와 동일한 구조로 출력 → 기존 채점 도구 재사용
OUT="$DATA/repro/seedrep_cand/$MODEL/$TARGET/$DEPTH"
mkdir -p "$OUT"
for r in $(seq 0 $((REPS-1))); do
  d="$OUT/seed$r"
  if DONE "$d"; then say "이미 있음 skip 반복$r"; continue; fi
  mkdir -p "$d"
  cp "$SRC/$INP" "$d/$INP"                      # a3m은 절대경로 참조 → 기준 실행 것을 그대로 사용
  say "반복 $r 실행 ..."
  if RUN "$d" && DONE "$d"; then say "  OK 반복 $r"
  else say "  !! 실패 반복 $r → $d/run.log (tail: $(tail -1 "$d/run.log" 2>/dev/null))"; fi
done
say "완료. 채점:  python eval_dump_seedrep.py --data $DATA/repro --only $TARGET"
say "  → 'seed0/1/2'는 반복 실행이다(조성 아님). 값이 같으면 결정적, 갈리면 실행 잡음 존재."
