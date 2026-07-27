#!/usr/bin/env bash
# ══════════════════════════════════════════════════════════════════════════════
# [핵심 실험] 조성 × 반복 — "조성이 성공 확률을 정하는가"를 올바른 단위로 측정.
#
# 왜 필요한가(2026-07-27):
#   이전 설계는 조성마다 **1회씩만** 돌리고 자세 40개를 독립 표본처럼 셌다. 그런데
#   ①한 실행 안 자세들은 서로 거의 같고(상관됨) ②같은 입력을 다시 돌리면 결과가 바뀐다
#   (8ulr seed0: 0.588 → 0.011, 입력 md5 동일). 따라서 유효 표본 = **실행 횟수**이지 자세 수가 아니다.
#
#   → 조성마다 여러 번 돌려 **조성별 성공률**을 낸다. 조성별로 성공률이 갈리면 조성이 변수이고,
#     전부 비슷하면 실행 잡음뿐이다. 어느 쪽이든 분모가 맞는 결론이 된다.
#
# 설계: 조성 8개(기존 seed a3m 재사용) × 반복 4회 = 32회. Protenix ~94초 → 약 50분.
# 출력: seedrep_cand와 같은 구조 → dump_seedrep_full.py로 그대로 채점.
#       폴더명 seed<조성>_r<반복>  (예: seed0_r2 = 조성0의 3번째 반복)
#
# 사용:
#   bash comp_x_reps.sh                                   # 8ulr_HL protenix, 조성 0~7 × 4회
#   COMPS="0 1 2 3" REPS=6 bash comp_x_reps.sh
#   TARGET=9y0a_AB MODEL=boltz DEPTH=d35 bash comp_x_reps.sh
#   SMOKE=1 bash comp_x_reps.sh                           # 1건만
# ══════════════════════════════════════════════════════════════════════════════
set -uo pipefail
cd "$(cd "$(dirname "$0")" && pwd)"
DATA="${DATA:-/mnt/data/admuser/msadepth}"
TARGET="${TARGET:-8ulr_HL}"; MODEL="${MODEL:-protenix}"; DEPTH="${DEPTH:-}"
COMPS="${COMPS:-0 1 2 3 4 5 6 7}"; REPS="${REPS:-4}"
SAMP="${SAMP:-5}"; SEED="${SEED:-0}"; SMOKE="${SMOKE:-0}"
LIST="${LIST:-sweep_targets.csv}"
say(){ echo "[$(date '+%m-%d %H:%M:%S')] $*"; }
export XDG_CACHE_HOME="${XDG_CACHE_HOME:-$DATA/.cache}" TORCH_EXTENSIONS_DIR="${TORCH_EXTENSIONS_DIR:-$DATA/.cache/torch_ext}"
export HF_HOME="${HF_HOME:-$DATA/.cache/hf}" PIP_CACHE_DIR="${PIP_CACHE_DIR:-$DATA/.cache/pip}"
mkdir -p "$XDG_CACHE_HOME" "$TORCH_EXTENSIONS_DIR" "$HF_HOME"
source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null

cj="targets/$TARGET/chains.json"
[ -f "$cj" ] || { say "!! chains.json 없음: $cj"; exit 1; }
ag=$(awk -F, -v t="$TARGET" 'NR>1 && $1==t{print $6; exit}' "$LIST")
[ -n "$ag" ] || { say "!! ag_chains 조회 실패"; exit 1; }
IFS='|' read -ra AGC <<< "$ag"

# 조성 a3m 위치: pipeline/seedrep_cand/<target>_<chain>/d*/seed<N>.a3m
c0="${AGC[0]}"
[ -n "$DEPTH" ] || DEPTH="$(basename "$(ls -d seedrep_cand/${TARGET}_${c0}/d* 2>/dev/null | head -1)")"
[ -n "$DEPTH" ] || { say "!! 조성 a3m 폴더를 못 찾음: seedrep_cand/${TARGET}_${c0}/d*"; exit 1; }
say "타깃 $TARGET · 모델 $MODEL · 깊이 $DEPTH · 조성 [$COMPS] × 반복 $REPS"

case "$MODEL" in
  protenix)
    conda activate "${PROT_ENV:-protenix}" 2>/dev/null
    command -v protenix >/dev/null || { say "!! protenix 없음"; exit 1; }
    export PROTENIX_ROOT_DIR="${PROTENIX_ROOT_DIR:-/mnt/data/admuser/protenix_weights}" LAYERNORM_TYPE=torch
    PROT_MODEL="${PROT_MODEL:-protenix_base_default_v1.0.0}"; EXT="json"
    RUN(){ ( cd "$2" && protenix pred -i "$1" -o results -n "$PROT_MODEL" -s "$SEED" -e "$SAMP" \
             --trimul_kernel torch --triatt_kernel torch --enable_fusion False >"$3" 2>&1 ); }
    DONE(){ find "$1/results" -name '*sample*.cif' 2>/dev/null | grep -q .; } ;;
  boltz)
    conda activate boltz 2>/dev/null
    command -v boltz >/dev/null || { say "!! boltz 없음"; exit 1; }
    export BOLTZ_CACHE="${BOLTZ_CACHE:-$DATA/boltz_cache}"; mkdir -p "$BOLTZ_CACHE"; EXT="yaml"
    RUN(){ ( cd "$2" && boltz predict "$1" --out_dir results --cache "$BOLTZ_CACHE" --no_kernels --diffusion_samples "$SAMP" >"$3" 2>&1 ); }
    DONE(){ find "$1/results" -name '*_model_*.cif' 2>/dev/null | grep -q .; } ;;
  *) say "!! MODEL은 protenix|boltz"; exit 1;;
esac

OUT="$DATA/compreps/seedrep_cand/$MODEL/$TARGET/$DEPTH"
n_ok=0; n_fail=0
for c in $COMPS; do
  map=""; okk=1
  for ch in "${AGC[@]}"; do
    f="seedrep_cand/${TARGET}_${ch}/${DEPTH}/seed${c}.a3m"
    [ -f "$f" ] || { say "  !! 조성$c: $ch 사슬 a3m 없음($f) — 건너뜀"; okk=0; break; }
    map="${map:+$map,}$ch=$(readlink -f "$f")"
  done
  [ "$okk" = 1 ] || continue
  for r in $(seq 0 $((REPS-1))); do
    out="$OUT/seed${c}_r${r}"
    if DONE "$out"; then n_ok=$((n_ok+1)); continue; fi
    mkdir -p "$out"; inp="$out/input.$EXT"
    python make_input.py --cofolder "$MODEL" --chains "$cj" --ag-a3m "$map" --dir "$out" --out "$inp" >"$out/mk.log" 2>&1 \
      || { say "  !! make_input 실패 조성$c 반복$r"; n_fail=$((n_fail+1)); continue; }
    say "조성 $c · 반복 $r 실행 ..."
    if RUN "$inp" "$out" "$out/run.log"; then :; fi
    if DONE "$out"; then say "  OK 조성$c 반복$r"; n_ok=$((n_ok+1))
    else say "  !! 산출물 없음 조성$c 반복$r → $out/run.log"; n_fail=$((n_fail+1)); fi
    [ "$SMOKE" = 1 ] && { say "SMOKE=1 → 종료. 확인: $out/results"; exit 0; }
  done
done
say "완료: 성공 $n_ok · 실패 $n_fail"
say "채점:  python dump_seedrep_full.py --data $DATA/compreps --only $TARGET --csv-out results/compreps_poses.csv"
