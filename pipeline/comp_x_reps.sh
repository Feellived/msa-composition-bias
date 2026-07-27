#!/usr/bin/env bash
# ══════════════════════════════════════════════════════════════════════════════
# [핵심 실험] 같은 조건을 여러 번 돌려 **성공률**을 잰다.
#
# 왜: 이 모델들은 확산(diffusion) 기반이라 같은 입력에도 매번 답이 조금씩 다르다.
#     (8ulr에서 입력 md5가 같은데 0.588 → 0.011 관측) 그래서 한 번 돌린 값으로 판단할 수 없고,
#     **여러 번 돌려 몇 번 성공하는지**를 세야 한다. 유효 표본 = 실행 횟수(자세 수가 아님:
#     한 실행 안의 자세 5개는 서로 거의 같다).
#
# 두 조건을 같은 방식으로 재서 비교한다:
#   · 원래 MSA      : COMPS="full"          (rung0을 그대로 사용)
#   · 얕은 깊이 MSA : COMPS="0 1 2 3 ..."   (같은 개수로 다르게 뽑은 목록들)
#
# 출력은 seedrep_cand와 같은 구조 → dump_seedrep_full.py로 그대로 채점.
#   폴더명 seed<조성>_r<반복>   (예: seed2_r3, seedfull_r7)
#
# 사용:
#   # 8ulr 본실험(조성 8개 × 4회)
#   bash comp_x_reps.sh
#   # 원래 MSA 10회(대조)
#   COMPS="full" REPS=10 bash comp_x_reps.sh
#   # 다른 복합체 예비검정(조성 5개 1회씩 + 원래 5회) — RUNG 주면 조성 a3m을 자동 생성
#   RUNG=3 TARGET=9azr_HL COMPS="0 1 2 3 4" REPS=1 bash comp_x_reps.sh
#   RUNG=3 TARGET=9azr_HL COMPS="full"      REPS=5 bash comp_x_reps.sh
# ══════════════════════════════════════════════════════════════════════════════
set -uo pipefail
cd "$(cd "$(dirname "$0")" && pwd)"
DATA="${DATA:-/mnt/data/admuser/msadepth}"
TARGET="${TARGET:-8ulr_HL}"; MODEL="${MODEL:-protenix}"
RUNG="${RUNG:-}"                 # 지정 시 이 rung 깊이로 조성 a3m 자동 생성(없을 때만)
REPLICAS="${REPLICAS:-8}"        # 자동 생성할 조성 개수
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
[ -n "$ag" ] || { say "!! ag_chains 조회 실패: $TARGET"; exit 1; }
IFS='|' read -ra AGC <<< "$ag"

# ── 사슬별 깊이 결정 (+ 필요하면 조성 a3m 생성). 사슬마다 서열 수가 달라 개별 관리 ──
DEPTHS=()
for i in "${!AGC[@]}"; do
  ch="${AGC[$i]}"
  if [ -n "$RUNG" ]; then
    f0="$DATA/ladders/$TARGET/$ch/rung0.a3m"; fr="$DATA/ladders/$TARGET/$ch/rung${RUNG}.a3m"
    [ -f "$f0" ] && [ -f "$fr" ] || { say "!! $TARGET/$ch rung0/rung$RUNG a3m 없음"; exit 1; }
    d=$(grep -c '^>' "$fr"); DEPTHS[$i]="d$d"
    if ! ls "seedrep_cand/${TARGET}_${ch}/d${d}"/seed[0-9]*.a3m >/dev/null 2>&1; then
      say "조성 생성: $TARGET/$ch — 깊이 ${d}서열 × ${REPLICAS}가지"
      python seed_replicate.py --a3m "$f0" --depths "$d" --replicas "$REPLICAS" \
             --outdir "seedrep_cand/${TARGET}_${ch}" >/dev/null || { say "!! seed_replicate 실패"; exit 1; }
    fi
  else
    DEPTHS[$i]="$(basename "$(ls -d "seedrep_cand/${TARGET}_${ch}"/d* 2>/dev/null | head -1)")"
    [ -n "${DEPTHS[$i]}" ] || { say "!! 조성 폴더 없음. RUNG=<칸번호>를 지정하면 자동 생성됨"; exit 1; }
  fi
  # 'full'을 하나의 조성처럼 다루기 위해 rung0을 seedfull.a3m으로 준비
  case " $COMPS " in
    *" full "*) cp -n "$DATA/ladders/$TARGET/$ch/rung0.a3m" \
                      "seedrep_cand/${TARGET}_${ch}/${DEPTHS[$i]}/seedfull.a3m" 2>/dev/null ;;
  esac
done
DEPTH="${DEPTHS[0]}"
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
  for i in "${!AGC[@]}"; do
    ch="${AGC[$i]}"; f="seedrep_cand/${TARGET}_${ch}/${DEPTHS[$i]}/seed${c}.a3m"
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
    RUN "$inp" "$out" "$out/run.log"
    if DONE "$out"; then say "  OK 조성$c 반복$r"; n_ok=$((n_ok+1))
    else say "  !! 산출물 없음 조성$c 반복$r → $out/run.log"; n_fail=$((n_fail+1)); fi
    [ "$SMOKE" = 1 ] && { say "SMOKE=1 → 종료. 확인: $out/results"; exit 0; }
  done
done
say "완료: 성공 $n_ok · 실패 $n_fail"
say "채점: python dump_seedrep_full.py --data $DATA/compreps --only $TARGET --csv-out results/compreps_${TARGET}.csv"
