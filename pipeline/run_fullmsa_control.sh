#!/usr/bin/env bash
# ══════════════════════════════════════════════════════════════════════════════
# [결정적 통제] full MSA를 같은 자세 예산으로 뽑아보기 — "얕은 깊이가 연 것인가, 그냥 무작위성인가"
#
# 왜 필요한가:
#   seed 복제는 얕은 깊이에서 조성 8가지 × 자세 5개 = 자세 40개를 만들었다.
#   반면 full MSA는 조성이 하나뿐이라 자세 5개만 봤다. 이 상태로는
#   "얕은 깊이가 정답을 열었다"와 "그냥 많이 뽑아서 걸렸다"가 구분되지 않는다.
#   → full MSA에서도 자세 40개를 뽑아, 같은 예산에서 정답이 나오는지 본다.
#
#   full 40개 중 성공 없음  →  얕은 깊이가 실제로 정답을 여는 것 (가설 지지)
#   full 40개 중에도 성공   →  깊이가 아니라 표본 수 문제 (가설 기각)
#
# ⚠️ 조성을 못 바꾸므로(전부 = 1가지) 자세 표본 수로 예산을 맞춘다.
# ⚠️ 재랭커 아님 — MSA깊이 가설(Arm B) 통제.
#
# 사용(tmux 권장):
#   cd ~/projects/bk21-msa-depth-bias/pipeline && git pull
#   SMOKE=1 bash run_fullmsa_control.sh
#   bash run_fullmsa_control.sh 2>&1 | tee fullctl_$(date +%m%d).log
#   # 채점: python dump_seedrep_full.py --data-sub fullmsa_ctl   (또는 score 스크립트)
# ══════════════════════════════════════════════════════════════════════════════
set -uo pipefail
cd "$(cd "$(dirname "$0")" && pwd)"
DATA="${DATA:-/mnt/data/admuser/msadepth}"
LADDIR="$DATA/ladders"
LIST="${LIST:-sweep_targets.csv}"
CAND="${CAND:-seedrep_cand.csv}"
OUTROOT="${OUTROOT:-$DATA/fullmsa_ctl}"
NPOSE="${NPOSE:-40}"          # seed복제(8 seed × 5 자세)와 예산 맞춤
SEED="${SEED:-0}"; SMOKE="${SMOKE:-0}"
say(){ echo "[$(date '+%m-%d %H:%M:%S')] $*"; }
export XDG_CACHE_HOME="${XDG_CACHE_HOME:-$DATA/.cache}" TORCH_EXTENSIONS_DIR="${TORCH_EXTENSIONS_DIR:-$DATA/.cache/torch_ext}"
export HF_HOME="${HF_HOME:-$DATA/.cache/hf}" PIP_CACHE_DIR="${PIP_CACHE_DIR:-$DATA/.cache/pip}"
mkdir -p "$XDG_CACHE_HOME" "$TORCH_EXTENSIONS_DIR" "$HF_HOME"
source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null

setup_model(){
  case "$1" in
    boltz)
      conda activate boltz 2>/dev/null
      command -v boltz >/dev/null || { say "!! boltz 없음"; return 1; }
      export BOLTZ_CACHE="${BOLTZ_CACHE:-$DATA/boltz_cache}"; mkdir -p "$BOLTZ_CACHE"
      EXT="yaml"
      RUN(){ ( cd "$2" && boltz predict "$1" --out_dir results --cache "$BOLTZ_CACHE" --no_kernels --diffusion_samples "$NPOSE" >"$3" 2>&1 ); }
      DONE(){ find "$1/results" -name '*_model_*.cif' 2>/dev/null | grep -q .; } ;;
    protenix)
      conda activate "${PROT_ENV:-protenix}" 2>/dev/null
      command -v protenix >/dev/null || { say "!! protenix 없음"; return 1; }
      export PROTENIX_ROOT_DIR="${PROTENIX_ROOT_DIR:-/mnt/data/admuser/protenix_weights}" LAYERNORM_TYPE=torch
      PROT_MODEL="${PROT_MODEL:-protenix_base_default_v1.0.0}"; EXT="json"
      RUN(){ ( cd "$2" && protenix pred -i "$1" -o results -n "$PROT_MODEL" -s "$SEED" -e "$NPOSE" \
               --trimul_kernel torch --triatt_kernel torch --enable_fusion False >"$3" 2>&1 ); }
      DONE(){ find "$1/results" -name '*sample*.cif' 2>/dev/null | grep -q .; } ;;
    *) return 1;;
  esac
  return 0
}

ag_of(){ awk -F, -v t="$1" 'NR>1 && $1==t{print $6; exit}' "$LIST"; }

while IFS=, read -r target model peak_rung replicas obs_dq obs_rec; do
  [ -z "${target:-}" ] && continue; [ "$target" = "target" ] && continue
  ag=$(ag_of "$target"); [ -n "$ag" ] || { say "skip $target (ag_chains 조회 실패)"; continue; }
  cj="targets/$target/chains.json"; [ -f "$cj" ] || { say "skip $target (chains.json 없음)"; continue; }
  IFS='|' read -ra AGC <<< "$ag"
  setup_model "$model" || { say "skip $target ($model 준비 실패)"; continue; }

  # 각 항원 사슬의 full(rung0) a3m을 그대로 사용 = 조성 고정
  map=""; ok=1
  for c in "${AGC[@]}"; do
    full="$LADDIR/$target/$c/rung0.a3m"
    [ -f "$full" ] || { say "  !! $target/$c rung0.a3m 없음"; ok=0; break; }
    map="${map:+$map,}$c=$full"
  done
  [ "$ok" = 1 ] || continue

  out="$OUTROOT/$model/$target/full_n${NPOSE}"
  if DONE "$out"; then say "이미 있음 skip $target"; continue; fi
  mkdir -p "$out"; inp="$out/input.$EXT"
  python make_input.py --cofolder "$model" --chains "$cj" --ag-a3m "$map" --dir "$out" --out "$inp" >"$out/mk.log" 2>&1 \
    || { say "  !! make_input 실패 $target (로그 $out/mk.log)"; continue; }
  say "run(full MSA, 자세 ${NPOSE}개) $model $target ..."
  if RUN "$inp" "$out" "$out/run.log" && DONE "$out"; then
    n=$(find "$out/results" -name '*.cif' | wc -l); say "  OK $target — 자세 $n개"
  else
    say "  !! 실패 $target → $out/run.log (tail: $(tail -1 "$out/run.log" 2>/dev/null))"
  fi
  [ "$SMOKE" = 1 ] && { say "SMOKE=1 → 1건 후 종료. 확인: $out/results"; exit 0; }
done < "$CAND"
say "완료. 채점하면 'full MSA에서도 40개 중 성공이 나오나'를 볼 수 있다."
