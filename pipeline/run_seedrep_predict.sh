#!/usr/bin/env bash
# Track A Stage 0 예측 — seedrep a3m을 co-folder(Protenix/Chai)로 예측.
# run_sweep.sh의 '검증된' env·RUN·DONE을 그대로 복사, 반복만 rung → seedrep(depth×seed)로 교체.
#   bash run_seedrep_predict.sh protenix
#   SMOKE=1 bash run_seedrep_predict.sh protenix   # 1건 후 성공(0)/실패(1)로 종료 = 오케스트레이터 게이트용
#   ANCHORS="8wpy_AB" bash run_seedrep_predict.sh chai
# 입력 = pipeline/seedrep/<t>_<c>/d<depth>/seed<s>.a3m (run_track_a_seedrep.sh 산출)
# 출력 = $DATA/seedrep_pred/<COF>/<t>/d<depth>/seed<s>/results/*.cif
set -uo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"; cd "$HERE"
COF="${1:?protenix|chai}"
DATA="${DATA:-/mnt/data/admuser/msadepth}"
SEEDREP="${SEEDREP:-$HERE/seedrep}"
OUTROOT="${OUTROOT:-$DATA/seedrep_pred}"
LIST="${LIST:-$HERE/sweep_targets.csv}"
DIVERSE="${DIVERSE:-$HERE/../runs_diverse}"
SAMP="${SAMP:-5}"; SEED="${SEED:-0}"; SMOKE="${SMOKE:-0}"
ANCHORS="${ANCHORS:-8wpy_AB 8k3k_D 8k46_I 9y0a_AB}"
say(){ echo "[$(date '+%m-%d %H:%M:%S')] $*"; }
export XDG_CACHE_HOME="${XDG_CACHE_HOME:-$DATA/.cache}" TORCH_EXTENSIONS_DIR="${TORCH_EXTENSIONS_DIR:-$DATA/.cache/torch_ext}"
export HF_HOME="${HF_HOME:-$DATA/.cache/hf}" PIP_CACHE_DIR="${PIP_CACHE_DIR:-$DATA/.cache/pip}"
mkdir -p "$XDG_CACHE_HOME" "$TORCH_EXTENSIONS_DIR" "$HF_HOME"
source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null

case "$COF" in
  protenix)
    conda activate "${PROT_ENV:-protenix}" 2>/dev/null
    export PROTENIX_ROOT_DIR="${PROTENIX_ROOT_DIR:-/mnt/data/admuser/protenix_weights}" LAYERNORM_TYPE=torch
    command -v protenix >/dev/null || { say "!! protenix 없음"; exit 1; }
    PROT_MODEL="${PROT_MODEL:-protenix_base_default_v1.0.0}"; EXT="json"   # leakage-free(컷 2021-09-30)
    RUN(){ ( cd "$2" && protenix pred -i "$1" -o results -n "$PROT_MODEL" -s "$SEED" -e "$SAMP" \
             --trimul_kernel torch --triatt_kernel torch --enable_fusion False >"$3" 2>&1 ); }
    DONE(){ find "$1/results" -name '*sample*.cif' 2>/dev/null | grep -q .; } ;;
  chai)
    conda activate "${CHAI_ENV:-chai}" 2>/dev/null
    command -v chai-lab >/dev/null || { say "!! chai-lab 없음 (conda activate ${CHAI_ENV:-chai})"; exit 1; }
    export CHAI_DOWNLOADS_DIR="${CHAI_DOWNLOADS_DIR:-$DATA/chai_downloads}"; mkdir -p "$CHAI_DOWNLOADS_DIR"; EXT="fasta"
    RUN(){ rm -rf "$2/results"; ( cd "$2" && chai-lab fold --msa-directory "$2/msa" "$1" results >"$3" 2>&1 ); }
    DONE(){ find "$1/results" -name '*.cif' 2>/dev/null | grep -q .; } ;;
  *) say "지원: protenix|chai"; exit 1;;
esac

say "=== seedrep predict $COF | samples $SAMP | out $OUTROOT/$COF (SMOKE=$SMOKE) ==="
n_run=0; n_ok=0; n_fail=0
while IFS=, read -r target pdb group ab dirtype ag_chains label; do
  [ -z "$target" ] && continue
  case " $ANCHORS " in *" $target "*) : ;; *) continue;; esac
  if [ "$dirtype" = "diverse" ]; then tdir="$DIVERSE/$pdb"; else tdir="$HERE/targets/$target"; fi
  cj="$tdir/chains.json"; [ -f "$cj" ] || { say "skip $target (chains.json 없음: $cj)"; continue; }
  IFS='|' read -ra AGC <<< "$ag_chains"
  c0="${AGC[0]}"; base="$SEEDREP/${target}_${c0}"
  [ -d "$base" ] || { say "skip $target (seedrep 없음: $base — seed_replicate 먼저)"; continue; }
  for ddir in "$base"/d*/; do
    [ -d "$ddir" ] || continue; depth="$(basename "$ddir")"
    for a3m0 in "$ddir"seed*.a3m; do
      [ -e "$a3m0" ] || continue; s="$(basename "$a3m0" .a3m)"
      map=""; ok=1
      for c in "${AGC[@]}"; do
        f="$SEEDREP/${target}_${c}/${depth}/${s}.a3m"
        [ -f "$f" ] || { ok=0; break; }; map="${map:+$map,}$c=$f"
      done
      [ "$ok" = 1 ] || { say "  skip $target $depth $s (항원사슬 a3m 불완전)"; continue; }
      out="$OUTROOT/$COF/$target/$depth/$s"
      if DONE "$out"; then n_ok=$((n_ok+1)); [ "$SMOKE" = 1 ] && { say "이미 됨 → smoke OK"; exit 0; }; continue; fi
      mkdir -p "$out"; inp="$out/input.$EXT"
      python "$HERE/make_input.py" --cofolder "$COF" --chains "$cj" --ag-a3m "$map" --dir "$out" --out "$inp" >"$out/makeinput.log" 2>&1 \
        || { say "  !! make_input 실패 $target $depth $s (로그 $out/makeinput.log)"; n_fail=$((n_fail+1)); [ "$SMOKE" = 1 ] && exit 1; continue; }
      say "run $COF $target $depth $s (ag=$ag_chains) ..."
      if RUN "$inp" "$out" "$out/run.log" && DONE "$out"; then say "  OK $target $depth $s"; n_ok=$((n_ok+1)); rc=0
      else say "  !! 실패 $target $depth $s → $out/run.log (tail: $(tail -1 "$out/run.log" 2>/dev/null))"; n_fail=$((n_fail+1)); rc=1; fi
      n_run=$((n_run+1))
      [ "$SMOKE" = 1 ] && { say "SMOKE 종료 rc=$rc"; exit $rc; }
    done
  done
done < <(tail -n +2 "$LIST")
say "완료 $COF: run=$n_run ok=$n_ok fail=$n_fail → $OUTROOT/$COF/"
