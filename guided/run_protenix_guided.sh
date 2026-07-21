#!/usr/bin/env bash
# HADDOCK 패치 → Protenix pocket-guided 재실행 (guided-Boltz의 Protenix판, §실험 4).
# 커맨드는 run_protenix_msa_depth.sh와 동일(protenix pred -i -o -n -s -e). self-healing.
# 사용: PROT_ENV=protenix OUTDIR=runs_guided_prot GEN_ARGS="--topn 20" \
#         tmux new -s protguided; bash scripts/run_protenix_guided.sh "8XSI 9SBB 8SIS 9ML8 8SDF 8SIT 9ML9" 5
#   ⚠️ 스모크 먼저: SMOKE=1 bash scripts/run_protenix_guided.sh 8XSI
set -uo pipefail
cd "$(dirname "$0")/.." || exit 1
source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null
TARGETS="${1:-8XSI 9SBB 8SIS 9ML8 8SDF 8SIT 9ML9}"; SAMP="${2:-5}"
SMOKE="${SMOKE:-0}"; OUTDIR="${OUTDIR:-runs_guided_prot}"; GEN_ARGS="${GEN_ARGS:-}"
PROT_ENV="${PROT_ENV:-protenix}"; PROT_MODEL="${PROT_MODEL:-protenix_base_20250630_v1.0.0}"; SEED="${SEED:-101}"
say(){ echo "[$(date '+%m-%d %H:%M:%S')] $*"; }
conda activate "$PROT_ENV" 2>/dev/null
export PROTENIX_ROOT_DIR=/mnt/data/admuser/protenix_weights LAYERNORM_TYPE=torch
command -v protenix >/dev/null || { say "!! protenix 없음 (conda activate $PROT_ENV?)"; exit 1; }

# 1) HADDOCK 패치 → guided JSON 생성
MP=8; [ "$SMOKE" = "1" ] && { MP=1; SAMP=1; }
python scripts/haddock_to_protenix_pocket.py --targets "$TARGETS" --max-patches "$MP" --outdir "$OUTDIR" $GEN_ARGS \
  || { say "!! JSON 생성 실패"; exit 1; }

# 2) 각 patch JSON 실행
for T in $TARGETS; do
  t=$(echo "$T" | tr 'A-Z' 'a-z'); d="$OUTDIR/$t"
  [ -d "$d" ] || { say "skip $T (JSON 없음)"; continue; }
  for j in "$d"/protenix_${t}_patch*.json; do
    [ -e "$j" ] || continue
    pk=$(basename "$j" .json | sed 's/.*_//')      # patchN
    od="results_$pk"
    if find "$d/$od" -name '*sample*.cif' 2>/dev/null | grep -q .; then say "skip $T $pk (있음)"; continue; fi
    say "protenix guided $T $pk (samples=$SAMP) ..."
    ( cd "$d" && protenix pred -i "$(basename "$j")" -o "$od" -n "$PROT_MODEL" -s "$SEED" -e "$SAMP" \
        --trimul_kernel torch --triatt_kernel torch > "protenix_$pk.log" 2>&1 ) \
      && say "  OK $T $pk" || { say "  !! $T $pk 실패 → $d/protenix_$pk.log (tail:)"; tail -5 "$d/protenix_$pk.log"; }
    [ "$SMOKE" = "1" ] && { say "SMOKE=1 → 1개만 실행 후 종료. 확인: find $d/$od -name '*sample*.cif'"; exit 0; }
  done
done
say "완료 → $OUTDIR/<t>/results_patchN/. 분석: python scripts/dockq_protenix_guided.py --guided-dir $OUTDIR --targets \"$TARGETS\""
