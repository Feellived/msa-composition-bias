#!/usr/bin/env bash
# HADDOCK 패치 → tFold-Ag ICF-guided 재실행 (guided-Boltz/Protenix의 tFold판, §실험 4).
# ICF 손구성(--icf .pt, --model_version ppi = epitope-only). self-healing.
# 사용: TFOLD_ENV=tfold OUTDIR=runs_guided_tfold GEN_ARGS="--topn 20" \
#         tmux new -s tfguided; bash scripts/run_tfold_guided.sh "8XSI 9SBB 8SIS 9ML8 8SDF 8SIT 9ML9"
#   ⚠️ 스모크 먼저: SMOKE=1 bash scripts/run_tfold_guided.sh 8XSI
#   ⚠️ tFold JSON 모드 icf_path 버그 회피 = FASTA+--icf 경로. predict.py에 --icf 있는지 스모크로 확인.
set -uo pipefail
cd "$(dirname "$0")/.." || exit 1
ROOT="$(pwd)"
source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null
TARGETS="${1:-8XSI 9SBB 8SIS 9ML8 8SDF 8SIT 9ML9}"
SMOKE="${SMOKE:-0}"; OUTDIR="${OUTDIR:-runs_guided_tfold}"; GEN_ARGS="${GEN_ARGS:-}"
TFOLD_ENV="${TFOLD_ENV:-tfold}"; TFOLD_DIR="${TFOLD_DIR:-$ROOT/tfold}"; SRC="${SRC:-runs_rbd}"
SEEDS="${SEEDS:-42 43 44 45 46}"; MODEL_VERSION="${MODEL_VERSION:-ppi}"
say(){ echo "[$(date '+%m-%d %H:%M:%S')] $*"; }
conda activate "$TFOLD_ENV" 2>/dev/null
python -c "import tfold, torch" 2>/dev/null || say "!! tfold/torch import 실패(env $TFOLD_ENV) — 계속 시도"

# 1) HADDOCK 패치 → ICF .pt (torch+biopython 필요)
MP=8; [ "$SMOKE" = "1" ] && { MP=1; SEEDS="42"; }
python scripts/haddock_to_tfold_icf.py --targets "$TARGETS" --src "$SRC" --max-patches "$MP" --outdir "$OUTDIR" $GEN_ARGS \
  || { say "!! ICF 생성 실패"; exit 1; }

# 2) 각 patch ICF 실행 (predict.py --fasta --msa --icf --model_version ppi)
for T in $TARGETS; do
  t=$(echo "$T"|tr 'A-Z' 'a-z'); d="$OUTDIR/$t"
  [ -d "$d" ] || { say "skip $T (ICF 없음)"; continue; }
  fasta="$ROOT/$d/tfold_${t}.fasta"
  a3m="${AGA3M:-$ROOT/$SRC/$t/msa_$t/A.a3m}"        # 항원 full a3m (depth-sweep 아님)
  [ -f "$a3m" ] || { say "skip $T (항원 a3m 없음: $a3m)"; continue; }
  for icf in "$ROOT/$d"/tfold_${t}_patch*.pt; do     # ⚠️ 절대경로(predict를 tfold/에서 cd해 돌리므로)
    [ -e "$icf" ] || continue
    pk=$(basename "$icf" .pt | sed 's/.*_//')        # patchN
    od="$ROOT/$d/out_tfold_$pk"; mkdir -p "$od"
    for s in $SEEDS; do
      outpdb="$od/${t}_seed${s}.pdb"
      [ -s "$outpdb" ] && { say "skip $T $pk seed$s (있음)"; continue; }
      say "tfold-guided $T $pk seed$s ..."
      ( cd "$TFOLD_DIR" && python projects/tfold_ag/predict.py \
          --fasta "$fasta" --msa "$a3m" --icf "$icf" --output "$outpdb" \
          --model_version "$MODEL_VERSION" --seed "$s" > "$od/${t}_seed${s}.log" 2>&1 ) \
        && say "  OK $T $pk seed$s" || { say "  !! $T $pk seed$s 실패 → $od/${t}_seed${s}.log (tail:)"; tail -3 "$od/${t}_seed${s}.log"; }
      [ "$SMOKE" = "1" ] && { say "SMOKE=1 → 1건 후 종료. 확인: $outpdb (성공하면 --icf 경로 OK)"; exit 0; }
    done
  done
done
say "완료 → $OUTDIR/<t>/out_tfold_patchN/. 분석: python scripts/dockq_tfold_guided.py --guided-dir $OUTDIR --targets \"$TARGETS\""
