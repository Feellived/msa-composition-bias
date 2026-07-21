#!/usr/bin/env bash
# HADDOCK 클러스터 patch → Boltz pocket-guided 재실행 (staged refine 실험, §실험 4).
# blind baseline = runs_msad_143. self-healing(이미 있는 patch skip).
# 사용: bash scripts/run_boltz_guided.sh "8XSI 9SBB 8SIS 9ML8 8SDF 8SIT 9ML9" [samples=5]
#   ⚠️ 먼저 스모크(1타깃 1패치 1샘플)로 pocket YAML 포맷 검증:
#      SMOKE=1 bash scripts/run_boltz_guided.sh 8XSI
set -uo pipefail
cd "$(dirname "$0")/.." || exit 1
source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null
conda activate boltz 2>/dev/null
TARGETS="${1:-8XSI 9SBB 8SIS 9ML8 8SDF 8SIT 9ML9}"; SAMP="${2:-5}"
SMOKE="${SMOKE:-0}"
OUTDIR="${OUTDIR:-runs_guided}"          # 패치 YAML·결과 폴더(변형 실험이면 딴 이름)
GEN_ARGS="${GEN_ARGS:-}"                  # haddock_to_boltz_pocket.py 추가옵션(예: "--topn 20 --union --soft")
say(){ echo "[$(date '+%m-%d %H:%M:%S')] $*"; }
command -v boltz >/dev/null || { say "!! boltz 없음 (conda activate boltz?)"; exit 1; }

# 1) HADDOCK 클러스터 → guided YAML 생성
MP=8; [ "$SMOKE" = "1" ] && { MP=1; SAMP=1; }
python scripts/haddock_to_boltz_pocket.py --targets "$TARGETS" --max-patches "$MP" --outdir "$OUTDIR" $GEN_ARGS || { say "!! YAML 생성 실패"; exit 1; }

# 2) 각 patch YAML 실행
for T in $TARGETS; do
  t=$(echo "$T" | tr 'A-Z' 'a-z'); d="$OUTDIR/$t"
  [ -d "$d" ] || { say "skip $T (YAML 없음)"; continue; }
  for y in "$d"/boltz_${t}_patch*.yaml; do
    [ -e "$y" ] || continue
    pk=$(basename "$y" .yaml | sed 's/.*_//')      # patchN
    od="results_$pk"
    if find "$d/$od" -name '*_model_*.cif' 2>/dev/null | grep -q .; then say "skip $T $pk (있음)"; continue; fi
    say "boltz guided $T $pk (samples=$SAMP) ..."
    ( cd "$d" && boltz predict "$(basename "$y")" --out_dir "$od" --no_kernels --diffusion_samples "$SAMP" > "boltz_$pk.log" 2>&1 ) \
      && say "  OK $T $pk" || { say "  !! $T $pk 실패 → $d/boltz_$pk.log (tail:)"; tail -5 "$d/boltz_$pk.log"; }
    [ "$SMOKE" = "1" ] && { say "SMOKE=1 → 1개만 실행하고 종료. 포즈 확인: find $d/$od -name '*_model_*.cif'"; exit 0; }
  done
done
say "완료 → runs_guided/<t>/results_patchN/. 분석: python scripts/dockq_boltz_guided.py --targets \"$TARGETS\""
