#!/usr/bin/env bash
# ══════════════════════════════════════════════════════════════════════════════
# [D] 가짜 쌍 예측 러너 — prep_decoy_pairs.py 가 만든 targets/decoy_* 를 protenix 로 돌린다.
#
# 정답이 없는 쌍이라 MSA 를 새로 만들 필요가 없다. 가짜 쌍의 항원은 원래 타깃과 완전히
# 같으므로(잔기까지 동일), 그 항원의 **기존 Neff80 사다리(rung)를 조성 대신 그대로
# 재사용**한다. 조성을 재추첨하는 run_seed_replicate.py 단계를 건너뛰어 시간이 크게 준다.
#
# rung 3개(낮음·중간·높음)를 '조성 3가지'로 삼아 각각 1회 예측한다 — 정답이 없어
# DockQ/결합자리 채점을 안 하므로 반복(seed) 대신 rung 다양성으로 변주를 준다.
# 채점은 예측이 끝난 뒤 analyze_collect_iptm.py 로 ipTM 만 모은다(정답 불필요).
#
# ⚠️ GPU 를 오래 쓴다 — 서버에 영향이 갈 수 있는 작업이라 CLAUDE.md 규칙대로
#    돌리기 전에 반드시 사용자 확인을 거쳤다는 전제.
#
# 사용(tmux 권장):
#   cd ~/projects/bk21-msa-depth-bias/pipeline && git pull
#   SMOKE=1 bash run_decoy_predict.sh          # 가짜쌍 1개 × rung 1개만 스모크
#   bash run_decoy_predict.sh 2>&1 | tee logs/decoy_predict_$(date +%m%d).log
#   완료 후: python analyze_collect_iptm.py --data "$DATA/decoy/protenix" --out results/decoy_iptm.csv
# ══════════════════════════════════════════════════════════════════════════════
set -uo pipefail
cd "$(cd "$(dirname "$0")" && pwd)"
DATA="${DATA:-/mnt/data/admuser/msadepth}"
LADDIR="$DATA/ladders"
MANIFEST="${MANIFEST:-results/decoy_manifest.csv}"
OUTROOT="${OUTROOT:-$DATA/decoy/protenix}"    # ⚠️ analyze_collect_iptm.py 가 경로에서 리터럴 "/protenix/<타깃>/<깊이>/<실행>/results/" 를 찾는다
RUNGS="${RUNGS:-0 4 8}"          # 항원 사다리에서 재사용할 rung 번호 (낮음·중간·높음)
SAMP="${SAMP:-5}"; SEED="${SEED:-0}"; SMOKE="${SMOKE:-0}"
say(){ echo "[$(date '+%m-%d %H:%M:%S')] $*"; }

[ -f "$MANIFEST" ] || { say "!! 없음: $MANIFEST — prep_decoy_pairs.py --write 를 먼저 실행할 것"; exit 1; }

# run_sweep.sh 와 동일한 캐시 보호 + protenix 실행 관례
export XDG_CACHE_HOME="${XDG_CACHE_HOME:-$DATA/.cache}" TORCH_EXTENSIONS_DIR="${TORCH_EXTENSIONS_DIR:-$DATA/.cache/torch_ext}"
export HF_HOME="${HF_HOME:-$DATA/.cache/hf}" PIP_CACHE_DIR="${PIP_CACHE_DIR:-$DATA/.cache/pip}"
mkdir -p "$XDG_CACHE_HOME" "$TORCH_EXTENSIONS_DIR" "$HF_HOME"
source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null
conda activate "${PROT_ENV:-protenix}" 2>/dev/null
export PROTENIX_ROOT_DIR="${PROTENIX_ROOT_DIR:-/mnt/data/admuser/protenix_weights}" LAYERNORM_TYPE=torch
command -v protenix >/dev/null || { say "!! protenix 없음 (conda activate ${PROT_ENV:-protenix})"; exit 1; }
PROT_MODEL="${PROT_MODEL:-protenix_base_default_v1.0.0}"     # ⚠️ leakage-free 컷오프. base_20250630 쓰지 말 것

RUN(){ ( cd "$2" && protenix pred -i "$1" -o results -n "$PROT_MODEL" -s "$SEED" -e "$SAMP" \
         --trimul_kernel torch --triatt_kernel torch --enable_fusion False >"$3" 2>&1 ); }
DONE(){ find "$1/results" -name '*sample*.cif' 2>/dev/null | grep -q .; }

total=0; ok=0; skip=0; fail=0
while IFS=, read -r decoy ag ag_grp ab ab_grp; do
  [ -z "$decoy" ] && continue
  cj="targets/$decoy/chains.json"
  if [ ! -f "$cj" ]; then say "! $decoy: chains.json 없음 — prep_decoy_pairs.py --write 확인"; continue; fi
  agchain=$(python3 -c "import json; d=json.load(open('$cj')); print(d['antigen'][0])") || { say "! $decoy: chains.json 못 읽음"; continue; }

  for r in $RUNGS; do
    a3m="$LADDIR/$ag/$agchain/rung${r}.a3m"
    if [ ! -f "$a3m" ]; then say "! $decoy rung$r: a3m 없음 ($a3m) — 건너뜀"; skip=$((skip+1)); continue; fi
    # ⚠️ analyze_collect_iptm.py 의 정규식이 <타깃>/<깊이>/<실행>/results/ 3단을 기대한다 —
    #    <decoy>/rung<r>/r0/results/ 로 맞춘다(타깃=decoy, 깊이=rung<r>, 실행=r0).
    od="$OUTROOT/$decoy/rung${r}/r0"
    total=$((total+1))
    if DONE "$od"; then say "= $decoy rung$r 이미 완료 — skip"; ok=$((ok+1)); continue; fi
    mkdir -p "$od"
    python make_input.py --cofolder protenix --chains "$cj" --ag-a3m "$agchain=$a3m" --dir "$od" --out "$od/input.json" \
      || { say "! $decoy rung$r: make_input 실패"; fail=$((fail+1)); continue; }
    say "▶ $decoy rung$r (항원=$ag 항체=$ab)"
    RUN "$od/input.json" "$od" "$od/run.log"
    if DONE "$od"; then say "✓ $decoy rung$r 완료"; ok=$((ok+1)); else say "✗ $decoy rung$r 실패 — $od/run.log 확인"; fail=$((fail+1)); fi
    [ "$SMOKE" = "1" ] && { say "[SMOKE] 1건만 실행하고 종료"; exit 0; }
  done
done < <(tail -n +2 "$MANIFEST")

say "끝. 총 $total 건 시도 (완료 $ok · 스킵 $skip · 실패 $fail)"
say "다음: python analyze_collect_iptm.py --data \"$DATA/decoy/protenix\" --out results/decoy_iptm.csv"
say "     (⚠️ 정답이 없다 — DockQ·결합자리 채점 스크립트는 돌리지 말 것)"
