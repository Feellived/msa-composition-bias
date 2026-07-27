#!/usr/bin/env bash
# ══════════════════════════════════════════════════════════════════════════════
# [사다리 생성 — 새 타깃만] sweep_targets.csv를 훑어 사다리가 아직 없는 항원 사슬에만
# build_ladder.py 를 돌린다.
#
# ⚠️ 왜 "새 타깃만"인가: 기존 49종의 사다리를 다시 만들면 각 칸의 **조성(뽑힌 서열 목록)**이
#    바뀌어 지금까지의 결과가 전부 무효가 된다. 그래서 `rung1.a3m`이 이미 있으면 건너뛴다.
#    (rung0=full은 import_rbd가 먼저 만들어 두므로 존재 여부 판정에 쓰면 안 된다.)
#
# 전제: $DATA/ladders/<타깃>/<사슬>/rung0.a3m 이 있어야 한다(= 전체 MSA).
#       import_rbd.py --stage msa 가 만들거나, gen_msa.sh 산출물을 놓아둔 상태.
#
# 사용:
#   bash build_ladders_new.sh                 # dry-run: 무엇을 만들지만 출력
#   bash build_ladders_new.sh --apply
#   ONLY="8sis_HL 8p5m_GL" bash build_ladders_new.sh --apply
#   RUNGS=12 bash build_ladders_new.sh --apply
# ══════════════════════════════════════════════════════════════════════════════
set -uo pipefail
cd "$(cd "$(dirname "$0")" && pwd)"
DATA="${DATA:-/mnt/data/admuser/msadepth}"
LIST="${LIST:-sweep_targets.csv}"
RUNGS="${RUNGS:-12}"          # run_sweep.sh 기본값과 반드시 같아야 한다(rung0~rung11)
MINROWS="${MINROWS:-1}"       # 최심 칸 = 단일서열
SEED="${SEED:-0}"
ONLY="${ONLY:-}"
APPLY=0; [ "${1:-}" = "--apply" ] && APPLY=1

say(){ echo "[$(date '+%m-%d %H:%M:%S')] $*"; }
[ $APPLY -eq 1 ] && say "=== 실제 생성 (rungs=$RUNGS, min-rows=$MINROWS, seed=$SEED) ===" \
                 || say "=== dry-run — 아무것도 만들지 않음 ==="

n_new=0; n_skip=0; n_miss=0
while IFS=, read -r target pdb group ab dirtype ag_chains label; do
  [ -z "$target" ] && continue
  [ "$target" = "target" ] && continue
  if [ -n "$ONLY" ]; then
    case " $ONLY " in *" $target "*) ;; *) continue;; esac
  fi
  IFS='|' read -ra AGC <<< "$ag_chains"
  for c in "${AGC[@]}"; do
    [ -z "$c" ] && continue
    d="$DATA/ladders/$target/$c"
    if [ ! -f "$d/rung0.a3m" ]; then
      say "  없음  $target/$c — rung0.a3m 이 없다(전체 MSA 미준비)"; n_miss=$((n_miss+1)); continue
    fi
    if [ -f "$d/rung1.a3m" ]; then
      n_skip=$((n_skip+1)); continue     # 이미 사다리 있음 → 절대 덮지 않는다
    fi
    nseq=$(grep -c '^>' "$d/rung0.a3m")
    say "  생성  $target/$c  (전체 $nseq 서열 → $RUNGS 칸)"
    n_new=$((n_new+1))
    [ $APPLY -eq 0 ] && continue
    cp -n "$d/rung0.a3m" "$d/full_backup.a3m"          # 만일을 대비한 원본 사본
    python build_ladder.py --a3m "$d/full_backup.a3m" --outdir "$d" \
           --rungs "$RUNGS" --min-rows "$MINROWS" --seed "$SEED" || say "  !! 실패 $target/$c"
  done
done < "$LIST"

echo ""
say "생성 대상 $n_new · 이미 있어 건너뜀 $n_skip · rung0 없음 $n_miss"
if [ $APPLY -eq 0 ]; then
  echo "→ 실제로 만들려면:  bash build_ladders_new.sh --apply"
else
  echo "→ 다음: python check_msa_match.py   (사다리 칸까지 확인하려면 --rung rung5 등)"
  echo "→ 그다음: bash run_sweep.sh protenix <시간>"
fi
