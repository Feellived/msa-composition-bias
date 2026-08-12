#!/usr/bin/env bash
# ══════════════════════════════════════════════════════════════════════════════
# [일괄 실행] 후보 명단(candidates.csv)의 복합체를 전부 같은 절차로 돌린다.
#
# 복합체마다: ① 줄인 MSA — 조성 여러 가지 × 반복   ② 원래 MSA — 반복
#            (make_composition_reps.sh가 조성 a3m 생성·seedfull 준비까지 알아서 함)
#
# 명단 형식(쉼표 구분, comps는 공백으로 구분):
#   target,model,rung,comps,reps,full_reps,note
#
# 사용:
#   bash run_candidates.sh                       # 명단 전부
#   ONLY="8k3k_D 8tx3_FK" bash run_candidates.sh # 일부만
#   ANALYZE=1 bash run_candidates.sh             # 끝나고 분석까지
#   DRY=1 bash run_candidates.sh                 # 뭘 돌릴지만 보기(실행 안 함)
# ══════════════════════════════════════════════════════════════════════════════
set -uo pipefail
cd "$(cd "$(dirname "$0")" && pwd)"
CAND="${CAND:-candidates.csv}"; ONLY="${ONLY:-}"; ANALYZE="${ANALYZE:-0}"; DRY="${DRY:-0}"
say(){ echo "[$(date '+%m-%d %H:%M:%S')] $*"; }
[ -f "$CAND" ] || { say "!! 명단 없음: $CAND"; exit 1; }

done_list=""
while IFS=, read -r target model rung comps reps full_reps note; do
  [ -z "${target:-}" ] && continue
  case "$target" in target|\#*) continue;; esac
  if [ -n "$ONLY" ]; then case " $ONLY " in *" $target "*) ;; *) continue;; esac; fi
  echo ""
  say "████ $target ($model, rung$rung) — $note"
  if [ "$DRY" = 1 ]; then
    echo "   ① RUNG=$rung TARGET=$target MODEL=$model COMPS=\"$comps\" REPS=$reps"
    echo "   ② RUNG=$rung TARGET=$target MODEL=$model COMPS=\"full\" REPS=$full_reps"
    done_list="$done_list $target"; continue
  fi
  RUNG="$rung" TARGET="$target" MODEL="$model" COMPS="$comps" REPS="$reps" bash make_composition_reps.sh \
    || say "  !! ① 실패 $target (계속 진행)"
  RUNG="$rung" TARGET="$target" MODEL="$model" COMPS="full" REPS="$full_reps" bash make_composition_reps.sh \
    || say "  !! ② 실패 $target (계속 진행)"
  done_list="$done_list $target"
done < "$CAND"

say "실행 완료:$done_list"
if [ "$ANALYZE" = 1 ] && [ "$DRY" != 1 ] && [ -n "$done_list" ]; then
  say "분석 시작 (DockQ 환경이어야 함)"
  bash run_analyze_target.sh $done_list
else
  say "분석은 따로:  bash run_analyze_target.sh$done_list"
fi
