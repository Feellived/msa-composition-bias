#!/usr/bin/env bash
# ══════════════════════════════════════════════════════════════════════════════
# [표준 분석] 반복실행(comp_x_reps) 결과를 한 타깃에 대해 정해진 순서로 전부 분석.
#
#   ① 자세 단위 채점        dump_seedrep_full.py  → results/compreps_<타깃>.csv
#   ② 성공률·통계           score_compreps.py     → 원래 vs 얕은, 조성 간 이질성
#   ③ 기제(결합자리 이동)   epitope_cluster.py    → 성공/실패 자리 겹침·크기
#
# ⚠️ ②의 이질성 검정은 **조성당 반복 2회 이상**일 때만 나온다(예비검정은 1회라 안 나옴).
#    예비검정에서 유망하면 반복을 늘려 다시 돌릴 것.
#
# 사용(DockQ env):
#   bash analyze_target.sh 8ulr_HL
#   bash analyze_target.sh 9azr_HL 8k5g_HL 8q7s_C 8ume_HL      # 여러 개
# ══════════════════════════════════════════════════════════════════════════════
set -uo pipefail
cd "$(cd "$(dirname "$0")" && pwd)"
DATA="${DATA:-/mnt/data/admuser/msadepth}"
[ $# -ge 1 ] || { echo "사용: bash analyze_target.sh <타깃> [타깃...]"; exit 1; }
for T in "$@"; do
  echo ""
  echo "████████████████ $T ████████████████"
  csv="results/compreps_${T}.csv"
  python dump_seedrep_full.py --data "$DATA/compreps" --only "$T" --csv-out "$csv" > "results/dump_${T}.txt" 2>&1
  if [ ! -s "$csv" ]; then
    echo "  !! 자료 없음 (예측이 아직 없거나 실패) — results/dump_${T}.txt 확인"; continue
  fi
  echo "── ② 성공률·통계 ──"
  python score_compreps.py --csv "$csv"
  echo "── ③ 기제(결합자리) ──"
  python epitope_cluster.py --csv "$csv" --data "$DATA/compreps"
done
echo ""
echo "요약표: results/compreps_summary.csv · 그림용: results/epitope_cluster_<타깃>.csv"
