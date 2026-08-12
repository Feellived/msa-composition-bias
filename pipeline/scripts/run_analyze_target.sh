#!/usr/bin/env bash
# ══════════════════════════════════════════════════════════════════════════════
# [표준 분석] 반복실행(comp_x_reps) 결과를 한 타깃에 대해 정해진 순서로 전부 분석.
#
#   ① 자세 단위 채점        eval_dump_seedrep.py  → results/compreps_<타깃>.csv
#   ② 성공률·통계           eval_compreps.py     → 세 지표로 각각(원래 vs 얕은 + 조성 간 이질성)
#        ②-a DockQ(자세 정확도)  ②-b 진짜 자리 겹침  ②-c 흔한 자리 겹침(낮을수록 편향 이탈)
#   ③ 기제(결합자리 이동)   analyze_epitope_cluster.py    → 성공/실패 자리 겹침·크기
#   ④ 핵심 통계             analyze_site_reproducibility.py → (조성 내)/(조성 간) 결합자리 겹침 + 뒤섞기 검정
#                                                     + 서로 구별되는 자리 후보 개수와 각 후보의 진짜 자리 덮음
#
# ⚠️ ②의 이질성 검정은 **조성당 반복 2회 이상**일 때만 나온다(예비검정은 1회라 안 나옴).
#    예비검정에서 유망하면 반복을 늘려 다시 돌릴 것.
#
# 사용(DockQ env):
#   bash run_analyze_target.sh 8ulr_HL
#   bash run_analyze_target.sh 9azr_HL 8k5g_HL 8q7s_C 8ume_HL      # 여러 개
# ══════════════════════════════════════════════════════════════════════════════
set -uo pipefail
cd "$(cd "$(dirname "$0")" && pwd)"
DATA="${DATA:-/mnt/data/msadepth}"
[ $# -ge 1 ] || { echo "사용: bash run_analyze_target.sh <타깃> [타깃...]"; exit 1; }
for T in "$@"; do
  echo ""
  echo "████████████████ $T ████████████████"
  csv="results/compreps_${T}.csv"
  log="results/analyze_${T}.txt"
  : > "$log"                      # 화면에 나오는 통계를 파일로도 남긴다(터미널을 잃어도 보존)
  # 자세 단위 채점(DockQ 160회)이 제일 비싸다 — 이미 있으면 다시 계산하지 않는다.
  if [ -s "$csv" ] && [ "${REDO:-0}" != "1" ]; then
    echo "  (자세 단위 채점 결과가 이미 있다 — 재계산 생략. 강제로 다시 하려면 REDO=1)" | tee -a "$log"
  else
    # DEPTH=d90 로 깊이 폴더를 못박을 수 있다. 안 주면 maintest.csv 의 n_rows 로 고른다.
    # 깊이 폴더가 여럿인데 어느 것인지 못 정하면 dump 가 종료코드 4로 멈춘다(섞이는 것보다 낫다).
    python eval_dump_seedrep.py --data "$DATA/compreps" --only "$T" \
           ${DEPTH:+--depth "$DEPTH"} --csv-out "$csv" > "results/dump_${T}.txt" 2>&1
    rc=$?
    if [ "$rc" = "5" ]; then
      echo "  !! DockQ 값이 하나도 안 나왔다 — DockQ 환경 밖에서 돌린 것 같다." | tee -a "$log"
      echo "     conda activate DockQ 후 다시 실행할 것 (recall·overrep 만으로는 ③④가 안 된다)." | tee -a "$log"
      continue
    fi
    if [ "$rc" = "4" ]; then
      echo "  !! 깊이 폴더가 여러 개라 멈췄다 — 설계가 다른 실행이 섞이는 것을 막은 것." | tee -a "$log"
      grep -E "깊이 폴더|--depth" "results/dump_${T}.txt" | sed 's/^/     /' | tee -a "$log"
      echo "     예: DEPTH=d90 bash run_analyze_target.sh $T" | tee -a "$log"; continue
    fi
  fi
  if [ ! -s "$csv" ]; then
    echo "  !! 자료 없음 (예측이 아직 없거나 실패) — results/dump_${T}.txt 확인" | tee -a "$log"; continue
  fi
  # ⚠️ --out 을 지정하지 않으면 세 지표가 전부 results/compreps_summary.csv 한 파일에
  #    덮어써져, 타깃 29개 x 지표 3종 = 87번 중 마지막 하나만 남는다(2026-07-29 발견).
  echo "── ②-a 자세 정확도(DockQ, 높을수록 좋음) ──" | tee -a "$log"
  python eval_compreps.py --csv "$csv" --label dockq \
         --out "results/summary_${T}_dockq.csv" 2>&1 | tee -a "$log"
  echo "── ②-b 진짜 결합자리 겹침(높을수록 좋음) ──" | tee -a "$log"
  python eval_compreps.py --csv "$csv" --label recall --succ-th 0.4 \
         --out "results/summary_${T}_recall.csv" 2>&1 | tee -a "$log"
  echo "── ②-c 흔한 자리 겹침(낮을수록 좋음 = 편향 이탈) ──" | tee -a "$log"
  python eval_compreps.py --csv "$csv" --label overrep --lower-better --succ-th 0.3 \
         --out "results/summary_${T}_overrep.csv" 2>&1 | tee -a "$log"
  echo "── ③ 기제(결합자리) ──" | tee -a "$log"
  python analyze_epitope_cluster.py --csv "$csv" --data "$DATA/compreps" 2>&1 | tee -a "$log"
  echo "── ④ 핵심: 조성이 자리를 정하나 + 후보 몇 개 ──" | tee -a "$log"
  # --dump-sites: 후보 자리의 잔기 목록을 results/sites_<타깃>.json 으로. 유도 재도킹 데모 입력.
  python analyze_site_reproducibility.py --csv "$csv" --data "$DATA/compreps" --dump-sites 2>&1 | tee -a "$log"
done
echo ""
echo "타깃마다 남는 것:"
echo "  results/compreps_<타깃>.csv      자세 단위 원자료(실행 x 자세)"
echo "  results/summary_<타깃>_<지표>.csv 지표별 요약(이질성 p, 성공 수, 검정 p)"
echo "  results/epitope_cluster_<타깃>.csv 결합자리 무리"
echo "  results/site_repro_<타깃>.csv     조성이 자리를 정하나(핵심)"
echo "  results/analyze_<타깃>.txt        위 화면 출력 전체"
echo ""
echo "전부 끝나면 한 파일로 모으기:  python analyze_collect_results.py"
