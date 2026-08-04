#!/usr/bin/env bash
# ══════════════════════════════════════════════════════════════════════════════
# [지표 통일 ①] sites_<타깃>.json 을 지금의 ms_site_repro_*.csv 와 같은 기준으로
# 30종 전부 다시 쓴다.
#
# 왜 필요한가 — pipeline/results/sites_*.json 은 후보를 "합집합"으로 묶던 옛 방식으로
# 만들어진 뒤 갱신되지 않았다. report/data/ms_site_repro_*.csv 는 "투표" 방식(현재
# 기본값, --merge-frac 0.75)으로 다시 돌린 최신판이다. 두 자료가 30종 중 18종에서
# 후보 개수 자체가 달라, 4.6절 선택기 비교(F1 등)가 4.4·4.5절과 다른 후보를 보고 있었다.
#
# 이 스크립트는 site_reproducibility.py 를 --dump-sites 를 붙여 그대로 다시 돌릴 뿐,
# 새 로직은 없다 — analyze_target.sh ④단계와 완전히 같은 호출이다.
#
# 사용 (conda activate boltz · pipeline/ 에서):
#   bash rerun_sites_all.sh 2>&1 | tee results/rerun_sites_all.log
# ══════════════════════════════════════════════════════════════════════════════
set -uo pipefail
DATA="${DATA:-/mnt/data/admuser/msadepth}"

TARGETS="8k3k_D 8k46_I 8k5g_HL 8k5h_HL 8kep_HL 8q7s_C 8q7s_H 8q7s_O 8siq_HL 8sis_HL \
8sit_HL 8t4a_PR 8t4d_OQ 8tp5_HL 8u44_ST 8ulr_HL 8ume_HL 8xsi_HL 9azr_HL 9azt_HL \
9azv_HL 9b7g_QP 9bdg_FI 9kkj_HL 9ml9_HL 9mqr_DE 9w43_AB 9y0a_AB 9yc6_HL 9zdu_HL"

ok=0; fail=0
for T in $TARGETS; do
  csv="results/compreps_${T}.csv"
  if [ ! -f "$csv" ]; then
    echo "!! $T: $csv 없음 — 건너뜀"; fail=$((fail+1)); continue
  fi
  echo "== $T =="
  python -u site_reproducibility.py --csv "$csv" --data "$DATA/compreps" \
         --dump-sites --out "results/ms_site_repro_${T}.csv" \
    && ok=$((ok+1)) || fail=$((fail+1))
done
echo ""
echo "완료 $ok · 실패 $fail (총 30종)"
echo "→ results/sites_<타깃>.json 30개, results/ms_site_repro_<타깃>.csv 30개 갱신됨"
