#!/usr/bin/env bash
# ══════════════════════════════════════════════════════════════════════════════
# [지표 통일 ②] analyze_seed_vs_comp.py 를 30종 전부에 대해 돌린다.
#
# 왜 27종이었나 — analyze_seed_vs_comp.py 는 타깃 폴더 밑에 깊이(depth) 하위폴더가 정확히
# 1개일 때만 자동으로 그 폴더를 쓰고, 여러 개면 "--depth-dir 필요"라며 건너뛴다.
# 본 검정에 쓴 깊이는 이미 results/compreps_<타깃>.csv 의 depth 열에 남아 있으므로,
# 그 값을 --depth-dir 로 못박아 넘기면 이 스킵이 사라진다. 로직은 그대로다 — 실행
# 방식만 "타깃 하나 + 정확한 깊이 지정"으로 30번 나눠 돌린다.
#
# 각 타깃을 results/seedcomp_by_target/<타깃>.csv 에 따로 쓰고, 끝나면
# analyze_merge_seedcomp.py 로 하나의 표와 30종 기준 요약(부호검정 포함)을 만든다.
#
# 사용 (conda activate boltz · pipeline/ 에서):
#   bash run_rerun_seedcomp.sh 2>&1 | tee results/rerun_seedcomp_all.log
#   python -u analyze_merge_seedcomp.py
# ══════════════════════════════════════════════════════════════════════════════
set -uo pipefail
DATA="${DATA:-/mnt/data/admuser/msadepth}"
mkdir -p results/seedcomp_by_target

# "타깃:깊이" — results/compreps_<타깃>.csv 의 depth 열에서 그대로 읽은 값.
MAP="8k3k_D:d57 8k46_I:d23 8k5g_HL:d90 8k5h_HL:d15 8kep_HL:d37 8q7s_C:d86 \
8q7s_H:d55 8q7s_O:d35 8siq_HL:d670 8sis_HL:d4169 8sit_HL:d4035 8t4a_PR:d1824 \
8t4d_OQ:d150 8tp5_HL:d307 8u44_ST:d413 8ulr_HL:d1746 8ume_HL:d653 8xsi_HL:d8894 \
9azr_HL:d552 9azt_HL:d4 9azv_HL:d293 9b7g_QP:d625 9bdg_FI:d173 9kkj_HL:d2220 \
9ml9_HL:d8769 9mqr_DE:d3538 9w43_AB:d51 9y0a_AB:d289 9yc6_HL:d611 9zdu_HL:d11017"

ok=0; fail=0
for pair in $MAP; do
  T="${pair%%:*}"; D="${pair##*:}"
  echo "== $T (깊이 $D) =="
  python -u analyze_seed_vs_comp.py --only "$T" --depth-dir "$D" \
         --data "$DATA/compreps/seedrep_cand" \
         --out "results/seedcomp_by_target/${T}.csv" \
    && ok=$((ok+1)) || fail=$((fail+1))
done
echo ""
echo "완료 $ok · 실패 $fail (총 30종)"
echo "→ results/seedcomp_by_target/*.csv 30개. 이어서: python -u analyze_merge_seedcomp.py"
