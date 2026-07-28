#!/usr/bin/env bash
# ══════════════════════════════════════════════════════════════════════════════
# [정리] protenix 예측 폴더를 원래 상태로 되돌린다 — 2026-07-27~28 사고 수습.
#
# 무슨 일이 있었나:
#   · 2026-07-27 15:05 — a3m 오염 수습 과정에서 $DATA/protenix 를 통째로
#     $DATA/protenix_stale_0727_1505 로 옮김. 그런데 **Protenix는 그 사고의 영향을
#     받지 않았다**(make_input.py가 2026-07-22부터 clean_a3m 적용). 즉 불필요한 이동이었다.
#   · 2026-07-27 밤 — run_sweep.sh 는 "출력이 있으면 건너뛴다"로 동작하는데 폴더가 비어
#     있으니 44종을 처음부터 다시 돌렸고, 12시간 예산을 다 쓰고 정지. RBD 10종은 목록
#     맨 뒤라 차례가 오지 않았다. 결과: 쓸모없는 재계산 26종.
#
# 하는 일:
#   ① protenix/ 안에서 **stale 쪽에도 같은 이름이 있는 것만** 삭제(= 확실한 중복 재계산).
#      한쪽에만 있는 것은 절대 건드리지 않는다.
#   ② stale/ 에 남은 것을 protenix/ 로 되돌린다(= 원래 자리, pose_features.csv와 일치).
#   ③ 무슨 일이 있었는지 $DATA/README_protenix.txt 에 남긴다.
#
# 되돌리고 나면 run_sweep.sh 가 44종을 정상적으로 건너뛰므로 목록을 따로 만들 필요가 없다.
#
# ⚠️ 기본은 dry-run. 실제 수행은 --apply. $DATA 밖은 건드리지 않는다.
#
# 사용:
#   bash fix_protenix_dirs.sh              # 무엇을 지우고 옮길지만 출력
#   bash fix_protenix_dirs.sh --apply
# ══════════════════════════════════════════════════════════════════════════════
set -uo pipefail
DATA="${DATA:-/mnt/data/admuser/msadepth}"
NEW="$DATA/protenix"
OLD="${OLD:-$DATA/protenix_stale_0727_1505}"
APPLY=0; [ "${1:-}" = "--apply" ] && APPLY=1

[ -d "$OLD" ] || { echo "!! 원본 보관 폴더 없음: $OLD"; exit 1; }
mkdir -p "$NEW"
echo "$([ $APPLY -eq 1 ] && echo '[실제 수행]' || echo '[dry-run — 아무것도 지우거나 옮기지 않음]')"
echo "  현재 자리 : $NEW"
echo "  원본 보관 : $OLD"
echo ""

dup=(); uniq_new=(); back=()
for p in "$NEW"/*/; do
  [ -d "$p" ] || continue
  b="$(basename "$p")"
  if [ -d "$OLD/$b" ]; then dup+=("$b"); else uniq_new+=("$b"); fi
done
for p in "$OLD"/*/; do
  [ -d "$p" ] || continue
  back+=("$(basename "$p")")
done

echo "① 지울 것 — 어젯밤 재계산분(원본에도 같은 이름이 있는 것) ${#dup[@]}개"
printf '     %s\n' "${dup[@]:-(없음)}"
echo ""
if [ ${#uniq_new[@]} -gt 0 ]; then
  echo "⚠️ 현재 자리에만 있는 것 ${#uniq_new[@]}개 — 원본에 없으므로 **건드리지 않음**"
  printf '     %s\n' "${uniq_new[@]}"
  echo ""
fi
echo "② 되돌릴 것 — 원본 보관 → 현재 자리 ${#back[@]}개"
echo "     ${back[*]:-(없음)}"
echo ""

if [ $APPLY -eq 0 ]; then
  echo "→ 실제로 수행하려면:  bash fix_protenix_dirs.sh --apply"
  exit 0
fi

for b in "${dup[@]:-}"; do
  [ -z "$b" ] && continue
  rm -rf "${NEW:?}/$b" && echo "  삭제 $b"
done
for b in "${back[@]:-}"; do
  [ -z "$b" ] && continue
  if [ -e "$NEW/$b" ]; then echo "  !! 이미 있음, 건너뜀 $b"; continue; fi
  mv "$OLD/$b" "$NEW/$b" && echo "  복귀 $b"
done
rmdir "$OLD" 2>/dev/null && echo "  빈 보관 폴더 제거: $OLD"

cat > "$DATA/README_protenix.txt" <<'EOF'
protenix 예측 폴더 이력 (2026-07-28 정리)

2026-07-27 15:05  a3m 질의행 오염 수습 중 protenix/ 를 protenix_stale_0727_1505/ 로 이동.
                  → 불필요한 이동이었다. Protenix는 그 사고의 영향을 받지 않았다
                    (make_input.py 가 2026-07-22부터 clean_a3m 을 protenix·chai 에 적용).
                    실제로 MSA를 잃은 것은 boltz 뿐이다.

2026-07-27 밤     run_sweep.sh protenix 12 실행. 출력 폴더가 비어 있어 self-heal skip 이
                  작동하지 않았고, sweep_targets.csv 를 위에서부터 다시 돌려 44종 중 26종을
                  재계산한 뒤 12시간 예산 소진으로 정지. RBD 10종(목록 맨 뒤)은 미실행.

2026-07-28        fix_protenix_dirs.sh --apply 로 정리:
                  · 어젯밤 재계산분(원본에 같은 이름이 있는 것) 삭제
                  · protenix_stale_0727_1505/ 의 원본을 protenix/ 로 복귀

정본(source of truth) = pipeline/results/pose_features.csv 의 protenix 2640행.
이 파일은 복귀시킨 원본 예측으로 채점된 것이며, 2026-07-27 거울상 검정도 이 데이터로 냈다.
⚠️ pose_features.py 를 --rescore 로 돌리지 말 것(이미 채점된 자세를 다시 계산할 이유가 없다).

교훈: 출력 폴더를 옮기면 run_sweep.sh 의 "이미 있으면 건너뛰기"가 무력화되어
      목록 앞쪽부터 전부 재실행된다. 새 타깃만 돌리려면 LIST= 로 목록을 제한할 것.
EOF
echo ""
echo "→ $DATA/README_protenix.txt 에 이력 기록"
echo "→ 이제 run_sweep.sh 가 44종을 건너뛴다:  bash run_sweep.sh protenix 14"
