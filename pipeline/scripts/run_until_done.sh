#!/usr/bin/env bash
# ══════════════════════════════════════════════════════════════════════════════
# [밤새 이어달리기] run_maintest.sh 를 남은 타깃이 없어질 때까지 반복 호출한다.
#
# 왜 필요한가: run_maintest.sh 는 타깃 하나를 시작하기 전에 남은 예산(HOURS)을 보고
# 모자라면 멈춘다(타깃 중간에 죽이지 않으려는 안전장치). 그래서 예산이 짧으면 아침에
# 와서 보면 중간에 서 있다. 이 래퍼가 예산이 끝날 때마다 다시 불러 이어서 가게 한다.
# 완료된 타깃은 run_maintest.sh 가 알아서 건너뛰므로 같은 일을 두 번 하지 않는다.
#
# 안전장치
#   · 이미 돌고 있는 본 검정이 있으면 그것이 끝날 때까지 기다렸다가 시작한다
#     (GPU를 두 개가 동시에 잡으면 둘 다 느려지고 결과 폴더가 꼬인다).
#   · 진척 없음 감시 — 한 바퀴 돌았는데 남은 실행 수가 줄지 않으면(무한 반복 = 매번 같은
#     타깃에서 실패) 2회 연속에서 멈추고 끝낸다. 밤새 헛도는 것을 막는다.
#   · 최대 반복 횟수(MAX_ITERS)로도 막는다.
#   · 채점(run_analyze_target.sh)은 DockQ 환경이 필요하므로 conda run 으로 부른다.
#     AUTO_ANALYZE=0 으로 끌 수 있다.
#
# 사용:
#   tmux new -s overnight
#   bash run_until_done.sh                  # 예산 6시간짜리로 나눠 끝까지
#   HOURS=8 MAX_ITERS=12 bash run_until_done.sh
#   AUTO_ANALYZE=0 bash run_until_done.sh   # 예측만 하고 채점은 손으로
#
# 로그는 $DATA/logs/until_done_<시각>.log 에도 남는다.
# ══════════════════════════════════════════════════════════════════════════════
set -uo pipefail
cd "$(cd "$(dirname "$0")" && pwd)"
DATA="${DATA:-/mnt/data/msadepth}"
HOURS="${HOURS:-6}"
MAX_ITERS="${MAX_ITERS:-20}"
AUTO_ANALYZE="${AUTO_ANALYZE:-1}"
DOCKQ_ENV="${DOCKQ_ENV:-DockQ}"
WAIT_SEC="${WAIT_SEC:-300}"

if [ -z "${UD_TEED:-}" ]; then
  UDLOG="${UDLOG:-$DATA/logs/until_done_$(date +%m%d_%H%M%S).log}"
  mkdir -p "$(dirname "$UDLOG")" 2>/dev/null || UDLOG="/tmp/$(basename "$UDLOG")"
  UD_TEED=1 UDLOG="$UDLOG" bash "$0" "$@" 2>&1 | tee -a "$UDLOG"
  rc=${PIPESTATUS[0]}
  echo "[로그] $UDLOG"
  exit "$rc"
fi
say(){ echo "[$(date '+%m-%d %H:%M:%S')] [이어달리기] $*"; }

# ── 남은 실행 수 (dry-run 의 마지막 집계 줄에서 읽는다) ────────────────────────
remaining(){
  bash run_maintest.sh 2>/dev/null | sed -n 's/.*남은 실행 약 \([0-9]*\)회.*/\1/p' | tail -1
}
all_done(){ bash run_maintest.sh 2>/dev/null | grep -q "전부 완료 상태다"; }

# ── 다른 본 검정이 돌고 있으면 기다린다 ───────────────────────────────────────
busy(){ pgrep -f "run_maintest.sh --apply" >/dev/null || pgrep -f "make_composition_reps.sh" >/dev/null; }
if busy; then
  say "이미 돌고 있는 본 검정이 있다 — 끝나면 이어서 시작한다 (${WAIT_SEC}초마다 확인)"
  while busy; do sleep "$WAIT_SEC"; done
  say "앞의 실행이 끝났다. 이어서 간다."
fi

# ── 본체 ──────────────────────────────────────────────────────────────────────
stall=0
for ((i=1; i<=MAX_ITERS; i++)); do
  if all_done; then say "남은 타깃 없음 — 예측 단계 완료"; break; fi
  before=$(remaining); before=${before:-0}
  say "$i 바퀴 · 남은 실행 약 ${before}회 · 이번 예산 ${HOURS}시간"
  HOURS="$HOURS" bash run_maintest.sh --apply
  after=$(remaining); after=${after:-0}
  say "$i 바퀴 끝 · 남은 실행 ${before} → ${after}"
  if [ "$after" -ge "$before" ]; then
    stall=$((stall+1))
    say "⚠️ 진척 없음 ${stall}회 — 같은 곳에서 막혔을 수 있다"
    if [ "$stall" -ge 2 ]; then
      say "!! 두 바퀴 연속 진척이 없어 멈춘다. 로그에서 마지막 타깃의 실패 메시지를 확인할 것."
      exit 3
    fi
  else
    stall=0
  fi
done

if ! all_done; then
  say "!! 최대 반복(${MAX_ITERS})까지 왔는데 남은 타깃이 있다. 같은 명령을 다시 실행하면 이어서 간다."
  exit 4
fi

# ── 채점 ──────────────────────────────────────────────────────────────────────
if [ "$AUTO_ANALYZE" != "1" ]; then
  say "예측 완료. 채점은 직접: conda activate $DOCKQ_ENV && bash run_analyze_target.sh <타깃들>"
  exit 0
fi
mapfile -t TGTS < <(python3 - maintest.csv <<'PY'
import csv, sys
for r in csv.DictReader(open(sys.argv[1])):
    if r.get("status") == "run":
        print(r["target"])
PY
)
say "예측 완료 → 채점 시작 (${#TGTS[@]}개, DockQ 환경 '$DOCKQ_ENV', 이미 채점된 것은 건너뜀)"
fail=0
for t in "${TGTS[@]}"; do
  conda run -n "$DOCKQ_ENV" --no-capture-output bash run_analyze_target.sh "$t" || { fail=$((fail+1)); say "  ! $t 채점 실패"; }
done
say "채점 끝 (실패 ${fail}개). 다음: conda run -n $DOCKQ_ENV python analyze_collect_results.py"
[ "$fail" -eq 0 ] || exit 5
