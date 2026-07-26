#!/usr/bin/env bash
# ══════════════════════════════════════════════════════════════════════════════
# Track A 밤샘 오케스트레이터 (집 가서 돌리기용) — Stage 0(게이트)까지만 무인.
#   ① Phase0(pose_features.py) 끝나길 대기  → ② seed_replicate(CPU)
#   ③ GPU 비길(tFold 등 종료) 대기          → ④ SMOKE 게이트(깨졌으면 즉시 정지)
#   ⑤ 전체 Protenix + Chai 예측             → ⑥ 정지(채점·nested/LOCO/template은 게이트 결과 보고 결정)
#
# ⚠️ 왜 Stage 0까지만: 적대검증 판정 = 신호가 seed 게이트를 넘어야만 nested/LOCO로.
#    게이트 결과를 안 보고 다 돌리면 GPU 낭비 + forking-paths. 그래서 여기서 멈추고 아침에 판단.
# ⚠️ SMOKE 게이트: 새 러너가 깨졌으면 전체 몇 시간 태우기 전에 1건에서 멈춤(밤 낭비 방지).
#
# 사용(tmux 권장 — 세션 끊겨도 유지):
#   cd ~/projects/bk21-msa-depth-bias/pipeline && git pull
#   tmux new -s trackA
#   bash run_track_a_overnight.sh 2>&1 | tee overnight_$(date +%m%d).log
#   # Ctrl-b d 로 detach 후 집에. 아침에 overnight_*.log 마지막부분 붙여주세요.
# ══════════════════════════════════════════════════════════════════════════════
set -uo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"; cd "$HERE"
say(){ echo "[$(date '+%m-%d %H:%M:%S')] $*"; }

say "======== Track A 밤샘 시작 ========"

# ① Phase 0(pose_features.py) 종료 대기
if pgrep -f pose_features.py >/dev/null 2>&1; then
  say "① Phase0(pose_features.py) 실행 중 — 종료 대기(2분 간격)..."
  while pgrep -f pose_features.py >/dev/null 2>&1; do sleep 120; done
fi
say "① Phase0 종료 확인."

# ② seed_replicate (CPU) — a3m 생성. 이미 있으면 seed_replicate가 알아서 건너뜀/재생성(무해).
say "② seed_replicate(CPU) 실행..."
bash run_track_a_seedrep.sh || { say "!! seed_replicate 실패 — 중단(로그 확인)"; exit 1; }

# ③ GPU 유휴 대기(tFold 등 종료). nvidia-smi 컴퓨트앱 비면 통과.
if command -v nvidia-smi >/dev/null 2>&1; then
  say "③ GPU 비길 대기(tFold 등)..."
  while nvidia-smi --query-compute-apps=pid --format=csv,noheader 2>/dev/null | grep -q .; do
    say "   GPU 사용 중 — 5분 후 재확인"; sleep 300
  done
  say "③ GPU 유휴 확인."
else
  say "③ ⚠️ nvidia-smi 없음 — GPU 상태 확인 불가, 그대로 진행."
fi

# ④ SMOKE 게이트 — 러너 정상인지 1건으로 확인. 실패 시 정지(밤 낭비 방지).
say "④ SMOKE — Protenix 1건..."
if SMOKE=1 bash run_seedrep_predict.sh protenix; then
  say "④ Protenix smoke OK"
else
  say "!! Protenix smoke 실패 — 밤샘 정지. 아침에 run.log 보고 고칩니다(tFold·Phase0는 이미 완료)."
  exit 1
fi
CHAI_OK=1
say "④ SMOKE — Chai 1건..."
if SMOKE=1 bash run_seedrep_predict.sh chai; then say "④ Chai smoke OK"; else say "!! Chai smoke 실패 — Protenix만 전체 진행"; CHAI_OK=0; fi

# ⑤ 전체 예측
say "⑤ 전체 Protenix 예측..."
bash run_seedrep_predict.sh protenix || say "!! Protenix 전체 일부 실패(개별 run.log 확인)"
if [ "$CHAI_OK" = 1 ]; then
  say "⑤ 전체 Chai 예측..."
  bash run_seedrep_predict.sh chai || say "!! Chai 전체 일부 실패"
fi

# ⑥ 정지 — 게이트까지만.
DATA="${DATA:-/mnt/data/admuser/msadepth}"
np=$(find "$DATA/seedrep_pred" -path '*/results/*sample*.cif' 2>/dev/null | wc -l)
nc=$(find "$DATA/seedrep_pred/chai" -path '*/results/*.cif' 2>/dev/null | wc -l)
say "======== Track A Stage 0 완료 ========"
say "  Protenix pose(cif): $np · Chai pose(cif): $nc"
say "  ⚠️ nested/LOCO/template = 자동 실행 안 함(게이트 결과 보고 결정)."
say "  아침에 붙여주세요: 이 로그 마지막 20줄 + 위 pose 개수. epitope-recall 채점(같은 깊이 5seed 분산)은 이어서 진행."
