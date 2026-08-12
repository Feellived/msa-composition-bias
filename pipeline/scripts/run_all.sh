#!/usr/bin/env bash
# Boltz → Protenix 순차 실행 (총 예산 공유). run_sweep.sh를 두 번 부르되 deadline을 나눠 씀.
# Protenix 기본모델=protenix_base_default_v1.0.0 (학습컷오프 2021-09-30=leakage-free). 바꾸려면 PROT_MODEL 지정.
#   bash run_all.sh 22                                  # 총 22시간(밤샘~주말), Boltz 먼저 → 남는 시간 Protenix
#   PROT_ENV=protenix PROTENIX_ROOT_DIR=/mnt/data/admuser/protenix_weights bash run_all.sh 40
set -uo pipefail
# 스크립트는 pipeline/scripts/ 에 있고, 기준 디렉토리는 상위 pipeline/ 이다
HERE="$(cd "$(dirname "$0")/.." && pwd)"; cd "$HERE" || exit 1
HOURS="${1:-22}"
export PROT_MODEL="${PROT_MODEL:-protenix_base_default_v1.0.0}"   # 컷오프 2021-09-30=leakage-free(공식 권장·AF3 능가). 20250630=2025컷오프=leaky 금지. v2는 다운로드 403이라 보류
say(){ echo "[$(date '+%m-%d %H:%M:%S')] [run_all] $*"; }
START=$(date +%s); DEADLINE=$(( START + HOURS*3600 ))

# 1) Boltz — 남은 예산 전부
rem_h(){ echo $(( ( DEADLINE - $(date +%s) + 3599 ) / 3600 )); }   # 남은 시간(시, 올림)
H=$(rem_h)
say "=== [1/2] Boltz | 남은 예산 ${H}h ==="
bash "$HERE/run_sweep.sh" boltz "$H"

# 2) Protenix — Boltz가 끝나고 남은 예산으로
H=$(rem_h)
if [ "$H" -le 0 ]; then
  say "예산 소진 — Protenix 미시작. 재실행하면 self-heal로 Boltz 잔여 → Protenix 순서로 이어감."
  exit 0
fi
say "=== [2/2] Protenix | 남은 예산 ${H}h (PROT_MODEL=$PROT_MODEL) ==="
bash "$HERE/run_sweep.sh" protenix "$H"
say "완료 스캔 끝. 채점: python eval_dockq_sweep.py --models boltz protenix (DockQ env)."
