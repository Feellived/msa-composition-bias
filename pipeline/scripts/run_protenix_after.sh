#!/usr/bin/env bash
# Boltz(run_sweep boltz)가 끝나면 → Protenix 스모크 1건 → 성공(cif 생성)해야 전체 Protenix sweep.
# GPU 하나를 Boltz와 안 겹치게: Boltz 프로세스가 사라질 때까지 대기 후 시작.
#   bash run_protenix_after.sh [전체HOURS=22]
# (별도 tmux 창에서 실행 권장. 이미 수동 스모크로 cif가 있으면 재스모크 생략하고 바로 전체.)
set -uo pipefail
# 스크립트는 pipeline/scripts/ 에 있고, 기준 디렉토리는 상위 pipeline/ 이다
HERE="$(cd "$(dirname "$0")/.." && pwd)"; cd "$HERE" || exit 1
HOURS="${1:-22}"
DATA="${DATA:-/mnt/data/admuser/msadepth}"
say(){ echo "[$(date '+%m-%d %H:%M:%S')] [prot-after] $*"; }
has_cif(){ find "$DATA/protenix" -name '*sample*.cif' 2>/dev/null | grep -q .; }

say "Boltz(run_sweep boltz) 종료 대기 (5분 간격 확인)..."
while pgrep -f 'run_sweep.sh boltz' >/dev/null; do sleep 300; done
say "Boltz 종료 감지."

if has_cif; then
  say "이미 Protenix 출력 있음(수동 스모크됨) → 재스모크 생략, 바로 전체 sweep."
else
  say "Protenix 스모크 1건..."
  SMOKE=1 bash run_sweep.sh protenix 1
  if ! has_cif; then
    say "!! 스모크 실패(cif 없음) — 전체 실행 중단. 로그 확인: \$DATA/protenix/<타깃>/rung0/run.log"
    exit 1
  fi
  say "스모크 성공(cif 생성 확인)."
fi

say "전체 Protenix sweep 시작 (${HOURS}h)..."
bash run_sweep.sh protenix "$HOURS"
say "Protenix sweep 종료. 채점: python eval_dockq_sweep.py --models boltz protenix"
