#!/usr/bin/env bash
# ══════════════════════════════════════════════════════════════════════════════
# Track A · Stage 0 (게이트) — seed_replicate: 앵커 항원 MSA를 같은 깊이×다른 seed로.
#
# 왜 이게 0단계인가 (적대검증 만장일치, 2026-07-26):
#   존재증명(case study)에서 어떤 '이탈/rescue'도 best-of-N·seed 잡음 밴드를 넘어야만 유효.
#   seed_replicate = 개입이 아니라 그 잡음 밴드를 재는 협상불가 게이트(2026-07-23 linchpin).
#   이게 nested/LOCO/깊이스캔보다 논리적으로 위 — 신호가 실재하는지부터 확인.
#
# 이 스크립트 = CPU만(GPU 불필요). 각 앵커의 '항원' a3m을 같은 깊이(행수)×다른 seed 5개로 재추첨.
#   대상 = 항원 MSA (편향을 나르는 사슬; ladders/<t>/<항원사슬>/). 항체 MSA 아님.
#   깊이 선택 = neff.tsv에서 full·single 제외하고 얕은~중간 3개(사전등록 고정, 적응형 금지=forking-paths).
#
# 이후(GPU, tFold 후): make_input.py로 각 seed a3m → Protenix/Chai 입력 → 예측 → epitope-recall 채점.
#   판독 = 같은 깊이 5 seed의 예측 에피토프가 (a)서로 비슷=깊이(개수)가 원인·조성 무관
#          (b)draw마다 딴판=특정 서열(조성)이 원인 → Exp2 nested·Exp3 LOCO로.
#   ⚠️채점은 DockQ 아니라 epitope-recall(예측이 진짜 쪽으로 이동)이 정직한 종점.
#
# 사용:
#   cd ~/projects/bk21-msa-depth-bias/pipeline && git pull
#   bash run_track_a_seedrep.sh                     # 기본 앵커 4개
#   ANCHORS="8wpy_AB" bash run_track_a_seedrep.sh   # 한 앵커만
#   REPLICAS=5 OUT=seedrep bash run_track_a_seedrep.sh
# ══════════════════════════════════════════════════════════════════════════════
set -euo pipefail

DATA="${DATA:-/mnt/data/admuser/msadepth}"
LADDIR="$DATA/ladders"
OUT="${OUT:-seedrep}"                 # pipeline/seedrep/<target>_<chain>/
REPLICAS="${REPLICAS:-5}"
# 앵커 = 살아있는 신호(에피토프-이탈)가 실재하는 사례 + 9y0a_AB(DockQ 교차확인, 이미 게이트 통과)
ANCHORS="${ANCHORS:-8wpy_AB 8k3k_D 8k46_I 9y0a_AB}"

pick_depths(){  # neff.tsv(rung n_rows neff80) → full·1 제외 얕은~중간 3개 행수, 쉼표구분
  python3 - "$1" <<'PY'
import sys
rows=[int(l.split()[1]) for l in open(sys.argv[1]).read().splitlines()[1:] if l.split()]
if not rows: print(""); sys.exit()
vals=sorted({r for r in rows if 1 < r < max(rows)}, reverse=True)
if not vals: print(""); sys.exit()
if len(vals)<=3: print(",".join(map(str,vals)))
else:
    idx=sorted({int(len(vals)*f) for f in (0.2,0.5,0.8)})
    print(",".join(str(vals[i]) for i in idx))
PY
}

echo "== Track A Stage 0 seed_replicate | DATA=$DATA | replicas=$REPLICAS =="
for t in $ANCHORS; do
  [ -d "$LADDIR/$t" ] || { echo "skip $t (ladders 없음: $LADDIR/$t)"; continue; }
  for chdir in "$LADDIR/$t"/*/; do
    [ -d "$chdir" ] || continue
    ch=$(basename "$chdir")
    full="$chdir/rung0.a3m"; neff="$chdir/neff.tsv"
    [ -f "$full" ] && [ -f "$neff" ] || { echo "  skip $t/$ch (rung0.a3m 또는 neff.tsv 없음)"; continue; }
    depths=$(pick_depths "$neff")
    [ -n "$depths" ] || { echo "  skip $t/$ch (유효 깊이 없음 — MSA 너무 얕음)"; continue; }
    fullrows=$(sed -n '2p' "$neff" | awk '{print $2}')
    echo "  → $t/$ch  full=${fullrows}행  seedrep 깊이=$depths"
    python seed_replicate.py --a3m "$full" --depths "$depths" \
        --replicas "$REPLICAS" --outdir "$OUT/${t}_${ch}"
  done
done
echo "완료 → pipeline/$OUT/<target>_<chain>/d<depth>/seed{0..N}.a3m + neff.tsv"
echo "다음(GPU, tFold 후): 각 seed a3m → make_input.py → Protenix/Chai 예측 → epitope-recall 채점."
