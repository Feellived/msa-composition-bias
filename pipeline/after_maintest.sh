#!/usr/bin/env bash
# ══════════════════════════════════════════════════════════════════════════════
# [본 검정 뒤처리] 예측이 끝나기를 기다렸다가 채점 → 후보 자리 → (선택)데모 훑기까지 한 번에.
#
# run_until_done.sh 와 다른 점:
#   · 예측을 다시 돌리지 않는다(예산이 충분할 때 쓴다). 끝나기를 기다리기만 한다.
#   · ⭐ 채점할 때 **설계값과 실행 수가 다른 타깃은 강제로 다시 채점한다**.
#     analyze_target.sh 는 CSV 가 있으면 건너뛰므로, 옛 실험의 낡은 CSV 가 남아 있으면
#     새 예측이 끝나도 옛 통계가 그대로 살아남는다(9azr_HL 이 실제로 그 상태였다).
#     설계값 = maintest.csv 의 n_comp × n_reps + n_full.
#   · 끝나면 데모가 성립하는 복합체(고른 후보 ≠ 원래 MSA 자리)를 표로 뽑아 준다.
#
# 사용:
#   bash after_maintest.sh                     # 무엇을 할지만 출력(dry-run)
#   bash after_maintest.sh --apply             # 대기 → 채점 → 후보 자리
#   SCREEN=1 bash after_maintest.sh --apply    # 위 + boltz 'ours' 훑기까지
#   NOWAIT=1 bash after_maintest.sh --apply    # 기다리지 않고 지금 끝난 것만 처리
#
# env: DATA · CSV(maintest.csv) · WAIT_SEC · SCREEN · CD(consensus_docking 경로) · NOWAIT
# 로그는 $DATA/logs/after_maintest_<시각>.log 에도 남는다.
# ══════════════════════════════════════════════════════════════════════════════
set -uo pipefail
cd "$(cd "$(dirname "$0")" && pwd)"
DATA="${DATA:-/mnt/data/admuser/msadepth}"
CSV="${CSV:-maintest.csv}"
WAIT_SEC="${WAIT_SEC:-300}"
SCREEN="${SCREEN:-0}"
NOWAIT="${NOWAIT:-0}"
CD="${CD:-$HOME/projects/bk21-antibody-ml/consensus_docking}"
# ⚠️ 채점에는 DockQ 도구가 필요하다. tmux 가 물려받은 env 에 없으면 dockq 열만 조용히 비고
#    자리겹침만 나온다(2026-07-30 실제로 그랬다). 그래서 env 를 명시해 부른다.
DOCKQ_ENV="${DOCKQ_ENV:-boltz}"
# ⚠️ 설계가 일괄 규격과 다른 복합체는 실행 수 대조에서 빼야 한다.
#    8ulr_HL = 확정 실험(조성 8 × 반복 5 + 예산 맞춘 통제 → 42회). maintest.csv 에는
#    일괄 설계(32)가 적혀 있어 '낡음'으로 오탐된다. 재채점해도 같은 자료를 다시 읽을 뿐이지만
#    확정 결과 파일을 건드릴 이유가 없다.
SPECIAL="${SPECIAL:-8ulr_HL}"
APPLY=0
for a in "$@"; do case "$a" in --apply) APPLY=1 ;; *) echo "!! 모르는 인자: $a"; exit 1 ;; esac; done

if [ -z "${AM_TEED:-}" ]; then
  AMLOG="${AMLOG:-$DATA/logs/after_maintest_$(date +%m%d_%H%M%S).log}"
  mkdir -p "$(dirname "$AMLOG")" 2>/dev/null || AMLOG="/tmp/$(basename "$AMLOG")"
  echo "[로그] $AMLOG   ← 다른 창에서:  tail -f $AMLOG"
  AM_TEED=1 AMLOG="$AMLOG" bash "$0" "$@" 2>&1 | tee -a "$AMLOG"
  rc=${PIPESTATUS[0]}; echo "[로그] $AMLOG"; exit "$rc"
fi
say(){ echo "[$(date '+%m-%d %H:%M:%S')] [뒤처리] $*"; }
[ -f "$CSV" ] || { say "!! $CSV 없음"; exit 1; }
# ── 채점 환경 ────────────────────────────────────────────────────────────────
# dockq_sweep.dockq() 는 subprocess.run(["DockQ", ...]) 로 PATH 의 실행파일을 부르고,
# 없으면 except 로 None 을 돌려준다 → dockq 열만 조용히 비고 자리겹침만 나온다
# (2026-07-30 실제로 그랬다. 종료코드 5 가드가 CSV 를 안 써서 손실은 없었다).
# conda run 은 PATH 처리에 변덕이 있어 comp_x_reps.sh 와 같은 방식(source + activate)을 쓴다.
source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null
conda activate "$DOCKQ_ENV" 2>/dev/null
if ! command -v DockQ >/dev/null 2>&1; then
  say "!! env '$DOCKQ_ENV' 에서 DockQ 실행파일을 못 찾았다 — 채점하면 dockq 열이 통째로 빈다."
  say "   확인: conda activate $DOCKQ_ENV && which DockQ"
  say "   다른 env 면: DOCKQ_ENV=<이름> bash after_maintest.sh --apply"
  exit 1
fi
say "채점 env = $DOCKQ_ENV · DockQ = $(command -v DockQ)"

# ── 1) 예측이 끝나기를 기다린다 ───────────────────────────────────────────────
busy(){ pgrep -f "comp_x_reps.sh" >/dev/null || pgrep -f "run_maintest.sh --apply" >/dev/null; }
if [ "$NOWAIT" != "1" ] && busy; then
  say "본 검정이 돌고 있다 — 끝나면 이어서 한다 (${WAIT_SEC}초마다 확인)"
  [ $APPLY -eq 1 ] || { say "(dry-run 이라 기다리지 않고 지금 상태로 계속)"; }
  if [ $APPLY -eq 1 ]; then
    while busy; do sleep "$WAIT_SEC"; done
    say "예측 종료 감지 — 30초 뒤 시작(마지막 파일 쓰기 여유)"; sleep 30
  fi
fi

# ── 2) 채점 — 설계값과 실행 수가 다르면 강제 재채점 ───────────────────────────
say "채점 대상 판정 (설계값 = 조성 수 × 반복 수 + 원래 MSA 횟수)"
PLAN=$(SPECIAL="$SPECIAL" python3 - "$CSV" <<'PY'
import csv, os, sys
special = set(os.environ.get("SPECIAL", "").split())
for r in csv.DictReader(open(sys.argv[1])):
    if r.get("status") != "run":
        continue
    t = r["target"]
    want = int(r.get("n_comp") or 6) * int(r.get("n_reps") or 4) + int(r.get("n_full") or 8)
    p = f"results/compreps_{t}.csv"
    have, deps = 0, 0
    if os.path.exists(p):
        rows = list(csv.DictReader(open(p)))
        have = len({x["seed"] for x in rows})
        deps = len({x["depth"] for x in rows})
    # ⚠️ compreps CSV 만 보면 안 된다. analyze_target.sh 는 채점(dump) → 통계 → 결합자리 순으로
    #    돌고, 중간에 끊기면 compreps CSV 는 온전한데 하위 산출물이 없는 상태가 된다.
    #    그 상태를 'skip' 으로 판정하면 그 복합체는 영원히 채점되지 않고 조용히 빠진다.
    down = [f"results/summary_{t}_recall.csv", f"results/summary_{t}_dockq.csv",
            f"results/site_repro_{t}.csv", f"results/sites_{t}.json"]
    missing_down = [os.path.basename(x) for x in down if not os.path.exists(x)]
    if not os.path.exists(p):          act = "new"
    elif t in special:                 act = "special"   # 설계가 달라 실행 수 대조를 안 한다
    elif have != want or deps != 1:    act = "redo"
    elif missing_down:                 act = "redo"      # 하위 산출물이 빠졌다(중단된 흔적)
    else:                              act = "skip"
    print(f"{t}\t{act}\t{have}\t{want}\t{deps}\t{','.join(missing_down)}")
PY
)
echo "$PLAN" | awk -F'\t' 'BEGIN{printf "  %-12s %-8s %s\n","타깃","조치","실행/설계·깊이수"}
  {n=($2=="redo"?($6!=""?"   ← 하위 산출물 없음("$6") 재채점":"   ← 낡음, 강제 재채점"):($2=="special"?"   ← 설계가 달라 그대로 둠":""));
   printf "  %-12s %-8s %s/%s·%s%s\n",$1,$2,$3,$4,$5,n}'

if [ $APPLY -eq 0 ]; then
  say "dry-run 끝. 실제로 하려면 --apply"
  exit 0
fi

TODO=$(echo "$PLAN" | awk -F'\t' '$2=="new"||$2=="redo"' | wc -l | tr -d ' ')
say "채점할 복합체 ${TODO}개 (하나에 5~8분 걸린다 — DockQ 를 자세 160개에 대해 돌린다)"
T0=$(date +%s); I=0; NEW=0; REDONE=0; FAIL=0
while IFS=$'\t' read -r t act have want deps miss; do
  [ -n "$t" ] || continue
  case "$act" in skip|special) continue ;; esac
  I=$((I+1)); TS=$(date +%s)
  if [ "$act" = "new" ]; then
    say "[$I/$TODO] $t 채점 시작 (새로) · 경과 $(( (TS-T0)/60 ))분"
    bash analyze_target.sh "$t" </dev/null && NEW=$((NEW+1)) || FAIL=$((FAIL+1))
  else
    say "[$I/$TODO] $t 재채점 시작 (실행 $have/설계 $want · 깊이 $deps${miss:+ · 없던 파일 $miss}) · 경과 $(( (TS-T0)/60 ))분"
    REDO=1 bash analyze_target.sh "$t" </dev/null && REDONE=$((REDONE+1)) || FAIL=$((FAIL+1))
  fi
  say "[$I/$TODO] $t 끝 ($(( ($(date +%s)-TS)/60 ))분 걸림)"
done <<< "$PLAN"
say "채점 끝 — 새로 $NEW · 재채점 $REDONE · 실패 $FAIL · 총 $(( ($(date +%s)-T0)/60 ))분"

# ── 3) 데모가 성립하는 복합체 목록 ────────────────────────────────────────────
say "데모 성립 여부 (고른 후보 ≠ 원래 MSA 자리여야 대조가 성립)"
OKLIST=$(python3 - <<'PY'
import json, glob
ok = []
for f in sorted(glob.glob("results/sites_*.json")):
    d = json.load(open(f)); cs = d["candidates"]
    full = next((c for c in cs if c["from_full_msa"]), None)
    pick = max(cs, key=lambda c: c["n_comp"])
    good = len(cs) > 1 and full and full["cand"] != pick["cand"]
    print(("✅" if good else "  "), f"{d['target']:11}", f"후보{len(cs)}",
          f"고른것 {pick['cand']}(조성{pick['n_comp']}·{len(pick['residues'])}잔기)",
          f"원래MSA {full['cand'] if full else '-'}")
    if good: ok.append(d["target"])
open("results/demo_candidates.txt", "w").write(" ".join(ok) + "\n")
PY
)
echo "$OKLIST"
TGTS=$(cat results/demo_candidates.txt 2>/dev/null)
say "데모 후보: ${TGTS:-없음}  (results/demo_candidates.txt)"

# ── 4) (선택) boltz 훑기 — 반드시 noconstraint 와 함께 ────────────────────────
#   ⚠️ 'ours' 만 훑으면 안 된다. 2026-07-29 에 8sit_HL 이 ours 에서 DockQ 0.787 로 대박처럼
#      보였는데, 제약 없는 팔이 0.818 이었다 = boltz 가 원래 혼자 푸는 복합체였다.
#      제약 없는 팔보다 나아야 우리 기여다.
if [ "$SCREEN" = "1" ] && [ -n "$TGTS" ]; then
  if [ -d "$CD" ]; then
    say "boltz 훑기 시작(noconstraint + ours) (타깃 $(echo $TGTS | wc -w)개)"
    ( cd "$CD" && ARMS="noconstraint ours" MODELS=boltz bash scripts/run_demo_guided.sh "$TGTS" ) || say "! 훑기 중 오류"
    ( cd "$CD" && python scripts/dockq_demo.py --targets "$TGTS" --models boltz ) || say "! 채점 중 오류"
  else
    say "! $CD 없음 — 훑기 건너뜀"
  fi
else
  say "훑기는 하지 않았다. 하려면:"
  say "  cd $CD && ARMS=\"noconstraint ours\" MODELS=boltz bash scripts/run_demo_guided.sh \"$TGTS\""
fi
say "완료."
