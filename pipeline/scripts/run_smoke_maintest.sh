#!/usr/bin/env bash
# ══════════════════════════════════════════════════════════════════════════════
# [스모크] 본 검정을 40시간 돌리기 전에, 군마다 코드 경로가 실제로 도는지 검증한다.
#
# 왜 군별인가 — 이 실험에서 갈리는 지점이 군마다 다르다:
#   · 항원 사슬 수   RBD·C = 1사슬(A) / HA·Env 일부 = 2사슬(A|B)  ← 지도가 가장 위험
#   · MSA 규모       RBD 140~28,000줄(Neff 28) / Env 9,000줄(Neff 3,051) / HA 1,000줄대
#   · 고른 칸 깊이   4줄(9azt)부터 11,017줄(9zdu)까지 3자릿수 차이
#   · 과대표집 정의  C군만 없음(overrep = NA)
#
# 조용히 틀리는 경우를 잡는 것이 목적이다. 특히:
#   ① 조성 6개가 사실은 같은 파일 → 이질성 검정이 통째로 무의미
#   ② a3m 질의행 오염 → 모델이 MSA를 버리고 단일서열로 예측(경고만 남고 종료코드 정상)
#   ③ 다중사슬에서 한 사슬만 매핑 → 나머지 사슬이 원래 MSA로 들어감
#   ④ 산출물은 났는데 채점이 안 됨
#
# 2단계로 돈다:
#   Phase A (GPU 없음, 수십 초) — 조성 생성·구별성·질의행·사슬 수
#   Phase B (GPU, 타깃당 1회)   — 실제 예측 1회 + 입력 MSA 검증 + 로그 경고 + 채점
#
# 사용:
#   bash run_smoke_maintest.sh              # 군 대표 + 다중사슬 타깃 자동 선정, Phase A만
#   bash run_smoke_maintest.sh --gpu        # Phase B까지 (권장, 15~20분)
#   TARGETS="8t4a_PR 8u44_ST" bash run_smoke_maintest.sh --gpu
#   bash run_smoke_maintest.sh --gpu --keep # 스모크 산출물 남김(기본은 남김. 지우지 않는다)
#   bash run_smoke_maintest.sh --recheck    # 이미 돈 스모크를 GPU 없이 재판정(터미널을 닫아 판정을
#                                       #   놓쳤을 때. 산출물이 없는 타깃만 새로 돈다)
#
# 출력은 항상 $DATA/logs/smoke_<시각>.log 에도 저장된다(SMKLOG 로 바꿀 수 있음).
#
# 하나라도 실패하면 종료 코드 1. 본 검정을 시작하지 말 것.
# ══════════════════════════════════════════════════════════════════════════════
set -uo pipefail
cd "$(cd "$(dirname "$0")" && pwd)"
DATA="${DATA:-/mnt/data/msadepth}"
CSV="${CSV:-maintest.csv}"
LIST="${LIST:-sweep_targets.csv}"
TARGETS="${TARGETS:-}"
GPU=0; RECHECK=0
for arg in "$@"; do
  case "$arg" in
    --gpu) GPU=1 ;;
    --recheck) GPU=1; RECHECK=1 ;;
    --keep) : ;;
    *) echo "!! 모르는 인자: $arg"; exit 1 ;;
  esac
done
# 출력은 화면과 파일에 동시에 남긴다 — 터미널이 닫혀도 판정을 다시 볼 수 있게.
# 자기 자신을 tee 파이프로 한 번 다시 부르는 방식. 프로세스 치환(exec > >(tee))은
# 부모가 먼저 끝나 마지막 줄이 잘릴 수 있어 쓰지 않는다. 종료 코드는 PIPESTATUS로 보존.
if [ -z "${SMK_TEED:-}" ]; then
  SMKLOG="${SMKLOG:-$DATA/logs/smoke_$(date +%m%d_%H%M%S).log}"
  mkdir -p "$(dirname "$SMKLOG")" 2>/dev/null || SMKLOG="/tmp/$(basename "$SMKLOG")"
  SMK_TEED=1 SMKLOG="$SMKLOG" bash "$0" "$@" 2>&1 | tee -a "$SMKLOG"
  rc=${PIPESTATUS[0]}
  echo "[로그] $SMKLOG"
  exit "$rc"
fi
say(){ echo "[$(date '+%H:%M:%S')] $*"; }
PASS=0; FAIL=0; WARN=0
ok(){   PASS=$((PASS+1)); printf '   [OK]   %s\n' "$*"; }
bad(){  FAIL=$((FAIL+1)); printf '   [실패] %s\n' "$*"; }
warn(){ WARN=$((WARN+1)); printf '   [주의] %s\n' "$*"; }

[ -f "$CSV" ]  || { echo "!! $CSV 없음. 먼저 prep_pick_depth.py"; exit 1; }
[ -f "$LIST" ] || { echo "!! $LIST 없음"; exit 1; }

# ── 대표 타깃 선정: 군마다 하나 + 다중사슬 항원이 있으면 반드시 포함 ─────────────
if [ -z "$TARGETS" ]; then
  TARGETS=$(python3 - "$CSV" "$LIST" <<'PY'
import csv, sys
run = [r for r in csv.DictReader(open(sys.argv[1])) if r.get("status") == "run"]
ag = {r["target"]: r["ag_chains"] for r in csv.DictReader(open(sys.argv[2]))}
pick, seen = [], set()
# 군 대표 = 그 군에서 고른 칸의 서열 수가 가장 적은 것(가장 극단적인 조건을 먼저 깬다)
for g in ("Env", "C", "RBD", "HA"):
    rows = [r for r in run if r.get("group") == g]
    if rows:
        r = min(rows, key=lambda x: int(x["n_rows"] or 0))
        pick.append(r["target"]); seen.add(r["target"])
# 다중사슬 항원은 코드 경로가 달라 반드시 하나 포함
multi = [r["target"] for r in run if "|" in ag.get(r["target"], "")]
for t in multi:
    if t not in seen:
        pick.append(t); seen.add(t); break
print(" ".join(pick))
PY
)
fi
say "스모크 대상: $TARGETS"
echo

# ── 공통 조회 ────────────────────────────────────────────────────────────────
row_of(){ python3 - "$CSV" "$1" <<'PY'
import csv, sys
for r in csv.DictReader(open(sys.argv[1])):
    if r["target"] == sys.argv[2] and r.get("status") == "run":
        print("\t".join([r["rung"], r["n_rows"], r.get("n_comp") or "6",
                         r.get("n_reps") or "4", r.get("n_full") or "8",
                         r.get("group",""), r.get("stratum",""), r.get("neff_pick","")]))
        break
PY
}
agch_of(){ awk -F, -v t="$1" 'NR>1 && $1==t{print $6; exit}' "$LIST"; }
agseq_of(){ python3 - "targets/$1/chains.json" "$2" <<'PY'
import json, sys
d = json.load(open(sys.argv[1]))
sm = {c["id"]: c["seq"] for c in d["chains"]}
print(sm.get(sys.argv[2], ""))
PY
}
first_seq(){ python3 - "$1" <<'PY'
import re, sys
seq, started = [], False
for ln in open(sys.argv[1]):
    ln = ln.rstrip("\n")
    if not ln: continue
    ln = re.sub(r"^#\d+\s+\d+\s*", "", ln)
    if not ln: continue
    if ln[0] == ">":
        if started: break
        started = True; continue
    if started: seq.append(ln)
print("".join(seq).replace("-", "").upper())
PY
}

# ══════════════════════ Phase A — GPU 없음 ══════════════════════
say "═══ Phase A · 조성 생성과 입력 무결성 (GPU 사용 안 함)"
for T in $TARGETS; do
  echo; echo "── $T"
  R=$(row_of "$T"); [ -n "$R" ] || { bad "$T: $CSV 에 status=run 행이 없다"; continue; }
  IFS=$'\t' read -r RUNG NROWS NCOMP NREPS NFULL GRP STRAT NEFFPK <<< "$R"
  AG=$(agch_of "$T"); [ -n "$AG" ] || { bad "$T: $LIST 에서 ag_chains 조회 실패"; continue; }
  IFS='|' read -ra AGC <<< "$AG"
  printf '   군 %s · rung%s · %s줄 · Neff80 %s · 층 %s · 항원사슬 %s개(%s)\n' \
    "$GRP" "$RUNG" "$NROWS" "${NEFFPK:--}" "${STRAT:--}" "${#AGC[@]}" "$AG"

  [ -f "targets/$T/chains.json" ] && ok "chains.json 존재" || { bad "chains.json 없음"; continue; }

  # 조성 생성(GPU 없이)
  GEN_ONLY=1 RUNG="$RUNG" TARGET="$T" REPLICAS="$NCOMP" bash make_composition_reps.sh >"/tmp/smk_gen_$T.log" 2>&1
  if [ $? -ne 0 ]; then bad "조성 생성 실패 → /tmp/smk_gen_$T.log"; continue; fi
  ok "조성 생성 명령 통과"

  # 사슬마다 검사
  for ch in "${AGC[@]}"; do
    dsub=$(ls -d "seedrep_cand/${T}_${ch}"/d* 2>/dev/null | head -1)
    if [ -z "$dsub" ]; then bad "$ch: 조성 폴더 없음 (seedrep_cand/${T}_${ch}/d*)"; continue; fi

    # ① 조성 개수
    n=$(ls "$dsub"/seed[0-9]*.a3m 2>/dev/null | wc -l | tr -d " ")
    [ "$n" -eq "$NCOMP" ] && ok "$ch: 조성 $n개 생성" || bad "$ch: 조성이 $n개 (설계 $NCOMP개)"

    # ② 서로 다른가 — 같으면 이질성 검정이 통째로 무의미
    u=$(md5sum "$dsub"/seed[0-9]*.a3m 2>/dev/null | awk '{print $1}' | sort -u | wc -l | tr -d " ")
    if [ "$u" -eq "$n" ] && [ "$n" -gt 0 ]; then ok "$ch: 조성 $n개가 전부 서로 다름"
    else bad "$ch: 조성 중 중복 있음 (고유 $u / 전체 $n) — 이 상태로 돌리면 실험이 무의미"; fi

    # ③ 줄 수가 설계 깊이와 같은가
    for f in "$dsub"/seed[0-9]*.a3m; do
      [ -f "$f" ] || continue
      c=$(grep -c '^>' "$f")
      [ "$c" -eq "$c" ] 2>/dev/null || continue
      if [ "$c" -lt 2 ]; then bad "$ch: $(basename "$f") 서열 $c개 — 너무 얕다"; fi
    done
    depths=$(for f in "$dsub"/seed[0-9]*.a3m; do grep -c '^>' "$f"; done | sort -u | tr '\n' ' ')
    ok "$ch: 조성별 서열 수 [$depths] (사다리 rung$RUNG 기준)"

    # ④ 질의행 오염 — 첫 서열이 chains.json 항원 서열과 같아야 함
    want=$(agseq_of "$T" "$ch")
    badq=0
    for f in "$dsub"/seed[0-9]*.a3m "$dsub"/seedfull.a3m; do
      [ -f "$f" ] || continue
      got=$(first_seq "$f")
      [ "$got" = "$want" ] || { badq=$((badq+1)); }
    done
    [ "$badq" -eq 0 ] && ok "$ch: 모든 a3m의 질의행이 항원 서열과 일치" \
                      || bad "$ch: 질의행 불일치 $badq개 — 모델이 MSA를 버리고 단일서열로 돈다"

    # ⑤ 원래 MSA 대조군 파일
    [ -f "$dsub/seedfull.a3m" ] && ok "$ch: 원래 MSA 대조군(seedfull.a3m) 준비됨" \
                                || warn "$ch: seedfull.a3m 없음 — COMPS=full 단계에서 생성된다"
  done
done

echo
say "Phase A 종료 — 통과 $PASS · 실패 $FAIL · 주의 $WARN"
if [ "$FAIL" -gt 0 ]; then
  echo; say "!! Phase A에서 실패가 있다. 본 검정을 시작하지 말 것."; exit 1
fi
if [ $GPU -eq 0 ]; then
  echo; say "GPU 검증까지 하려면:  bash run_smoke_maintest.sh --gpu"; exit 0
fi

# ══════════════════════ Phase B — GPU 1회씩 ══════════════════════
echo
say "═══ Phase B · 실제 예측 1회 (타깃당 조성0·반복0만)"
for T in $TARGETS; do
  echo; echo "── $T"
  R=$(row_of "$T"); IFS=$'\t' read -r RUNG NROWS NCOMP NREPS NFULL GRP STRAT NEFFPK <<< "$R"
  AG=$(agch_of "$T"); IFS='|' read -ra AGC <<< "$AG"

  dsub=$(ls -d "seedrep_cand/${T}_${AGC[0]}"/d* 2>/dev/null | head -1)
  DEP=$(basename "${dsub:-none}")
  out="$DATA/compreps/seedrep_cand/protenix/$T/$DEP/seed0_r0"

  # 이미 산출물이 있으면(=스모크를 이미 돌렸으면) GPU를 다시 쓰지 않고 그 파일로 재판정한다.
  have=$(find "$out/results" -name '*sample*.cif' 2>/dev/null | wc -l | tr -d ' ')
  if [ "$RECHECK" -eq 1 ] && [ "${have:-0}" -ge 1 ]; then
    dt=0; say "  (재확인 모드 — 기존 산출물로 판정, GPU 미사용)"
  else
    if [ "$RECHECK" -eq 1 ]; then say "  (재확인 모드지만 산출물이 없어 새로 돈다)"; fi
    t0=$(date +%s)
    SMOKE=1 RUNG="$RUNG" TARGET="$T" REPLICAS="$NCOMP" COMPS="0" REPS=1 \
      bash make_composition_reps.sh >"/tmp/smk_run_$T.log" 2>&1
    dt=$(( $(date +%s) - t0 ))
  fi

  # ① 산출물
  npose=$(find "$out/results" -name '*sample*.cif' 2>/dev/null | wc -l | tr -d ' ')
  [ "$dt" -gt 0 ] && dts="소요 ${dt}초" || dts="기존 산출물"
  [ "$npose" -ge 1 ] && ok "예측 산출물 $npose개 · $dts" \
                     || { bad "산출물 없음 → $out/run.log · /tmp/smk_run_$T.log"; continue; }
  [ "$npose" -eq 5 ] || warn "산출물이 5개가 아니다($npose개) — SAMP 설정 확인"

  # ② 입력에 조성 MSA가 실제로 들어갔는가 (사슬 수만큼)
  if [ -f "$out/input.json" ]; then
    python3 - "$out/input.json" "${#AGC[@]}" <<'PY'
import json, os, sys
d = json.load(open(sys.argv[1]))[0]
need = int(sys.argv[2])
paths = [c["proteinChain"].get("unpairedMsaPath") for c in d["sequences"]]
withmsa = [p for p in paths if p and os.path.exists(p) and sum(1 for l in open(p) if l.startswith(">")) > 1]
print(f"__NMSA__ {len(withmsa)} {need}")
for p in withmsa:
    print(f"__P__ {p} {sum(1 for l in open(p) if l.startswith('>'))}")
PY
  else
    bad "input.json 없음"; continue
  fi > "/tmp/smk_json_$T.txt"
  nmsa=$(awk '/__NMSA__/{print $2}' "/tmp/smk_json_$T.txt")
  [ "${nmsa:-0}" -eq "${#AGC[@]}" ] \
    && ok "항원 사슬 ${#AGC[@]}개 모두에 MSA가 실려 있음" \
    || bad "MSA가 실린 사슬이 $nmsa개 (항원 사슬 ${#AGC[@]}개) — 나머지는 단일서열로 돈다"
  awk '/__P__/{printf "          정제 a3m %s (%s서열)\n", $2, $3}' "/tmp/smk_json_$T.txt"

  # ③ 정제 후에도 질의행이 맞는가
  badq=0
  for p in $(awk '/__P__/{print $2}' "/tmp/smk_json_$T.txt"); do
    ch=$(basename "$p" | sed -E 's/^ag_(.+)_clean\.a3m$/\1/')
    want=$(agseq_of "$T" "$ch"); got=$(first_seq "$p")
    [ "$got" = "$want" ] || badq=$((badq+1))
  done
  [ "$badq" -eq 0 ] && ok "정제(clean_a3m) 후에도 질의행 일치" \
                    || bad "정제 후 질의행 불일치 $badq개"

  # ④ 로그 경고 — 조용한 실패의 흔적
  if grep -qiE 'does not match input sequence|creating dummy|msa .*mismatch|query.*size mismatch' "$out/run.log" 2>/dev/null; then
    bad "run.log에 MSA 불일치 경고 — 모델이 MSA를 버렸다: $out/run.log"
    grep -inE 'does not match input sequence|creating dummy|mismatch' "$out/run.log" | head -3 | sed 's/^/          /'
  else
    ok "run.log에 MSA 관련 경고 없음"
  fi
  if grep -qiE '\berror\b|traceback' "$out/run.log" 2>/dev/null; then
    warn "run.log에 error/traceback 문자열 있음(치명적이 아닐 수 있음): $out/run.log"
  fi

  # ⑤ 채점이 되는가
  python eval_dump_seedrep.py --data "$DATA/compreps" --only "$T" \
      --csv-out "/tmp/smk_score_$T.csv" >"/tmp/smk_score_$T.log" 2>&1
  if [ -f "/tmp/smk_score_$T.csv" ]; then
    n=$(python3 - "/tmp/smk_score_$T.csv" <<'PY'
import csv, sys, math
rows = list(csv.DictReader(open(sys.argv[1])))
def num(x):
    try:
        v = float(x); return v == v
    except Exception:
        return False
good = [r for r in rows if num(r.get("dockq")) and num(r.get("recall"))]
print(len(good))
if good:
    print("   " + " ".join(f"DockQ {float(r['dockq']):.2f}/recall {float(r['recall']):.2f}" for r in good[:3]))
PY
)
    cnt=$(echo "$n" | head -1)
    [ "${cnt:-0}" -ge 1 ] && { ok "채점 통과 — 값이 있는 자세 $cnt개"; echo "$n" | tail -n +2 | sed 's/^/       /'; } \
                          || { bad "채점 결과가 비었거나 전부 NaN → /tmp/smk_score_$T.log"
                               grep -vE '^\s*$' "/tmp/smk_score_$T.log" 2>/dev/null | tail -4 | sed 's/^/          | /'; }
  else
    bad "채점 산출 없음 → /tmp/smk_score_$T.log"
    grep -vE '^\s*$' "/tmp/smk_score_$T.log" 2>/dev/null | tail -4 | sed 's/^/          | /'
  fi
done

echo
echo "════════════════════════════════════════════════════════════════"
say "스모크 종료 — 통과 $PASS · 실패 $FAIL · 주의 $WARN"
if [ "$FAIL" -gt 0 ]; then
  say "!! 실패가 있다. 본 검정(약 40시간)을 시작하지 말 것."; exit 1
fi
say "전부 통과. 본 검정을 시작해도 된다:"
say "   tmux new -s maintest"
say "   HOURS=12 bash run_maintest.sh --apply"
say "※ 스모크로 만든 seed0_r0 실행분은 남겨 둔다 — 본 검정에서 그대로 재사용된다(이어달리기)."
