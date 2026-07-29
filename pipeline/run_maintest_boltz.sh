#!/usr/bin/env bash
# ══════════════════════════════════════════════════════════════════════════════
# [선택 실행] 지정한 복합체만 **boltz** 로 본 검정과 같은 설계로 돌린다.
#
# ⚠️ 이것은 연구 주장용이 아니다. run_maintest.sh(protenix, 사전 확정 명단 29종 전수)가
#    주 결과이고, 이 스크립트는 **데모용으로 고른 복합체**를 boltz 로도 돌려 보기 위한 것이다.
#    결과를 빈도·재현율로 보고하려면 명단을 결과와 무관하게 정해야 한다(6.4 규율).
#    골라서 돌린 것은 발표·문서에서 "골랐다"를 반드시 함께 적을 것.
#
# 설계는 maintest.csv 를 그대로 따른다 — 깊이(rung)·조성 수·반복 수·원래 MSA 횟수를
# 여기서 바꾸지 않는다. 모델만 boltz 로 바꾼다.
#
# ⭐ 조성 a3m 은 **다시 만들지 않는다.** comp_x_reps.sh 가 없을 때만 만들고, protenix 실행 때
#    이미 만들어 뒀다. 즉 boltz 는 **같은 조성 목록**을 받는다 → 모델 비교가 깨끗해진다.
#    출력은 seedrep_cand/<모델>/ 아래로 갈리므로 protenix 결과와 섞이지 않는다.
#
# 사용:
#   bash run_maintest_boltz.sh                      # 무엇을 돌릴지만 출력(dry-run)
#   bash run_maintest_boltz.sh --apply 8k3k_D
#   bash run_maintest_boltz.sh --apply "8k3k_D 8sit_HL 8t4d_OQ"
#   HOURS=8 bash run_maintest_boltz.sh --apply "..."
#
# 끝나면 채점(모델을 boltz 로 지정해야 한다):
#   python dump_seedrep_full.py --data $DATA/compreps --only <타깃> \
#          --csv-out results/compreps_<타깃>_boltz.csv
#   ※ dump 는 maintest.csv 의 model(protenix)을 보므로, boltz 결과를 채점하려면
#     --cand 로 model=boltz 인 한 줄짜리 CSV 를 주는 것이 가장 확실하다. 아래 안내가 나온다.
# ══════════════════════════════════════════════════════════════════════════════
set -uo pipefail
cd "$(cd "$(dirname "$0")" && pwd)"
DATA="${DATA:-/mnt/data/admuser/msadepth}"
CSV="${CSV:-maintest.csv}"
HOURS="${HOURS:-8}"
APPLY=0; TARGETS=""
for arg in "$@"; do
  case "$arg" in
    --apply) APPLY=1 ;;
    -*) echo "!! 모르는 인자: $arg"; exit 1 ;;
    *) TARGETS="$TARGETS $arg" ;;
  esac
done
TARGETS="$(echo "$TARGETS" | xargs)"
say(){ echo "[$(date '+%m-%d %H:%M:%S')] $*"; }
[ -f "$CSV" ] || { say "!! $CSV 없음"; exit 1; }
[ -n "$TARGETS" ] || { say "사용: bash run_maintest_boltz.sh [--apply] \"<타깃> [타깃...]\""; exit 1; }

START=$(date +%s); BUDGET=$(python3 -c "print(int(float('$HOURS')*3600))")
left(){ echo $(( BUDGET - ($(date +%s) - START) )); }

# ── 항원 사슬 명단 준비 ───────────────────────────────────────────────────────
# comp_x_reps.sh 는 항원 사슬을 $LIST(기본 sweep_targets.csv)의 ag_chains 열에서 읽고,
# 없으면 그 자리에서 종료한다. 그런데 세트 3 복합체(8sit_HL·8siq_HL·8sis_HL·8xsi_HL)는
# sweep_targets.csv 에 없다 → 그대로 부르면 즉시 죽는다.
# 그래서 필요한 행만 모은 임시 명단을 만들어 넘긴다. 없는 행은 targets/<타깃>/chains.json 의
# role=="antigen" 에서 직접 만든다(같은 정보의 다른 출처라 값이 달라질 여지가 없다).
LIST_SRC="${LIST:-sweep_targets.csv}"
TMPLIST="$(mktemp -t maintest_boltz_list.XXXXXX.csv)"
trap 'rm -f "$TMPLIST"' EXIT
MISSING=$(LIST_SRC="$LIST_SRC" TMPLIST="$TMPLIST" TARGETS="$TARGETS" python3 <<'PY'
import csv, json, os
src, out, tg = os.environ["LIST_SRC"], os.environ["TMPLIST"], os.environ["TARGETS"].split()
have = {}
if os.path.exists(src):
    for r in csv.DictReader(open(src)):
        have[r["target"]] = r
cols = ["target", "pdb", "group", "ab", "dirtype", "ag_chains", "label"]
miss = []
with open(out, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore"); w.writeheader()
    for t in tg:
        if t in have:
            w.writerow(have[t]); continue
        cj = os.path.join("targets", t, "chains.json")
        if not os.path.exists(cj):
            miss.append(t); continue
        d = json.load(open(cj))
        ag = [str(c["id"]) for c in d.get("chains", []) if c.get("role") == "antigen"]
        if not ag:
            miss.append(t); continue
        w.writerow(dict(target=t, pdb=t.split("_")[0], group="", ab="",
                        dirtype="targets", ag_chains="|".join(ag), label=""))
print(" ".join(miss))
PY
)
[ -z "$MISSING" ] || say "!! 항원 사슬을 못 찾은 복합체(건너뜀): $MISSING"

say "$([ $APPLY -eq 1 ] && echo '[실제 실행]' || echo '[dry-run — 아무것도 실행하지 않음]') 모델 boltz · 예산 ${HOURS}시간"
echo
printf '%-13s %-7s %-8s %-14s %s\n' 타깃 칸 서열수 설계 상태
printf -- '-%.0s' {1..64}; echo

for T in $TARGETS; do
  ROW=$(python3 - "$CSV" "$T" <<'PY'
import csv, sys
for r in csv.DictReader(open(sys.argv[1])):
    if r.get("target") == sys.argv[2] and r.get("status") == "run":
        print("\t".join([str(r["rung"]), str(r.get("n_rows","")), str(r.get("n_comp") or 6),
                         str(r.get("n_reps") or 4), str(r.get("n_full") or 8)])); break
PY
)
  if [ -z "$ROW" ]; then
    printf '%-13s %s\n' "$T" "!! maintest.csv 에 status=run 인 행이 없다 — 건너뜀"; continue
  fi
  IFS=$'\t' read -r RUNG NROWS NCOMP NREPS NFULL <<< "$ROW"
  grep -q "^$T," "$TMPLIST" || { printf '%-13s %s\n' "$T" "!! 항원 사슬 미확인 — 건너뜀"; continue; }
  WANT=$(( NCOMP * NREPS + NFULL ))
  # 이미 끝난 실행 수(boltz 쪽만)
  BASE="$DATA/compreps/seedrep_cand/boltz/$T"
  HAVE=0
  if [ -d "$BASE" ]; then
    for d in "$BASE"/d*/seed*_r*/; do
      [ -d "$d" ] || continue
      find "$d/results" -name '*_model_*.cif' 2>/dev/null | grep -q . && HAVE=$((HAVE+1))
    done
  fi
  if [ "$HAVE" -ge "$WANT" ]; then ST="완료 ($HAVE/$WANT) — 건너뜀"
  elif [ "$HAVE" -gt 0 ]; then ST="이어서 ($HAVE/$WANT)"
  else ST="새로 (${WANT}회)"; fi
  printf '%-13s %-7s %-8s %-14s %s\n' "$T" "rung$RUNG" "$NROWS" "${NCOMP}×${NREPS}+${NFULL}" "$ST"
  [ $APPLY -eq 1 ] || continue
  [ "$HAVE" -ge "$WANT" ] && continue

  REM=$(left)
  if [ "$REM" -lt 600 ]; then say "  예산이 10분 미만 남았다 — 여기서 멈춘다"; break; fi
  COMPLIST=$(python3 -c "print(' '.join(str(i) for i in range($NCOMP)))")

  say "  ① 조성 $NCOMP가지 × 반복 $NREPS회 (boltz)"
  MODEL=boltz LIST="$TMPLIST" RUNG="$RUNG" TARGET="$T" REPLICAS="$NCOMP" COMPS="$COMPLIST" REPS="$NREPS" \
    bash comp_x_reps.sh </dev/null || say "  ! 조성 실행에서 오류 — 로그 확인"

  say "  ② 원래 MSA $NFULL회 (boltz)"
  MODEL=boltz LIST="$TMPLIST" RUNG="$RUNG" TARGET="$T" COMPS="full" REPS="$NFULL" \
    bash comp_x_reps.sh </dev/null || say "  ! 원래 MSA 실행에서 오류 — 로그 확인"
done

echo
if [ $APPLY -eq 1 ]; then
  say "끝. 채점은 모델을 boltz 로 지정해서 해야 한다 — 타깃마다:"
  for T in $TARGETS; do
    echo "  printf 'target,model,peak_rung,replicas,obs_dq,obs_rec\\n${T},boltz,0,8,,\\n' > cand_boltz_${T}.csv"
    echo "  python dump_seedrep_full.py --data \$DATA/compreps --only $T --cand-first \\"
    echo "         --cand cand_boltz_${T}.csv --csv-out results/compreps_${T}_boltz.csv"
  done
  say "그다음 site_reproducibility.py --csv results/compreps_<타깃>_boltz.csv --dump-sites 로 후보 자리를 뽑는다."
else
  say "실제로 돌리려면 --apply 를 붙일 것."
fi
