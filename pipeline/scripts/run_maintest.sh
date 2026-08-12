#!/usr/bin/env bash
# ══════════════════════════════════════════════════════════════════════════════
# [본 검정 실행기] prep_pick_depth.py가 만든 maintest.csv를 읽어 타깃마다
#   조성 N가지 × 반복 M회  +  원래 MSA F회
# 를 make_composition_reps.sh로 돌린다. 설계값(N·M·F)과 깊이는 CSV에 이미 적혀 있으므로
# 여기서 바꾸지 않는다 — 판정 기준을 사후에 옮겼다는 반박을 받지 않기 위함.
#
# 순서(2026-07-28 확정): Env → C → RBD → HA.  근거는 아래.
#   ① Env — 확정 사례 8ulr과 같은 항원 계열이고 조성 다양성 regime도 같다(Neff80 수천).
#      "8ulr이 재현되는가"를 가장 직접 묻는 자리라 먼저 본다.
#   ② C  — 과대표집되지 않은 항원. 효과가 과대표집 항원에 국한되는지(적용 범위)를 가른다.
#   ③ RBD — 조성 다양성 없음 층 전체(Neff80 rung0 ≈ 28). 빠르고, 여기까지 오면 한 층이
#      통째로 완결되며 세트 3의 편향 없는 빈도(M/10)도 확보된다.
#   ④ HA — 항원이 커서 가장 느리다.
#   ORDER 환경변수로 뒤집을 수 있다(예: ORDER="RBD C Env HA").
#
# 안전장치
#   · 기본 dry-run. 실제 실행은 --apply (GPU를 하루 단위로 점유하므로).
#   · HOURS 예산 — 타깃 하나를 시작하기 전에 남은 시간을 확인하고, 모자라면 멈춘다.
#     (타깃 중간에 죽이지 않는다. make_composition_reps.sh 자체가 이어달리기라 다음에 이어서 감)
#   · prep_a3m_check_match.py 게이트 — MSA 질의행이 어긋난 채로 조용히 단일서열 예측이
#     되는 사고(2026-07-27 boltz)를 막는다. --skip-gate 로만 건너뛸 수 있다.
#   · 이미 끝난 타깃은 건너뛴다(출력 폴더의 실행 수를 세어 판단).
#
# 사용:
#   bash run_maintest.sh                      # 무엇을 돌릴지만 출력
#   bash run_maintest.sh --apply              # 실제 실행 (기본 예산 12시간)
#   HOURS=6 bash run_maintest.sh --apply
#   ONLY="8sis_HL 9zdu_HL" bash run_maintest.sh --apply
#   bash run_maintest.sh --apply --skip-gate  # 게이트 생략(권장하지 않음)
#
# 출력은 항상 $DATA/logs/maintest_<시각>.log 에도 저장된다(MTLOG 로 바꿀 수 있음).
# tmux가 죽어도 이 파일로 어디까지 갔는지 확인할 수 있고, 같은 명령을 다시 넣으면
# 완료된 타깃을 건너뛰고 이어서 간다.
#
# 끝나면 타깃마다:  bash run_analyze_target.sh <타깃>
# ══════════════════════════════════════════════════════════════════════════════
set -uo pipefail
cd "$(cd "$(dirname "$0")" && pwd)"
DATA="${DATA:-/mnt/data/msadepth}"
CSV="${CSV:-maintest.csv}"
LIST="${LIST:-sweep_targets.csv}"
HOURS="${HOURS:-12}"
ONLY="${ONLY:-}"
ORDER="${ORDER:-Env C RBD HA}"
APPLY=0; GATE=1
# 실행 1회에 걸리는 시간(초). 첫 타깃은 이 값으로 예산을 판단하고, 이후에는 실측으로 갱신.
EST="${EST:-180}"
OVERRUN="${OVERRUN:-0}"   # 1이면 예산이 모자라도 다음 타깃을 시작한다(예전 동작)
for arg in "$@"; do
  case "$arg" in
    --apply) APPLY=1 ;;
    --skip-gate) GATE=0 ;;
    *) echo "!! 모르는 인자: $arg"; exit 1 ;;
  esac
done
# 출력을 화면과 파일에 동시에 남긴다. 40시간짜리라 tmux가 죽으면 화면 기록이 통째로
# 날아간다(계산 결과는 디스크에 남지만, 어느 타깃이 왜 건너뛰어졌는지가 사라진다).
# 자기 자신을 tee 파이프로 다시 부르는 방식 — 마지막 줄까지 확실히 저장되고 종료 코드도 보존.
if [ -z "${MT_TEED:-}" ]; then
  MTLOG="${MTLOG:-$DATA/logs/maintest_$(date +%m%d_%H%M%S).log}"
  mkdir -p "$(dirname "$MTLOG")" 2>/dev/null || MTLOG="/tmp/$(basename "$MTLOG")"
  MT_TEED=1 MTLOG="$MTLOG" bash "$0" "$@" 2>&1 | tee -a "$MTLOG"
  rc=${PIPESTATUS[0]}
  echo "[로그] $MTLOG"
  exit "$rc"
fi
say(){ echo "[$(date '+%m-%d %H:%M:%S')] $*"; }

[ -f "$CSV" ] || { say "!! $CSV 없음. 먼저 python prep_pick_depth.py 를 돌릴 것"; exit 1; }

START=$(date +%s)
BUDGET=$(python3 -c "print(int(float('$HOURS')*3600))")
left(){ echo $(( BUDGET - ($(date +%s) - START) )); }

# ── CSV에서 status=run 인 행만, RBD → 음성대조 → 나머지 순으로 정렬 ──────────────
ROWS=()
while IFS= read -r __ln; do ROWS+=("$__ln"); done < <(python3 - "$CSV" "$ORDER" <<'PY'
import csv, sys
order = sys.argv[2].split()
rows = [r for r in csv.DictReader(open(sys.argv[1])) if r.get("status") == "run"]
def key(r):
    g = r.get("group", "")
    return (order.index(g) if g in order else len(order), r["target"])
for r in sorted(rows, key=key):
    print("\t".join([r["target"], r.get("group",""), r.get("model","protenix"),
                     str(r["rung"]), str(r.get("n_rows","")),
                     str(r.get("n_comp") or 6), str(r.get("n_reps") or 4),
                     str(r.get("n_full") or 8), r.get("stratum", ""), r.get("neff_pick", "")]))
PY
)
[ "${#ROWS[@]}" -gt 0 ] || { say "!! $CSV 에 status=run 인 행이 없다"; exit 1; }

# ── 이미 끝난 실행 수 세기 ────────────────────────────────────────────────────
depth_of(){   # $1=target $2=칸번호 → 깊이 폴더 이름 d<서열수>
  # ⚠️ maintest.csv 의 서열 수(neff.tsv 유래)와 실제 폴더 이름이 어긋날 수 있다.
  #    폴더는 make_composition_reps.sh 가 a3m 의 '>' 줄을 직접 세어 만든다(8sit_HL 에서 4034 대 4035
  #    로 하나 차이가 났고, 그 바람에 다 끝난 32회를 0회로 세었다).
  #    → 여기서도 생성기와 똑같은 방법으로 구해 출처를 하나로 맞춘다.
  local t="$1" r="$2" ch f
  ch=$(python3 - "$LIST" "$t" <<'DEPTHPY'
import csv, sys
for row in csv.DictReader(open(sys.argv[1])):
    if row.get("target") == sys.argv[2]:
        print((row.get("ag_chains") or "A").split("|")[0]); break
DEPTHPY
)
  [ -n "$ch" ] || ch="A"
  f="$DATA/ladders/$t/$ch/rung${r}.a3m"
  [ -f "$f" ] || { echo ""; return; }
  echo "d$(grep -c '^>' "$f" | tr -d ' ')"
}

done_runs(){   # $1=target $2=model $3=깊이폴더(d<서열수>) → 그 깊이의 산출물 있는 실행 수
  # ⚠️ 예전에는 "$base"/*/ 로 모든 깊이를 셌다. 같은 타깃 아래 이전 실험의 다른 깊이
  #    폴더가 남아 있으면 진행률이 부풀고, 그 수가 목표에 닿으면 이번 실행을 통째로
  #    건너뛴다(2026-07-28 발견). 이번 검정의 깊이 폴더 d<서열수>만 센다.
  local t="$1" m="$2" dep="$3" base n=0 d
  base="$DATA/compreps/seedrep_cand/$m/$t"
  [ -n "$dep" ] && [ -d "$base/$dep" ] || { echo 0; return; }
  for d in "$base/$dep"/seed*_r*/; do
    [ -d "$d" ] || continue
    find "$d/results" -name '*sample*.cif' -o -name '*_model_*.cif' 2>/dev/null | grep -q . && n=$((n+1))
  done
  echo "$n"
}

echo
if [ $APPLY -eq 1 ]; then say "[실제 실행] 예산 ${HOURS}시간"; else say "[dry-run — 아무것도 실행하지 않음]"; fi
printf '%-13s %-4s %-8s %-7s %-8s %-8s %-11s %-13s %s\n' \
  타깃 군 모델 칸 서열수 Neff80 층 "설계" 상태
printf -- '-%.0s' {1..108}; echo

PLAN=(); TOTAL=0
for row in "${ROWS[@]}"; do
  IFS=$'\t' read -r t grp model rung nrows ncomp nreps nfull stratum neffpk <<< "$row"
  if [ -n "$ONLY" ] && [[ " $ONLY " != *" $t "* ]]; then continue; fi
  want=$(( ncomp * nreps + nfull ))
  dep=$(depth_of "$t" "$rung")
  if [ -n "$dep" ] && [ "$dep" != "d${nrows}" ]; then
    say "  ※ $t: 계획표 서열 수 ${nrows} 와 실제 깊이 폴더 ${dep} 가 다르다 — 실제 폴더를 따른다."
  fi
  have=$(done_runs "$t" "$model" "$dep")
  if [ "$have" -ge "$want" ]; then st="완료 ($have/$want) — 건너뜀"
  elif [ "$have" -gt 0 ]; then st="이어서 ($have/$want)"; PLAN+=("$row"); TOTAL=$((TOTAL+want-have))
  else st="새로 ($want회)"; PLAN+=("$row"); TOTAL=$((TOTAL+want)); fi
  case "$stratum" in rich) sl="다양성있음";; poor) sl="다양성없음";; *) sl="-";; esac
  printf '%-13s %-4s %-8s %-7s %-8s %-8s %-11s %-13s %s\n' \
    "$t" "$grp" "$model" "rung$rung" "$nrows" "${neffpk:--}" "$sl" "${ncomp}×${nreps}+${nfull}" "$st"
done
echo
say "돌릴 타깃 ${#PLAN[@]}개 · 남은 실행 약 ${TOTAL}회 · 순서 [$ORDER]"
say "※ 층은 성적이 아니라 입력의 Neff80으로 나눈 것이다. 층별 빈도를 따로 내되,"
say "  이 자료에서는 층이 항원 계열과 거의 겹치므로 '조성 다양성 때문'과 '항원 계열 때문'을 분리할 수 없다."
[ "${#PLAN[@]}" -gt 0 ] || { say "전부 완료 상태다. run_analyze_target.sh 로 넘어갈 것."; exit 0; }

if [ $APPLY -eq 0 ]; then
  echo
  say "실제로 돌리려면:  bash run_maintest.sh --apply"
  say "끝난 뒤 분석:     for t in <타깃들>; do bash run_analyze_target.sh \$t; done"
  exit 0
fi

# ── 실행 ─────────────────────────────────────────────────────────────────────
n_done=0; n_skip=0
for row in "${PLAN[@]}"; do
  IFS=$'\t' read -r t grp model rung nrows ncomp nreps nfull stratum neffpk <<< "$row"
  rem=$(left)
  # 예산은 타깃 경계에서만 본다(타깃을 중간에 끊지 않는다). 그래서 "10분 넘게 남았으면
  # 시작" 같은 기준으로는 77분짜리 타깃을 38분 남기고 시작해 예산을 크게 넘긴다.
  # → 이 타깃에 실제로 필요한 시간을 어림해서, 모자라면 시작하지 않는다.
  dep=$(depth_of "$t" "$rung")
  have0=$(done_runs "$t" "$model" "$dep")
  want0=$(( ncomp * nreps + nfull ))
  left_runs=$(( want0 - have0 ))
  [ "$left_runs" -gt 0 ] || left_runs=0
  need=$(( left_runs * EST ))
  if [ "$rem" -le 600 ] || { [ "$need" -gt "$rem" ] && [ "$OVERRUN" != "1" ]; }; then
    if [ "$need" -gt "$rem" ] && [ "$rem" -gt 600 ]; then
      say "다음 타깃 $t 은 약 $((need/60))분 필요한데 예산이 $((rem/60))분 남았다 → 시작하지 않는다."
      say "  (타깃을 중간에 끊지 않기 위함. 예산을 넘겨서라도 돌리려면 OVERRUN=1)"
    else
      say "예산 소진 — 남은 타깃은 다음 실행에서 이어간다(이어달리기 지원)."
    fi
    break
  fi
  t_start=$(date +%s)
  say "───── $t (군 $grp · $model · rung$rung · ${nrows}서열 · Neff80 ${neffpk:--} · ${stratum:--})"
  say "  남은 예산 $((rem/60))분 · 이 타깃 예상 $((need/60))분 (남은 실행 ${left_runs}회 × 1회 $((EST/60))분 $((EST%60))초)"

  if [ $GATE -eq 1 ]; then
    # 판정은 낱말 검색이 아니라 결과 표(CSV)의 verdict 열로 한다.
    # prep_a3m_check_match.py 는 요약줄 "정상 N · 머리말오염 N · 서열자체다름 N" 을 항상 찍으므로,
    # '오염'·'다름' 같은 낱말을 grep 하면 정상인 타깃까지 전부 실패로 잡힌다(2026-07-28에 겪음).
    python prep_a3m_check_match.py --only "$t" --out "/tmp/gate_$t.csv" >"/tmp/gate_$t.log" 2>&1
    verdict=$(python3 - "/tmp/gate_$t.csv" <<'GATEPY'
import csv, os, sys
p = sys.argv[1]
if not os.path.exists(p):
    print("NOCSV"); raise SystemExit
rows = list(csv.DictReader(open(p)))
if not rows:
    print("EMPTY"); raise SystemExit
bad = [r for r in rows if r.get("verdict") != "OK"]
if bad:
    print("BAD " + " ".join(f"{r['chain']}:{r['verdict']}" for r in bad))
else:
    print(f"OK {len(rows)}")
GATEPY
    )
    case "$verdict" in
      OK\ *)  say "  MSA 게이트 통과 (항원 사슬 ${verdict#OK }개 질의행 일치)" ;;
      BAD\ *) say "  !! MSA 질의행 이상 → 건너뜀: ${verdict#BAD }  · /tmp/gate_$t.log"
              n_skip=$((n_skip+1)); continue ;;
      *) say "  !! MSA 게이트를 판정할 수 없음($verdict) → 건너뜀. /tmp/gate_$t.log 확인"
         n_skip=$((n_skip+1)); continue ;;
    esac
  fi

  comps=$(seq -s' ' 0 $((ncomp-1)))
  say "  ① 조성 $ncomp가지 × 반복 $nreps회"
  RUNG="$rung" TARGET="$t" MODEL="$model" REPLICAS="$ncomp" \
    COMPS="$comps" REPS="$nreps" bash make_composition_reps.sh || say "  !! 조성 단계에서 오류(로그 확인)"

  say "  ② 원래 MSA $nfull회"
  RUNG="$rung" TARGET="$t" MODEL="$model" \
    COMPS="full" REPS="$nfull" bash make_composition_reps.sh || say "  !! 원래 MSA 단계에서 오류(로그 확인)"

  have=$(done_runs "$t" "$model" "$dep"); want=$(( ncomp * nreps + nfull ))
  # 실행 1회 소요 시간을 실측해 다음 타깃의 예산 판단에 쓴다(군마다 항원 크기가 달라
  # 1회가 2분에서 4분까지 벌어진다).
  did=$(( have - have0 ))
  if [ "$did" -gt 0 ]; then
    EST=$(( ($(date +%s) - t_start) / did ))
    [ "$EST" -ge 30 ] || EST=30
    say "  실행 1회 약 $((EST/60))분 $((EST%60))초 (다음 타깃 예산 판단에 사용)"
  fi
  say "  $t 완료 — 실행 $have/$want"
  n_done=$((n_done+1))
done

echo
say "═══ 종료: 타깃 $n_done개 진행 · 게이트 실패로 건너뜀 $n_skip개 · 경과 $(( ($(date +%s)-START)/60 ))분"
say "다음: 타깃마다  bash run_analyze_target.sh <타깃>"
say "      판정 기준은 인수인계서 Ⅱ 6.4절(결과를 보기 전에 확정한 것)을 따른다."
