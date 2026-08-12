#!/usr/bin/env bash
# ══════════════════════════════════════════════════════════════════════════════
# [선택기 교체 재도킹] 후보를 ncomp 가 아니라 ncomp×AbEpi-max 로 고른 뒤 네 팔을 돌린다.
#
# 왜: 2026-08-01 크기-정규화 지표로 다시 재니 선택기 순위가 바뀌었다.
#     F1 천장대비  ncomp 0.68 · AbEpi-max 0.74 · ncomp×AbEpi-max 0.76 · 무작위 0.61
#     (덮음만 보면 '가장 큰 후보'가 1등으로 나오는데 그건 지표의 인공물이다)
#
# 사용:
#   bash run_reselect.sh                 # 무엇을 할지만 출력
#   bash run_reselect.sh --apply         # 실제 실행
#   SEL=abepi_max bash run_reselect.sh --apply
#
# env: MSAD(msa-depth pipeline) · CD(consensus_docking) · SEL · ARMS · MODELS
# ══════════════════════════════════════════════════════════════════════════════
set -uo pipefail
MSAD="${MSAD:-$HOME/projects/bk21-msa-depth-bias/pipeline}"
CD="${CD:-$HOME/projects/bk21-antibody-ml/pipeline}"
SEL="${SEL:-ncomp_x_abemax}"
ARMS="${ARMS:-noconstraint fullmsa sizematch ours}"
MODELS="${MODELS:-boltz}"
APPLY=0; for a in "$@"; do [ "$a" = "--apply" ] && APPLY=1; done

cd "$MSAD" || { echo "!! $MSAD 없음"; exit 1; }
ABE="$CD/results/abepiscore_all.csv"
[ -f "$ABE" ] || { echo "!! $ABE 없음 — eval_abepitope.py 먼저"; exit 1; }

echo "[1/3] 선택기 '$SEL' 로 후보 고르기"
python select_eval_selectors.py --sites results --abepi "$ABE" \
       --pick "$SEL" --out results/selected_sites.csv || exit 1

echo ""
echo "[2/3] ncomp 와 달라지는 타깃 (여기만 다시 돌리면 된다)"
python - <<'PY'
import csv, json, os
sel = {r["target"]: r for r in csv.DictReader(open("results/selected_sites.csv"))}
chg = []
for t, r in sorted(sel.items()):
    d = json.load(open(f"results/sites_{t}.json"))
    nc = max(d["candidates"], key=lambda c: (c["n_comp"], -len(c["residues"])))["cand"]
    if int(r["cand"]) != nc:
        chg.append((t, nc, r["cand"], r["true_covered"], r["f1"]))
print(f"  {'타깃':<11}{'ncomp':>6}{'새선택':>7}{'덮음':>8}{'F1':>7}")
for t, a, b, cov, f in chg:
    print(f"  {t:<11}{a:>6}{b:>7}{float(cov):>8.3f}{float(f):>7.3f}")
print(f"\n  바뀌는 타깃 {len(chg)}종 / 전체 {len(sel)}종")
open("results/reselect_targets.txt", "w").write(" ".join(t for t, *_ in chg) + "\n")
PY
TGTS=$(cat results/reselect_targets.txt 2>/dev/null)
echo "  → results/reselect_targets.txt : ${TGTS:-없음}"

if [ $APPLY -eq 0 ]; then echo -e "\n(dry-run. 실제로 하려면 --apply)"; exit 0; fi
[ -n "$TGTS" ] || { echo "바뀌는 타깃이 없다 — 끝"; exit 0; }

echo ""
echo "[3/3] 유도 재도킹 (팔: $ARMS · 모델: $MODELS)"
cd "$CD" || exit 1
for t in $TGTS; do
  CAND=$(python - "$t" <<'PY'
import csv,sys
for r in csv.DictReader(open(f"{__import__('os').environ['MSAD']}/results/selected_sites.csv")):
    if r["target"]==sys.argv[1]: print(r["cand"]); break
PY
)
  echo "── $t (후보 $CAND) ──"
  for arm in $ARMS; do
    python scripts/make_pocket_from_sites.py --sites "$MSAD/results/sites_${t}.json" \
        --chains "$MSAD/targets/${t}/chains.json" --arm "$arm" --cand "$CAND" \
        --emit boltz,protenix || echo "  ! $arm 입력 생성 실패"
  done
done
ARMS="$ARMS" MODELS="$MODELS" bash scripts/run_demo_guided.sh "$TGTS"
python scripts/eval_demo.py --targets "$TGTS" --models "$MODELS" \
       --out results/demo_dockq_reselect.csv
echo "완료 → $CD/results/demo_dockq_reselect.csv"
