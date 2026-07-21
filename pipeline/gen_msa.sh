#!/usr/bin/env bash
# 항원 MSA 생성 + depth 사다리 빌드 (GPU 불필요, CPU/네트워크). 재실행 self-heal.
# A/B: targets/<id>/ 항원 사슬별 colabfold MSA → build_ladder. C: runs_diverse 기존 a3m 재활용.
# 출력 사다리: $DATA/ladders/<target>/<chain>/rung{0..R}.a3m + neff.tsv
#   bash gen_msa.sh            # 전체
#   ONLY=8q7s_O bash gen_msa.sh  # 한 타깃만(smoke)
# ⚠️ MSA_CMD를 서버 colabfold에 맞게 확인(기본 colabfold_batch --msa-only).
set -uo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"; cd "$HERE" || exit 1
DATA="${DATA:-/mnt/data/admuser/msadepth}"; LADDIR="$DATA/ladders"; RUNGS="${RUNGS:-6}"
DIVERSE="${DIVERSE:-$HERE/../runs_diverse}"; LIST="${LIST:-$HERE/sweep_targets.csv}"; ONLY="${ONLY:-}"
# 항원 서열 → a3m. 서버 colabfold에 맞게 조정. 입력=단일서열 fasta, 출력=<outdir>/*.a3m
MSA_CMD="${MSA_CMD:-colabfold_batch --msa-only}"
say(){ echo "[$(date '+%m-%d %H:%M:%S')] $*"; }
source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null; conda activate "${MSA_ENV:-colabfold}" 2>/dev/null

gen_full_a3m(){ # $1=fasta $2=outdir → echo a3m경로 (self-heal)
  local fa="$1" od="$2"; mkdir -p "$od"
  local a3m; a3m=$(find "$od" -maxdepth 1 -name '*.a3m' 2>/dev/null | head -1)
  [ -n "$a3m" ] && { echo "$a3m"; return; }
  $MSA_CMD "$fa" "$od" >"$od/msa.log" 2>&1 || return 1
  find "$od" -maxdepth 1 -name '*.a3m' 2>/dev/null | head -1
}

while IFS=, read -r target pdb group ab dirtype ag_chains label; do
  [ -z "$target" ] && continue
  [ -n "$ONLY" ] && [ "$target" != "$ONLY" ] && continue
  IFS='|' read -ra AGC <<< "$ag_chains"
  if [ "$dirtype" = "diverse" ]; then
    # C: 기존 a3m 재활용(runs_diverse/<pdb>/ 아무 a3m)
    src=$(find "$DIVERSE/$pdb" -name '*.a3m' 2>/dev/null | head -1)
    [ -z "$src" ] && { say "skip $target (diverse a3m 없음: $DIVERSE/$pdb)"; continue; }
    od="$LADDIR/$target/A"
    [ -f "$od/rung$((RUNGS-1)).a3m" ] && { say "skip $target (사다리 있음)"; continue; }
    python "$HERE/build_ladder.py" --a3m "$src" --outdir "$od" --rungs "$RUNGS" && say "ladder(C) $target ← $(basename "$src")"
  else
    cj="$HERE/targets/$target/chains.json"
    [ -f "$cj" ] || { say "skip $target (chains.json 없음)"; continue; }
    for c in "${AGC[@]}"; do
      od="$LADDIR/$target/$c"
      [ -f "$od/rung$((RUNGS-1)).a3m" ] && { say "skip $target/$c (사다리 있음)"; continue; }
      # 이 항원 사슬만 fasta 추출
      fa="$HERE/targets/$target/ag_${c}.fasta"
      python -c "
import json,sys
d=json.load(open('$cj')); sm={x['id']:x['seq'] for x in d['chains']}
open('$fa','w').write('>${target}_${c}\n'+sm['$c']+'\n')"
      msad="$HERE/targets/$target/msa_$c"
      a3m=$(gen_full_a3m "$fa" "$msad") || { say "  !! $target/$c MSA 실패 → $msad/msa.log"; continue; }
      python "$HERE/build_ladder.py" --a3m "$a3m" --outdir "$od" --rungs "$RUNGS" && say "ladder $target/$c ← $(basename "$a3m")"
    done
  fi
done < <(tail -n +2 "$LIST")
say "gen_msa 완료 → $LADDIR/<target>/<chain>/. 다음: bash run_sweep.sh boltz 11"
