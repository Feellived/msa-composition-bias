#!/usr/bin/env bash
# ══════════════════════════════════════════════════════════════════════════════
# seed-복제 케이스스터디 러너 — cross-check가 찾은 depth-rescue가 진짜인지 통제.
#
# 각 후보 (target, model, peak_rung)에 대해:
#   ① peak rung의 a3m에서 서열 개수(n_rows)를 세고 = 그 rung이 성공했던 '깊이'
#   ② 그 깊이로 FULL a3m을 N seed 재추첨(seed_replicate) = 개수 고정, 조성만 흔들기
#   ③ 그 모델(boltz/protenix)로 각 seed 예측
#   → 채점(score_seedrep_cand.py)에서 "그 깊이서 관측했던 peak(예 8ulr 0.66)이 N seed에서 재현되나
#     vs best-of-5 단발 운"을 판정.
#
# ⚠️ 재랭커 아님 — MSA깊이 가설 케이스스터디(Arm B).
# ⚠️ 후보 전부 항원 사슬 단일(A)이라 단일-항원사슬 가정. RUN은 run_sweep.sh 검증본 그대로.
#
# 사용(tmux 권장):
#   cd ~/projects/bk21-msa-depth-bias/pipeline && git pull
#   SMOKE=1 bash run_seedrep_cand.sh          # 후보1개×seed1개 스모크(모델·경로 검증)
#   bash run_seedrep_cand.sh 2>&1 | tee seedrep_cand_$(date +%m%d).log
#   # 채점: python score_seedrep_cand.py  (DockQ env)
# ══════════════════════════════════════════════════════════════════════════════
set -uo pipefail
cd "$(cd "$(dirname "$0")" && pwd)"
DATA="${DATA:-/mnt/data/admuser/msadepth}"
LADDIR="$DATA/ladders"
LIST="${LIST:-sweep_targets.csv}"          # ag_chains 조회
CAND="${CAND:-seedrep_cand.csv}"           # target,model,peak_rung,replicas,obs_dq,obs_rec
OUTROOT="${OUTROOT:-$DATA/seedrep_cand}"
SAMP="${SAMP:-5}"; SEED="${SEED:-0}"; SMOKE="${SMOKE:-0}"
say(){ echo "[$(date '+%m-%d %H:%M:%S')] $*"; }
# 캐시 홈 보호(run_sweep와 동일)
export XDG_CACHE_HOME="${XDG_CACHE_HOME:-$DATA/.cache}" TORCH_EXTENSIONS_DIR="${TORCH_EXTENSIONS_DIR:-$DATA/.cache/torch_ext}"
export HF_HOME="${HF_HOME:-$DATA/.cache/hf}" PIP_CACHE_DIR="${PIP_CACHE_DIR:-$DATA/.cache/pip}"
mkdir -p "$XDG_CACHE_HOME" "$TORCH_EXTENSIONS_DIR" "$HF_HOME"
source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null

# 모델별 env·RUN·DONE·EXT (run_sweep.sh 검증본 그대로)
setup_model(){
  case "$1" in
    boltz)
      conda activate boltz 2>/dev/null
      command -v boltz >/dev/null || { say "!! boltz 없음"; return 1; }
      export BOLTZ_CACHE="${BOLTZ_CACHE:-$DATA/boltz_cache}"; mkdir -p "$BOLTZ_CACHE"
      EXT="yaml"
      RUN(){ ( cd "$2" && boltz predict "$1" --out_dir results --cache "$BOLTZ_CACHE" --no_kernels --diffusion_samples "$SAMP" >"$3" 2>&1 ); }
      DONE(){ find "$1/results" -name '*_model_*.cif' 2>/dev/null | grep -q .; } ;;
    protenix)
      conda activate "${PROT_ENV:-protenix}" 2>/dev/null
      command -v protenix >/dev/null || { say "!! protenix 없음"; return 1; }
      export PROTENIX_ROOT_DIR="${PROTENIX_ROOT_DIR:-/mnt/data/admuser/protenix_weights}" LAYERNORM_TYPE=torch
      PROT_MODEL="${PROT_MODEL:-protenix_base_default_v1.0.0}"; EXT="json"   # leakage-free
      RUN(){ ( cd "$2" && protenix pred -i "$1" -o results -n "$PROT_MODEL" -s "$SEED" -e "$SAMP" \
               --trimul_kernel torch --triatt_kernel torch --enable_fusion False >"$3" 2>&1 ); }
      DONE(){ find "$1/results" -name '*sample*.cif' 2>/dev/null | grep -q .; } ;;
    chai)
      conda activate "${CHAI_ENV:-chai}" 2>/dev/null
      command -v chai-lab >/dev/null || { say "!! chai-lab 없음"; return 1; }
      export CHAI_DOWNLOADS_DIR="${CHAI_DOWNLOADS_DIR:-$DATA/chai_downloads}"; mkdir -p "$CHAI_DOWNLOADS_DIR"; EXT="fasta"
      RUN(){ rm -rf "$2/results"; ( cd "$2" && chai-lab fold --msa-directory "$2/msa" "$1" results >"$3" 2>&1 ); }
      DONE(){ find "$1/results" -name '*.cif' 2>/dev/null | grep -q .; } ;;
    *) return 1;;
  esac
  # numpy 프리플라이트(seed_replicate = numpy + neff_ladder)
  python -c 'import numpy' 2>/dev/null || { say "!! $1 env에 numpy 없음 — seed_replicate 불가"; return 1; }
  return 0
}

ag_of(){ awk -F, -v t="$1" 'NR>1 && $1==t{print $6; exit}' "$LIST"; }

while IFS=, read -r target model peak_rung replicas obs_dq obs_rec; do
  [ -z "${target:-}" ] && continue; [ "$target" = "target" ] && continue
  replicas="${replicas:-8}"
  ag=$(ag_of "$target"); [ -n "$ag" ] || { say "skip $target (ag_chains 조회 실패: $LIST)"; continue; }
  cj="targets/$target/chains.json"; [ -f "$cj" ] || { say "skip $target (chains.json 없음)"; continue; }
  IFS='|' read -ra AGC <<< "$ag"
  setup_model "$model" || { say "skip $target ($model 준비 실패)"; continue; }

  # ①·② 각 항원 사슬을 peak 깊이 × N seed 재추첨
  # ⚠️ 사슬마다 서열 개수가 다르다(예: 8txu_HL A=413, B=579). 그래서 사슬별 깊이를 따로 기억한다.
  #    (예전 버그: 첫 사슬의 깊이 폴더 이름으로 모든 사슬을 찾아 다중사슬 항원이 통째로 skip 됨)
  DEPTHS=()   # AGC와 같은 순서(인덱스 배열 = bash 3.2에서도 동작)
  ok=1
  for i in "${!AGC[@]}"; do
    c="${AGC[$i]}"
    full="$LADDIR/$target/$c/rung0.a3m"; peaka="$LADDIR/$target/$c/rung${peak_rung}.a3m"
    [ -f "$full" ] && [ -f "$peaka" ] || { say "  !! $target/$c rung0/rung$peak_rung a3m 없음"; ok=0; break; }
    depth=$(grep -c '^>' "$peaka" 2>/dev/null || echo 0)
    [ "$depth" -ge 1 ] || { say "  !! $target/$c rung$peak_rung 서열 0"; ok=0; break; }
    DEPTHS[$i]="$depth"
    say "$target/$c: peak rung$peak_rung = ${depth}서열 × ${replicas} seed 재추첨"
    python seed_replicate.py --a3m "$full" --depths "$depth" --replicas "$replicas" \
        --outdir "seedrep_cand/${target}_${c}" >/dev/null \
      || { say "  !! seed_replicate 실패 $target/$c"; ok=0; break; }
  done
  [ "$ok" = 1 ] || continue

  # ③ seed 축으로 순회(사슬마다 깊이는 달라도 seed 번호는 공통) → 각 사슬 자기 깊이 폴더에서 같은 seed를 꺼내 map 구성
  c0="${AGC[0]}"; d0="d${DEPTHS[0]}"; base="seedrep_cand/${target}_${c0}/$d0"
  for a3m0 in "$base"/seed*.a3m; do
    [ -e "$a3m0" ] || continue; s="$(basename "$a3m0" .a3m)"
    map=""; okk=1
    for i in "${!AGC[@]}"; do
      c="${AGC[$i]}"
      f="seedrep_cand/${target}_${c}/d${DEPTHS[$i]}/${s}.a3m"
      [ -f "$f" ] || { say "  !! $target $s: $c 사슬 a3m 없음($f)"; okk=0; break; }
      map="${map:+$map,}$c=$f"
    done
    [ "$okk" = 1 ] || continue
    out="$OUTROOT/$model/$target/$d0/$s"
    if DONE "$out"; then continue; fi
    mkdir -p "$out"; inp="$out/input.$EXT"
    python make_input.py --cofolder "$model" --chains "$cj" --ag-a3m "$map" --dir "$out" --out "$inp" >"$out/mk.log" 2>&1 \
      || { say "  !! make_input 실패 $target $d0 $s (로그 $out/mk.log)"; continue; }
    say "run $model $target $d0 $s ..."
    if RUN "$inp" "$out" "$out/run.log" && DONE "$out"; then say "  OK $target $d0 $s"
    else say "  !! 실패 $target $d0 $s → $out/run.log (tail: $(tail -1 "$out/run.log" 2>/dev/null))"; fi
    [ "$SMOKE" = 1 ] && { say "SMOKE=1 → 1건 후 종료. 확인: $out/results"; exit 0; }
  done
done < "$CAND"
say "완료. 채점: python score_seedrep_cand.py  (DockQ env)"
