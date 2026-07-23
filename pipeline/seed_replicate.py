#!/usr/bin/env python3
"""[통제 Exp1 — Seed 복제] 같은 깊이(서열 개수)를 서로 다른 seed로 M번 독립 subsample.
목적: depth-response가 '깊이(개수)' 때문인지 '조성(어느 서열이 뽑혔나)' 때문인지 가르는 1차 필터.
  - 이후 co-folder로 M개를 각각 예측 → DockQ/recall 분포.
  - 분산 작음(다 비슷) = 깊이(개수)가 원인, 조성 무관 → 그 타깃 조성실험 불필요.
  - 분산 큼(draw마다 다름)   = 특정 서열(조성)이 원인 → Exp2(nested)·Exp3(LOCO/AOCI)로.
현재 사다리는 rung마다 '독립 랜덤 추첨'이라 개수 vs 조성이 섞여 있음 — 이 스크립트가 '개수 고정, 조성만 흔들기'.

사용(pipeline/에서, GPU 불필요):
  python seed_replicate.py --a3m <full.a3m> --depths 18,30 --replicas 10 --outdir seedrep/<target>
출력: seedrep/<target>/d<depth>/seed{s}.a3m ×M + neff.tsv(depth, seed, n_rows, neff80)
   → 이후 make_input.py로 각 seed a3m을 co-folder 입력으로 변환 후 예측(GPU).
"""
import argparse, os
import numpy as np
import neff_ladder as NL

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--a3m", required=True)
    ap.add_argument("--depths", required=True, help="쉼표구분 서열 개수(전이 깊이). 예: 18,30,64")
    ap.add_argument("--replicas", type=int, default=10, help="seed 복제 수 M")
    ap.add_argument("--start-seed", type=int, default=0)
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()
    headers, seqs = NL.read_raw(a.a3m)
    N = len(seqs)
    depths = [int(x) for x in a.depths.split(",") if x.strip()]
    os.makedirs(a.outdir, exist_ok=True)
    rows = []
    for depth in depths:
        d = min(depth, N)
        dd = os.path.join(a.outdir, f"d{depth}"); os.makedirs(dd, exist_ok=True)
        for s in range(a.start_seed, a.start_seed + a.replicas):
            rng = np.random.default_rng(s)                       # seed=조성만 바꾸는 축
            if d >= N:
                idx = list(range(N))                             # full이면 seed 무관(경고성)
            else:
                rest = rng.choice(np.arange(1, N), size=d - 1, replace=False)
                idx = [0] + sorted(rest.tolist())                # query(row0) 항상 포함
            out = os.path.join(dd, f"seed{s}.a3m")
            NL.write_a3m(out, headers, seqs, idx)
            nf = NL.neff80(NL.read_a3m_match_columns(out))
            rows.append((depth, s, len(idx), round(nf, 3)))
        print(f"[seedrep] depth={depth} ×{a.replicas} seeds (N_full={N})")
    with open(os.path.join(a.outdir, "neff.tsv"), "w") as f:
        f.write("depth\tseed\tn_rows\tneff80\n")
        for depth, s, n, nf in rows: f.write(f"{depth}\t{s}\t{n}\t{nf}\n")
    print(f"→ {a.outdir}/neff.tsv ({len(rows)} 파일). 다음: make_input으로 seed a3m → co-folder 입력(GPU 예측).")

if __name__ == "__main__":
    main()
