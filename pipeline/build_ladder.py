#!/usr/bin/env python3
"""full a3m → rung-인덱스 사다리(rung0=full … rung{R-1}=single-seq). 사슬 간 정렬 위해 인덱스로.
halving 스케줄: 행수 [N, N/2, N/4, …, 1]. 각 rung의 실측 Neff80 기록.
사용: python build_ladder.py --a3m msa/A.a3m --outdir ladder/A --rungs 6 [--seed 0]
출력: ladder/A/rung{k}.a3m + ladder/A/neff.tsv (rung, n_rows, neff80)
"""
import argparse, os
import neff_ladder as NL   # read_raw, read_a3m_match_columns, neff80, write_a3m
import numpy as np

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--a3m", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--rungs", type=int, default=6)
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    headers, seqs = NL.read_raw(a.a3m)
    N = len(seqs)
    # halving 행수(중복 제거, 내림차순), 마지막은 1(single-seq=query만)
    counts = []
    c = N
    for _ in range(a.rungs - 1):
        counts.append(max(1, int(round(c)))); c /= 2.0
    counts.append(1)
    counts = sorted(set(counts), reverse=True)
    while len(counts) < a.rungs: counts.append(1)   # N이 작아 겹치면 1로 채움(뒤에서 dedup)
    counts = counts[:a.rungs]
    rng = np.random.default_rng(a.seed)
    rows = []
    for k, cnt in enumerate(counts):
        cnt = min(cnt, N)
        if cnt >= N: idx = list(range(N))
        else:
            rest = rng.choice(np.arange(1, N), size=cnt - 1, replace=False)
            idx = [0] + sorted(rest.tolist())
        out = os.path.join(a.outdir, f"rung{k}.a3m")
        NL.write_a3m(out, headers, seqs, idx)
        nf = NL.neff80(NL.read_a3m_match_columns(out))
        rows.append((k, cnt, round(nf, 3)))
    with open(os.path.join(a.outdir, "neff.tsv"), "w") as f:
        f.write("rung\tn_rows\tneff80\n")
        for k, cnt, nf in rows: f.write(f"{k}\t{cnt}\t{nf}\n")
    print(f"[ladder] {a.a3m} N={N} → {len(rows)} rungs " +
          " ".join(f"r{k}:{cnt}({nf})" for k, cnt, nf in rows))

if __name__ == "__main__":
    main()
