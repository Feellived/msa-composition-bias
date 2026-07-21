#!/usr/bin/env python3
"""full a3m → rung-인덱스 사다리(rung0=full … 최심=min_rows 서열). 사슬 간 정렬 위해 인덱스로.
스케줄: geomspace(N → min_rows), 로그균등 → 얕은 구간까지 촘촘.
⚠️ 순수 single-seq(1행)는 Boltz 데이터로더가 폭주(수십GB)·stall → 최심 rung을 min_rows(기본 4)로 바닥.
각 rung의 실측 Neff80 기록.
사용: python build_ladder.py --a3m msa/A.a3m --outdir ladder/A --rungs 6 [--min-rows 4] [--seed 0]
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
    ap.add_argument("--min-rows", type=int, default=4,
                    help="최심 rung 서열수(1=순수 single-seq는 Boltz 데이터로더 폭주 → 금지, 기본 4)")
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    headers, seqs = NL.read_raw(a.a3m)
    N = len(seqs)
    # geomspace(N → floor) 로그균등. ⚠️ 순수 single-seq(1) 금지: 최소 min_rows 서열 바닥(Boltz 폭주 회피).
    FLOOR = max(1, min(a.min_rows, N))
    if N <= FLOOR:
        counts = [N] * a.rungs                              # depth range 없음(작은 MSA) — 전부 full
    else:
        counts = [int(round(x)) for x in np.geomspace(N, FLOOR, a.rungs)]
    counts = sorted(set(counts), reverse=True)
    while len(counts) < a.rungs: counts.append(FLOOR)       # 반올림 중복 시 floor로 채움
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
