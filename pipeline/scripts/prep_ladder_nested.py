#!/usr/bin/env python3
"""[통제 Exp2 — Nested 사다리] rung_{k+1} ⊂ rung_k (재추첨 없이 '빼기만').
목적: 깊이 전이를 '그 단계에서 빠진 특정 서열군'에 귀속. 독립추첨 사다리(build_ladder)는 rung간
  서열집합이 달라 rung→rung 비교로 인과서열을 못 짚음. nested는 중첩이라 '빠진 서열 = 전이 원인 후보'.
방법: 비-query 서열을 한 번만 셔플 → rung_k = query + (셔플 순서 앞 cnt_k−1개). counts 감소 → 자동 중첩.
예상: 특정 단계에서 DockQ 급락 → 그 단계 dropped 서열 = 인과 후보 / 완만하면 개수 자체(깊이)가 중요.

사용(pipeline/에서, GPU 불필요):
  python prep_ladder_nested.py --a3m <full.a3m> --outdir nested/<target> --rungs 12 [--min-rows 1] [--seed 0]
출력: nested/<target>/rung{k}.a3m + neff.tsv + membership.tsv(rung별 포함 서열 index)
      + dropped.tsv(전이 k-1→k 마다 빠진 서열 index·header) ← 인과 후보
"""
import argparse, os
import numpy as np
import neff_ladder as NL

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--a3m", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--rungs", type=int, default=12)
    ap.add_argument("--min-rows", type=int, default=1)
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    headers, seqs = NL.read_raw(a.a3m)
    N = len(seqs)
    FLOOR = max(1, min(a.min_rows, N))
    if N <= FLOOR:
        counts = [N] * a.rungs
    else:
        counts = [int(round(x)) for x in np.geomspace(N, FLOOR, a.rungs)]
    counts = sorted(set(counts), reverse=True)
    while len(counts) < a.rungs: counts.append(FLOOR)
    counts = counts[:a.rungs]

    rng = np.random.default_rng(a.seed)
    perm = rng.permutation(np.arange(1, N))          # ★ 비-query 순서를 '한 번만' 정함 → 중첩 보장
    neff_rows, memb_rows, drop_rows = [], [], []
    prev = None
    for k, cnt in enumerate(counts):
        cnt = min(cnt, N)
        idx = list(range(N)) if cnt >= N else [0] + sorted(perm[:cnt - 1].tolist())
        out = os.path.join(a.outdir, f"rung{k}.a3m")
        NL.write_a3m(out, headers, seqs, idx)
        nf = NL.neff80(NL.read_a3m_match_columns(out))
        neff_rows.append((k, cnt, round(nf, 3)))
        cur = set(idx)
        memb_rows.append((k, cnt, " ".join(map(str, sorted(cur)))))
        if prev is not None:
            dropped = sorted(prev - cur)             # 이 단계에서 빠진 서열 = 전이 원인 후보
            for di in dropped:
                drop_rows.append((k - 1, k, di, headers[di][:60]))
        prev = cur

    with open(os.path.join(a.outdir, "neff.tsv"), "w") as f:
        f.write("rung\tn_rows\tneff80\n")
        for k, cnt, nf in neff_rows: f.write(f"{k}\t{cnt}\t{nf}\n")
    with open(os.path.join(a.outdir, "membership.tsv"), "w") as f:
        f.write("rung\tn_rows\tindices\n")
        for k, cnt, ix in memb_rows: f.write(f"{k}\t{cnt}\t{ix}\n")
    with open(os.path.join(a.outdir, "dropped.tsv"), "w") as f:
        f.write("from_rung\tto_rung\tdropped_idx\theader\n")
        for fr, to, di, h in drop_rows: f.write(f"{fr}\t{to}\t{di}\t{h}\n")
    print(f"[nested] {a.a3m} N={N} → {len(neff_rows)} rungs (중첩 보장) " +
          " ".join(f"r{k}:{c}({nf})" for k, c, nf in neff_rows))
    print(f"→ {a.outdir}/ (neff.tsv·membership.tsv·dropped.tsv). 전이 나는 단계의 dropped 서열 = 인과 후보.")

if __name__ == "__main__":
    main()
