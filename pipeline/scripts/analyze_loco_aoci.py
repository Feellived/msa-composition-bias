#!/usr/bin/env python3
"""[통제 Exp3 — LOCO/AOCI] full MSA를 %동일성 그리디 클러스터 → cluster별 leave-out / add-in a3m 생성.
목적: 어느 '진화 클러스터(동족 서열군)'가 co-folder의 배치를 좌우하는지 직접 인과 검정(클러스터 단위 = compute 유계).
  - LOCO(leave-one-cluster-out): full − cluster_i → 예측. 큰 하락 = cluster_i가 정답 배치에 '필요'.
  - AOCI(add-one-cluster-in):   query + cluster_i → 예측. 큰 rescue = cluster_i만으로 정답 유도.
  - 교차확인: 인과 cluster 제거가 깊이 전이 재현? 그 추가가 single-seq rescue?
query(row0)는 모든 파일에 항상 포함(클러스터링 대상에서 제외).

사용(pipeline/에서, GPU 불필요):
  python analyze_loco_aoci.py --a3m <full.a3m> --outdir loco/<target> [--ident 0.62] [--min-clust 3] [--max-clusters 15] [--max-seqs 4000]
출력: loco/<target>/{clusters.tsv, loco_c{i}.a3m, aoci_c{i}.a3m, neff.tsv}
   ident 정의 = neff80과 동일(eq / non-gap 합집합). 큰 MSA는 --max-seqs로 subsample 후 클러스터(명시).
"""
import argparse, os
import numpy as np
import neff_ladder as NL

GAP = ord('-')

def greedy_cluster(arr, ident):
    """arr (M,L) uint8(=비-query 서열의 match columns). 그리디: 미할당 첫 서열=centroid, ident 이상=같은 클러스터."""
    M = arr.shape[0]
    nongap = (arr != GAP)
    unassigned = np.ones(M, dtype=bool)
    labels = -np.ones(M, dtype=int)
    cid = 0
    order = np.arange(M)
    while unassigned.any():
        i = order[unassigned][0]                                   # 미할당 첫 서열 = centroid
        ci, cn = arr[i], nongap[i]
        both = cn & nongap                                         # (M,L)
        eq = (both & (arr == ci)).sum(1)
        uni = (cn | nongap).sum(1)
        idn = eq / np.maximum(uni, 1)
        members = unassigned & (idn >= ident)
        members[i] = True
        labels[members] = cid; unassigned[members] = False; cid += 1
    return labels

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--a3m", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--ident", type=float, default=0.62, help="클러스터 동일성 임계")
    ap.add_argument("--min-clust", type=int, default=3, help="이 크기 미만 클러스터는 LOCO/AOCI 생략")
    ap.add_argument("--max-clusters", type=int, default=15, help="큰 클러스터 상위 N개만(compute 유계)")
    ap.add_argument("--max-seqs", type=int, default=4000, help="이보다 크면 클러스터링용 subsample(query 포함)")
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    headers, seqs = NL.read_raw(a.a3m)
    N = len(seqs)
    arr = NL.read_a3m_match_columns(a.a3m)                          # (N,L), row0=query
    # 대용량이면 비-query에서 subsample (query=0 항상 유지)
    idx_all = np.arange(N)
    if N > a.max_seqs:
        rng = np.random.default_rng(a.seed)
        keep = np.concatenate([[0], rng.choice(np.arange(1, N), size=a.max_seqs - 1, replace=False)])
        keep.sort(); idx_all = keep
        print(f"[loco] N={N} > max_seqs={a.max_seqs} → 클러스터링용 {len(keep)} subsample")
    nonq = idx_all[idx_all != 0]                                   # 비-query 원본 index
    labels = greedy_cluster(arr[nonq], a.ident)
    # 클러스터별 원본 index 모으기
    clusters = {}
    for lab, gi in zip(labels, nonq):
        clusters.setdefault(lab, []).append(int(gi))
    big = sorted(clusters.items(), key=lambda kv: -len(kv[1]))
    big = [(l, m) for l, m in big if len(m) >= a.min_clust][:a.max_clusters]

    crows, neff_rows = [], []
    all_nonq = set(nonq.tolist())
    for rank_i, (lab, mem) in enumerate(big):
        mem_set = set(mem)
        # LOCO: query + (모든 비-query − 이 클러스터)
        loco_idx = [0] + sorted(all_nonq - mem_set)
        lp = os.path.join(a.outdir, f"loco_c{rank_i}.a3m"); NL.write_a3m(lp, headers, seqs, loco_idx)
        # AOCI: query + 이 클러스터만
        aoci_idx = [0] + sorted(mem_set)
        apth = os.path.join(a.outdir, f"aoci_c{rank_i}.a3m"); NL.write_a3m(apth, headers, seqs, aoci_idx)
        nf_l = NL.neff80(NL.read_a3m_match_columns(lp)); nf_a = NL.neff80(NL.read_a3m_match_columns(apth))
        neff_rows.append((rank_i, "loco", len(loco_idx), round(nf_l, 3)))
        neff_rows.append((rank_i, "aoci", len(aoci_idx), round(nf_a, 3)))
        # 대표 헤더(taxon 힌트)
        crows.append((rank_i, len(mem), headers[mem[0]][:70]))
        print(f"  c{rank_i}: size={len(mem):4d}  loco(N={len(loco_idx)},Neff{nf_l:.1f}) aoci(N={len(aoci_idx)},Neff{nf_a:.1f})  ~{headers[mem[0]][:40]}")

    with open(os.path.join(a.outdir, "clusters.tsv"), "w") as f:
        f.write("cluster\tsize\trep_header\n")
        for ci, sz, h in crows: f.write(f"{ci}\t{sz}\t{h}\n")
    with open(os.path.join(a.outdir, "neff.tsv"), "w") as f:
        f.write("cluster\tkind\tn_rows\tneff80\n")
        for ci, kind, n, nf in neff_rows: f.write(f"{ci}\t{kind}\t{n}\t{nf}\n")
    tot = sum(len(m) for _, m in clusters.items())
    print(f"[loco] 클러스터 {len(clusters)}개(임계 {a.ident}), LOCO/AOCI 대상 상위 {len(big)}개 → {a.outdir}/")
    print("  다음: 각 loco_/aoci_ a3m을 make_input → co-folder 예측 → DockQ/recall (LOCO 큰 하락=필요 / AOCI 큰 rescue=유도)")

if __name__ == "__main__":
    main()
