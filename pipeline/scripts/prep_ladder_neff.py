#!/usr/bin/env python3
"""항원별 Neff80 depth-사다리 빌더 (143 하드코드 대체).
Neff80 = θ=0.8로 가중한 per-residue Neff의 median (AF 방식; NEFFy/AF2 supp).
각 항원의 full a3m에서:
  ① full Neff80 측정 → ② 로그간격으로 행수 subsample(query 항상 포함) → ③ 각 subsample의 실측 Neff80 라벨
  → depth 축을 절대행수(143)가 아니라 항원별 full→single 상대 + 실측 Neff80로.
사용: python prep_ladder_neff.py --a3m runs/<t>/msa_<t>/A.a3m --outdir ladders/<t> --rungs 8
"""
import argparse, os, numpy as np
import re

def read_a3m_match_columns(path):
    """a3m -> (N,L) uint8. 소문자/'.'(insertion) 제거, 대문자+'-'(match)만."""
    names, seqs, cur = [], [], []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line: continue
            if line[0] == ">":
                if cur: seqs.append("".join(cur)); cur = []
                names.append(line[1:])
            else: cur.append(line)
    if cur: seqs.append("".join(cur))
    keep = lambda s: "".join(c for c in s if c.isupper() or c == "-")
    seqs = [keep(s) for s in seqs]
    L = len(seqs[0])
    seqs = [s for s in seqs if len(s) == L]        # 길이 안 맞는 행 방어적 제거
    arr = np.frombuffer("".join(seqs).encode("ascii"), dtype=np.uint8).reshape(len(seqs), L)
    return arr

def neff80(arr, theta=0.8):
    """AF식 per-residue Neff의 median. identity 분모='둘 중 하나라도 non-gap'(합집합)."""
    GAP = ord('-'); nongap = (arr != GAP); N, L = arr.shape
    neigh = np.ones(N, dtype=np.float64)
    for i in range(N):
        both = nongap[i] & nongap
        eq = (both & (arr[i] == arr)).sum(1)
        uni = (nongap[i] | nongap).sum(1)
        ident = eq / np.maximum(uni, 1)
        neigh[i] = (ident >= theta).sum()
    w = 1.0 / neigh
    neff_col = (nongap * w[:, None]).sum(0)
    return float(np.median(neff_col))

def write_a3m(path, headers, seqs, idx):
    with open(path, "w") as f:
        for i in idx:
            f.write(f">{headers[i]}\n{seqs[i]}\n")

def read_raw(path):
    """a3m 원본을 (헤더, 서열)로 읽는다.

    ⚠️ 2026-07-27 버그 수정: ColabFold a3m은 첫 줄이 메타 주석(`#<길이>\t<개수>`)이다.
    이전 판은 (a) 주석 줄을 서열로 취급했고 (b) 첫 `>`에서 h가 None이라 cur을 비우지
    않아, 주석이 **질의 서열 앞에 그대로 붙었다**(예: `#440\t1NLWVT...`).
    그 결과 boltz는 "MSA does not match input sequence"로 MSA를 통째로 버리고,
    Protenix는 질의행만 밀린 정렬을 그대로 썼다. 아래 두 줄이 그 수정이다.
    """
    headers, seqs, cur, h = [], [], [], None
    for line in open(path):
        line = line.rstrip("\n")
        if not line: continue
        if line[0] == "#":                    # 메타 주석(`#<길이>\t<개수>`)
            rest = re.sub(r"^#\d+\t\d+", "", line)
            if not rest: continue             # 단독 주석 줄 → 버림(원본 a3m의 정상 형태)
            line = rest                       # 서열이 뒤에 붙어 있으면(옛 손상 파일) 서열만 살림
        if line[0] == ">":
            if h is not None: headers.append(h); seqs.append("".join(cur))
            cur = []                          # h가 None일 때도 반드시 초기화
            h = line[1:]
        else: cur.append(line)
    if h is not None: headers.append(h); seqs.append("".join(cur))
    return headers, seqs

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--a3m", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--rungs", type=int, default=8, help="full~single 사이 로그 눈금 수")
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    headers, seqs = read_raw(a.a3m)
    N = len(seqs)
    arr_full = read_a3m_match_columns(a.a3m)
    nf_full = neff80(arr_full)
    print(f"[{a.a3m}] N={N} rows  full Neff80={nf_full:.2f}")
    # 로그간격 행수: N → 1 (query=row0 항상 포함)
    counts = sorted(set(int(round(c)) for c in np.geomspace(N, 1, a.rungs)), reverse=True)
    if 1 not in counts: counts.append(1)
    rng = np.random.default_rng(a.seed)
    rows = []
    for c in counts:
        if c >= N: idx = list(range(N))
        else:
            rest = rng.choice(np.arange(1, N), size=c-1, replace=False)
            idx = [0] + sorted(rest.tolist())
        out = os.path.join(a.outdir, f"depth_{c}.a3m")
        write_a3m(out, headers, seqs, idx)
        nf = neff80(read_a3m_match_columns(out))
        rows.append((c, round(nf, 3), out))
        print(f"  rows={c:5d}  Neff80={nf:7.3f}  → {os.path.basename(out)}")
    with open(os.path.join(a.outdir, "ladder.tsv"), "w") as f:
        f.write("n_rows\tneff80\tpath\n")
        for c, nf, p in rows: f.write(f"{c}\t{nf}\t{p}\n")
    print(f"→ {a.outdir}/ladder.tsv ({len(rows)} rungs; 축=실측 Neff80)")

if __name__ == "__main__":
    main()
