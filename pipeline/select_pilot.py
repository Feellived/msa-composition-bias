#!/usr/bin/env python3
"""labeled 매니페스트 → 서열 dedup → 항원별 A/B 균형 파일럿 선별.
- 구조에서 VH(+VL) 관측서열 추출 → (항원,AB,H+L서열) 정확 dedup(결정학 사본·동일항체 제거).
- 불량 flag 제외(no-RBD-contact·HA-no-map·error·low-classmatch·ambiguous·borderline·mixed).
- 항원별 A/B 각 N개, 에피토프 클래스·해상도 다양성 우선 선별.
사용: python select_pilot.py --n 10 --out locked_pilot_overrep.csv
"""
import argparse, csv, os, warnings
warnings.filterwarnings("ignore")
from collections import defaultdict, OrderedDict
from Bio.PDB import MMCIFParser
import classify_epitope as C

BAD = ("no-RBD-contact","HA-no-map-review","error","low-classmatch","ambiguous-lowoverlap","borderline")

def hlseq(model, r):
    hs = "".join(C.chain_seq(C.get_chain_residues(model, c))[0] for c in C.split_chains(r["Hchain"]))
    ls = "".join(C.chain_seq(C.get_chain_residues(model, c))[0] for c in C.split_chains(r["Lchain"]))
    return hs + "|" + ls

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--labeled", default="manifest_labeled.csv")
    ap.add_argument("--struct", default="structures")
    ap.add_argument("--n", type=int, default=10)
    ap.add_argument("--out", default="locked_pilot_overrep.csv")
    a = ap.parse_args()
    parser = MMCIFParser(QUIET=True)
    rows = list(csv.DictReader(open(a.labeled)))
    # 서열 추출 + dedup
    seen = {}   # (antigen,AB,hlseq) -> row(best res)
    uniq_ab = defaultdict(set)
    kept = []
    for r in rows:
        if r["AB"] not in ("A","B"): continue
        if any(b in (r.get("flags") or "") for b in BAD): continue
        if (r.get("epitope_class") or "") in ("head/stem-mixed","other/variable"): continue
        cif = os.path.join(a.struct, f"{r['pdb'].lower()}.cif")
        if not os.path.exists(cif): continue
        try:
            model = parser.get_structure(r["pdb"], cif)[0]
            hl = hlseq(model, r)
        except Exception: continue
        if len(hl.replace("|","")) < 100: continue      # 불완전 항체 제외
        key = (r["antigen"], r["AB"], hl)
        uniq_ab[(r["antigen"], r["AB"])].add(hl)
        r["_res"] = float(r.get("res", 9) or 9)
        if key not in seen or r["_res"] < seen[key]["_res"]:
            seen[key] = r
    print("=== 서열 dedup 후 unique 항체 수 ===")
    for g in ["RBD","HA","Env"]:
        print(f"  {g}: A={len(uniq_ab[(g,'A')])}  B={len(uniq_ab[(g,'B')])}")
    # 항원별 A/B 선별: 에피토프클래스 라운드로빈 + 해상도 우선
    pool = defaultdict(list)
    for r in seen.values(): pool[(r["antigen"], r["AB"])].append(r)
    picked = []
    for (g, ab), lst in pool.items():
        byc = defaultdict(list)
        for r in sorted(lst, key=lambda x: x["_res"]): byc[r["epitope_class"]].append(r)
        order = []
        while any(byc.values()) and len(order) < a.n:
            for c in list(byc.keys()):
                if byc[c]: order.append(byc[c].pop(0))
                if len(order) >= a.n: break
        picked += order[:a.n]
    # 출력
    cols = ["pdb","antigen","AB","epitope_class","Hchain","Lchain","antigen_chain","n_epi","res"]
    with open(a.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
        w.writeheader(); w.writerows(picked)
    print(f"\n=== 선별(항원별 A/B 각 {a.n}) → {a.out} ===")
    cnt = defaultdict(lambda: defaultdict(int))
    for r in picked: cnt[r["antigen"]][r["AB"]] += 1
    for g in ["RBD","HA","Env"]:
        print(f"  {g}: A={cnt[g]['A']} B={cnt[g]['B']}")
        for ab in ("A","B"):
            cc = defaultdict(int)
            for r in picked:
                if r["antigen"]==g and r["AB"]==ab: cc[r["epitope_class"]]+=1
            print(f"     {ab}: "+", ".join(f"{k}×{v}" for k,v in cc.items()))
    print(f"  총 {len(picked)} 복합체")

if __name__ == "__main__":
    main()
