#!/usr/bin/env python3
"""locked_pilot.csv + structures/ → 복합체별 생성 입력 준비.
각 복합체: 항체 H(+L) + 항원 서열 추출, 역할 배정, chains.json + 항원 fasta(MSA용) 작성.
항원 = 항체가 가장 많이 접촉하는 단백질 사슬 1개(trimer면 한 protomer). 큰 항원은 길이 flag.
사용: python prep_targets.py --struct structures --outdir targets [--csv locked_pilot.csv]
출력 targets/<id>/: chains.json, antigen.fasta, native.cif(심링크)
"""
import argparse, csv, json, os, warnings
warnings.filterwarnings("ignore")
from collections import defaultdict
from Bio.PDB import MMCIFParser, NeighborSearch
import classify_epitope as C

def contacted_antigen_chains(model, ab_ids, ag_ids, min_contacts=3):
    """항체가 접촉하는 항원 단백질 사슬 전부(접촉순). HA1+HA2·gp120+gp41 같은 다사슬 항원 유지."""
    ab_atoms = [a for cid in ab_ids if cid in model for r in model[cid] for a in r if not r.id[0].strip()]
    if not ab_atoms: return []
    ns = NeighborSearch(ab_atoms)
    hits = []
    for cid in ag_ids:
        res = C.get_chain_residues(model, cid)
        if len(res) < 25: continue
        n = sum(1 for r in res if any(ns.search(a.coord, 5.0) for a in r))
        if n >= min_contacts: hits.append((cid, n))
    return [c for c, _ in sorted(hits, key=lambda x: -x[1])]

def chain_seq_cropped(model, cid, lo=None, hi=None):
    """사슬 서열(옵션: auth 잔기번호 lo~hi 크롭). 반환 (seq, kept_ids)."""
    res = [r for r in C.get_chain_residues(model, cid) if (lo is None or lo <= r.id[1] <= hi)]
    return "".join(C.one(r.resname) for r in res), [r.id[1] for r in res]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", default="locked_pilot.csv")
    ap.add_argument("--struct", default="structures")
    ap.add_argument("--outdir", default="targets")
    a = ap.parse_args()
    parser = MMCIFParser(QUIET=True)
    os.makedirs(a.outdir, exist_ok=True)
    rows = list(csv.DictReader(open(a.csv)))
    lens = []
    ok = 0
    idx = defaultdict(int)
    summary = []
    for r in rows:
        pid = r["pdb"].strip().lower()
        cif = os.path.join(a.struct, f"{pid}.cif")
        # 한 PDB 여러 항체 → id에 항체사슬 붙여 유일화
        tid = f"{pid}_{r['Hchain'].strip()}{r['Lchain'].strip()}".replace(" ", "").replace("|", "")
        rec = {"target": tid, "pdb": pid, "antigen_grp": r["antigen"], "AB": r["AB"], "label": r.get("label", "")}
        if not os.path.exists(cif):
            rec["status"] = "no-struct"; summary.append(rec); continue
        try:
            model = parser.get_structure(pid, cif)[0]
            H = C.split_chains(r["Hchain"]); L = C.split_chains(r["Lchain"])
            AG = C.split_chains(r["antigen_chain"])
            agcs = contacted_antigen_chains(model, H + L, AG)
            if not agcs: rec["status"] = "no-antigen-chain"; summary.append(rec); continue
            # RBD 그룹: 큰 사슬(full spike)만 RBD 도메인(319-541 Wuhan)으로 크롭
            crop = (319, 541) if r["antigen"] == "RBD" else None
            ag_chains, ag_total = [], 0
            for k, cid in enumerate(agcs[:2]):           # 최대 2 사슬(HA1+HA2, gp120+gp41)
                lo = hi = None
                if crop:
                    full = C.chain_seq(C.get_chain_residues(model, cid))[0]
                    if len(full) > 400: lo, hi = crop    # RBD 구성체(작음)는 그대로
                seq, ids = chain_seq_cropped(model, cid, lo, hi)
                if len(seq) < 25: continue
                aid = chr(ord("A") + k)
                ag_chains.append({"id": aid, "role": "antigen", "seq": seq, "src": cid,
                                  "crop": [lo, hi] if lo else None})
                ag_total += len(seq)
            if not ag_chains: rec["status"] = "antigen-empty-after-crop"; summary.append(rec); continue
            h_seq = "".join(C.chain_seq(C.get_chain_residues(model, c))[0] for c in H)
            l_seq = "".join(C.chain_seq(C.get_chain_residues(model, c))[0] for c in L)
            if len(h_seq) < 80: rec["status"] = "short-heavy"; summary.append(rec); continue
            td = os.path.join(a.outdir, tid); os.makedirs(td, exist_ok=True)
            hid = chr(ord("A") + len(ag_chains))
            chains = ag_chains + [{"id": hid, "role": "heavy", "seq": h_seq}]
            lid = None
            if l_seq:
                lid = chr(ord("A") + len(ag_chains) + 1)
                chains.append({"id": lid, "role": "light", "seq": l_seq})
            ag_ids = [c["id"] for c in ag_chains]
            cj = {"pdb_id": pid, "target": tid, "antigen_grp": r["antigen"], "AB": r["AB"],
                  "label": r.get("label", ""), "antigen": ag_ids, "antibody": [hid] + ([lid] if lid else []),
                  "chains": chains, "src_chains": {"antigen": agcs[:2], "H": H, "L": L}}
            json.dump(cj, open(os.path.join(td, "chains.json"), "w"), indent=2)
            with open(os.path.join(td, "antigen.fasta"), "w") as f:
                for c in ag_chains: f.write(f">{tid}_{c['id']}\n{c['seq']}\n")
            src = os.path.abspath(cif); ln = os.path.join(td, "native.cif")
            if not os.path.exists(ln): os.symlink(src, ln)
            lens.append(ag_total); ok += 1
            rec.update(status="ok", ag_chains=ag_ids, ag_len=ag_total, n_ag=len(ag_chains),
                       h_len=len(h_seq), l_len=len(l_seq), nano=(not l_seq))
        except Exception as e:
            rec["status"] = f"error:{type(e).__name__}"
        summary.append(rec)
    json.dump(summary, open(os.path.join(a.outdir, "_prep_summary.json"), "w"), indent=2)
    print(f"=== prep: {ok}/{len(rows)} 성공 → {a.outdir}/<id>/ ===")
    st = defaultdict(int)
    for s in summary: st[s["status"]] += 1
    print("  status:", dict(st))
    if lens:
        lens.sort()
        import statistics as S
        print(f"  항원 길이: min={lens[0]} median={int(S.median(lens))} max={lens[-1]}  (>400={sum(l>400 for l in lens)}개 = 크롭 검토)")
    # 항원그룹×AB 성공 카운트
    cc = defaultdict(lambda: defaultdict(int))
    for s in summary:
        if s.get("status") == "ok": cc[s["antigen_grp"]][s["AB"]] += 1
    for g in ["RBD", "HA", "Env", "C"]:
        print(f"  {g}: "+", ".join(f"{k}={v}" for k, v in cc[g].items()))

if __name__ == "__main__":
    main()
