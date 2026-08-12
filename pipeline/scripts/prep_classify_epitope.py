#!/usr/bin/env python3
"""각 복합체의 항원 접촉잔기(에피토프) → 에피토프 클래스 → A(지배)/B(비지배) 라벨.
설계: 분류기는 '객관적 클래스'를 뽑고, A/B는 epitope_defs의 매핑으로 분리(순환성 통제).
 - RBD: Wuhan auth 넘버링 직접. RBM(437-508+417) overlap 비율 → A(RBM) / B(class3·4).
 - HA:  사슬 정체 기반. HA2(fusion peptide GLFGAIAGF 함유/짧은 사슬) 접촉비율 → stem(B) / head(A).
        HA0 단일사슬이면 참조정렬 fallback + review 플래그.
 - Env: HXB2(P04578) 참조정렬 → bnAb class footprint overlap → class → A/B(ENV_DOMINANT_MODE).
출력: manifest_labeled.csv (pdb, antigen, AB, epitope_class, n_epi, frac, flags, epitope_residues).
의존: biopython. 실행 env 예) conda activate boltz (or DockQ env).
사용: python prep_classify_epitope.py --manifest manifest_candidates.csv --struct structures --out manifest_labeled.csv
"""
import argparse, csv, os, sys, urllib.request, warnings
from collections import defaultdict
import lib_epitope_defs as E
warnings.filterwarnings("ignore")
from Bio.PDB import MMCIFParser, NeighborSearch
from Bio.Align import PairwiseAligner, substitution_matrices

CUTOFF = 4.5
FP_MOTIF = "GLFGAIAGF"     # HA2 N-말단 fusion peptide (group1/2 근사 공통)
_3to1 = {  # 표준 20 + 흔한 변형
 'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H',
 'ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W',
 'TYR':'Y','VAL':'V','MSE':'M','SEC':'U','PYL':'O','HSD':'H','HSE':'H','HSP':'H','CSO':'C','SEP':'S','TPO':'T','PTR':'Y'}
def one(rn): return _3to1.get(rn.upper(), 'X')

_ALN = PairwiseAligner()
_ALN.substitution_matrix = substitution_matrices.load("BLOSUM62")
_ALN.open_gap_score = -11; _ALN.extend_gap_score = -1
_ALN.mode = "global"; _ALN.target_end_gap_score = 0; _ALN.query_end_gap_score = 0

def split_chains(field):
    return [c.strip() for c in str(field).replace("|", " ").split() if c.strip() and c.strip().upper() != "NA"]

def get_chain_residues(model, cid):
    """(auth_id) 순 aa 잔기 리스트 반환."""
    if cid not in model: return []
    out = []
    for res in model[cid]:
        if res.id[0].strip(): continue                 # HETATM(hetflag) skip
        if res.resname.upper() not in _3to1: continue
        out.append(res)
    return out

def chain_seq(residues):
    return "".join(one(r.resname) for r in residues), [r.id[1] for r in residues]

def epitope(model, ab_ids, ag_ids):
    """항원 접촉잔기: 항체 원자 CUTOFF 이내 항원 잔기. 반환 {ag_chain: {auth_id: resname}}."""
    ab_atoms = [a for cid in ab_ids if cid in model for r in model[cid] for a in r if not r.id[0].strip()]
    if not ab_atoms: return {}
    ns = NeighborSearch(ab_atoms)
    epi = defaultdict(dict)
    for cid in ag_ids:
        for r in get_chain_residues(model, cid):
            hit = any(ns.search(a.coord, CUTOFF) for a in r)
            if hit: epi[cid][r.id[1]] = r.resname
    return epi

# ── 참조 서열 fetch(캐시) ──
def fetch_ref(acc, cache="refs"):
    os.makedirs(cache, exist_ok=True)
    p = os.path.join(cache, f"{acc}.fasta")
    if not os.path.exists(p):
        urllib.request.urlretrieve(f"https://rest.uniprot.org/uniprotkb/{acc}.fasta", p)
    seq = "".join(l.strip() for l in open(p) if not l.startswith(">"))
    return seq

def map_to_ref(qseq, qids, ref):
    """query(auth_id) → ref position(1-based) 매핑 dict. 정렬 실패시 {}."""
    if not qseq or not ref: return {}
    try: aln = _ALN.align(ref, qseq)[0]
    except Exception: return {}
    m = {}; ri = qi = 0
    # aligned indices
    a_ref, a_q = aln.aligned  # blocks
    for (rs, re_), (qs, qe) in zip(a_ref, a_q):
        for k in range(re_ - rs):
            refpos = rs + k + 1           # 1-based ref
            qidx = qs + k
            if qidx < len(qids): m[qids[qidx]] = refpos
    return m

# ── 분류기 ──
def classify_rbd(epi):
    res = {rid for cid in epi for rid in epi[cid] if E.RBD_RANGE[0] <= rid <= E.RBD_RANGE[1]}
    if not res: return {"AB":"?","class":"none","n_epi":0,"frac":0.0,"flags":"no-RBD-contact","epi":res}
    f_rbm = len(res & E.RBD_RBM) / len(res)
    f_c3  = len(res & E.RBD_CLASS3) / len(res); f_c4 = len(res & E.RBD_CLASS4) / len(res)
    if f_rbm >= 0.5: cls, ab = ("class1/2(RBM)", "A")
    elif f_c4 >= 0.15: cls, ab = ("class4", "B")
    elif f_c3 >= 0.15: cls, ab = ("class3", "B")
    else: cls, ab = ("off/other", "B")
    flags = "" if (f_rbm>=0.5 or f_c3>=0.15 or f_c4>=0.15) else "ambiguous-lowoverlap"
    return {"AB":ab,"class":cls,"n_epi":len(res),"frac":round(f_rbm,2),"flags":flags,"epi":res}

def classify_ha(model, epi, ag_ids, ha_ref):
    # 참조(P03437) 정렬 → 접촉잔기를 참조 위치로 → head globular 범위 내 비율.
    # subtype·HA1only·단일사슬 HA0 모두 처리(사슬 정체에 의존 안 함).
    refpos = []
    for cid in ag_ids:
        residues = get_chain_residues(model, cid)
        qseq, qids = chain_seq(residues)
        if len(qseq) < 50: continue
        m = map_to_ref(qseq, qids, ha_ref)
        for rid in epi.get(cid, {}):
            if rid in m: refpos.append(m[rid])
    if not refpos:
        return {"AB":"?","class":"none","n_epi":0,"frac":0.0,"flags":"HA-no-map-review","epi":set()}
    lo, hi = E.HA_HEAD_REFRANGE
    f_head = sum(lo <= p <= hi for p in refpos) / len(refpos)
    if f_head >= 0.6: cls, ab = ("head", "A")
    elif f_head <= 0.35: cls, ab = ("stem", "B")
    else: cls, ab = ("head/stem-mixed", "?")
    return {"AB":ab,"class":cls,"n_epi":len(refpos),"frac":round(f_head,2),
            "flags":("borderline" if 0.35<f_head<0.6 else ""),"epi":set(refpos)}

def classify_env(model, epi, ag_ids, hxb2):
    # 각 항원사슬을 HXB2에 정렬 → epitope를 HXB2 위치로
    refpos = set()
    for cid in ag_ids:
        residues = get_chain_residues(model, cid)
        qseq, qids = chain_seq(residues)
        m = map_to_ref(qseq, qids, hxb2)
        for rid in epi.get(cid, {}):
            if rid in m: refpos.add(m[rid])
    if not refpos:
        return {"AB":"?","class":"none","n_epi":0,"frac":0.0,"flags":"no-map","epi":refpos}
    best, bo = None, 0
    for cls, s in E.ENV_CLASSES.items():
        ov = len(refpos & s)
        if ov > bo: bo, best = ov, cls
    if best is None or bo < 3:
        return {"AB":"B","class":"other/variable","n_epi":len(refpos),"frac":0.0,
                "flags":"low-classmatch","epi":refpos}
    ab = E.env_class_to_ab(best)
    return {"AB":ab,"class":best,"n_epi":len(refpos),"frac":round(bo/len(refpos),2),"flags":"","epi":refpos}

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", default="manifest_candidates.csv")
    ap.add_argument("--struct", default="structures")
    ap.add_argument("--out", default="manifest_labeled.csv")
    a = ap.parse_args()
    parser = MMCIFParser(QUIET=True)
    hxb2 = haref = None
    rows = list(csv.DictReader(open(a.manifest)))
    out = []
    counts = defaultdict(lambda: defaultdict(int))
    for i, r in enumerate(rows, 1):
        pid = r["pdb"].strip().lower(); antg = r["antigen"]
        cif = os.path.join(a.struct, f"{pid}.cif")
        rec = {"pdb":pid,"antigen":antg,"Hchain":r["Hchain"],"Lchain":r["Lchain"],
               "antigen_chain":r["antigen_chain"],"AB":"?","epitope_class":"","n_epi":0,
               "frac":0.0,"flags":"","epitope_residues":""}
        if not os.path.exists(cif):
            rec["flags"]="no-struct"; out.append(rec); continue
        try:
            model = parser.get_structure(pid, cif)[0]
            ab_ids = split_chains(r["Hchain"]) + split_chains(r["Lchain"])
            ag_ids = split_chains(r["antigen_chain"])
            epi = epitope(model, ab_ids, ag_ids)
            if antg == "RBD": res = classify_rbd(epi)
            elif antg == "HA":
                if haref is None: haref = fetch_ref(E.UNIPROT_REF["HA"])
                res = classify_ha(model, epi, ag_ids, haref)
            elif antg == "Env":
                if hxb2 is None: hxb2 = fetch_ref(E.UNIPROT_REF["Env"])
                res = classify_env(model, epi, ag_ids, hxb2)
            else: res = {"AB":"?","class":"?","n_epi":0,"frac":0,"flags":"unknown-antigen","epi":set()}
            rec.update(AB=res["AB"], epitope_class=res["class"], n_epi=res["n_epi"],
                       frac=res["frac"], flags=res["flags"],
                       epitope_residues=" ".join(map(str, sorted(res["epi"]))))
            counts[antg][res["AB"]] += 1
        except Exception as e:
            rec["flags"] = f"error:{type(e).__name__}"
        out.append(rec)
        if i % 25 == 0: print(f"  ... {i}/{len(rows)}")
    with open(a.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(out[0].keys())); w.writeheader(); w.writerows(out)
    print(f"\n=== A/B 요약 (Env dominant mode={E.ENV_DOMINANT_MODE}) ===")
    for antg in ["RBD","HA","Env"]:
        c = counts[antg]; print(f"  {antg:4s}: A={c['A']:3d}  B={c['B']:3d}  ?={c['?']:3d}")
    print(f"→ {a.out}")
    print("  ⚠️ HA 'HA2-unresolved-review'·RBD 'ambiguous'·Env 'low-classmatch' 행은 수동 확인.")

if __name__ == "__main__":
    main()
