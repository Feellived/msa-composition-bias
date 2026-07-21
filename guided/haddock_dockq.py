#!/usr/bin/env python3
"""HADDOCK 포즈 DockQ: recall 상위 K개 포즈의 DockQ(near-native 검증) + CSV 저장.
   항체 H+L 병합. 사용(boltz, DockQ 설치): python scripts/haddock_dockq.py [--stage 4_emref|1_rigidbody] [--topk 15]"""
import argparse, csv, glob, gzip, json, os, re, subprocess, tempfile
import numpy as np
from Bio.PDB import PDBParser, PDBIO
from Bio.PDB.Structure import Structure
from Bio.PDB.Model import Model as BModel
from Bio.PDB.Chain import Chain as BChain
import build_msafree_summary as B

def load_gz(p):
    h = gzip.open(p, "rt") if p.endswith(".gz") else open(p)
    return PDBParser(QUIET=True).get_structure("x", h)[0]

def write_merged(model, ag_id, ab_ids, out_pdb):
    s = Structure("m"); m = BModel(0); s.add(m); A = BChain("A"); Bc = BChain("B")
    i = 1
    for r in model[ag_id]:
        if "CA" not in r: continue
        rr = r.copy(); rr.id = (" ", i, " "); A.add(rr); i += 1
    j = 1
    for aid in ab_ids:
        for r in model[aid]:
            if "CA" not in r: continue
            rr = r.copy(); rr.id = (" ", j, " "); Bc.add(rr); j += 1
    m.add(A); m.add(Bc); io = PDBIO(); io.set_structure(s); io.save(out_pdb)

def dockq(model_pdb, native_pdb):
    try:
        r = subprocess.run(["DockQ", model_pdb, native_pdb], capture_output=True, text=True, timeout=180)
    except Exception:
        return None
    out = (r.stdout or "") + (r.stderr or "")
    vals = re.findall(r"\bDockQ\b[:\s]+([0-9]*\.[0-9]+)", out)
    return max(float(v) for v in vals) if vals else None

def id_and_recall(pm, ni):
    pid, _, ps, pca, _ = B.best(pm, ni["ag_seq"])
    abids = [ch.id for ch in pm if ch.id != pid]
    pab = []
    for aid in abids:
        _, ca_arr, _ = B.seqca(pm[aid])
        if len(ca_arr): pab.append(ca_arr)
    pab = np.vstack(pab) if pab else np.zeros((0, 3))
    mp = B.amap(ps, ni["ag_seq"]); pred = set()
    for i in B.cset(pca, pab):
        if i in mp and mp[i] < len(ni["auth"]): pred.add(ni["auth"][mp[i]])
    rec = len(pred & ni["true"]) / len(ni["true"]) if ni["true"] else 0
    return rec, pid, abids

def native_merged(t, td):
    rd = os.path.join("runs_rbd", t.lower()); d = json.load(open(os.path.join(rd, "chains.json")))
    ag = str(d["antigen"]); ab = d["antibody"]; ab = [ab] if isinstance(ab, str) else [str(x) for x in ab]
    nat = next((os.path.join(rd, f) for f in os.listdir(rd) if t.lower() in f.lower() and f.endswith((".pdb", ".cif"))
                and not any(x in f.lower() for x in ("boltz", "chai", "protenix", "colabfold", "tfold"))), None)
    if not nat: return None
    nm = B.load(nat); ids = [c.id for c in nm]; sm = {c["id"]: c["seq"] for c in d["chains"]}
    ag_id = ag if ag in ids else B.best(nm, sm[ag])[0]
    ab_ids = [i for i in ab if i in ids]
    if not ab_ids: return None
    out = os.path.join(td, "nat.pdb"); write_merged(nm, ag_id, ab_ids, out); return out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", default="4_emref"); ap.add_argument("--topk", type=int, default=15)
    ap.add_argument("--targets", default=" ".join(B.TARGETS)); ap.add_argument("--out", default=None)
    a = ap.parse_args()
    skey = "emref" if "emref" in a.stage else ("rigidbody" if "rigidbody" in a.stage else a.stage)
    out = a.out or f"results/haddock_dockq_{skey}.csv"
    print(f"[stage={a.stage}, DockQ on top{a.topk} by recall]")
    print(f"{'target':7}{'class':9}{'best_DockQ':11}{'(recall)':10}{'#DockQ>=.23'}")
    print("-"*52)
    rows = []
    for t in a.targets.split():
        t = t.upper(); ni = B.native_info(t)
        if ni is None: print(f"{t:7} native 없음"); continue
        poses = sorted(glob.glob(os.path.join("haddock", t.lower(), "run", a.stage, "*.pdb.gz")))
        scored = []
        for p in poses:
            try:
                rec, pid, abids = id_and_recall(load_gz(p), ni); scored.append((rec, p, pid, abids))
            except Exception: continue
        scored.sort(reverse=True)
        with tempfile.TemporaryDirectory() as td:
            natm = native_merged(t, td)
            if natm is None: print(f"{t:7} native merge 실패"); continue
            bq = None; brec = None; nhit = 0; ndone = 0
            for rec, p, pid, abids in scored[:a.topk]:
                try:
                    mp = os.path.join(td, "m.pdb"); write_merged(load_gz(p), pid, abids, mp)
                    q = dockq(mp, natm); ndone += 1
                    if q is not None:
                        if q >= 0.23: nhit += 1
                        if bq is None or q > bq: bq = q; brec = rec
                except Exception: continue
            bs = f"{bq:.3f}" if bq is not None else "  -  "
            rs = f"({brec:.2f})" if brec is not None else "-"
            print(f"{t:7}{B.EPICLASS.get(t,''):9}{bs:>10} {rs:>9}   {nhit}/{ndone}")
            rows.append({"target": t, "epitope_class": B.EPICLASS.get(t, ""), "stage": skey,
                         "best_dockq": f"{bq:.3f}" if bq is not None else "", "best_recall": f"{brec:.2f}" if brec is not None else "",
                         "n_hit": nhit, "n_scored": ndone})
    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["target", "epitope_class", "stage", "best_dockq", "best_recall", "n_hit", "n_scored"])
        w.writeheader(); w.writerows(rows)
    print(f"\n→ {out} 저장 ({len(rows)}행). 스코어카드: python scripts/haddock_scorecard.py")

if __name__ == "__main__": main()
