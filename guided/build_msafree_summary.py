#!/usr/bin/env python3
"""MSA-free 실험 통합 summary CSV (long-format, source of truth).
co-folder(runs_rbd) + MSA-depth(runs_msad_*) + HADDOCK(haddock/*/run) → 한 표.
각 (target,experiment,model,condition): epitope_recall, rbm_frac(=RBM 겹침), n_poses.
native/true는 항상 runs_rbd/<t>에서. DockQ는 --merged로 병합(선택).
사용(boltz env): python build_msafree_summary.py --out results/msafree_summary.csv"""
import argparse, csv, glob, json, os
from collections import Counter
import numpy as np
from Bio.PDB import MMCIFParser, PDBParser
from Bio.Align import PairwiseAligner
from Bio.Data.IUPACData import protein_letters_3to1
T3TO1={k.upper():v for k,v in protein_letters_3to1.items()}; EPI=8.0
_al=PairwiseAligner();_al.mode="global";_al.match_score=1;_al.mismatch_score=-1;_al.open_gap_score=-3;_al.extend_gap_score=-0.5
TARGETS="8P5M 8SDF 8SIQ 8SIS 8SIT 8XSI 9ML8 9ML9 9SBB 9ZDU".split()
EPICLASS={"8SIT":"core","9ML9":"core","9ZDU":"core","8SIQ":"core","8SIS":"core_other",
          "8SDF":"cryptic","8XSI":"cryptic","9ML8":"cryptic","9SBB":"other","8P5M":"RBM_edge"}
def is_rbm(auth): return 437<=auth<=508
def load(p): return (MMCIFParser(QUIET=True) if p.endswith(".cif") else PDBParser(QUIET=True)).get_structure("x",p)[0]
def seqca(ch):
    s,c,rr=[],[],[]
    for r in ch:
        if "CA" in r:
            aa=T3TO1.get(r.get_resname().upper())
            if aa: s.append(aa); c.append(r["CA"].get_coord()); rr.append(r)
    return "".join(s),(np.array(c) if c else np.zeros((0,3))),rr
def best(m,tgt,exclude=()):
    b=(None,-1e9,"",np.zeros((0,3)),[])
    for ch in m:
        if ch.id in exclude: continue
        s,ca,rr=seqca(ch)
        if len(s)<5: continue
        sc=_al.score(s,tgt)/max(len(s),len(tgt),1)
        if sc>b[1]: b=(ch.id,sc,s,ca,rr)
    return b
def amap(sf,st):
    if not sf or not st: return {}
    aln=_al.align(sf,st)[0]; m={}
    for (a0,a1),(b0,b1) in zip(aln.aligned[0],aln.aligned[1]):
        for k in range(a1-a0): m[a0+k]=b0+k
    return m
def cset(ag,ab):
    if len(ag)==0 or len(ab)==0: return set()
    d=np.sqrt(((ag[:,None,:]-ab[None,:,:])**2).sum(-1)); return set(np.where(d.min(1)<=EPI)[0])

def native_info(t):
    """runs_rbd/<t>: 항원서열·native 잔기(auth)·true 에피토프(auth집합)."""
    rd=os.path.join("runs_rbd",t.lower()); d=json.load(open(os.path.join(rd,"chains.json")))
    ag=str(d["antigen"]); ab=d["antibody"]; ab=[ab] if isinstance(ab,str) else [str(x) for x in ab]
    sm={c["id"]:c["seq"] for c in d.get("chains",[])}; ag_seq=sm.get(ag,""); ab_seqs=[sm.get(i,"") for i in ab]
    nat=[f for f in os.listdir(rd) if t.lower() in f.lower() and f.endswith((".pdb",".cif"))
         and not any(m in f.lower() for m in ("boltz","chai","protenix","colabfold","tfold"))]
    if not nat or not ag_seq: return None
    m=load(os.path.join(rd,nat[0])); ids=[c.id for c in m]
    if ag in ids: _,nca,nrr=seqca(m[ag])
    else: _,_,_,nca,nrr=best(m,ag_seq)
    nab=np.vstack([seqca(m[i])[1] for i in ab if i in ids]) if any(i in ids for i in ab) else np.zeros((0,3))
    auth=[r.id[1] for r in nrr]  # native auth 번호(순서=서열 index)
    true={auth[i] for i in cset(nca,nab)}
    return {"ag_seq":ag_seq,"ab_seqs":ab_seqs,"auth":auth,"true":true}

def pred_epitope(run_dir, t, ag_seq, ab_seqs, auth, glob_pat):
    """pose들 majority-vote 예측 에피토프 → native auth 집합 + n_poses."""
    rd=os.path.join(run_dir, t.lower())
    votes=Counter(); n=0
    for pose in sorted(glob.glob(os.path.join(rd,glob_pat),recursive=True)):
        try:
            pm=load(pose); n+=1
            pid,_,ps,pca,_=best(pm,ag_seq); used={pid}; pab=[]
            for s in ab_seqs:
                r=best(pm,s,exclude=used)
                if r[0]: used.add(r[0]); pab.append(r[3])
            pab=np.vstack(pab) if pab else np.zeros((0,3))
            mp=amap(ps,ag_seq)  # pose 항원 idx → native 항원 서열 idx
            for i in cset(pca,pab):
                if i in mp and mp[i]<len(auth): votes[auth[mp[i]]]+=1
        except Exception: continue
    if n==0: return None,0
    thr=max(1,(n+1)//2); pred={a for a,v in votes.items() if v>=thr}
    return pred,n

# (experiment, model, condition, run_dir, glob)
SOURCES=[
  ("cofolder","Boltz-2","server","runs_rbd","*boltz/**/*_model_*.cif"),
  ("cofolder","Chai-1","server","runs_rbd","out_chai/pred.model_idx_*.cif"),
  ("cofolder","Protenix-base","server","runs_rbd","out_protenix/**/*_sample_*.cif"),
  ("cofolder","AF2-Multimer","server","runs_rbd","*colabfold/*_unrelaxed_rank_*.pdb"),
  ("cofolder","tFold-Ag","server","runs_rbd","out_tfold/*seed*.pdb"),
  ("msa_depth","Boltz-2","d143","runs_msad_143","*boltz/**/*_model_*.cif"),
  ("msa_depth","Boltz-2","d8","runs_msad_8","*boltz/**/*_model_*.cif"),
  ("msa_depth","Boltz-2","d1","runs_msad_1","*boltz/**/*_model_*.cif"),
  ("msa_depth","Protenix-base","d143","runs_msad_143","out_protenix/**/*sample*.cif"),
  ("msa_depth","Protenix-base","d8","runs_msad_8","out_protenix/**/*sample*.cif"),
  ("msa_depth","Protenix-base","d1","runs_msad_1","out_protenix/**/*sample*.cif"),
  ("msa_depth","Chai-1","d1","runs_msad_1","out_chai/**/*.cif"),
  ("msa_depth","AF2-Multimer","d1","runs_msad_1","out_colabfold/*rank*.pdb"),
  ("msa_depth","tFold","d8","runs_msad_8","out_tfold/*seed*.pdb"),
  ("msa_depth","tFold","d1","runs_msad_1","out_tfold/*seed*.pdb"),
  ("haddock","HADDOCK","ab-initio","haddock","run/*seletopclusts*/*.pdb"),
]
def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--out",default="results/msafree_summary.csv"); a=ap.parse_args()
    os.makedirs(os.path.dirname(a.out) or ".",exist_ok=True)
    nat={t:native_info(t) for t in TARGETS}
    rows=[]
    for exp,model,cond,rundir,gl in SOURCES:
        for t in TARGETS:
            ni=nat.get(t)
            if not ni: continue
            pred,n=pred_epitope(rundir,t,ni["ag_seq"],ni["ab_seqs"],ni["auth"],gl)
            if pred is None: continue
            true=ni["true"]
            recall=len(pred&true)/len(true) if true else None
            rbm=len([a for a in pred if is_rbm(a)])/len(pred) if pred else None
            rows.append({"target":t,"epitope_class":EPICLASS.get(t,""),"experiment":exp,"model":model,
                         "condition":cond,"epitope_recall":f"{recall:.3f}" if recall is not None else "",
                         "rbm_frac":f"{rbm:.3f}" if rbm is not None else "","n_poses":n,"dockq":""})
    with open(a.out,"w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=["target","epitope_class","experiment","model","condition","epitope_recall","rbm_frac","n_poses","dockq"])
        w.writeheader(); w.writerows(rows)
    print(f"{len(rows)}행 → {a.out}")
    print("experiment별:", Counter(r["experiment"] for r in rows))
if __name__=="__main__": main()
