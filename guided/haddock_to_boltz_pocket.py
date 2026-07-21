#!/usr/bin/env python3
"""HADDOCK ab-initio 클러스터 → Boltz pocket-guided YAML (staged refine 실험, §실험 4).
각 클러스터 대표 pose의 항원 접촉 잔기(=후보 에피토프 패치)를 Boltz pocket 제약으로 주입.
YAML = 항원 full a3m(runs_msad_143의 A_d143.a3m) + 항체 empty + pocket(binder=항체 chain, contacts=항원 에피토프).
비교 baseline = runs_msad_143 blind Boltz (같은 MSA, 제약만 차이 → 순수 '힌트 효과' 격리).
사용(boltz env): python scripts/haddock_to_boltz_pocket.py --targets "8XSI 9SBB 8SIS 9ML8 8SDF 8SIT 9ML9" [--max-patches 8]"""
import argparse, csv, glob, gzip, json, os, re
import numpy as np
from Bio.PDB import PDBParser
import build_msafree_summary as B

def load_gz(p):
    h = gzip.open(p, "rt") if p.endswith(".gz") else open(p)
    return PDBParser(QUIET=True).get_structure("x", h)[0]

def topscore_poses(t, n):
    """클러스터링 전 — HADDOCK score 상위 N개 emref pose 경로(점수순). clustfcc.tsv에서 순위."""
    rundir = os.path.join("haddock", t.lower(), "run"); emdir = os.path.join(rundir, "4_emref")
    tsv = None
    for f in glob.glob(os.path.join(rundir, "*clustfcc*", "*.tsv")) + glob.glob(os.path.join(rundir, "**", "*.tsv"), recursive=True):
        try:
            hdr = next(csv.reader(open(f), delimiter="\t"))
            if {"model_name", "score"} <= set(hdr): tsv = f; break
        except Exception: pass
    if not tsv: return []
    rows = sorted(csv.DictReader(open(tsv), delimiter="\t"), key=lambda r: float(r["score"]))
    out = []
    for r in rows[:n]:
        p = os.path.join(emdir, r["model_name"] + ".gz")
        if not os.path.exists(p): p = os.path.join(emdir, r["model_name"])
        if os.path.exists(p): out.append(p)
    return out

def cluster_reps(t):
    """haddock/<t>/run/*seletopclusts*/ 각 클러스터 대표(model_1) pose 경로."""
    base = os.path.join("haddock", t.lower(), "run")
    sd = None
    for cand in sorted(glob.glob(os.path.join(base, "*seletopclusts*"))):
        sd = cand  # 마지막(=최종 seletopclusts)
    if not sd: return []
    reps = sorted(glob.glob(os.path.join(sd, "cluster_*_model_1.pdb.gz"))) + \
           sorted(glob.glob(os.path.join(sd, "cluster_*_model_1.pdb")))
    if not reps:  # model_1 없으면 클러스터별 아무 model 하나
        seen = {}
        for p in sorted(glob.glob(os.path.join(sd, "cluster_*.pdb.gz"))) + sorted(glob.glob(os.path.join(sd, "cluster_*.pdb"))):
            m = re.search(r"cluster_(\d+)_model", os.path.basename(p))
            seen.setdefault(m.group(1) if m else p, p)
        reps = list(seen.values())
    return reps

def patch_from_pose(pm, ni):
    """pose 항원 접촉 잔기 → Boltz 항원 서열 1-based 위치 집합. (best()로 서열매칭 → 사슬라벨 무관)"""
    pid, _, ps, pca, _ = B.best(pm, ni["ag_seq"])
    pab = []
    for ch in pm:
        if ch.id == pid: continue
        _, ca, _ = B.seqca(ch)
        if len(ca): pab.append(ca)
    pab = np.vstack(pab) if pab else np.zeros((0, 3))
    mp = B.amap(ps, ni["ag_seq"])   # pose 항원 idx → native 항원 서열 idx(0-based)
    return {mp[i] + 1 for i in B.cset(pca, pab) if i in mp}   # 1-based (Boltz 서열 위치)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--targets", default="8XSI 9SBB 8SIS 9ML8 8SDF 8SIT 9ML9")
    ap.add_argument("--src", default="runs_rbd")
    ap.add_argument("--max-patches", type=int, default=8)
    ap.add_argument("--max-distance", type=float, default=6.0)
    ap.add_argument("--min-contacts", type=int, default=4)
    ap.add_argument("--soft", action="store_true", help="force:false(소프트 조건화). 기본=force:true(패치로 실제 steering; Boltz #418 회피)")
    ap.add_argument("--topn", type=int, default=0, help=">0이면 클러스터 대신 HADDOCK score 상위 N개 emref pose에서 patch 추출(클러스터링 우회)")
    ap.add_argument("--union", action="store_true", help="추출한 patch들을 하나로 합쳐 단일 pocket 제약(넓은 영역 → --soft 권장)")
    ap.add_argument("--outdir", default="runs_guided")
    a = ap.parse_args()
    for t in a.targets.split():
        t = t.upper(); ni = B.native_info(t)
        if ni is None: print(f"{t}: native 없음 skip"); continue
        cjp = os.path.join(a.src, t.lower(), "chains.json")
        if not os.path.exists(cjp): print(f"{t}: chains.json 없음 skip"); continue
        cj = json.load(open(cjp))
        ag = str(cj["antigen"]); ab = cj["antibody"]; ab = [ab] if isinstance(ab, str) else [str(x) for x in ab]
        binder = str(ab[0])   # 첫 항체 사슬(보통 heavy) = binder
        a3m = os.path.abspath(os.path.join("runs_msad_143", t.lower(), "A_d143.a3m"))
        if not os.path.exists(a3m):
            a3m = os.path.abspath(os.path.join(a.src, t.lower(), f"msa_{t.lower()}", "A.a3m"))
        if not os.path.exists(a3m): print(f"{t}: a3m 없음 skip"); continue
        reps = topscore_poses(t, a.topn) if a.topn > 0 else cluster_reps(t)
        src = f"top{a.topn} score pose" if a.topn > 0 else "클러스터"
        if not reps: print(f"{t}: {src} 없음 skip"); continue
        maxp = a.topn if a.topn > 0 else a.max_patches
        patches = []
        for p in reps:
            try: pos = patch_from_pose(load_gz(p), ni)
            except Exception: continue
            if len(pos) < a.min_contacts: continue
            if any(len(pos & q) / max(1, len(pos | q)) > 0.8 for q in patches): continue  # dedup(Jaccard>0.8)
            patches.append(pos)
            if len(patches) >= maxp: break
        if a.union and patches:   # top-N 접촉 전체를 한 pocket으로
            patches = [set().union(*patches)]
        outd = os.path.join(a.outdir, t.lower()); os.makedirs(outd, exist_ok=True)
        for k, pos in enumerate(patches):
            contacts = "\n".join(f"        - [{ag}, {r}]" for r in sorted(pos))
            L = ["version: 1", "sequences:"]
            for c in cj["chains"]:
                L += ["  - protein:", f"      id: {c['id']}", f'      sequence: "{c["seq"]}"',
                      f"      msa: {a3m if c['id'] == ag else 'empty'}"]
            L += ["constraints:", "  - pocket:", f"      binder: {binder}", "      contacts:", contacts,
                  f"      max_distance: {a.max_distance}", f"      force: {'false' if a.soft else 'true'}"]
            open(os.path.join(outd, f"boltz_{t.lower()}_patch{k}.yaml"), "w").write("\n".join(L) + "\n")
        nres = len(patches[0]) if (a.union and patches) else "-"
        print(f"{t}: {src} {len(reps)} → 패치 {len(patches)}개 YAML"
              f"{f' (union {nres}잔기)' if a.union else ''} (binder={binder}, 항원={ag}, force={'false' if a.soft else 'true'})")

if __name__ == "__main__": main()
