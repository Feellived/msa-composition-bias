#!/usr/bin/env python3
"""HADDOCK 패치 → Protenix pocket-guided JSON (guided-Boltz의 Protenix판, §실험 4).
각 후보 에피토프 패치(HADDOCK 클러스터 또는 score 상위 N pose의 항원 접촉 잔기)를
Protenix constraint.pocket 으로 주입. 포맷은 make_protenix_msa_depth.py와 동일:
  항원 = entity 1 (unpairedMsaPath = A_d143.a3m), 항체 각 사슬 = entity 2.. (query-only a3m = single-seq).
pocket: binder_chain = 항체 heavy(entity 2), contact_residues = 항원(entity 1) 1-based 잔기, max_distance.
⚠️ Protenix pocket은 soft constraint(force 옵션 없음) — Boltz force:true보다 부드러운 steering.
비교 baseline = Protenix blind d143(runs_msad_143/<t>/out_protenix) 또는 server(summary.csv).
사용(protenix env): python scripts/haddock_to_protenix_pocket.py --targets "8XSI ..." [--topn 20] [--union]"""
import argparse, json, os
import build_msafree_summary as B
from haddock_to_boltz_pocket import load_gz, cluster_reps, topscore_poses, patch_from_pose


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--targets", default="8XSI 9SBB 8SIS 9ML8 8SDF 8SIT 9ML9")
    ap.add_argument("--src", default="runs_rbd")
    ap.add_argument("--max-patches", type=int, default=8)
    ap.add_argument("--max-distance", type=float, default=6.0)
    ap.add_argument("--min-contacts", type=int, default=4)
    ap.add_argument("--topn", type=int, default=0, help=">0이면 클러스터 대신 HADDOCK score 상위 N pose에서 patch 추출")
    ap.add_argument("--union", action="store_true", help="patch들을 하나로 합쳐 단일 pocket")
    ap.add_argument("--outdir", default="runs_guided_prot")
    a = ap.parse_args()
    for t in a.targets.split():
        t = t.upper(); ni = B.native_info(t)
        if ni is None: print(f"{t}: native 없음 skip"); continue
        cjp = os.path.join(a.src, t.lower(), "chains.json")
        if not os.path.exists(cjp): print(f"{t}: chains.json 없음 skip"); continue
        cj = json.load(open(cjp))
        ag = str(cj["antigen"]); ab = cj["antibody"]; ab = [ab] if isinstance(ab, str) else [str(x) for x in ab]
        sm = {c["id"]: c["seq"] for c in cj["chains"]}
        a3m = os.path.abspath(os.path.join("runs_msad_143", t.lower(), "A_d143.a3m"))
        if not os.path.exists(a3m):
            a3m = os.path.abspath(os.path.join(a.src, t.lower(), f"msa_{t.lower()}", "A.a3m"))
        if not os.path.exists(a3m): print(f"{t}: 항원 a3m 없음 skip"); continue
        reps = topscore_poses(t, a.topn) if a.topn > 0 else cluster_reps(t)
        src_lbl = f"top{a.topn} score pose" if a.topn > 0 else "클러스터"
        if not reps: print(f"{t}: {src_lbl} 없음 skip"); continue
        maxp = a.topn if a.topn > 0 else a.max_patches
        patches = []
        for p in reps:
            try: pos = patch_from_pose(load_gz(p), ni)
            except Exception: continue
            if len(pos) < a.min_contacts: continue
            if any(len(pos & q) / max(1, len(pos | q)) > 0.8 for q in patches): continue  # dedup Jaccard>0.8
            patches.append(pos)
            if len(patches) >= maxp: break
        if a.union and patches:
            patches = [set().union(*patches)]
        outd = os.path.join(a.outdir, t.lower()); os.makedirs(outd, exist_ok=True)
        # 항원=entity1, 항체 각 사슬=entity2.. (make_protenix_msa_depth와 동일 순서). binder=ab[0]=entity2.
        for k, pos in enumerate(patches):
            seqs = [{"proteinChain": {"sequence": sm[ag], "count": 1, "unpairedMsaPath": a3m}}]
            for i in ab:
                qa = os.path.abspath(os.path.join(outd, f"ab_{i}_query.a3m"))
                open(qa, "w").write(f">{i}\n{sm[i]}\n")   # 항체 single-seq
                seqs.append({"proteinChain": {"sequence": sm[i], "count": 1, "unpairedMsaPath": qa}})
            contacts = [{"entity": 1, "copy": 1, "position": int(r)} for r in sorted(pos)]
            job = {"name": t.lower(), "sequences": seqs,
                   "constraint": {"pocket": {
                       "binder_chain": {"entity": 2, "copy": 1},   # 첫 항체 사슬(heavy) = binder
                       "contact_residues": contacts,
                       "max_distance": a.max_distance}}}
            json.dump([job], open(os.path.join(outd, f"protenix_{t.lower()}_patch{k}.json"), "w"), indent=2)
        nres = len(patches[0]) if (a.union and patches) else "-"
        print(f"{t}: {src_lbl} {len(reps)} → 패치 {len(patches)}개 JSON"
              f"{f' (union {nres}잔기)' if a.union else ''} (binder=entity2={ab[0]}, 항원=entity1={ag}, soft pocket)")


if __name__ == "__main__": main()
