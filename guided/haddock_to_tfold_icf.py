#!/usr/bin/env python3
"""HADDOCK 패치 → tFold-Ag ICF(.pt) guided (guided-Boltz/Protenix의 tFold판, §실험 4).
tFold는 pocket 대신 ICF(inter-chain feature) 텐서로 계면 지정:
  calc_ppi_sites 포맷 = 1D int8, 길이 L_H+L_L+L_A, ordering [H,L,A], 잔기가 상대 CA 10Å내면 1.
  epitope 모드 = 항체부분 0, 항원(에피토프)만 1  → predict.py --icf X.pt --model_version ppi.
HADDOCK 패치(항원 접촉잔기, patch_from_pose=1-based 서열위치)를 icf[L_ab + (pos-1)]=1 로 직접 구성.
⚠️ torch + biopython 필요(tfold env 권장). ⚠️ tFold JSON 모드는 icf_path 버그 → FASTA+--icf 경로 사용.
사용(tfold env): python scripts/haddock_to_tfold_icf.py --targets "8XSI 9SBB ..." [--topn 20] [--union]"""
import argparse, json, os
import torch
import build_msafree_summary as B
from haddock_to_boltz_pocket import load_gz, cluster_reps, topscore_poses, patch_from_pose


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--targets", default="8XSI 9SBB 8SIS 9ML8 8SDF 8SIT 9ML9")
    ap.add_argument("--src", default="runs_rbd")
    ap.add_argument("--max-patches", type=int, default=8)
    ap.add_argument("--min-contacts", type=int, default=4)
    ap.add_argument("--topn", type=int, default=0, help=">0이면 클러스터 대신 HADDOCK score 상위 N pose")
    ap.add_argument("--union", action="store_true", help="patch들을 하나로 합쳐 단일 ICF")
    ap.add_argument("--outdir", default="runs_guided_tfold")
    a = ap.parse_args()
    for t in a.targets.split():
        t = t.upper(); ni = B.native_info(t)
        if ni is None: print(f"{t}: native 없음 skip"); continue
        cjp = os.path.join(a.src, t.lower(), "chains.json")
        if not os.path.exists(cjp): print(f"{t}: chains.json 없음 skip"); continue
        cj = json.load(open(cjp))
        ag = str(cj["antigen"]); ab = cj["antibody"]; ab = [ab] if isinstance(ab, str) else [str(x) for x in ab]
        sm = {c["id"]: c["seq"] for c in cj["chains"]}
        # tFold 순서 H=ab[0], L=ab[1](Fab), A=ag
        L_H = len(sm[ab[0]]); L_L = len(sm[ab[1]]) if len(ab) >= 2 else 0; L_A = len(sm[ag])
        L_ab = L_H + L_L
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
        # tfold fasta(>H,>L,>A) — make_tfold와 동일 매핑
        recs = [("H", sm[ab[0]])] + ([("L", sm[ab[1]])] if L_L else []) + [("A", sm[ag])]
        with open(os.path.join(outd, f"tfold_{t.lower()}.fasta"), "w") as f:
            for cid, seq in recs: f.write(f">{cid}\n{seq}\n")
        for k, pos in enumerate(patches):
            icf = torch.zeros(L_ab + L_A, dtype=torch.int8)
            nmark = 0
            for r in sorted(pos):
                i = int(r) - 1                       # 1-based 서열위치 → 0-based 항원 인덱스
                if 0 <= i < L_A: icf[L_ab + i] = 1; nmark += 1
            torch.save(icf, os.path.join(outd, f"tfold_{t.lower()}_patch{k}.pt"))
        nres = len(patches[0]) if (a.union and patches) else "-"
        print(f"{t}: {src_lbl} {len(reps)} → ICF {len(patches)}개(.pt) "
              f"[L_ab={L_ab} L_A={L_A}]{f' union {nres}잔기' if a.union else ''} (epitope-only, model_version=ppi)")


if __name__ == "__main__": main()
