#!/usr/bin/env python3
"""manifest_candidates.csv의 PDB를 RCSB에서 mmCIF로 내려받는다(서버 실행).
확장ID(pdb_0000xxxx)는 4char로 이미 매핑됨(build_manifest). 재실행 self-heal(있으면 skip).
사용: python prep_fetch_structures.py [--manifest manifest_candidates.csv] [--outdir structures]
"""
import argparse, csv, os, time, urllib.request, urllib.error

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", default="manifest_candidates.csv")
    ap.add_argument("--outdir", default="structures")
    ap.add_argument("--sleep", type=float, default=0.2)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    ids = []
    with open(a.manifest) as f:
        for r in csv.DictReader(f):
            ids.append(r["pdb"].strip())
    ids = sorted(set(ids))
    print(f"[fetch] {len(ids)} distinct PDB → {a.outdir}/")
    ok = miss = skip = 0
    for i, pid in enumerate(ids, 1):
        out = os.path.join(a.outdir, f"{pid.lower()}.cif")
        if os.path.exists(out) and os.path.getsize(out) > 1000:
            skip += 1; continue
        url = f"https://files.rcsb.org/download/{pid}.cif"
        try:
            urllib.request.urlretrieve(url, out)
            ok += 1
        except Exception as e:
            miss += 1; print(f"  !! {pid} 실패: {e}")
            if os.path.exists(out): os.remove(out)
        if i % 25 == 0: print(f"  ... {i}/{len(ids)} (ok={ok} skip={skip} miss={miss})")
        time.sleep(a.sleep)
    print(f"[done] ok={ok} skip={skip} miss={miss} → {a.outdir}/")

if __name__ == "__main__":
    main()
