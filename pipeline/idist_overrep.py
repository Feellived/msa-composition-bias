#!/usr/bin/env python3
"""iDist 계면 과대표집 — 각 테스트 항체-항원 계면이 '컷 이전 SAbDab'에 얼마나 흔한가 = 과대표집 점수.
가설: 우세부위(A) 계면은 컷 이전 PDB에 near-duplicate가 많고(과대표집 큼), 비우세(B)는 적다.
     이 점수가 우리 47복합체의 '깊이-취약성'과 상관되면 = 위치편향 인과 축을 데이터로 확정.

도구 = PPIRef의 iDist (Bushuiev et al., ICLR 2024, "Revealing data leakage in PPI benchmarks").
  설치(순수 파이썬, 외부 바이너리 불필요): pip install git+https://github.com/anton-bushuiev/PPIRef.git
  ※ 방법 novelty 없음 — iDist/iAlign은 기성 선점 도구. 여기선 '항체-항원 계면 과대표집 정량' 응용·재현으로만 인용.

흐름:
  (1) native → 2체인 PDB 병합(항원=A, 항체 H+L=B; dockq_sweep.merged_pdb 로직 재사용 → 내부 VH-VL 계면 제외)
  (2) PPIExtractor(kind='heavy', radius=6.0)로 항체-항원 계면 추출  ← radius=6.0 ↔ near-dup 임계 0.04 반드시 짝
  (3) IDist로 레퍼런스 계면 임베딩 → 각 테스트 계면 임베딩 → 유클리드 거리로 이웃수(near-dup + 반경들) 산출
매칭 = 항원군(RBD/HA/Env) 단위. 테스트=우리 A/B 복합체, 레퍼런스=같은 항원군의 컷 이전 SAbDab 항체.

사용(ppiref env, 서버):
  # 0) 레퍼런스 구조 준비: python idist_ref_manifest.py  →  python fetch_structures.py --manifest ref_manifest.csv --outdir ref_structures
  python idist_overrep.py --stage extract           # 테스트+레퍼런스 계면 추출(재사용 캐시)
  python idist_overrep.py --stage score             # 임베딩→이웃수→ results/overrep_idist.csv
  (또는 --stage all 로 한 번에)
"""
import argparse, csv, glob, json, os, re, sys, tempfile
import numpy as np
from Bio.PDB import MMCIFParser, PDBParser, PDBIO
from Bio.PDB.Structure import Structure
from Bio.PDB.Model import Model as BModel
from Bio.PDB.Chain import Chain as BChain

RADIUS = 6.0
NEAR_DUP = 0.04                 # 6Å 계면 전용(공식 보정값). 10Å면 0.03.
# 실측 거리 스케일(진단: 동일 항원군 계면 간 거리 min 0.03 ~ max 0.14) — 0.15·0.30은 포화(전부 근접)라 무의미.
# → near-dup(0.04) 주변을 촘촘히 + min_dist·최근접평균(연속값)을 주 과대표집 신호로.
RADII = [0.02, 0.03, 0.04, 0.05, 0.06]
KNN = 5                         # 최근접 k개 평균거리 = best-of-N 잡음에 덜 민감한 연속 과대표집 신호

# ---------- 구조 병합 (dockq_sweep.merged_pdb 재사용) ----------
def load(p): return (MMCIFParser(QUIET=True) if p.endswith(".cif") else PDBParser(QUIET=True)).get_structure("x", p)[0]

def merged_pdb(ag_groups, ab_groups, out):
    s = Structure("m"); m = BModel(0); s.add(m); A = BChain("A"); Bc = BChain("B")
    i = 1
    for grp in ag_groups:
        for r in grp:
            if "CA" not in r: continue
            rr = r.copy(); rr.id = (" ", i, " "); A.add(rr); i += 1
    j = 1
    for grp in ab_groups:
        for r in grp:
            if "CA" not in r: continue
            rr = r.copy(); rr.id = (" ", j, " "); Bc.add(rr); j += 1
    if i == 1 or j == 1: return False
    m.add(A); m.add(Bc); io = PDBIO(); io.set_structure(s); io.save(out); return True

def chain_res(model, cid):
    cid = str(cid)
    if cid not in model: return []
    return [r for r in model[cid] if "CA" in r]

def split_chains(s):
    return [x for x in re.split(r"[|,;\s]+", str(s)) if x and x.upper() not in ("NA", "NAN", "")]

def merge_test(cif, cj_path, out):
    cj = json.load(open(cj_path)); src = cj.get("src_chains", {})
    m = load(cif)
    crops = {}
    for c in cj["chains"]:
        if c["role"] == "antigen": crops[str(c.get("src"))] = c.get("crop")
    ag = []
    for cid in src.get("antigen", []):
        rr = chain_res(m, cid)
        if not rr: continue
        cr = crops.get(str(cid))
        if cr: rr = [r for r in rr if cr[0] <= r.id[1] <= cr[1]]
        ag.append(rr)
    ab = [chain_res(m, c) for c in list(src.get("H", [])) + list(src.get("L", []))]
    ab = [x for x in ab if x]
    if not any(ag) or not ab: return False
    return merged_pdb(ag, ab, out)

def merge_ref(cif, ag_chain, hchain, lchain, out):
    m = load(cif)
    ag = [chain_res(m, c) for c in split_chains(ag_chain)]; ag = [x for x in ag if x]
    ab = [chain_res(m, c) for c in split_chains(hchain) + split_chains(lchain)]; ab = [x for x in ab if x]
    if not any(ag) or not ab: return False
    return merged_pdb(ag, ab, out)

# ---------- PPIRef 계면 추출 ----------
def extract_iface(pdb, out_dir):
    base = os.path.splitext(os.path.basename(pdb))[0]
    pid = base.replace("_", "-")   # ⚠️ PPIExtractor는 pdb_id의 '_'를 '-'로 치환해 저장(chain 구분자 '_'와 충돌 방지)
    def find():
        for pat in (f"{pid}_A_B.pdb", f"{pid}_B_A.pdb", f"{base}_A_B.pdb", f"{base}_B_A.pdb"):
            h = glob.glob(os.path.join(out_dir, "**", pat), recursive=True)
            if h: return h[0]
        return None
    got = find()
    if got: return got                       # 이미 추출됨 → PPIExtractor 재호출 스킵(재실행 빠름)
    from ppiref.extraction import PPIExtractor
    ex = PPIExtractor(out_dir=out_dir, kind="heavy", radius=RADIUS, bsa=False)
    try: ex.extract(pdb, partners=["A", "B"])
    except Exception as e: print(f"    extract 실패 {os.path.basename(pdb)}: {e}"); return None
    return find()

# ---------- stage: 병합 + 계면 추출 ----------
def do_extract(a):
    os.makedirs(a.work, exist_ok=True)
    merged_dir = os.path.join(a.work, "merged"); iface_dir = os.path.join(a.work, "iface")
    os.makedirs(merged_dir, exist_ok=True); os.makedirs(iface_dir, exist_ok=True)
    rows = list(csv.DictReader(open(a.list)))
    man = {r["pdb"]: r for r in csv.DictReader(open(a.ref_manifest))} if os.path.exists(a.ref_manifest) else {}
    idx = []  # (kind, id, family, ab, iface_path)

    # 테스트
    for r in rows:
        tid = r["target"]; fam = r["group"]; ab = r["ab"]
        cif = os.path.join(a.targets_dir, tid, "native.cif"); cj = os.path.join(a.targets_dir, tid, "chains.json")
        if not (os.path.exists(cif) and os.path.exists(cj)): continue
        mp = os.path.join(merged_dir, f"T_{tid}.pdb")
        if not os.path.exists(mp) and not merge_test(cif, cj, mp): print(f"  merge 실패(test) {tid}"); continue
        ip = extract_iface(mp, iface_dir)
        if ip: idx.append(("test", tid, fam, ab, ip)); print(f"  test  {tid:12} {fam}/{ab}  iface={os.path.basename(ip)}")

    # 레퍼런스
    for pdb, r in man.items():
        cif = os.path.join(a.ref_struct, f"{pdb.lower()}.cif")   # fetch_structures가 소문자로 저장
        if not os.path.exists(cif): continue
        mp = os.path.join(merged_dir, f"R_{pdb}.pdb")
        if not os.path.exists(mp) and not merge_ref(cif, r.get("antigen_chain", ""), r.get("Hchain", ""), r.get("Lchain", ""), mp):
            continue
        ip = extract_iface(mp, iface_dir)
        if ip: idx.append(("ref", pdb, r.get("antigen", ""), "ref", ip))
    print(f"  레퍼런스 계면 {sum(1 for x in idx if x[0]=='ref')}개 추출")

    with open(os.path.join(a.work, "iface_index.csv"), "w", newline="") as f:
        w = csv.writer(f); w.writerow(["kind", "id", "family", "ab", "iface"]); w.writerows(idx)
    print(f"→ {a.work}/iface_index.csv ({len(idx)}개)")

# ---------- stage: 임베딩 + 이웃수 ----------
def embed_refs(fam, ref, cache, reembed):
    """항원군 레퍼런스 계면 임베딩 → (ref_mat[valid], ref_ids[valid]). NaN 행 제거 + npz 캐시(재실행 시 재임베딩 생략)."""
    from ppiref.comparison import IDist
    from pathlib import Path
    if os.path.exists(cache) and not reembed:
        z = np.load(cache, allow_pickle=False)
        print(f"  [{fam}] 캐시 로드: valid ref {len(z['mat'])} ({os.path.basename(cache)})")
        return z["mat"], list(z["ids"])
    idist = IDist(near_duplicate_threshold=NEAR_DUP)
    idist.embed_parallel([Path(r["iface"]) for r in ref])
    E = idist.get_embeddings()                          # index=계면 id, 값=임베딩 벡터(dim=20)
    ids = list(E.index); mat = np.asarray(E.values, dtype=float)
    valid = np.isfinite(mat).all(axis=1)                # ⚠️ embed_parallel은 실패 계면을 NaN 행으로 반환 → 제거
    mat = mat[valid]; ids = [i for i, v in zip(ids, valid) if v]
    dropped = int((~valid).sum())
    np.savez(cache, mat=mat, ids=np.array(ids))         # 문자열 배열(pickle 불필요)
    print(f"  [{fam}] 임베딩 {len(valid)}개 중 NaN {dropped}개 제거 → valid {len(mat)}, 캐시 저장")
    return mat, ids

def do_score(a):
    from ppiref.comparison import IDist
    from pathlib import Path                       # IDist.embed는 Path 객체를 요구(ppi.stem)
    idx = list(csv.DictReader(open(os.path.join(a.work, "iface_index.csv"))))
    fams = sorted({r["family"] for r in idx if r["kind"] == "test"})
    idist_q = IDist(near_duplicate_threshold=NEAR_DUP)   # 테스트 계면 임베딩 전용
    out = []
    for fam in fams:
        ref = [r for r in idx if r["kind"] == "ref" and r["family"] == fam]
        test = [r for r in idx if r["kind"] == "test" and r["family"] == fam]
        if not ref:
            print(f"  [{fam}] 레퍼런스 0 → 스킵(과대표집 미측정)"); continue
        ref_mat, ref_ids = embed_refs(fam, ref, os.path.join(a.work, f"refemb_{fam}.npz"), a.reembed)
        if len(ref_mat) == 0:
            print(f"  [{fam}] valid ref 0 → 스킵"); continue
        print(f"  [{fam}] 테스트 {len(test)}개 채점")
        for t in test:
            q = np.asarray(idist_q.embed(Path(t["iface"]), store=False), dtype=float).ravel()
            if q.shape[0] != ref_mat.shape[1] or not np.isfinite(q).all():
                print(f"    스킵 {t['id']} (dim 불일치/NaN)"); continue
            d = np.linalg.norm(ref_mat - q, axis=1)
            k = min(KNN, len(d)); knn = float(np.sort(d)[:k].mean())
            counts = {f"n_{rad}": int((d <= rad).sum()) for rad in RADII}
            j = int(np.argmin(d))
            out.append(dict(target=t["id"], family=fam, ab=t["ab"], n_ref=len(ref_mat),
                            min_dist=round(float(d[j]), 4), nearest=ref_ids[j], mean_knn=round(knn, 4),
                            frac_ndup=round(counts["n_0.04"] / len(ref_mat), 4), **counts))
            print(f"    {t['id']:12} {t['ab']}  min_d={d[j]:.3f}(→{ref_ids[j]})  knn5={knn:.3f}  "
                  f"ndup0.04={counts['n_0.04']}({out[-1]['frac_ndup']:.0%})")
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    cols = ["target", "family", "ab", "n_ref", "min_dist", "nearest", "mean_knn", "frac_ndup"] + [f"n_{r}" for r in RADII]
    with open(a.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols); w.writeheader(); w.writerows(out)
    print(f"\n→ {a.out} ({len(out)}행). 주 신호 = min_dist(작을수록 과대표집)·mean_knn·frac_ndup. "
          f"다음: idist_analyze.py 로 A/B 대비 + 깊이-취약성 상관.")

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--stage", choices=["extract", "score", "all"], default="all")
    p.add_argument("--list", default="sweep_targets.csv")
    p.add_argument("--targets-dir", default="targets")
    p.add_argument("--ref-manifest", dest="ref_manifest", default="ref_manifest.csv")
    p.add_argument("--ref-struct", dest="ref_struct", default="ref_structures")
    p.add_argument("--work", default=os.path.join(os.environ.get("DATA", "/mnt/data/admuser/msadepth"), "idist"))
    p.add_argument("--out", default="results/overrep_idist.csv")
    p.add_argument("--reembed", action="store_true", help="refemb 캐시 무시하고 레퍼런스 재임베딩")
    a = p.parse_args()
    if a.stage in ("extract", "all"): do_extract(a)
    if a.stage in ("score", "all"): do_score(a)

if __name__ == "__main__":
    main()
