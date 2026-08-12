#!/usr/bin/env python3
"""SAbDab2 summary → RBD/HA/Env 후보 매니페스트(post-2023-06 clean).
정밀 필터 = antigen_name 키워드 + organism/antigen_species. PDB 확장ID→4char 매핑.
paired(H+L) vs nanobody(H only) 세분. Fab-fiducial(G단백질/GPCR) 제외.
"""
import pandas as pd, re
df = pd.read_csv("sabdab2_summary.csv", low_memory=False)
df["res"] = pd.to_numeric(df["resolution"], errors="coerce")
df["dt"] = pd.to_datetime(df["date"], errors="coerce")
at = df["antigen_type"].fillna("").str.upper()
df["is_prot"] = at.str.split("|").map(lambda xs: any(s.strip()=="PROTEIN" for s in xs))
def s(x): return str(x).strip()
df["H"]=df["Hchain"].map(s); df["L"]=df["Lchain"].map(s); df["AG"]=df["antigen_chain"].map(s)
has_h = ~df["H"].str.upper().isin(["","NA","NAN"])
has_ag = ~df["AG"].str.upper().isin(["","NA","NAN"])
df["nb"] = df["L"].str.upper().isin(["","NA","NAN"])   # nanobody = no light chain
agn = df["antigen_name"].fillna("").str.lower()
org = (df["organism"].fillna("")+" "+df["antigen_species"].fillna("")).str.lower()

FID = agn.str.contains(r"guanine nucleotide|g protein|gpcr|beta.?arrestin|megabody|fiducial", regex=True)

def grp(row_agn, row_org):
    if re.search(r"hemagglutinin|haemagglutinin", row_agn): return "HA"
    if re.search(r"gp120|gp160|envelope glyco", row_agn) and re.search(r"immunodef|hiv|hxb", row_org): return "Env"
    if re.search(r"gp120|gp160", row_agn): return "Env"  # gp120/160 = HIV 사실상 고유
    if re.search(r"spike|receptor.?binding|rbd|\bs1\b|s2 subunit", row_agn) and re.search(r"coronavirus 2|sars.?cov.?2|2697049|severe acute resp", row_org): return "RBD"
    return None
df["grp"] = [grp(a,o) for a,o in zip(agn,org)]

CUT = pd.Timestamp("2023-06-01")
sel = df[df["is_prot"] & has_h & has_ag & ~FID & df["grp"].notna()
         & (df["dt"]>=CUT) & (df["res"]<=3.5)].copy()

def pdb4(x):
    m = re.match(r"pdb_0000([a-z0-9]{4})$", str(x))
    return m.group(1).upper() if m else str(x)
sel["pdb4"] = sel["PDB"].map(pdb4)

# distinct 구조 단위로 dedup(한 PDB에 여러 Fv면 첫 행). 한 PDB 여러 항체는 뒤 step서 사슬별 처리.
mani = sel.sort_values("dt").drop_duplicates(subset=["PDB","H","L","AG"])
out = mani[["pdb4","PDB","grp","H","L","AG","nb","date","res","antigen_name","organism"]].rename(
    columns={"pdb4":"pdb","PDB":"pdb_full","grp":"antigen","H":"Hchain","L":"Lchain","AG":"antigen_chain"})
out.to_csv("manifest_candidates.csv", index=False)

print("=== 후보 매니페스트 (post-2023-06, res≤3.5) ===")
for g in ["RBD","HA","Env"]:
    sub = out[out.antigen==g]
    npdb = sub["pdb"].nunique(); npair = sub[~sub.nb]["pdb"].nunique(); nnb = sub[sub.nb]["pdb"].nunique()
    print(f"  {g:4s}: rows(Fv)={len(sub):4d}  distinct PDB={npdb:4d}  (paired={npair}, nanobody={nnb})")
print(f"  전체 rows={len(out)}  distinct PDB={out['pdb'].nunique()}")
print(f"\n  organism 샘플(정밀필터 검증):")
for g in ["RBD","HA","Env"]:
    ex = out[out.antigen==g]["organism"].value_counts().head(3)
    print(f"   [{g}] "+" | ".join(f"{k[:35]}×{v}" for k,v in ex.items()))
print(f"\n→ manifest_candidates.csv ({len(out)} rows)")
