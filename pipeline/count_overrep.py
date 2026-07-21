#!/usr/bin/env python3
"""SAbDab2 summary → 과대표집 항원별 표본 풀 카운트.
- 행=Fv 페어링이므로 '구조 수'는 distinct PDB로 센다.
- 과대표집 = 항원(정규화 name/키워드 버킷)당 distinct PDB 수.
- leakage 프록시 = SABDABdepo_date(deposition, YYYYMMDD). 최종 lock은 RCSB release date로 재확인 필요.
"""
import pandas as pd, re, collections

df = pd.read_csv("sabdab2_summary.csv", low_memory=False)
print(f"[load] rows={len(df)}  cols={len(df.columns)}")

# --- 정제 ---
df["res"] = pd.to_numeric(df["resolution"], errors="coerce")
def depo(x):
    try: return int(float(x))
    except: return None
df["depo"] = df["SABDABdepo_date"].map(depo)          # YYYYMMDD int
at = df["antigen_type"].fillna("").str.upper()
is_prot = at.str.split("|").map(lambda xs: any(s.strip()=="PROTEIN" for s in xs))
has_ab  = df["Hchain"].notna() & (df["Hchain"].astype(str).str.strip()!="") & (df["Hchain"].astype(str).str.strip().str.upper()!="NA")
has_ag  = df["antigen_chain"].notna() & (df["antigen_chain"].astype(str).str.strip()!="") & (df["antigen_chain"].astype(str).str.strip().str.upper()!="NA")
base = df[is_prot & has_ab & has_ag].copy()
print(f"[filter] 단백질항원 + Ab + Ag사슬: rows={len(base)}  distinct PDB={base['PDB'].nunique()}")

def nname(s):
    s = str(s).lower().strip()
    s = re.sub(r"[^a-z0-9 ]"," ",s); s = re.sub(r"\s+"," ",s).strip()
    return s
base["agn"] = base["antigen_name"].map(nname)

# --- (1) 정규화 name 기준 top 과대표집 (distinct PDB) ---
g = base.groupby("agn")["PDB"].nunique().sort_values(ascending=False)
print("\n===== [A] 정규화 antigen_name 상위 25 (distinct PDB) =====")
for name, n in g.head(25).items():
    print(f"  {n:5d}  {name[:70]}")

# --- (2) 키워드 버킷(과대표집 후보 항원) ---
BUCKETS = {
 "SARS2 spike/RBD": r"(sars.?cov.?2|spike|receptor.?binding|rbd)",
 "influenza HA":    r"(hemagglutinin|haemagglutinin|influenza)",
 "HIV Env/gp120":   r"(gp120|gp160|envelope glycoprotein|hiv)",
 "lysozyme(HEL)":   r"(lysozyme)",
 "PD-1/PD-L1":      r"(pd.?1|pd.?l1|programmed cell death)",
 "TNF":             r"(tumou?r necrosis|tnf)",
 "VEGF":            r"(vegf|vascular endothelial)",
 "EGFR/HER2":       r"(egfr|epidermal growth factor receptor|her2|erbb2)",
 "RSV F":           r"(respiratory syncytial|rsv|fusion glycoprotein f)",
 "Ebola GP":        r"(ebola|glycoprotein gp)",
}
print("\n===== [B] 과대표집 후보 버킷 (distinct PDB; total / depo≥2023-06 / depo≥2023-01) =====")
for lbl, pat in BUCKETS.items():
    m = base[base["agn"].str.contains(pat, regex=True, na=False)]
    tot = m["PDB"].nunique()
    p2306 = m[m["depo"]>=20230601]["PDB"].nunique()
    p2301 = m[m["depo"]>=20230101]["PDB"].nunique()
    print(f"  {lbl:18s}  total={tot:5d}   post-2023-06={p2306:4d}   post-2023-01={p2301:4d}")

# --- (3) 전체 분포 요약: 과대표집 정의 감 잡기 ---
print("\n===== [C] 항원(정규화 name)당 distinct-PDB 분포 =====")
vc = g.values
import numpy as np
for q in [0.5,0.75,0.9,0.95,0.99]:
    print(f"  {int(q*100)}th pct = {np.quantile(vc,q):.0f} complexes/antigen")
print(f"  antigens with ≥10 PDB: {(g>=10).sum()}   ≥20: {(g>=20).sum()}   ≥50: {(g>=50).sum()}")
print(f"  총 distinct 항원(name) 수: {len(g)}")

# --- (4) post-2023-06 클린 풀 전체 규모 ---
clean = base[base["depo"]>=20230601]
print(f"\n===== [D] post-2023-06(deposition) 클린 풀 =====")
print(f"  rows={len(clean)}  distinct PDB={clean['PDB'].nunique()}")
print(f"  해상도≤3.5Å 만: distinct PDB={clean[clean['res']<=3.5]['PDB'].nunique()}")
# 나노바디(Lchain 없음) vs 통상항체
lc = clean["Lchain"].astype(str).str.strip().str.upper()
nb = clean[(lc=="")|(lc=="NA")]["PDB"].nunique()
print(f"  그 중 나노바디(VHH, Lchain 없음) distinct PDB≈{nb}")
