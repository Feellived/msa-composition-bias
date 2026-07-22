#!/usr/bin/env python3
"""iDist 레퍼런스 매니페스트 — 컷 이전(date < 2023-06-01) SAbDab 항체-항원 = 모델이 학습한 'prior'.
build_manifest.py의 항원군 분류(RBD/HA/Env)를 재사용하되 날짜 필터를 반전. 각 테스트 계면의
'컷 이전 near-duplicate 이웃수'를 세기 위한 참조 풀. 출력 = ref_manifest.csv → fetch_structures.py로 구조 다운.

⚠️ leakage 기준 = SAbDab 'date'(release/공개일). 최종 lock은 RCSB release date 재확인 권장(현재 1차 필터).
⚠️ 과대표집은 raw 이웃수라 레퍼런스 크기에 스케일 의존 → A/B '대비(비율)'로 해석. --max-per-family 로 상한(초과분 로그).

사용: python idist_ref_manifest.py [--cut 2023-06-01] [--max-per-family 1200]
"""
import argparse, re
import pandas as pd

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--summary", default="sabdab2_summary.csv")
    ap.add_argument("--cut", default="2023-06-01")
    ap.add_argument("--max-res", type=float, default=3.5)
    ap.add_argument("--max-per-family", type=int, default=1200)
    ap.add_argument("--out", default="ref_manifest.csv")
    a = ap.parse_args()

    df = pd.read_csv(a.summary, low_memory=False)
    df["res"] = pd.to_numeric(df["resolution"], errors="coerce")
    df["dt"] = pd.to_datetime(df["date"], errors="coerce")
    at = df["antigen_type"].fillna("").str.upper()
    df["is_prot"] = at.str.split("|").map(lambda xs: any(s.strip() == "PROTEIN" for s in xs))
    s = lambda x: str(x).strip()
    df["H"] = df["Hchain"].map(s); df["L"] = df["Lchain"].map(s); df["AG"] = df["antigen_chain"].map(s)
    has_h = ~df["H"].str.upper().isin(["", "NA", "NAN"])
    has_ag = ~df["AG"].str.upper().isin(["", "NA", "NAN"])
    df["nb"] = df["L"].str.upper().isin(["", "NA", "NAN"])
    agn = df["antigen_name"].fillna("").str.lower()
    org = (df["organism"].fillna("") + " " + df["antigen_species"].fillna("")).str.lower()
    FID = agn.str.contains(r"guanine nucleotide|g protein|gpcr|beta.?arrestin|megabody|fiducial", regex=True)

    def grp(ag, o):
        if re.search(r"hemagglutinin|haemagglutinin", ag): return "HA"
        if re.search(r"gp120|gp160|envelope glyco", ag) and re.search(r"immunodef|hiv|hxb", o): return "Env"
        if re.search(r"gp120|gp160", ag): return "Env"
        if re.search(r"spike|receptor.?binding|rbd|\bs1\b|s2 subunit", ag) and re.search(r"coronavirus 2|sars.?cov.?2|2697049|severe acute resp", o): return "RBD"
        return None
    df["grp"] = [grp(x, y) for x, y in zip(agn, org)]

    CUT = pd.Timestamp(a.cut)
    sel = df[df["is_prot"] & has_h & has_ag & ~FID & df["grp"].notna()
             & (df["dt"] < CUT) & (df["res"] <= a.max_res)].copy()          # ← 날짜 반전 = 컷 이전
    sel = sel.sort_values("dt").drop_duplicates(subset=["PDB", "H", "L", "AG"])

    def pdb4(x):
        m = re.match(r"pdb_0000([a-z0-9]{4})$", str(x)); return m.group(1).upper() if m else str(x)
    sel["pdb"] = sel["PDB"].map(pdb4)
    out = sel[["pdb", "grp", "H", "L", "AG", "nb", "date", "res", "antigen_name"]].rename(
        columns={"grp": "antigen", "H": "Hchain", "L": "Lchain", "AG": "antigen_chain"})

    # 항원군별 상한(초과분 로그 = 조용한 절단 금지)
    kept = []
    for g in ["RBD", "HA", "Env"]:
        sub = out[out.antigen == g].sort_values("date")
        if len(sub) > a.max_per_family:
            print(f"  ⚠️ {g}: {len(sub)}개 중 {a.max_per_family}개만 사용(오래된순), {len(sub)-a.max_per_family}개 드롭")
            sub = sub.iloc[:a.max_per_family]
        kept.append(sub)
    out = pd.concat(kept, ignore_index=True)
    out.to_csv(a.out, index=False)

    print(f"\n=== iDist 레퍼런스 (date < {a.cut}, res≤{a.max_res}) ===")
    for g in ["RBD", "HA", "Env"]:
        sub = out[out.antigen == g]
        print(f"  {g:4s}: {len(sub):4d} 복합체 (paired={ (~sub.nb).sum() }, nanobody={ sub.nb.sum() })")
    print(f"  전체 {len(out)} → {a.out}")
    print(f"  다음: python fetch_structures.py --manifest {a.out} --outdir ref_structures")

if __name__ == "__main__":
    main()