#!/usr/bin/env python3
"""통합 co-folder 입력생성기(다사슬 항원 대응). MSA-depth 실험용:
 항원 사슬 = subsampled a3m(깊이별), 항체 = single-seq(query only).
chains.json(prep_targets 포맷: antigen=[ids], antibody=[ids], chains=[{id,role,seq}]).
--ag-a3m "A=/abs/A_d64.a3m,B=/abs/B_d64.a3m" (항원 사슬별 그 깊이의 a3m).
사용:
 python make_input.py --cofolder boltz    --chains T/chains.json --ag-a3m "A=..,B=.." --out T/boltz_d64.yaml
 python make_input.py --cofolder protenix --chains T/chains.json --ag-a3m "A=.."      --dir T/prot_d64  --out T/prot_d64.json
 python make_input.py --cofolder chai     --chains T/chains.json --ag-a3m "A=..,B=.." --dir T/chai_d64  --out T/chai_d64.fasta
   chai: --out=FASTA, --dir/msa 에 항원 사슬별 aligned.pqt(파일명=서열 sha256) 생성. 항체는 pqt 없음 → chai가 single-seq.
"""
import argparse, json, os

def parse_map(s):
    m = {}
    for kv in s.split(","):
        if "=" in kv:
            k, v = kv.split("=", 1); m[k.strip()] = v.strip()
    return m

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cofolder", required=True, choices=["boltz", "protenix", "chai"])
    ap.add_argument("--chains", required=True)
    ap.add_argument("--ag-a3m", required=True, help='"A=/abs/A_dD.a3m,B=/abs/B_dD.a3m"')
    ap.add_argument("--out", required=True)
    ap.add_argument("--dir", help="protenix: query a3m 저장 위치")
    a = ap.parse_args()
    d = json.load(open(a.chains))
    ag = d["antigen"]; ag = [ag] if isinstance(ag, str) else [str(x) for x in ag]
    ab = d["antibody"]; ab = [ab] if isinstance(ab, str) else [str(x) for x in ab]
    sm = {c["id"]: c["seq"] for c in d["chains"]}
    a3m = parse_map(a.ag_a3m)
    for c in ag:
        if c not in a3m: raise SystemExit(f"항원 사슬 {c} a3m 경로 없음: --ag-a3m")

    if a.cofolder == "boltz":
        # 각 항원 사슬 = 자기 depth a3m, 항체 = empty(single-seq)
        L = ["version: 1", "sequences:"]
        for c in d["chains"]:
            msa = os.path.abspath(a3m[c["id"]]) if c["id"] in ag else "empty"
            L += ["  - protein:", f"      id: {c['id']}", f'      sequence: "{c["seq"]}"', f"      msa: {msa}"]
        open(a.out, "w").write("\n".join(L) + "\n")

    elif a.cofolder == "protenix":
        outdir = a.dir or os.path.dirname(a.out); os.makedirs(outdir, exist_ok=True)
        seqs = []
        for c in ag:
            seqs.append({"proteinChain": {"sequence": sm[c], "count": 1,
                                          "unpairedMsaPath": os.path.abspath(a3m[c])}})
        for i in ab:
            qa = os.path.join(outdir, f"ab_{i}_query.a3m"); open(qa, "w").write(f">{i}\n{sm[i]}\n")
            seqs.append({"proteinChain": {"sequence": sm[i], "count": 1, "unpairedMsaPath": os.path.abspath(qa)}})
        name = d.get("target", "job").lower()
        json.dump([{"name": name, "sequences": seqs}], open(a.out, "w"), indent=2)

    else:  # chai — FASTA(>protein|name=ID) + msa_directory에 항원 aligned.pqt(파일명=서열 sha256)
        from pathlib import Path
        # chai_lab 내부 변환기(설치된 버전 함수 사용) — chai env에서만 lazy import
        from chai_lab.data.parsing.msas.aligned_pqt import a3m_to_aligned_dataframe, expected_basename
        from chai_lab.data.parsing.msas.data_source import MSADataSource
        fa = []
        for c in d["chains"]:
            fa += [f">protein|name={c['id']}", c["seq"]]
        open(a.out, "w").write("\n".join(fa) + "\n")
        msadir = Path(a.dir or os.path.dirname(a.out)) / "msa"; msadir.mkdir(parents=True, exist_ok=True)
        for c in ag:                                        # 항원 사슬만 pqt(=그 rung의 depth). 항체는 pqt 없음→single-seq
            df = a3m_to_aligned_dataframe(os.path.abspath(a3m[c]),
                                          source_database=MSADataSource.UNIREF90, insert_pairing_key=False)
            df.to_parquet(msadir / expected_basename(sm[c]))    # expected_basename: sha256(uppercased seq)
    print(f"[{a.cofolder}] {d.get('target')} 항원={ag}(depth a3m) 항체={ab}(single-seq) -> {a.out}")

if __name__ == "__main__":
    main()
