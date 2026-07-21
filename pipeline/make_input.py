#!/usr/bin/env python3
"""통합 co-folder 입력생성기(다사슬 항원 대응). MSA-depth 실험용:
 항원 사슬 = subsampled a3m(깊이별), 항체 = single-seq(query only).
chains.json(prep_targets 포맷: antigen=[ids], antibody=[ids], chains=[{id,role,seq}]).
--ag-a3m "A=/abs/A_d64.a3m,B=/abs/B_d64.a3m" (항원 사슬별 그 깊이의 a3m).
사용:
 python make_input.py --cofolder boltz    --chains T/chains.json --ag-a3m "A=..,B=.." --out T/boltz_d64.yaml
 python make_input.py --cofolder protenix --chains T/chains.json --ag-a3m "A=.."      --dir T/prot_d64  --out T/prot_d64.json
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
    ap.add_argument("--cofolder", required=True, choices=["boltz", "protenix"])
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

    else:  # protenix
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
    print(f"[{a.cofolder}] {d.get('target')} 항원={ag}(depth a3m) 항체={ab}(single-seq) -> {a.out}")

if __name__ == "__main__":
    main()
