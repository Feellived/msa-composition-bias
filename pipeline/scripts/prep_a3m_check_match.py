#!/usr/bin/env python3
"""[무결성 점검] a3m의 질의서열이 실제 입력 서열과 같은가 — MSA가 진짜로 쓰였는지 확인.

왜 필요한가: boltz가 "MSA does not match input sequence, creating dummy" 경고와 함께
**항원 MSA를 버리고 더미로 대체**한 사례가 발견됨(8txu·9y0a·8y6a). 그러면 깊이·조성을
아무리 바꿔도 입력이 사실상 동일해져 실험 자체가 성립하지 않는다.
Protenix는 같은 상황에서 경고를 안 낼 수 있으므로 **로그가 아니라 서열을 직접 대조**한다.

각 타깃의 항원 사슬마다:
  chains.json 의 서열  vs  ladders/<target>/<chain>/rung0.a3m 의 첫 서열(질의행)
을 비교해 길이·동일성·불일치 위치 수를 보고한다.

판정:
  OK        = 완전히 같음 → MSA 정상 사용
  MISMATCH  = 다름 → 그 타깃의 MSA 실험은 무효(모델이 더미로 대체했을 수 있음)
  LEN       = 길이부터 다름 → 구조/구성물이 다른 서열로 a3m을 만든 것

사용(stdlib only):
  python prep_a3m_check_match.py                      # 전 타깃
  python prep_a3m_check_match.py --only 8ulr_HL 9azr_HL
  python prep_a3m_check_match.py --list sweep_targets.csv --data /mnt/data/msadepth
"""
import argparse, csv, json, os, re, sys

AA = set("ACDEFGHIKLMNPQRSTVWYXBZUO")


def read_a3m_query(path):
    """a3m 첫 레코드의 서열(정렬 문자 제거).

    ⚠️ ColabFold a3m은 첫 줄이 메타 주석(`#<길이>\t<개수>`)이므로 건너뛴다. 또한 옛 사다리
    파일은 그 주석이 질의 서열 앞에 붙어 있으므로(2026-07-27 버그), 서열 맨 앞의 군더더기도
    떼고 비교한다 — 그래야 "진짜 서열이 다른가"와 "머리말만 붙었나"를 구분할 수 있다.
    """
    seq, started = [], False
    with open(path) as f:
        for ln in f:
            ln = ln.rstrip("\n")
            if not ln:
                continue
            if ln.startswith("#"):
                rest = re.sub(r"^#\d+\t\d+", "", ln)
                if not rest:                   # 단독 주석 줄(원본 a3m의 정상 형태) → 무시
                    continue
                ln = ln                        # 서열이 뒤에 붙은 손상 형태 → 그대로 두고 아래서 감지
            if ln.startswith(">"):
                if started:
                    break
                started = True
                continue
            if started:
                seq.append(ln)
    raw = "".join(seq)
    clean = re.sub(r"^#\d+\s*\d*", "", raw)    # 서열 앞에 붙어버린 주석을 뗀 것
    norm = lambda s: "".join(c for c in s.upper() if c not in "-.")
    return norm(raw), norm(clean)


def find_chain_seqs(cj):
    """chains.json에서 {사슬ID: 서열} 과 항원 사슬 목록을 최대한 견고하게 추출."""
    seqs, ag = {}, []

    def looks_seq(v):
        return isinstance(v, str) and len(v) >= 20 and set(v.upper()) <= AA

    def walk(node, key_hint=""):
        if isinstance(node, dict):
            for k, v in node.items():
                if looks_seq(v) and len(str(k)) <= 4:
                    seqs[str(k)] = v.upper()
                elif isinstance(v, dict) and looks_seq(v.get("sequence", "")) and len(str(k)) <= 4:
                    seqs[str(k)] = v["sequence"].upper()
                else:
                    walk(v, str(k).lower())
        elif isinstance(node, list):
            # [{"id":"A","sequence":"..."} , ...] 형태
            for it in node:
                if isinstance(it, dict):
                    cid = it.get("id") or it.get("chain") or it.get("chain_id")
                    sq = it.get("sequence") or it.get("seq")
                    if cid and looks_seq(sq or ""):
                        seqs[str(cid)] = sq.upper()
                        if any(t in key_hint for t in ("antigen", "ag")):
                            ag.append(str(cid))
                    else:
                        walk(it, key_hint)
                elif isinstance(it, str) and len(it) <= 4 and any(
                        t in key_hint for t in ("antigen", "ag")):
                    ag.append(it)
                else:
                    walk(it, key_hint)

    walk(cj)
    # 항원 목록을 못 찾았으면 키 이름으로 재시도
    if not ag:
        for k, v in (cj.items() if isinstance(cj, dict) else []):
            if any(t in str(k).lower() for t in ("antigen", "ag_chain", "ag")):
                if isinstance(v, list):
                    ag += [str(x) for x in v if isinstance(x, (str, int)) and len(str(x)) <= 4]
                elif isinstance(v, dict):
                    ag += [str(x) for x in v.keys() if len(str(x)) <= 4]
    return seqs, sorted(set(ag))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--list", default="sweep_targets.csv")
    ap.add_argument("--targets-dir", default="targets")
    ap.add_argument("--data", default=os.environ.get("DATA", "/mnt/data/msadepth"))
    ap.add_argument("--rung", default="rung0", help="비교할 사다리 칸(기본 rung0=full)")
    ap.add_argument("--only", nargs="*", default=[])
    ap.add_argument("--out", default="results/msa_match_check.csv")
    a = ap.parse_args()

    rows = list(csv.DictReader(open(a.list)))
    laddir = os.path.join(a.data, "ladders")
    out, bad = [], []
    print(f"a3m 질의서열 vs chains.json 서열 대조 (사다리 {a.rung})\n")
    print(f"{'target':12}{'사슬':>5}{'chains.json':>12}{'a3m':>8}{'같은가':>9}  비고")
    print("-" * 78)

    for r in rows:
        t = r["target"]
        if a.only and t not in a.only:
            continue
        cjp = os.path.join(a.targets_dir, t, "chains.json")
        if not os.path.exists(cjp):
            print(f"{t:12}  chains.json 없음"); continue
        try:
            cj = json.load(open(cjp))
        except Exception as e:
            print(f"{t:12}  chains.json 읽기 실패: {e}"); continue
        seqs, ag = find_chain_seqs(cj)
        # sweep_targets.csv 6번째 열(ag_chains, 'A|B')을 우선 사용
        col = list(r.values())
        agc = (r.get("ag_chains") or (col[5] if len(col) > 5 else "") or "").strip()
        chains = [c for c in agc.split("|") if c] or ag
        if not chains:
            print(f"{t:12}  항원 사슬 판별 실패 (chains.json 키: {list(cj)[:6] if isinstance(cj,dict) else type(cj)})")
            continue
        for c in chains:
            ap3 = os.path.join(laddir, t, c, f"{a.rung}.a3m")
            if not os.path.exists(ap3):
                print(f"{t:12}{c:>5}{'-':>12}{'없음':>8}{'?':>9}  a3m 없음")
                continue
            q, qc = read_a3m_query(ap3)
            ref = seqs.get(c, "")
            if not ref:
                print(f"{t:12}{c:>5}{'?':>12}{len(q):>8}{'?':>9}  chains.json에 이 사슬 서열 없음")
                continue
            if q == ref:
                verdict, note = "OK", ""
            elif qc == ref:
                verdict = "머리말오염"
                note = f"앞 {len(q)-len(qc)}글자 군더더기 — 서열은 정확. prep_a3m_fix_query.py로 복구 가능"
                bad.append((t, c, note))
            else:
                verdict = "MISMATCH"
                if len(qc) != len(ref):
                    note = f"⚠️ 서열 자체가 다름 (길이 차 {abs(len(qc)-len(ref))})"
                else:
                    d = sum(1 for x, y in zip(qc, ref) if x != y)
                    note = f"⚠️ 서열 자체가 다름 ({d}곳, {100*(1-d/len(qc)):.1f}% 동일)"
                bad.append((t, c, note))
            print(f"{t:12}{c:>5}{len(ref):>12}{len(q):>8}{verdict:>11}  {note}")
            out.append(dict(target=t, chain=c, len_input=len(ref), len_a3m_raw=len(q),
                            len_a3m_clean=len(qc), verdict=verdict, note=note))

    print("\n" + "=" * 78)
    head = [r for r in out if r["verdict"] == "머리말오염"]
    real = [r for r in out if r["verdict"] == "MISMATCH"]
    okn = [r for r in out if r["verdict"] == "OK"]
    print(f"정상 {len(okn)} · 머리말오염 {len(head)} · 서열자체다름 {len(real)}")
    if head:
        ts = sorted({r["target"] for r in head})
        print(f"\n[머리말 오염] 타깃 {len(ts)}개 — 서열은 맞고 앞에 군더더기만 붙음(복구 가능):")
        print("   " + ", ".join(ts))
        print("   → python prep_a3m_fix_query.py  (dry-run 확인 후 --apply)")
        print("   ⚠️ 이 입력으로 만든 예측은 다시 돌려야 함(boltz는 MSA를 버렸고 Protenix는 밀린 정렬을 씀).")
    if real:
        ts = sorted({r["target"] for r in real})
        print(f"\n[서열 자체가 다름] 타깃 {len(ts)}개 — 복구 불가, a3m 재생성 필요:")
        for r in real:
            print(f"   {r['target']} / {r['chain']} — {r['note']}")
    if not head and not real:
        print("✅ 전부 정상 — MSA 질의행이 입력 서열과 일치.")

    if out:
        os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
        with open(a.out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(out[0].keys()))
            w.writeheader(); w.writerows(out)
        print(f"\n→ {a.out}")


if __name__ == "__main__":
    main()
