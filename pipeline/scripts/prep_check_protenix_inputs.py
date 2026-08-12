#!/usr/bin/env python3
"""protenix/chai 입력 사전검증 (CPU 병렬, GPU 불필요).
모든 target×ag_chain×rung의 사다리 a3m을 clean_a3m으로 정제 후 '쿼리 서열 == chains.json 서열'인지 확인
→ protenix 'MSA query/size mismatch'를 GPU 쓰기 전에 일괄 검출. (2026-07-22 protenix 전멸 재발 방지 프리플라이트)

사용(pipeline/에서): python prep_check_protenix_inputs.py [--rungs 12] [--workers 0=자동]
  DATA 환경변수로 사다리 위치 지정(기본 /mnt/data/msadepth). GPU·모델 실행 없음 = 수 초.
"""
import argparse, csv, json, os
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from make_input import clean_a3m_lines

DATA = os.environ.get("DATA", "/mnt/data/msadepth")
LADDIR = os.path.join(DATA, "ladders")

def query_of(a3m_path):
    """정제된 a3m에서 첫 '>'+다음 줄 = 쿼리 서열."""
    lines = clean_a3m_lines(a3m_path)
    for i, l in enumerate(lines):
        if l.startswith(">") and i + 1 < len(lines):
            return lines[i + 1]
    return None

def check_one(t):
    target, chain, seq, rung = t
    a3m = os.path.join(LADDIR, target, chain, f"rung{rung}.a3m")
    if not os.path.exists(a3m):
        return (target, chain, rung, "MISSING", len(seq), None)
    try:
        q = query_of(a3m)
    except Exception as e:
        return (target, chain, rung, f"ERR:{type(e).__name__}", len(seq), None)
    if q is None:
        return (target, chain, rung, "NO_QUERY", len(seq), None)
    if q == seq:
        return (target, chain, rung, "OK", len(seq), len(q))
    return (target, chain, rung, "MISMATCH", len(seq), len(q))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--list", default="sweep_targets.csv")
    ap.add_argument("--targets-dir", default="targets")
    ap.add_argument("--rungs", type=int, default=12)
    ap.add_argument("--workers", type=int, default=0, help="0=자동(코어수)")
    a = ap.parse_args()

    tasks = []
    for r in csv.DictReader(open(a.list)):
        tgt = r["target"]
        cj = os.path.join(a.targets_dir, tgt, "chains.json")
        if not os.path.exists(cj):
            print(f"  [skip] {tgt} chains.json 없음"); continue
        d = json.load(open(cj))
        ag = d["antigen"]; ag = [ag] if isinstance(ag, str) else [str(x) for x in ag]
        sm = {c["id"]: c["seq"] for c in d["chains"]}
        for c in ag:
            c = str(c)
            if c not in sm: continue
            for rung in range(a.rungs):
                tasks.append((tgt, c, sm[c], rung))
    print(f"검증 대상 = {len(tasks)}개 (target×ag_chain×rung)")

    workers = a.workers or (os.cpu_count() or 4)
    results = []
    with ProcessPoolExecutor(max_workers=workers) as ex:
        for res in ex.map(check_one, tasks, chunksize=8):
            results.append(res)

    cnt = Counter(r[3].split(":")[0] for r in results)
    ok = cnt.get("OK", 0)
    print(f"\n=== 결과: OK {ok}/{len(results)} (workers={workers}) ===")
    for k, v in sorted(cnt.items()):
        if k != "OK": print(f"  {k}: {v}")
    bad = [r for r in results if r[3] != "OK"]
    if bad:
        print(f"\n실패 상세(최대 40):")
        for target, chain, rung, status, ls, lq in bad[:40]:
            print(f"  {target} {chain} rung{rung}: {status} (seq_len={ls}, a3m_query_len={lq})")
        print(f"\n⚠️ 실패 {len(bad)}건 — protenix 돌리기 전 해결 필요.")
    else:
        print("\n✅ 전부 OK — protenix/chai 재실행 안전 (정제 a3m 쿼리 == 입력 서열).")

if __name__ == "__main__":
    main()
