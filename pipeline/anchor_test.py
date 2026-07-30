#!/usr/bin/env python3
"""[앵커 검정] 원래 MSA를 한 집단으로 포함한 조성 이질성 검정 — 실행 단위·순열.

8ulr_HL에 쓴 검정을 모든 타깃에 적용한다. 본 검정 설계가 이미
`조성 N × 반복 R + 원래 MSA(seedfull) × R'`이므로 GPU를 더 쓰지 않는다.

왜 이 검정인가(5.1에서 확정):
  · 실행 단위  — 한 실행의 자세 5개는 독립 표본이 아니다. 실행별 best로 성패를 정한다.
  · 통제 포함  — 원래 MSA(seedfull)를 9번째 조성으로 넣어 별도 2×2 비교를 만들지 않는다.
  · 정보 보존  — 조성들을 하나로 뭉개면(2×2) 조성 간 차이가 사라져 검정이 약해진다.

입력  : results/compreps_<target>.csv  (dump_seedrep_full.py 출력)
        컬럼 = target,model,depth,seed,pose,dockq,recall  (seed = seed<조성>_r<반복>)
출력  : results/anchor_tests.csv + 화면 표

사용(stdlib만, GPU 불필요):
  cd ~/projects/bk21-msa-depth-bias/pipeline && python anchor_test.py
  python anchor_test.py --metric recall --thr 0.4
  python anchor_test.py --only 8ulr_HL --perms 200000
"""
import argparse, csv, glob, os, random
from collections import defaultdict
from math import log


def g2(obs, ns):
    """우도비 통계량 — 집단별 성공률이 같다는 귀무가설에서의 이탈."""
    tot_o, tot_n = sum(obs), sum(ns)
    if tot_o in (0, tot_n):
        return 0.0
    p = tot_o / tot_n
    t = 0.0
    for o, m in zip(obs, ns):
        for c, e in ((o, m * p), (m - o, m * (1 - p))):
            if c > 0 and e > 0:
                t += 2 * c * log(c / e)
    return t


def perm_p(obs, ns, perms, seed=0):
    rng = random.Random(seed)
    stat = g2(obs, ns)
    if stat <= 0:
        return 1.0, stat
    pool = [1] * sum(obs) + [0] * (sum(ns) - sum(obs))
    hit = 0
    for _ in range(perms):
        rng.shuffle(pool)
        i, o = 0, []
        for m in ns:
            o.append(sum(pool[i:i + m])); i += m
        if g2(o, ns) >= stat - 1e-12:
            hit += 1
    return (hit + 1) / (perms + 1), stat


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default="results")
    ap.add_argument("--metric", default="dockq", choices=["dockq", "recall"])
    ap.add_argument("--thr", type=float, default=None,
                    help="성공 문턱 (기본: dockq 0.49 · recall 0.4)")
    ap.add_argument("--perms", type=int, default=50000)
    ap.add_argument("--only", default=None)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    thr = a.thr if a.thr is not None else (0.49 if a.metric == "dockq" else 0.4)

    files = sorted(glob.glob(os.path.join(a.dir, "compreps_*.csv")))
    if a.only:
        files = [f for f in files if a.only in f]
    if not files:
        raise SystemExit(f"!! {a.dir}/compreps_*.csv 없음")

    rows_out = []
    print(f"[{a.metric} ≥ {thr}] 실행 단위 · 순열 {a.perms:,}회 · 원래 MSA(seedfull) 포함\n")
    hdr = f"{'복합체':<12} {'조성':>4} {'실행':>5} {'성공':>5} {'조성성공률':>10} {'원래MSA':>9} {'통계량':>7} {'순열 p':>9}"
    print(hdr); print("-" * len(hdr))
    for f in files:
        runs = defaultdict(list)          # 실행(seedX_rY) -> [값…]
        tgt = None
        for r in csv.DictReader(open(f)):
            tgt = tgt or r["target"]
            v = r.get(a.metric)
            try:
                v = float(v)
            except (TypeError, ValueError):
                continue
            runs[r["seed"]].append(v)
        if not runs:
            continue
        groups = defaultdict(list)        # 조성 라벨 -> [실행별 성패]
        for run, vals in runs.items():
            comp = run.rsplit("_r", 1)[0]           # seed3_r2 -> seed3
            groups[comp].append(1 if max(vals) >= thr else 0)
        if len(groups) < 2:
            continue
        labels = sorted(groups, key=lambda k: (k == "seedfull", k))
        obs = [sum(groups[k]) for k in labels]
        ns = [len(groups[k]) for k in labels]
        p, stat = perm_p(obs, ns, a.perms)
        fi = labels.index("seedfull") if "seedfull" in labels else None
        comp_s = sum(o for i, o in enumerate(obs) if i != fi)
        comp_n = sum(m for i, m in enumerate(ns) if i != fi)
        full = f"{obs[fi]}/{ns[fi]}" if fi is not None else "없음"
        mark = " ★" if p < 0.05 and comp_s > 0 else ""
        print(f"{tgt:<12} {len(labels) - (1 if fi is not None else 0):>4} "
              f"{sum(ns):>5} {sum(obs):>5} {comp_s:>4}/{comp_n:<5} {full:>9} "
              f"{stat:>7.2f} {p:>9.5f}{mark}")
        rows_out.append({"target": tgt, "metric": a.metric, "thr": thr,
                         "n_comp": len(labels) - (1 if fi is not None else 0),
                         "n_runs": sum(ns), "comp_succ": comp_s, "comp_runs": comp_n,
                         "full_succ": obs[fi] if fi is not None else "",
                         "full_runs": ns[fi] if fi is not None else "",
                         "stat": round(stat, 3), "perm_p": round(p, 5),
                         "per_group": " ".join(f"{k}:{sum(groups[k])}/{len(groups[k])}"
                                               for k in labels)})
    out = a.out or os.path.join(a.dir, f"anchor_tests_{a.metric}.csv")
    if rows_out:
        with open(out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows_out[0]))
            w.writeheader(); w.writerows(rows_out)
        print(f"\n★ = 순열 p < 0.05 이고 조성 쪽 성공이 있음 (앵커 후보)")
        print(f"→ {out} ({len(rows_out)}행)")


if __name__ == "__main__":
    main()
