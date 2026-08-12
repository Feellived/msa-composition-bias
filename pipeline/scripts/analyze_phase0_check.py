#!/usr/bin/env python3
"""[Phase 0 확인] pose_features.csv 건전성 요약 — 재랭커(Phase 1) 들어가기 전 점검.

stdlib only (DockQ env 아니어도·로컬에서도 실행됨). lib_pose_features.py 산출 CSV를 읽어:
  · 규모(pose/target/model) · 커버리지((target,model)당 pose 수)
  · confidence 게이트(iptm/ptm/plddt 회수율) · pop-피처 가용성(group C는 pop set 없어 빔)
  · ⚠️ ipTM이 (target,model,rung) 그룹 내 상수인가(=그룹 내 pose 못 가름 → Phase 1 핵심 전제)
  · 타깃별 oracle DockQ(전 pose 최고)와 티어(0.23/0.49/0.80) · 모델별 DockQ 분포
사용:
  python analyze_phase0_check.py                         # results/pose_features.csv
  python analyze_phase0_check.py --csv 다른경로.csv
"""
import argparse, csv, math
from collections import defaultdict

NUM = ["dockq", "recall", "n_contact", "overrep", "true_rank", "pop_rank",
       "dcc_true", "dcc_pop", "iptm", "ptm", "plddt"]


def f(x):
    try:
        v = float(x)
        return None if math.isnan(v) else v
    except Exception:
        return None


def load(path):
    rows = []
    for r in csv.DictReader(open(path)):
        d = {k: r.get(k, "") for k in ("target", "group", "ab", "model", "rung", "pose", "label")}
        for k in NUM:
            d[k] = f(r.get(k, ""))
        rows.append(d)
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", default="results/pose_features.csv")
    a = ap.parse_args()
    rows = load(a.csv)
    n = len(rows)
    if n == 0:
        raise SystemExit(f"!! {a.csv} 비었음/없음")
    targets = sorted({r["target"] for r in rows})
    models = sorted({r["model"] for r in rows})
    groups = sorted({r["group"] for r in rows})
    print(f"== Phase 0 확인: {a.csv} ==")
    print(f"pose {n} · target {len(targets)} · model {models} · group {groups}")

    print("\n[confidence 게이트]  (100%면 Phase 1에서 ipTM 대비·신뢰도 피처 사용 가능)")
    for k in ("iptm", "ptm", "plddt"):
        c = sum(1 for r in rows if r[k] is not None)
        print(f"    {k:6}: {c}/{n} ({100*c/n:.0f}%)")
    print("[pop-피처 가용성]  (group C는 popular set 없어 비어있음 = 정상)")
    for k in ("overrep", "dcc_pop", "pop_rank"):
        c = sum(1 for r in rows if r[k] is not None)
        print(f"    {k:9}: {c}/{n} ({100*c/n:.0f}%)")

    # ⚠️ ipTM이 그룹 내 상수인가
    g = defaultdict(set)
    for r in rows:
        if r["iptm"] is not None:
            g[(r["target"], r["model"], r["rung"])].add(round(r["iptm"], 4))
    distinct = [len(v) for v in g.values()]
    if distinct:
        avg = sum(distinct) / len(distinct)
        print(f"\n[⚠️ ipTM 판별력]  (target,model,rung) 그룹 내 distinct iptm 평균 = {avg:.2f}")
        print(f"    → {'≈1.0 = 그룹 내 5 pose를 전혀 못 가름(ipTM은 그룹 내 무력 선택기)' if avg < 1.5 else '그룹 내 변동 있음'}")

    # 커버리지
    cov = defaultdict(int)
    for r in rows:
        cov[(r["target"], r["model"])] += 1
    cv = sorted(cov.values())
    print(f"\n[커버리지] (target,model)당 pose 수: min {cv[0]} · median {cv[len(cv)//2]} · max {cv[-1]}"
          f"  (기대 ≈ 12 rung × 5 = 60)")

    # 타깃별 oracle + 티어
    byt = defaultdict(list)
    for r in rows:
        if r["dockq"] is not None:
            byt[r["target"]].append(r["dockq"])
    tiers = {0.23: [], 0.49: [], 0.80: []}
    for t in targets:
        o = max(byt[t]) if byt[t] else None
        for th in tiers:
            if o is not None and o >= th:
                tiers[th].append(t)
    print(f"\n[타깃별 oracle DockQ(전 pose 최고) 티어]  (전체 {len(targets)} 타깃)")
    for th in (0.23, 0.49, 0.80):
        print(f"    ≥{th} (Acceptable/Medium/High): {len(tiers[th]):>3} 타깃")

    print("\n[모델별 DockQ 분포]")
    for m in models:
        v = sorted(r["dockq"] for r in rows if r["model"] == m and r["dockq"] is not None)
        if v:
            print(f"    {m:9} n={len(v):5}  mean {sum(v)/len(v):.3f}  median {v[len(v)//2]:.3f}  max {max(v):.3f}"
                  f"  (≥0.49: {sum(1 for x in v if x>=0.49)})")

    print("\n다음 → analyze_phase1_rerank.py : 무학습 재랭커 make-or-break (ipTM vs 배포가능 기하 피처).")


if __name__ == "__main__":
    main()
