#!/usr/bin/env python3
"""[선택기 천장] 후보를 **전부** 도킹해 보고, 우리가 고른 것이 얼마나 손해였는지 잰다.

지금 파이프라인은 후보 하나를 골라 그것만 제약으로 준다. 그래서 "다른 후보를 줬으면
어땠을까"를 몰랐다. 후보의 정답 덮음(자리 축)은 도킹 없이도 알지만, 그 자리를 줬을 때
DockQ(자세 축)가 얼마 나오는지는 돌려보지 않으면 알 수 없다.

이 스크립트가 답하는 것 두 가지:

  ① 선택기가 완벽했다면 얼마나 좋아졌을까  = 후보 중 최선(천장) − 우리가 고른 것
     이 값이 크면 병목은 **선택**이고, 0 에 가까우면 병목은 **후보 자체**다.

  ② 천장조차 제약 없음보다 낮은 복합체가 있는가
     있다면 그 복합체는 선택을 아무리 잘해도 안 된다. 즉 **제약을 주면 안 되는 복합체**다.
     (arm_verdict.py 의 '제약 자체가 해로움' 판정과 같은 것을 다른 길로 확인하는 셈이다.)

GPU 를 쓰지 않는다. 채점 CSV 만 읽는다.

사용:
  python -u allcand_ceiling.py \\
      --allcand ~/projects/bk21-antibody-ml/consensus_docking/results/allcand \\
      --dockq   ~/projects/bk21-antibody-ml/consensus_docking/results/demo_dockq.csv
"""
import argparse, csv, glob, os, re
import statistics as st


def f(x):
    """'0.372' → 0.372 · '  -  '·'' → None"""
    try:
        return float(str(x).strip())
    except (TypeError, ValueError):
        return None


def load_main(path, model):
    """원래 표에서 (타깃, 팔) → (DockQ, 자리겹침)"""
    dq, rc = {}, {}
    for r in csv.DictReader(open(path)):
        if r.get("model") != model:
            continue
        k = (r["target"], r["arm"])
        dq[k] = f(r.get("dockq_max")); rc[k] = f(r.get("recall_max"))
    return dq, rc


def load_cands(d, model):
    """후보 폴더별 CSV → 타깃 → [(후보번호, DockQ, 자리겹침, 자세수)]"""
    out = {}
    files = sorted(glob.glob(os.path.join(d, "c*.csv")),
                   key=lambda p: int(re.search(r"c(\d+)", os.path.basename(p)).group(1)))
    if not files:
        raise SystemExit(f"!! {d} 에 c*.csv 가 없다. 먼저 score_allcand.sh --apply")
    for p in files:
        k = int(re.search(r"c(\d+)", os.path.basename(p)).group(1))
        for r in csv.DictReader(open(p)):
            if r.get("model") != model or r.get("arm") != "ours":
                continue
            q = f(r.get("dockq_max"))
            if q is None:
                continue
            out.setdefault(r["target"], []).append(
                (k, q, f(r.get("recall_max")), int(r.get("n_pose") or 0)))
    return out, files


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--allcand", required=True, help="results/allcand 폴더")
    ap.add_argument("--dockq", required=True, help="원래 30종 채점 CSV")
    ap.add_argument("--model", default="boltz")
    ap.add_argument("--eps", type=float, default=0.02, help="이 미만 차이는 같은 것으로 본다")
    ap.add_argument("--out", default="results/allcand_ceiling.csv")
    a = ap.parse_args()

    cands, files = load_cands(a.allcand, a.model)
    dq, rc = load_main(a.dockq, a.model)
    print(f"후보 CSV {len(files)}개 · 타깃 {len(cands)}종 · 모델 {a.model}")

    rows = []
    for t in sorted(cands):
        cs = sorted(cands[t], key=lambda x: -x[1])
        top_k, top_q, top_r, _ = cs[0]
        ours = dq.get((t, "ours"))
        noc = dq.get((t, "noconstraint"))
        ours_r = rc.get((t, "ours")); noc_r = rc.get((t, "noconstraint"))
        # 우리 값이 후보들 사이에서 몇 등인가 (같은 값이면 더 좋은 쪽으로 세지 않는다)
        rank = (sum(1 for c in cs if ours is not None and c[1] > ours + a.eps) + 1
                if ours is not None else None)
        rows.append(dict(
            target=t, n_cand=len(cs), best_cand=top_k,
            dq_noconstraint=noc, dq_ours=ours, dq_ceiling=top_q,
            lost=(None if ours is None else round(top_q - ours, 3)),
            rank=rank,
            rc_noconstraint=noc_r, rc_ours=ours_r, rc_ceiling=top_r,
            n_pose_min=min(c[3] for c in cs)))

    W = [("target", "타깃", 12), ("n_cand", "후보수", 7), ("dq_noconstraint", "제약없음", 10),
         ("dq_ours", "우리선택", 10), ("dq_ceiling", "천장", 8), ("lost", "잃은양", 9),
         ("rank", "우리등수", 9), ("best_cand", "최선후보", 9)]
    print("\n" + "=" * 80)
    print("  후보를 전부 도킹했을 때 — 잃은 양이 큰 순")
    print("=" * 80)
    print("  " + "".join(f"{h:>{w}}" for _, h, w in W))
    for r in sorted(rows, key=lambda r: -(r["lost"] if r["lost"] is not None else -9)):
        print("  " + "".join(
            f"{(f'{r[c]:.3f}' if isinstance(r[c], float) else str(r[c] if r[c] is not None else '-')):>{w}}"
            for c, _, w in W))

    ok = [r for r in rows if None not in (r["dq_ours"], r["dq_noconstraint"])]
    lost = [r["lost"] for r in ok if r["lost"] is not None]
    print("\n" + "-" * 72)
    print("  ① 병목이 선택인가, 후보인가")
    print("-" * 72)
    if lost:
        print(f"  잃은 양(천장 − 우리 선택)   중앙값 {st.median(lost):+.3f} · 합 {sum(lost):+.3f}")
        first = [r for r in ok if r["rank"] == 1]
        print(f"  우리가 이미 최선을 골랐다   {len(first)}/{len(ok)}종")
        big = [r for r in ok if r["lost"] is not None and r["lost"] > 0.09]
        print(f"  크게 잃은 것(>0.09)        {len(big)}종   "
              f"{', '.join(r['target'] for r in big) or '-'}")

    print("\n" + "-" * 72)
    print("  ② 천장조차 제약 없음보다 낮은 복합체 = 어떤 후보를 줘도 안 된다")
    print("-" * 72)
    hopeless = [r for r in ok if r["dq_ceiling"] < r["dq_noconstraint"] - a.eps]
    winnable = [r for r in ok if r["dq_ceiling"] > r["dq_noconstraint"] + a.eps]
    actual = [r for r in ok if r["dq_ours"] > r["dq_noconstraint"] + a.eps]
    print(f"  천장이 제약없음보다 낮다   {len(hopeless):>2}종   "
          f"{', '.join(r['target'] for r in hopeless) or '-'}")
    print(f"  천장이 제약없음보다 높다   {len(winnable):>2}종   "
          f"{', '.join(r['target'] for r in winnable) or '-'}")
    print(f"  실제로 우리가 이긴 것      {len(actual):>2}종   "
          f"{', '.join(r['target'] for r in actual) or '-'}")
    if winnable:
        miss = [r for r in winnable if r not in actual]
        print(f"\n  ⭐ 이길 수 있었는데 놓친 것 {len(miss)}종: "
              f"{', '.join(r['target'] for r in miss) or '-'}")
        print("     → 이 종수만큼이 '선택을 고치면 얻을 수 있는 몫'이다.")
    if hopeless:
        print(f"\n  ⭐ {len(hopeless)}종은 후보를 전부 시도해도 제약 없음을 못 넘었다.")
        print("     선택의 문제가 아니라 **이 복합체에는 제약을 주면 안 된다**는 뜻이다.")

    thin = [r for r in rows if r["n_pose_min"] < 5]
    if thin:
        print(f"\n  ! 자세가 5개 미만인 후보가 있는 타깃 {len(thin)}종: "
              f"{', '.join(r['target'] for r in thin)}")

    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)
    print(f"\n→ {a.out}  ({len(rows)}종)")


if __name__ == "__main__":
    main()
