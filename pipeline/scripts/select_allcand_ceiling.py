#!/usr/bin/env python3
"""[선택기 천장] 후보를 **전부** 도킹해 보고, 우리가 고른 것이 얼마나 손해였는지 잰다.

지금 파이프라인은 후보 하나를 골라 그것만 제약으로 준다. 그래서 "다른 후보를 줬으면
어땠을까"를 몰랐다. 후보의 정답 덮음(자리 축)은 도킹 없이도 알지만, 그 자리를 줬을 때
DockQ(자세 축)가 얼마 나오는지는 돌려보지 않으면 알 수 없다.

이 스크립트가 답하는 것 두 가지:

  ① 선택기가 완벽했다면 얼마나 좋아졌을까  = 후보 중 최선(천장) − 우리가 고른 것
     이 값이 크면 병목은 **선택**이고, 0 에 가까우면 병목은 **후보 자체**다.

  ② 천장이 제약 없음을 넘는가 — **예산을 맞춰서**
     ⚠️ 천장은 후보 7개 × 자세 5개 = 자세 35개 중 최대값인데 제약 없음은 자세 5개 중
     최대값이다. 더 뽑았으니 최대값은 가만히 있어도 오른다. 그래서 같은 예산인
     '무작위 후보 하나'와 먼저 비교하고, 후보를 k 개만 봤을 때의 기대 최고값을 k 별로
     찍어 **표본의 몫**과 **선택의 몫**을 가른다. 특정 후보 하나 때문에 값이 껑충 뛰면
     그건 진짜 선택 실패이고, 완만히 오르면 그냥 많이 뽑아서다.

     이 ②는 eval_arm_verdict.py 의 '제약 자체가 해로움' 판정을 **뒤집을 수 있다.**
     arm_verdict 는 대체 자리를 둘(원래 MSA 자리·같은 크기 다른 자리)만 보는데,
     후보 목록 안에 좋은 자리가 따로 있으면 그 둘이 나빴을 뿐인 게 된다.

GPU 를 쓰지 않는다. 채점 CSV 만 읽는다.

사용:
  python -u select_allcand_ceiling.py \\
      --allcand ~/projects/epitope-guided-docking/pipeline/results/allcand \\
      --dockq   ~/projects/epitope-guided-docking/pipeline/results/demo_dockq.csv
"""
import argparse, csv, glob, os, re
import statistics as st


def f(x):
    """'0.372' → 0.372 · '  -  '·'' → None"""
    try:
        return float(str(x).strip())
    except (TypeError, ValueError):
        return None


def exp_max(vals, k):
    """후보 N 개 중 **k 개만** 뽑아 봤을 때 최고값의 기대치 (비복원).

    왜 필요한가 — '천장'은 후보 7개 × 자세 5개 = 자세 35개 중 최대값인데, 비교 대상인
    제약 없음은 자세 5개 중 최대값이다. 표본을 7배 더 뽑았으니 최대값은 가만히 있어도
    올라간다. k=1 이 **예산을 맞춘 비교**(후보 하나 = 자세 5개)이고, k 를 키우며 보면
    천장 중 얼마가 '선택을 잘해서'이고 얼마가 '많이 뽑아서'인지가 갈린다.

    오름차순 v[0..N-1] 에서 P(최대 = v[i]) = C(i, k-1) / C(N, k).
    """
    from math import comb
    v = sorted(vals); n = len(v)
    k = max(1, min(k, n))
    d = comb(n, k)
    return sum(x * comb(i, k - 1) for i, x in enumerate(v)) / d


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
        raise SystemExit(f"!! {d} 에 c*.csv 가 없다. 먼저 eval_all_candidates.sh --apply")
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
        vals = [c[1] for c in cs]
        rows.append(dict(
            target=t, n_cand=len(cs), best_cand=top_k,
            dq_noconstraint=noc, dq_ours=ours, dq_ceiling=top_q,
            dq_rand1=round(exp_max(vals, 1), 3),          # 예산을 맞춘 비교 기준
            lost=(None if ours is None else round(top_q - ours, 3)),
            rank=rank,
            rc_noconstraint=noc_r, rc_ours=ours_r, rc_ceiling=top_r,
            n_pose_min=min(c[3] for c in cs), _vals=vals))

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
        print("  ! 재도킹 대상이 '선택 실패로 분류된 타깃'뿐이면 이 두 줄은 그렇게 나올 수밖에")
        print("    없다(고른 집단이라서). 전체 집단의 선택기 성능으로 읽지 말 것.")
        big = [r for r in ok if r["lost"] is not None and r["lost"] > 0.09]
        print(f"  크게 잃은 것(>0.09)        {len(big)}종   "
              f"{', '.join(r['target'] for r in big) or '-'}")

    print("\n" + "-" * 72)
    print("  ② 천장이 제약없음을 넘는가 — 단, 예산을 맞춰야 한다")
    print("-" * 72)
    print("  천장은 후보를 전부 뒤진 최대값이라 자세를 후보 수만큼 더 뽑은 셈이다.")
    print("  같은 예산의 비교는 '무작위 후보 하나'(자세 5개) 대 '제약 없음'(자세 5개)이다.")
    fair_win = [r for r in ok if r["dq_rand1"] > r["dq_noconstraint"] + a.eps]
    fair_lose = [r for r in ok if r["dq_rand1"] < r["dq_noconstraint"] - a.eps]
    ceil_win = [r for r in ok if r["dq_ceiling"] > r["dq_noconstraint"] + a.eps]
    hopeless = [r for r in ok if r["dq_ceiling"] < r["dq_noconstraint"] - a.eps]
    print(f"\n  [예산 맞춤] 무작위 후보 하나가 제약없음보다")
    print(f"    높다  {len(fair_win):>2}종   {', '.join(r['target'] for r in fair_win) or '-'}")
    print(f"    낮다  {len(fair_lose):>2}종   {', '.join(r['target'] for r in fair_lose) or '-'}")
    print(f"\n  [예산 안 맞춤] 천장이 제약없음보다 높다 {len(ceil_win)}종 "
          f"· 천장조차 낮다 {len(hopeless)}종")
    if hopeless:
        print(f"    천장조차 낮은 것: {', '.join(r['target'] for r in hopeless)}")

    # ── 후보를 k 개만 봤을 때의 기대 최고값 — 어디부터가 '선택의 몫'인가 ──────
    print("\n" + "-" * 72)
    print("  ③ 후보를 k 개만 봤다면 (기대 최고값) — 표본의 몫과 선택의 몫을 가른다")
    print("-" * 72)
    kmax = min(7, max(r["n_cand"] for r in rows))
    head = "".join(f"{('k='+str(k)):>8}" for k in range(1, kmax + 1))
    print(f"  {'타깃':<12}{'제약없음':>10}{head}")
    for r in sorted(rows, key=lambda r: r["target"]):
        cells = "".join(
            f"{exp_max(r['_vals'], k):>8.3f}" if k <= r["n_cand"] else f"{'-':>8}"
            for k in range(1, kmax + 1))
        noc_s = f"{r['dq_noconstraint']:>10.3f}" if r["dq_noconstraint"] is not None else f"{'-':>10}"
        print(f"  {r['target']:<12}{noc_s}{cells}")
    print("\n  같은 줄에서 k 를 키울 때 값이 완만히 오르면 그건 **표본을 더 뽑아서**다.")
    print("  k=1 에서 이미 제약없음을 넘거나, 특정 후보 하나 때문에 값이 껑충 뛰면")
    print("  그건 **자리 자체가 옳아서**이고, 그때만 선택기를 고칠 값어치가 있다.")
    jump = [r for r in rows
            if r["n_cand"] >= 2 and r["dq_ceiling"] - exp_max(r["_vals"], r["n_cand"] - 1) > 0.09]
    if jump:
        print(f"\n  ⭐ 후보 하나가 판을 뒤집는 타깃 {len(jump)}종: "
              f"{', '.join(r['target'] for r in jump)}")
        print("     그 하나를 고르기만 하면 되는 것이라 **선택 실패가 맞다.**")

    thin = [r for r in rows if r["n_pose_min"] < 5]
    if thin:
        print(f"\n  ! 자세가 5개 미만인 후보가 있는 타깃 {len(thin)}종: "
              f"{', '.join(r['target'] for r in thin)}")

    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    slim = [{k: v for k, v in r.items() if not k.startswith("_")} for r in rows]
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(slim[0]))
        w.writeheader(); w.writerows(slim)
    print(f"\n→ {a.out}  ({len(rows)}종)")


if __name__ == "__main__":
    main()
