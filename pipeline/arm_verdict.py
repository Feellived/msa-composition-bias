#!/usr/bin/env python3
"""[재분류] 손해가 '우리가 자리를 잘못 골라서'인가, '제약 자체가 해로워서'인가.

지금까지는 우리 자리(ours)와 제약 없음(noconstraint)만 비교해 손해를 셌다. 그런데 대조 팔이
두 개 더 있다 — 원래 MSA 가 간 자리(fullmsa)와 같은 크기의 다른 자리(sizematch)다.
**그 둘도 같이 무너졌다면 우리 선택 탓이 아니다.** 그 복합체는 어떤 자리를 줘도 망가진다.

8sit_HL 이 그랬다. 제약 없이 0.808 인데 세 제약 조건이 전부 0.025~0.026 으로 무너졌다.
자리를 잘못 골라서가 아니라 **제약을 주는 행위 자체가 해로운 복합체**다.

판정 (네 팔이 모두 있는 타깃만):
  우리만 이득       ours 가 나머지 셋보다 모두 높다            ← 자리가 옳았다는 증거
  제약 자체가 해로움 세 제약이 모두 크게 떨어졌다               ← 우리 탓이 아니다
  우리만 해로움     ours 만 크게 떨어졌다                      ← 선택 문제
  차이 없음         그 밖

GPU 를 쓰지 않는다. 채점 CSV 하나만 읽는다.

사용:
  python -u arm_verdict.py --dockq ~/projects/bk21-antibody-ml/consensus_docking/results/demo_dockq.csv
"""
import argparse, csv, os
import statistics as st

ARMS = ["noconstraint", "fullmsa", "sizematch", "ours"]
NAME = {"noconstraint": "제약없음", "fullmsa": "원래MSA자리",
        "sizematch": "같은크기다른자리", "ours": "우리자리"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dockq", required=True)
    ap.add_argument("--model", default="boltz", help="이 모델의 행만 쓴다")
    ap.add_argument("--harm", type=float, default=0.09, help="이보다 더 떨어지면 '크게 떨어짐'")
    ap.add_argument("--eps", type=float, default=0.02, help="이 미만 차이는 같은 것으로 본다")
    ap.add_argument("--out", default="results/arm_verdict.csv")
    a = ap.parse_args()

    rows = list(csv.DictReader(open(a.dockq)))
    if not rows:
        raise SystemExit(f"!! 비어 있음: {a.dockq}")
    need = {"target", "arm", "model", "dockq_max", "recall_max"}
    miss = need - set(rows[0])
    if miss:
        raise SystemExit(f"!! 열이 없다: {sorted(miss)}\n   있는 열 = {list(rows[0])}")

    dq, rc = {}, {}
    for r in rows:
        if r["model"] != a.model:
            continue
        try:
            dq[(r["target"], r["arm"])] = float(r["dockq_max"])
            rc[(r["target"], r["arm"])] = float(r["recall_max"])
        except ValueError:
            pass
    tgts = sorted({t for t, _ in dq})
    full4 = [t for t in tgts if all((t, x) in dq for x in ARMS)]
    print(f"모델={a.model} · 타깃 {len(tgts)}종 · 네 팔이 모두 있는 것 {len(full4)}종")
    if len(tgts) > len(full4):
        miss_t = [t for t in tgts if t not in full4]
        print(f"  ! 네 팔이 없어 판정 불가 {len(miss_t)}종: {', '.join(miss_t)}")

    out = []
    for t in full4:
        d = {x: dq[(t, x)] for x in ARMS}
        base = d["noconstraint"]
        drop = {x: base - d[x] for x in ("fullmsa", "sizematch", "ours")}
        if all(d["ours"] > d[x] + a.eps for x in ("noconstraint", "fullmsa", "sizematch")):
            v = "우리만 이득"
        elif min(drop.values()) > a.harm:
            v = "제약 자체가 해로움"
        elif drop["ours"] > a.harm:
            v = "우리만 해로움"
        elif drop["ours"] < -a.eps:
            v = "우리가 나음"
        else:
            v = "차이 없음"
        out.append(dict(target=t, verdict=v, **{f"dq_{x}": round(d[x], 3) for x in ARMS},
                        drop_full=round(drop["fullmsa"], 3),
                        drop_size=round(drop["sizematch"], 3),
                        drop_ours=round(drop["ours"], 3),
                        rc_no=round(rc[(t, "noconstraint")], 3),
                        rc_ours=round(rc[(t, "ours")], 3)))

    out.sort(key=lambda r: r["drop_ours"], reverse=True)
    W = ["target", "verdict", "dq_noconstraint", "dq_fullmsa", "dq_sizematch", "dq_ours",
         "drop_full", "drop_size", "drop_ours"]
    H = {"target": "타깃", "verdict": "판정", "dq_noconstraint": "제약없음",
         "dq_fullmsa": "원래MSA", "dq_sizematch": "같은크기", "dq_ours": "우리자리",
         "drop_full": "원래하락", "drop_size": "크기하락", "drop_ours": "우리하락"}
    print("\n" + "=" * 108)
    print("  네 조건 비교 — 우리 하락이 큰 순")
    print("=" * 108)
    print("  " + "".join(f"{H[c]:>12}" if c != "verdict" else f"{H[c]:>18}" for c in W))
    for r in out:
        print("  " + "".join(f"{r[c]!s:>12}" if c != "verdict" else f"{r[c]!s:>18}" for c in W))

    print("\n" + "-" * 70)
    print("  판정별 종수")
    print("-" * 70)
    for v in ("우리만 이득", "우리가 나음", "차이 없음", "우리만 해로움", "제약 자체가 해로움"):
        s = [r for r in out if r["verdict"] == v]
        if s:
            print(f"  {v:<18} {len(s):>2}종   {', '.join(r['target'] for r in s)}")

    big = [r for r in out if r["drop_ours"] > a.harm]
    mine = [r for r in big if r["verdict"] == "우리만 해로움"]
    theirs = [r for r in big if r["verdict"] == "제약 자체가 해로움"]
    print("\n" + "-" * 70)
    print(f"  크게 망가진 {len(big)}종의 책임 소재  (하락 > {a.harm})")
    print("-" * 70)
    print(f"  우리 선택 탓        {len(mine):>2}종   {', '.join(r['target'] for r in mine) or '-'}")
    print(f"  제약 자체가 해로움  {len(theirs):>2}종   {', '.join(r['target'] for r in theirs) or '-'}")
    if theirs:
        print(f"\n  ⭐ 이 {len(theirs)}종은 원래 MSA 가 간 자리를 줘도, 같은 크기의 다른 자리를 줘도")
        print(f"     똑같이 무너졌다. 자리를 잘못 골라서가 아니라 **제약을 주면 안 되는 복합체**다.")

    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0]))
        w.writeheader(); w.writerows(out)
    print(f"\n→ {a.out}  ({len(out)}종)")


if __name__ == "__main__":
    main()
