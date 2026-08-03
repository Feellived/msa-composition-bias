#!/usr/bin/env python3
"""[⓪ 적용 기준] 정답을 안 보고도 "이 복합체에 우리 방법을 써도 되나"를 가를 수 있나.

지금 가진 신호는 **후보 개수** 하나뿐이고, 그것도 23종을 다 보고 정한 사후 규칙이다.
여기서는 **제약 없이 돌린 예측 그 자체**에서 신호를 더 캐낸다. 전부 정답이 필요 없다.

  ① iptm      — 모델이 스스로 매긴 계면 신뢰도 (Boltz 가 자세마다 남긴다)
  ② 자세일치   — 무제약 자세 5개가 **서로 같은 자리**를 가리키나 (쌍마다 겹침을 재 평균)
  ③ 자리개수   — 그 5개 자세가 몇 군데로 흩어지나 (겹침 0.5 이상이면 같은 자리로 묶음)
  ④ 후보개수   — 이미 쓰고 있는 것 (sites JSON 에서 읽는다)

그리고 이 신호들이 두 가지를 얼마나 맞히는지 잰다. **여기서만 정답을 쓴다(채점용).**
  A. 크게 망가질 것인가        — 우리 자리를 줬을 때 DockQ 가 0.09 넘게 떨어지나
  B. 무제약이 실패인가          — 제약 없이 돌린 DockQ 가 0.23 미만인가

⚠️ 신호가 A 를 잘 맞히면 **개입 여부를 정답 없이 정할 수 있다**는 뜻이다.
   B 를 잘 맞히면 ④단계(새 복합체 선별)에서 **실패할 것 같은 종을 미리 고를 수 있다**.
   둘 다 되면 같은 장치로 두 문제가 풀린다.

GPU 를 쓰지 않는다. 이미 만들어 둔 구조 파일과 confidence JSON 만 읽는다.

사용 (conda activate boltz · pipeline/ 에서):
  python -u gate_signals.py
  python -u gate_signals.py --runs ~/projects/bk21-antibody-ml/consensus_docking/runs_sites_guided_honest
"""
import argparse, csv, glob, json, os
import statistics as st
from collections import defaultdict

import pose_features as PF
import epitope_cluster as EC

AB = os.path.expanduser("~/projects/bk21-antibody-ml/consensus_docking")


# ── 채점용 도구 ─────────────────────────────────────────────────────────────
def auc(pos, neg):
    """순위 기반 판별력. 0.5 면 동전던지기, 1.0 이면 완전히 갈린다. 동점 처리 포함."""
    pos = [v for v in pos if v == v]
    neg = [v for v in neg if v == v]
    if not pos or not neg:
        return float("nan"), len(pos), len(neg)
    xs = sorted([(v, 1) for v in pos] + [(v, 0) for v in neg])
    r, i = [0.0] * len(xs), 0
    while i < len(xs):
        j = i
        while j + 1 < len(xs) and xs[j + 1][0] == xs[i][0]:
            j += 1
        avg = (i + j) / 2 + 1
        for k in range(i, j + 1):
            r[k] = avg
        i = j + 1
    rs = sum(r[k] for k in range(len(xs)) if xs[k][1] == 1)
    n1, n0 = len(pos), len(neg)
    return (rs - n1 * (n1 + 1) / 2) / (n1 * n0), n1, n0


def best_cut(vals, lab):
    """양성을 가장 잘 걸러내는 문턱 하나. (문턱, 방향, 정확도) — 어디까지나 참고용."""
    pairs = [(v, l) for v, l in zip(vals, lab) if v == v]
    if not pairs:
        return None
    best = None
    for thr in sorted({v for v, _ in pairs}):
        for sign in (+1, -1):
            hit = sum(1 for v, l in pairs
                      if ((v >= thr) if sign > 0 else (v < thr)) == bool(l))
            if best is None or hit > best[2]:
                best = (thr, sign, hit)
    thr, sign, hit = best
    return thr, ("이상" if sign > 0 else "미만"), hit / len(pairs)


# ── 신호 뽑기 ───────────────────────────────────────────────────────────────
def poses_of(runs, tgt, arm):
    """그 팔의 예측 자세 파일들. 경로 안에 팔 이름이 들어 있다."""
    pat = os.path.join(runs, tgt.lower(), "**", "*model*.cif")
    return sorted(p for p in glob.glob(pat, recursive=True)
                  if f"/{arm}/" in p or f"_{arm}_" in p or f"/{arm}_" in p)


def iptm_of(cif):
    """자세 파일 옆의 confidence JSON 에서 계면 신뢰도. 없으면 NaN."""
    d, base = os.path.dirname(cif), os.path.basename(cif)
    stem = base.rsplit(".", 1)[0]
    cands = [os.path.join(d, f"confidence_{stem}.json"),
             os.path.join(d, f"{stem}_confidence.json")]
    cands += sorted(glob.glob(os.path.join(d, "confidence*.json")))
    for c in cands:
        if os.path.exists(c):
            try:
                j = json.load(open(c))
            except Exception:
                continue
            for k in ("iptm", "ipTM", "interface_ptm", "complex_iplddt"):
                if k in j:
                    try:
                        return float(j[k])
                    except Exception:
                        pass
    return float("nan")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", default=os.path.join(AB, "runs_sites_guided_honest"))
    ap.add_argument("--diag", default="results/harm_diag.csv",
                    help="diagnose_harm.py 결과 (타깃·변화량·후보수)")
    ap.add_argument("--sites", default="results/honest")
    ap.add_argument("--targets-dir", default="targets")
    ap.add_argument("--arm", default="noconstraint", help="신호를 뽑을 팔")
    ap.add_argument("--cutoff", type=float, default=5.0)
    ap.add_argument("--link", type=float, default=0.5,
                    help="자세끼리 이만큼 겹치면 같은 자리로 본다")
    ap.add_argument("--harm", type=float, default=-0.09,
                    help="이보다 더 떨어지면 '크게 망가짐'")
    ap.add_argument("--out", default="results/gate_signals.csv")
    a = ap.parse_args()

    if not os.path.isdir(a.runs):
        raise SystemExit(f"!! 예측 폴더가 없다: {a.runs}\n"
                         f"   --runs 로 경로를 줄 것 "
                         f"(run_honest_guided.sh 의 OUTDIR)")
    diag = list(csv.DictReader(open(a.diag)))
    if not diag:
        raise SystemExit(f"!! {a.diag} 가 비었다 — diagnose_harm.py 를 먼저 돌릴 것")
    print(f"타깃 {len(diag)}종 · 팔={a.arm} · 예측폴더={a.runs}\n")

    rows, nofile = [], []
    for d in diag:
        tgt = d["target"]
        cifs = poses_of(a.runs, tgt, a.arm)
        if not cifs:
            nofile.append(tgt)
            continue
        cj = json.load(open(os.path.join(a.targets_dir, tgt, "chains.json")))

        eps, ip = [], []
        for c in cifs:
            ep, _ = EC.pred_epitope(cj, c, a.cutoff)
            if ep:
                eps.append(ep)
            v = iptm_of(c)
            if v == v:
                ip.append(v)
        if len(eps) < 2:
            nofile.append(tgt + "(자세 부족)")
            continue

        # ② 자세끼리 얼마나 같은 자리를 가리키나
        jj = [EC.jac(eps[i], eps[j]) for i in range(len(eps)) for j in range(i + 1, len(eps))]
        jj = [v for v in jj if v == v]
        agree = st.mean(jj) if jj else float("nan")

        # ③ 몇 군데로 흩어지나 (단일연결로 묶는다)
        par = list(range(len(eps)))
        def find(x):
            while par[x] != x:
                par[x] = par[par[x]]; x = par[x]
            return x
        for i in range(len(eps)):
            for j in range(i + 1, len(eps)):
                if EC.jac(eps[i], eps[j]) >= a.link:
                    ra, rb = find(i), find(j)
                    if ra != rb:
                        par[ra] = rb
        nsite = len({find(i) for i in range(len(eps))})

        rows.append(dict(
            target=tgt, n_pose=len(eps),
            iptm=round(st.mean(ip), 4) if ip else float("nan"),
            agree=round(agree, 4), nsite=nsite,
            n_cand=int(d["n_cand"]),
            delta=float(d["delta"]), dq_no=float(d["dq_no"]), why=d["why"]))

    if nofile:
        print(f"  ! 자세를 못 읽은 타깃 {len(nofile)}종: {', '.join(nofile[:8])}"
              f"{' …' if len(nofile) > 8 else ''}")
    if not rows:
        raise SystemExit("!! 아무 타깃도 못 읽었다 — --runs 경로와 --arm 이름을 확인할 것")

    SIG = [("iptm", "모델 신뢰도"), ("agree", "자세끼리 같은 자리를 가리키나"),
           ("nsite", "자세가 몇 군데로 흩어지나"), ("n_cand", "만들어진 후보 개수")]

    print("\n" + "=" * 96)
    print("  타깃별 신호 — 변화량 오름차순")
    print("=" * 96)
    W = ["target", "iptm", "agree", "nsite", "n_cand", "dq_no", "delta", "why"]
    H = {"target": "타깃", "iptm": "신뢰도", "agree": "자세일치", "nsite": "자리수",
         "n_cand": "후보수", "dq_no": "무제약", "delta": "변화", "why": "원인"}
    rows.sort(key=lambda r: r["delta"])
    print("  ".join(f"{H[c]:>8}" for c in W))
    for r in rows:
        print("  ".join(f"{r[c]!s:>8}" for c in W))

    def report(title, is_pos, note):
        print("\n" + "-" * 84)
        print(f"  {title}")
        print("-" * 84)
        lab = [1 if is_pos(r) else 0 for r in rows]
        npos = sum(lab)
        if npos in (0, len(rows)):
            print(f"  한쪽으로 몰려 잴 수 없다 (양성 {npos}/{len(rows)})"); return
        print(f"  양성 {npos}종 · 음성 {len(rows)-npos}종      {note}")
        print(f"  {'신호':<28}{'판별력':>7}   {'가장 잘 가르는 문턱':>22}")
        for key, name in SIG:
            vals = [r[key] for r in rows]
            A, n1, n0 = auc([v for v, l in zip(vals, lab) if l],
                            [v for v, l in zip(vals, lab) if not l])
            # 0.5 아래면 방향이 반대라는 뜻 — 크기만 보면 되므로 뒤집어 함께 적는다
            flag = "" if A != A else ("  (방향 반대)" if A < 0.5 else "")
            cut = best_cut(vals, lab)
            ctxt = "-" if not cut else f"{key} {cut[0]:g} {cut[1]}  정확도 {cut[2]:.0%}"
            print(f"  {name:<28}{A:>7.3f}   {ctxt:>22}{flag}")

    report("A. 크게 망가질 것인가 — 개입 여부를 정할 수 있나",
           lambda r: r["delta"] < a.harm,
           f"(변화 < {a.harm})")
    report("B. 무제약이 실패인가 — 새 복합체를 미리 고를 수 있나",
           lambda r: r["dq_no"] < 0.23,
           "(무제약 DockQ < 0.23)")

    print("\n  ※ 판별력 0.5 = 동전던지기 · 0.7 이상이면 쓸 만함 · 1.0 = 완전히 갈림")
    print("  ※ 문턱과 정확도는 같은 자료에서 고른 값이라 낙관적이다. 새 세트로 확인해야 한다.")

    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)
    print(f"\n→ {a.out}  ({len(rows)}종)")


if __name__ == "__main__":
    main()
