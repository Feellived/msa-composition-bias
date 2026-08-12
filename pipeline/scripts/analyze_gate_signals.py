#!/usr/bin/env python3
"""[⓪ 적용 기준] 정답을 안 보고도 "이 복합체에 우리 방법을 써도 되나"를 가를 수 있나.

지금 가진 신호는 **후보 개수** 하나뿐이고, 그것도 23종을 다 보고 정한 사후 규칙이다.
여기서는 **제약 없이 돌린 예측 그 자체**에서 신호를 더 캐낸다. 전부 정답이 필요 없다.

  ① 신뢰도    — 모델이 스스로 매긴 계면 점수 (Boltz 가 자세마다 남긴다)
  ② 자세일치   — 무제약 자세들이 **서로 같은 자리**를 가리키나 (쌍마다 겹침을 재 평균)
  ③ 자리수    — 그 자세들이 몇 군데로 흩어지나 (겹침 0.5 이상이면 같은 자리로 묶음)
  ④ 후보개수   — 이미 쓰고 있는 것 (analyze_harm.py 결과에서 읽는다)

그리고 이 신호들이 두 가지를 얼마나 맞히는지 잰다. **여기서만 정답을 쓴다(채점용).**
  A. 크게 망가질 것인가        — 우리 자리를 줬을 때 DockQ 가 정해진 값보다 더 떨어지나
  B. 무제약이 실패인가          — 제약 없이 돌린 DockQ 가 0.23 미만인가

⚠️ 표본이 23종뿐이라 판별력 하나만 보면 안 된다. 그래서 셋을 같이 낸다.
   · **판별력** = 방향과 무관한 크기(0.5=동전던지기, 1.0=완전히 갈림)
   · **95% 구간** = 타깃을 복원추출해 다시 잰 값의 범위. 좁아야 믿을 수 있다
   · **p** = 딱지를 무작위로 다시 붙였을 때 이만큼 갈릴 확률
⚠️ 정확도는 **기저율**(다수쪽만 찍었을 때의 정확도)과 반드시 같이 본다. 양성이 6/23 이면
   "아무 일도 안 난다"고만 답해도 74% 다.

GPU 를 쓰지 않는다. 이미 만들어 둔 구조 파일과 confidence JSON 만 읽는다.

사용 (conda activate boltz · pipeline/ 에서):
  python -u analyze_gate_signals.py --selftest      # 통계 부분만 자체 검증 (자료 없이 즉시)
  python -u analyze_gate_signals.py
  python -u analyze_gate_signals.py --runs ~/projects/bk21-antibody-ml/pipeline/runs_sites_guided_honest
"""
import argparse, csv, glob, json, os, random, re
import statistics as st

AB = os.path.expanduser("~/projects/bk21-antibody-ml/pipeline")


# ── 채점용 도구 ─────────────────────────────────────────────────────────────
def auc(pos, neg):
    """순위 기반 판별력. 0.5 면 동전던지기, 1.0 이면 완전히 갈린다. 동점 처리 포함."""
    pos = [v for v in pos if v == v]
    neg = [v for v in neg if v == v]
    if not pos or not neg:
        return float("nan")
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
    return (rs - n1 * (n1 + 1) / 2) / (n1 * n0)


def stats(vals, lab, nboot=2000, nperm=5000, seed=0):
    """판별력 + 95% 구간(복원추출) + p(뒤섞기). 방향은 따로 알려준다."""
    pairs = [(v, l) for v, l in zip(vals, lab) if v == v]
    P = [v for v, l in pairs if l]
    N = [v for v, l in pairs if not l]
    if not P or not N:
        return None
    A = auc(P, N)
    D = max(A, 1 - A)                       # 방향과 무관한 크기
    rng = random.Random(seed)

    bs = []
    for _ in range(nboot):                  # 타깃을 복원추출해 다시 잰다
        s = [pairs[rng.randrange(len(pairs))] for _ in pairs]
        p = [v for v, l in s if l]
        n = [v for v, l in s if not l]
        if not p or not n:
            continue
        v = auc(p, n)
        if v == v:                          # 전부 한쪽이면 NaN 이 나온다 — 버린다
            bs.append(max(v, 1 - v))
    bs.sort()
    lo = bs[int(0.025 * len(bs))] if len(bs) >= 40 else float("nan")
    hi = bs[int(0.975 * len(bs))] if len(bs) >= 40 else float("nan")

    vs = [v for v, _ in pairs]
    ls = [l for _, l in pairs]
    hit = 0
    for _ in range(nperm):                  # 딱지를 무작위로 다시 붙인다
        rng.shuffle(ls)
        p = [v for v, l in zip(vs, ls) if l]
        n = [v for v, l in zip(vs, ls) if not l]
        a = auc(p, n)
        if max(a, 1 - a) >= D - 1e-12:
            hit += 1
    return dict(disc=D, lo=lo, hi=hi, p=(hit + 1) / (nperm + 1), n=len(pairs),
                dirn="클수록 위험" if A > 0.5 else "작을수록 위험")


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
    """그 팔의 예측 자세 파일들. 경로 안에 팔 이름이 폴더로 들어 있다."""
    pat = os.path.join(runs, tgt.lower(), "**", "*model*.cif")
    return sorted(p for p in glob.glob(pat, recursive=True)
                  if os.sep + arm + os.sep in p)


def iptm_of(cif):
    """자세 **그 자세의** confidence JSON 에서 계면 신뢰도. 못 찾으면 NaN.

    ⚠️ 예전 판은 못 찾으면 폴더 안 아무 confidence*.json 이나 집었다. 그러면 자세마다
       같은 값이 박혀 평균이 그 값 하나가 되고, 그걸 알아챌 방법이 없다. 그래서
       **이름이 정확히 맞거나 자세 번호가 같은 것만** 쓴다.
    """
    d = os.path.dirname(cif)
    stem = os.path.basename(cif).rsplit(".", 1)[0]
    cands = [os.path.join(d, f"confidence_{stem}.json"),
             os.path.join(d, f"{stem}_confidence.json"),
             os.path.join(d, f"{stem}.json")]
    m = re.search(r"model_(\d+)", stem)
    if m:                                    # 자세 번호가 같은 것만 추가로 본다
        cands += sorted(glob.glob(os.path.join(d, f"*model_{m.group(1)}.json")))
    for c in cands:
        if not os.path.exists(c):
            continue
        try:
            j = json.load(open(c))
        except Exception:
            continue
        for k in ("iptm", "ipTM", "interface_ptm"):
            if k in j:
                try:
                    return float(j[k])
                except Exception:
                    pass
    return float("nan")


# ── 자체 검증 ───────────────────────────────────────────────────────────────
def selftest():
    """통계 부분만 자료 없이 확인한다. 자료가 준비되기 전에 돌려볼 수 있다."""
    ok = True

    def chk(name, got, want, tol=1e-9):
        nonlocal ok
        good = (got != got and want != want) or abs(got - want) <= tol
        ok &= good
        print(f"  {'통과' if good else '실패':<4} {name:<34} 나온값 {got!s:<10} 기대 {want}")

    print("[1] 판별력(AUC) 계산")
    chk("완전히 갈림", auc([3, 4, 5], [0, 1, 2]), 1.0)
    chk("방향이 반대", auc([0, 1, 2], [3, 4, 5]), 0.0)
    chk("반반", auc([1, 4], [2, 3]), 0.5)
    chk("전부 동점", auc([2, 2], [2, 2]), 0.5)
    chk("한쪽이 비었음", auc([], [1, 2]), float("nan"))

    print("\n[2] 문턱 찾기")
    c = best_cut([1, 2, 3, 9, 8, 7], [0, 0, 0, 1, 1, 1])
    chk("정확도", c[2], 1.0); print(f"       → 문턱 {c[0]} {c[1]}")
    c = best_cut([9, 8, 7, 1, 2, 3], [0, 0, 0, 1, 1, 1])
    print(f"  {'통과' if c[1] == '미만' else '실패'}  방향이 뒤집힌 경우도 잡는다 → {c[1]}")
    ok &= c[1] == "미만"

    print("\n[3] 구간과 p — 완전히 갈리는 신호")
    s = stats([9, 8, 7, 6, 1, 2, 3, 4], [1, 1, 1, 1, 0, 0, 0, 0], nboot=500, nperm=2000)
    print(f"       판별력 {s['disc']:.3f}  구간 [{s['lo']:.2f}, {s['hi']:.2f}]  p={s['p']:.4f}")
    good = s["disc"] == 1.0 and s["p"] < 0.05
    ok &= good; print(f"  {'통과' if good else '실패'}  판별력 1.0 이고 p 가 작아야 한다")

    print("\n[4] 구간과 p — 아무 정보도 없는 신호 (딱지와 무관)")
    rng = random.Random(7)
    v = [rng.random() for _ in range(23)]
    l = [1] * 6 + [0] * 17
    s = stats(v, l, nboot=500, nperm=2000)
    print(f"       판별력 {s['disc']:.3f}  구간 [{s['lo']:.2f}, {s['hi']:.2f}]  p={s['p']:.4f}")
    good = s["p"] > 0.05
    ok &= good; print(f"  {'통과' if good else '실패'}  p 가 커야 한다(유의하면 안 된다)")

    print("\n[5] 자세별 신뢰도 짝맞추기 — 아무 파일이나 집지 않는가")
    import tempfile
    with tempfile.TemporaryDirectory() as td:
        json.dump({"iptm": 0.11}, open(os.path.join(td, "confidence_x_model_0.json"), "w"))
        json.dump({"iptm": 0.99}, open(os.path.join(td, "confidence_x_model_1.json"), "w"))
        a = iptm_of(os.path.join(td, "x_model_0.cif"))
        b = iptm_of(os.path.join(td, "x_model_1.cif"))
        c = iptm_of(os.path.join(td, "x_model_7.cif"))       # 짝이 없는 자세
        chk("0번 자세", a, 0.11); chk("1번 자세", b, 0.99)
        chk("짝 없으면 NaN", c, float("nan"))

    print("\n" + ("전부 통과" if ok else "!! 실패한 항목이 있다"))
    return 0 if ok else 1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--selftest", action="store_true", help="통계 부분만 자료 없이 확인")
    ap.add_argument("--runs", default=os.path.join(AB, "runs_sites_guided_honest"))
    ap.add_argument("--diag", default="results/harm_diag.csv",
                    help="analyze_harm.py 결과 (타깃·변화량·후보수)")
    ap.add_argument("--targets-dir", default="targets")
    ap.add_argument("--arm", default="noconstraint", help="신호를 뽑을 팔")
    ap.add_argument("--cutoff", type=float, default=5.0)
    ap.add_argument("--link", type=float, default=0.5,
                    help="자세끼리 이만큼 겹치면 같은 자리로 본다")
    ap.add_argument("--harm", type=float, default=-0.09,
                    help="이보다 더 떨어지면 '크게 망가짐'")
    ap.add_argument("--fail-at", type=float, default=0.13,
                    help="무제약 DockQ 가 이 미만이면 '완전히 실패'로 본다. "
                         "이 구간에서 손해가 0 이었으므로, 이것을 정답 없이 맞히면 실전 규칙이 선다")
    ap.add_argument("--nboot", type=int, default=2000)
    ap.add_argument("--nperm", type=int, default=5000)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", default="results/gate_signals.csv")
    a = ap.parse_args()

    if a.selftest:
        raise SystemExit(selftest())

    import lib_pose_features as PF          # 자체 검증에는 필요 없으므로 여기서 불러온다
    import epitope_cluster as EC

    if not os.path.isdir(a.runs):
        raise SystemExit(f"!! 예측 폴더가 없다: {a.runs}\n"
                         f"   --runs 로 경로를 줄 것 (run_honest_guided.sh 의 OUTDIR)")
    if not os.path.exists(a.diag):
        raise SystemExit(f"!! {a.diag} 가 없다 — analyze_harm.py 를 먼저 돌릴 것")
    diag = list(csv.DictReader(open(a.diag)))
    need = {"target", "delta", "dq_no", "n_cand", "why"}
    miss = need - set(diag[0] if diag else {})
    if miss:
        raise SystemExit(f"!! {a.diag} 에 열이 없다: {sorted(miss)}\n"
                         f"   있는 열 = {list(diag[0]) if diag else []}")
    print(f"타깃 {len(diag)}종 · 팔={a.arm} · 예측폴더={a.runs}\n")

    rows, nofile, no_ip = [], [], 0
    for d in diag:
        tgt = d["target"]
        cifs = poses_of(a.runs, tgt, a.arm)
        if not cifs:
            nofile.append(tgt); continue
        cj = json.load(open(os.path.join(a.targets_dir, tgt, "chains.json")))

        eps, ip = [], []
        for c in cifs:
            ep, _ = EC.pred_epitope(cj, c, a.cutoff)
            if ep:
                eps.append(ep)
            v = iptm_of(c)
            (ip.append(v) if v == v else None)
        no_ip += len(cifs) - len(ip)
        if len(eps) < 2:
            nofile.append(tgt + "(자세 부족)"); continue

        jj = [EC.jac(eps[i], eps[j]) for i in range(len(eps)) for j in range(i + 1, len(eps))]
        jj = [v for v in jj if v == v]
        agree = st.mean(jj) if jj else float("nan")

        par = list(range(len(eps)))          # 겹치는 자세끼리 묶어 자리 수를 센다
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
            agree=round(agree, 4), nsite=nsite, n_cand=int(d["n_cand"]),
            delta=float(d["delta"]), dq_no=float(d["dq_no"]), why=d["why"]))

    if nofile:
        print(f"  ! 자세를 못 읽은 타깃 {len(nofile)}종: {', '.join(nofile[:8])}"
              f"{' …' if len(nofile) > 8 else ''}")
    if no_ip:
        print(f"  ! 짝맞는 confidence JSON 이 없던 자세 {no_ip}개 — 그 자세는 신뢰도에서 뺐다")
    if not rows:
        raise SystemExit("!! 아무 타깃도 못 읽었다 — --runs 경로와 --arm 이름을 확인할 것")
    nan_ip = sum(1 for r in rows if r["iptm"] != r["iptm"])
    if nan_ip:
        print(f"  ! 신뢰도를 전혀 못 구한 타깃 {nan_ip}종 — 그 신호는 그만큼 표본이 준다")

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
        print("\n" + "-" * 94)
        print(f"  {title}")
        print("-" * 94)
        lab = [1 if is_pos(r) else 0 for r in rows]
        npos = sum(lab)
        if npos in (0, len(rows)):
            print(f"  한쪽으로 몰려 잴 수 없다 (양성 {npos}/{len(rows)})"); return
        base = max(npos, len(rows) - npos) / len(rows)
        print(f"  양성 {npos}종 · 음성 {len(rows)-npos}종   {note}")
        print(f"  ⚠️ 기저율 {base:.0%} — 다수쪽만 찍어도 이만큼 맞는다. 정확도는 이보다 훨씬 커야 뜻이 있다")
        print(f"  {'신호':<26}{'판별력':>6}{'95% 구간':>16}{'p':>9}   {'방향':<12}{'문턱(정확도)'}")
        for key, name in SIG:
            vals = [r[key] for r in rows]
            s = stats(vals, lab, a.nboot, a.nperm, a.seed)
            cut = best_cut(vals, lab)
            ctxt = "-" if not cut else f"{cut[0]:g} {cut[1]} ({cut[2]:.0%})"
            if not s:
                print(f"  {name:<26}{'-':>6}{'-':>16}{'-':>9}   {'-':<12}{ctxt}")
                continue
            band = "[{:.2f}, {:.2f}]".format(s["lo"], s["hi"])
            print(f"  {name:<26}{s['disc']:>6.3f}{band:>16}"
                  f"{s['p']:>9.4f}   {s['dirn']:<12}{ctxt}")

    report("A. 크게 망가질 것인가 — 개입 여부를 정할 수 있나",
           lambda r: r["delta"] < a.harm, f"(변화 < {a.harm})")
    report(f"B. 무제약이 완전히 실패인가 — 개입해도 안전한 구간을 정답 없이 고를 수 있나",
           lambda r: r["dq_no"] < a.fail_at, f"(무제약 DockQ < {a.fail_at})")
    report("B'. 무제약이 합격선을 못 넘는가 (참고)",
           lambda r: r["dq_no"] < 0.23, "(무제약 DockQ < 0.23)")

    print("\n  ※ 판별력 0.5 = 동전던지기 · 1.0 = 완전히 갈림. 방향과 무관한 크기로 적었다")
    print("  ※ 95% 구간이 0.5 를 품으면 그 신호는 우연과 구별되지 않는다")
    print("  ※ 문턱과 정확도는 같은 자료에서 고른 값이라 낙관적이다. 새 세트로 확인해야 한다")

    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)
    print(f"\n→ {a.out}  ({len(rows)}종)")


if __name__ == "__main__":
    main()
