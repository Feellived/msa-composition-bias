#!/usr/bin/env python3
"""[③ 준비] 시드를 몇 번 돌리면 찾을 자리를 다 찾는가 — 포화 곡선.

연구실 경험칙은 **시드 12면 후보가 포화된다**는 것이다. 그 말이 우리 자료에서도
맞는지 먼저 확인한다. 맞으면 시드-대조 실험을 12회로 설계할 근거가 생기고,
안 맞으면 설계를 다시 잡아야 한다.

쓰는 자료는 **이미 있는 것**이다 — 원래 MSA 를 시드만 바꿔 여러 번 돌린 대조군
(fullmsa_ctl). GPU 를 쓰지 않는다.

재는 값 두 가지. 실행 수 n 을 1 부터 늘려 가며 무작위로 n 개를 뽑아 평균낸다.
  ① 자리 수    — 그 n 개가 찾아낸 **서로 구별되는 결합 자리**가 몇 개인가
  ② 정답 도달률 — 그 n 개 중 정답 자리를 충분히 덮는 실행이 **하나라도** 있을 확률
                  (여기서만 정답을 쓴다. 채점용이다)

①이 평평해지는 지점이 "더 돌려도 새 자리가 안 나온다"는 뜻이고, 그것이 포화다.

사용 (conda activate boltz · pipeline/ 에서):
  python -u seed_saturation.py --selftest        # 계산 부분만 자료 없이 확인
  python -u seed_saturation.py
  python -u seed_saturation.py --only "8ulr_HL 8k3k_D" --model protenix
"""
import argparse, csv, glob, json, os, random
import statistics as st


# ── 계산 ────────────────────────────────────────────────────────────────────
def n_sites(idxs, eps, jac, link):
    """뽑힌 실행들이 찾아낸 서로 구별되는 자리 수(겹침 link 이상이면 같은 자리)."""
    par = {i: i for i in idxs}
    def find(x):
        while par[x] != x:
            par[x] = par[par[x]]; x = par[x]
        return x
    for a_ in range(len(idxs)):
        for b_ in range(a_ + 1, len(idxs)):
            i, j = idxs[a_], idxs[b_]
            if jac(eps[i], eps[j]) >= link:
                ra, rb = find(i), find(j)
                if ra != rb:
                    par[ra] = rb
    return len({find(i) for i in idxs})


def curve(eps, hit, jac, link, nboot, seed=0):
    """n = 1..N 에서 (평균 자리 수, 정답 도달률). hit[i] = i번 실행이 정답을 덮었나."""
    rng = random.Random(seed)
    N = len(eps)
    out = []
    for n in range(1, N + 1):
        sites, reach = [], 0
        for _ in range(nboot):
            pick = rng.sample(range(N), n)
            sites.append(n_sites(pick, eps, jac, link))
            reach += any(hit[i] for i in pick)
        out.append((n, st.mean(sites), reach / nboot))
    return out


def plateau(vals, frac=0.95):
    """마지막 값의 frac 에 처음 닿는 n. 곡선이 평평해지기 시작하는 지점으로 읽는다."""
    if not vals:
        return None
    top = vals[-1]
    if top <= 0:
        return None
    for n, v in enumerate(vals, 1):
        if v >= frac * top:
            return n
    return len(vals)


def bar(v, top, width=34):
    return "█" * max(0, min(width, round(width * v / top))) if top > 0 else ""


# ── 자체 검증 ───────────────────────────────────────────────────────────────
def selftest():
    ok = True

    def jac(a, b):
        u = len(a | b)
        return len(a & b) / u if u else 0.0

    def chk(name, got, want):
        nonlocal ok
        good = got == want
        ok &= good
        print(f"  {'통과' if good else '실패':<4} {name:<40} 나온값 {got}  기대 {want}")

    print("[1] 자리 수 세기")
    A, B, C = {1, 2, 3}, {1, 2, 3}, {9, 8, 7}          # A·B 는 같은 자리, C 는 다른 자리
    eps = [A, B, C]
    chk("같은 자리 둘 + 다른 자리 하나", n_sites([0, 1, 2], eps, jac, 0.5), 2)
    chk("같은 자리끼리만", n_sites([0, 1], eps, jac, 0.5), 1)
    chk("하나만", n_sites([2], eps, jac, 0.5), 1)

    print("\n[2] 포화 지점 찾기")
    chk("3에서 천장에 닿는 곡선", plateau([1, 2, 3, 3, 3]), 3)
    chk("끝까지 오르는 곡선", plateau([1, 2, 3, 4, 5]), 5)

    print("\n[3] 곡선 — 자리가 하나뿐이면 아무리 늘려도 1")
    eps = [{1, 2, 3}] * 8
    c = curve(eps, [False] * 8, jac, 0.5, 200)
    good = all(abs(m - 1.0) < 1e-9 for _, m, _ in c)
    ok &= good; print(f"  {'통과' if good else '실패'}  전 구간 1.0 (마지막 {c[-1][1]:.2f})")

    print("\n[4] 곡선 — 자리가 전부 다르면 n 을 그대로 따라간다")
    eps = [{i} for i in range(8)]
    c = curve(eps, [False] * 8, jac, 0.5, 200)
    good = all(abs(m - n) < 1e-9 for n, m, _ in c)
    ok &= good; print(f"  {'통과' if good else '실패'}  n 과 같음 (마지막 {c[-1][1]:.2f} / 8)")

    print("\n[5] 정답 도달률 — 8개 중 1개만 정답이면 n=8 에서 1.0, n=1 에서 1/8")
    hit = [True] + [False] * 7
    c = curve([{i} for i in range(8)], hit, jac, 0.5, 4000, seed=1)
    good = c[-1][2] == 1.0 and abs(c[0][2] - 0.125) < 0.03
    ok &= good; print(f"  {'통과' if good else '실패'}  n=1 {c[0][2]:.3f} (기대 0.125) · n=8 {c[-1][2]:.3f}")

    print("\n" + ("전부 통과" if ok else "!! 실패한 항목이 있다"))
    return 0 if ok else 1


# ── 본 계산 ─────────────────────────────────────────────────────────────────
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--data", default=os.environ.get("DATA", "") + "/fullmsa_ctl",
                    help="원래 MSA 시드 대조군 폴더 ($DATA/fullmsa_ctl)")
    ap.add_argument("--model", default="protenix")
    ap.add_argument("--depth-dir", default="full_n40", help="타깃 아래 깊이 폴더 이름")
    ap.add_argument("--targets-dir", default="targets")
    ap.add_argument("--only", default="", help="공백으로 구분한 타깃 목록")
    ap.add_argument("--cutoff", type=float, default=5.0)
    ap.add_argument("--link", type=float, default=0.5,
                    help="실행끼리 이만큼 겹치면 같은 자리로 본다")
    ap.add_argument("--pose-frac", type=float, default=0.5,
                    help="한 실행 안에서 자세 몇 비율에 나와야 그 실행의 자리로 볼까")
    ap.add_argument("--recall-ok", type=float, default=0.4, help="정답 도달 기준")
    ap.add_argument("--nboot", type=int, default=400)
    ap.add_argument("--mark", type=int, default=12, help="눈여겨볼 실행 수(연구실 경험칙)")
    ap.add_argument("--out", default="results/seed_saturation.csv")
    a = ap.parse_args()

    if a.selftest:
        raise SystemExit(selftest())

    import pose_features as PF
    import epitope_cluster as EC

    if not os.path.isdir(a.data):
        raise SystemExit(f"!! 대조군 폴더가 없다: {a.data}\n"
                         f"   $DATA 가 비어 있지 않은지 확인하고 --data 로 직접 줄 것\n"
                         f"   예: --data /mnt/data/admuser/msadepth/fullmsa_ctl")
    root = os.path.join(a.data, a.model)
    if not os.path.isdir(root):
        raise SystemExit(f"!! 모델 폴더가 없다: {root}\n"
                         f"   있는 것 = {sorted(os.listdir(a.data))[:8]}")

    want = a.only.split() if a.only else None
    tgts = sorted(t for t in os.listdir(root) if os.path.isdir(os.path.join(root, t)))
    if want:
        tgts = [t for t in tgts if t in want]
    if not tgts:
        raise SystemExit(f"!! 타깃이 없다 — {root}")
    print(f"모델={a.model} · 타깃 {len(tgts)}종 · 대조군 폴더={root}\n")

    rows = []
    for tgt in tgts:
        base = os.path.join(root, tgt, a.depth_dir)
        if not os.path.isdir(base):
            sub = [d for d in os.listdir(os.path.join(root, tgt))
                   if os.path.isdir(os.path.join(root, tgt, d))]
            print(f"  ! {tgt}: {a.depth_dir} 없음 (있는 것 = {sub[:4]}) — 건너뜀"); continue
        runs = sorted(d for d in os.listdir(base) if os.path.isdir(os.path.join(base, d)))
        if len(runs) < 4:
            print(f"  ! {tgt}: 실행이 {len(runs)}개뿐 — 건너뜀"); continue

        cj = json.load(open(os.path.join(a.targets_dir, tgt, "chains.json")))
        tr = PF.native_true(cj, os.path.join(a.targets_dir, tgt, "native.cif"), a.cutoff)
        if tr is None:
            print(f"  ! {tgt}: 정답 결합자리 계산 실패 — 건너뜀"); continue
        true = set(tr[0])

        eps, hit, npose = [], [], []
        for r in runs:
            cifs = sorted(glob.glob(os.path.join(base, r, "results", "**", "*.cif"),
                                    recursive=True))
            ps = []
            for c in cifs:
                ep, _ = EC.pred_epitope(cj, c, a.cutoff)
                if ep:
                    ps.append(ep)
            if not ps:
                continue
            e = ps[0] if len(ps) == 1 else EC.consensus(ps, a.pose_frac)
            if not e:                       # 자세가 서로 너무 다르면 합의가 빈다
                e = set().union(*ps)
            eps.append(e)
            hit.append(len(e & true) / len(true) >= a.recall_ok)
            npose.append(len(ps))
        if len(eps) < 4:
            print(f"  ! {tgt}: 자리를 읽은 실행이 {len(eps)}개뿐 — 건너뜀"); continue

        c = curve(eps, hit, EC.jac, a.link, a.nboot)
        sat_s = plateau([m for _, m, _ in c])
        sat_r = plateau([r for _, _, r in c])
        N = len(eps)
        print("=" * 84)
        print(f"  {tgt}   실행 {N}개 · 자세 중앙값 {st.median(npose):.0f}개 · "
              f"정답 도달 실행 {sum(hit)}/{N}")
        print("=" * 84)
        top = max(m for _, m, _ in c)
        print(f"  {'실행수':>5}{'자리수':>7}{'정답도달률':>10}   자리 수")
        for n, m, r in c:
            mk = "  ← 경험칙" if n == a.mark else ""
            print(f"  {n:>5}{m:>7.2f}{r:>10.2f}   {bar(m, top)}{mk}")
        print(f"\n  자리 수가 천장의 95% 에 닿는 실행 수 = {sat_s}"
              f"   ·   정답 도달률은 {sat_r}")
        if sat_s and sat_s <= a.mark:
            print(f"  → {a.mark}회면 새 자리가 거의 안 나온다. 경험칙과 맞는다")
        else:
            print(f"  → {a.mark}회로는 아직 새 자리가 나온다. 설계를 다시 볼 것")

        for n, m, r in c:
            rows.append(dict(target=tgt, model=a.model, n_run=n,
                             n_site=round(m, 3), reach=round(r, 3),
                             n_total=N, n_hit=sum(hit),
                             sat_site=sat_s, sat_reach=sat_r))

    if not rows:
        raise SystemExit("!! 아무 타깃도 못 읽었다 — --data 와 --model 을 확인할 것")

    print("\n" + "-" * 72)
    print(f"  타깃별 포화 지점 (자리 수가 천장의 95% 에 닿는 실행 수)")
    print("-" * 72)
    per = {}
    for r in rows:
        per[r["target"]] = (r["sat_site"], r["sat_reach"], r["n_total"])
    for t, (s, rr, n) in sorted(per.items()):
        flag = "경험칙 안" if s and s <= a.mark else "더 필요"
        print(f"  {t:<12} 자리 {s!s:>4} · 정답도달 {rr!s:>4} · 총 실행 {n:>3}   {flag}")
    vals = [s for s, _, _ in per.values() if s]
    if vals:
        print(f"\n  중앙값 {st.median(vals):.0f}회 · 최대 {max(vals)}회"
              f"   (경험칙 {a.mark}회)")

    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)
    print(f"\n→ {a.out}  ({len(per)}종)")


if __name__ == "__main__":
    main()
