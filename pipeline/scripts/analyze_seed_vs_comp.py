#!/usr/bin/env python3
"""[③] 조성을 바꾸는 것이 시드를 바꾸는 것과 다른가 — 같은 실행 수에서 비교.

이 연구의 첫 단계(조성 재추첨)가 존재할 이유 자체를 재는 실험이다. 반박은 하나다.
"그냥 여러 번 돌린 것과 뭐가 다르냐."

같은 실행 수 n 에서 세 가지를 비교한다.
  A 조성군     서로 다른 조성 n 개에서 각 1회        ← 우리 방식
  B 시드군     **조성 하나를 고정**하고 시드만 n 회   ← 결정적 대조 (깊이도 조성도 같다)
  C 원래MSA군  원래 MSA 로 시드만 n 회               ← 실무에서의 기본값

B 가 핵심이다. C 는 깊이가 달라 "서열 수가 줄어서 그런 것 아니냐"를 막지 못한다.

재는 값은 analyze_seed_saturation.py 와 같다.
  ① 자리 수    — 그 n 회가 찾아낸 서로 구별되는 결합 자리
  ② 정답 도달률 — 그 n 회 중 정답 자리를 충분히 덮는 실행이 하나라도 있을 확률

⚠️ 시드군에서 어느 조성을 고정하느냐가 결과를 좌우한다. 좋은 조성을 고정하면 시드군이 이기고
   나쁜 조성을 고정하면 진다. 그래서 **한 조성으로 고정하지 않고 모든 조성에 대해 평균**낸다
   (매번 조성 하나를 무작위로 골라 그 안에서 n 회를 뽑는다). 이것이 "조성은 그대로 두고 시드만
   바꿨을 때 기대되는 성능"이다. `--fix-comp` 로 하나만 볼 수도 있지만 그러면 그 조성 하나의
   운에 결과가 걸린다.
⚠️ 조성마다 반복이 4회뿐이면 n 은 4 까지만 볼 수 있다. 12 까지 늘리려면 GPU 가 필요하다.

GPU 를 쓰지 않는다. 이미 있는 실행만 읽는다.

사용 (conda activate boltz · pipeline/ 에서):
  python -u analyze_seed_vs_comp.py --selftest
  export DATA=/mnt/data/msadepth && python -u analyze_seed_vs_comp.py
  python -u analyze_seed_vs_comp.py --only "8k3k_D 8k46_I 8tp5_HL"
"""
import argparse, csv, glob, json, os, random, re
import statistics as st

RUN_RE = re.compile(r"^(?P<comp>.+?)_r(?P<rep>\d+)$")     # seed3_r2 · seedfull_r7
FULL = "seedfull"


def split_run(name):
    m = RUN_RE.match(name)
    return (m.group("comp"), m.group("rep")) if m else (name, "0")


def n_sites(picks, jac, link):
    """뽑힌 자리 집합들이 몇 군데로 나뉘나(겹침 link 이상이면 같은 자리)."""
    par = list(range(len(picks)))
    def find(x):
        while par[x] != x:
            par[x] = par[par[x]]; x = par[x]
        return x
    for i in range(len(picks)):
        for j in range(i + 1, len(picks)):
            if jac(picks[i], picks[j]) >= link:
                a, b = find(i), find(j)
                if a != b:
                    par[a] = b
    return len({find(i) for i in range(len(picks))})


def arm_curve(groups, mode, jac, link, nmax, nboot, seed=0):
    """mode='across' 서로 다른 조성에서 하나씩 · 'within' 한 조성 안에서 여러 개.
    groups = {조성: [(자리집합, 정답덮었나), ...]}"""
    rng = random.Random(seed)
    out = []
    keys = sorted(groups)
    for n in range(1, nmax + 1):
        sites, reach, ok = [], 0, 0
        for _ in range(nboot):
            if mode == "across":
                if len(keys) < n:
                    continue
                sel = [rng.choice(groups[k]) for k in rng.sample(keys, n)]
            else:
                cand = [k for k in keys if len(groups[k]) >= n]
                if not cand:
                    continue
                k = rng.choice(cand)
                sel = rng.sample(groups[k], n)
            ok += 1
            sites.append(n_sites([e for e, _ in sel], jac, link))
            reach += any(h for _, h in sel)
        out.append((n, st.mean(sites) if sites else float("nan"),
                    reach / ok if ok else float("nan"), ok))
    return out


def selftest():
    ok = True

    def jac(a, b):
        u = len(a | b)
        return len(a & b) / u if u else 0.0

    def chk(name, cond):
        nonlocal ok
        ok &= cond
        print(f"  {'통과' if cond else '실패':<4} {name}")

    print("[1] 조성이 자리를 정하는 경우 — 조성군만 자리가 늘어야 한다")
    # 조성 4개, 각 조성은 자기 자리만 낸다. 조성 안에서는 늘 같은 자리.
    g = {f"seed{i}": [({i * 10 + k for k in range(3)}, i == 0) for _ in range(4)]
         for i in range(4)}
    ac = arm_curve(g, "across", jac, 0.5, 4, 300)
    wi = arm_curve(g, "within", jac, 0.5, 4, 300)
    chk(f"조성군 n=4 자리 {ac[-1][1]:.2f} (기대 4)", abs(ac[-1][1] - 4) < 1e-9)
    chk(f"시드군 n=4 자리 {wi[-1][1]:.2f} (기대 1)", abs(wi[-1][1] - 1) < 1e-9)
    chk(f"조성군 정답도달 {ac[-1][2]:.2f} > 시드군 {wi[-1][2]:.2f}", ac[-1][2] > wi[-1][2])

    print("\n[2] 조성이 아무 상관 없는 경우 — 두 팔이 같아야 한다")
    # 모든 실행이 무작위로 네 자리 중 하나. 조성 딱지는 의미가 없다.
    rng = random.Random(3)
    g = {f"seed{i}": [({rng.randrange(4) * 10 + k for k in range(3)}, False)
                      for _ in range(4)] for i in range(4)}
    ac = arm_curve(g, "across", jac, 0.5, 4, 2000, seed=1)
    wi = arm_curve(g, "within", jac, 0.5, 4, 2000, seed=1)
    d = abs(ac[-1][1] - wi[-1][1])
    chk(f"n=4 자리 차이 {d:.2f} 가 작다 (조성군 {ac[-1][1]:.2f} · 시드군 {wi[-1][1]:.2f})", d < 0.35)

    print("\n" + ("전부 통과" if ok else "!! 실패한 항목이 있다"))
    return 0 if ok else 1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--data", default=os.environ.get("DATA", "") + "/compreps/seedrep_cand")
    ap.add_argument("--model", default="protenix")
    ap.add_argument("--depth-dir", default="", help="비우면 타깃마다 자동으로 찾는다")
    ap.add_argument("--targets-dir", default="targets")
    ap.add_argument("--only", default="")
    ap.add_argument("--fix-comp", default="",
                    help="시드군에서 조성 하나만 보고 싶을 때. 비우면 **모든 조성에 대해 평균**낸다")
    ap.add_argument("--cutoff", type=float, default=5.0)
    ap.add_argument("--link", type=float, default=0.5)
    ap.add_argument("--pose-frac", type=float, default=0.5)
    ap.add_argument("--recall-ok", type=float, default=0.4)
    ap.add_argument("--nboot", type=int, default=600)
    ap.add_argument("--out", default="results/seed_vs_comp.csv")
    a = ap.parse_args()

    if a.selftest:
        raise SystemExit(selftest())

    import lib_pose_features as PF
    import epitope_cluster as EC

    root = os.path.join(a.data, a.model)
    if not os.path.isdir(root):
        raise SystemExit(f"!! 폴더가 없다: {root}\n   $DATA 확인 후 --data 로 줄 것")
    want = a.only.split() if a.only else None
    tgts = sorted(t for t in os.listdir(root) if os.path.isdir(os.path.join(root, t)))
    if want:
        tgts = [t for t in tgts if t in want]
    print(f"모델={a.model} · 타깃 {len(tgts)}종 · 폴더={root}\n")

    rows, summ = [], []
    for tgt in tgts:
        tdir = os.path.join(root, tgt)
        subs = sorted(d for d in os.listdir(tdir) if os.path.isdir(os.path.join(tdir, d)))
        depth = a.depth_dir or (subs[0] if len(subs) == 1 else "")
        if not depth:
            print(f"  ! {tgt}: 깊이 폴더가 {len(subs)}개 {subs[:4]} — --depth-dir 필요"); continue
        base = os.path.join(tdir, depth)

        cj = json.load(open(os.path.join(a.targets_dir, tgt, "chains.json")))
        tr = PF.native_true(cj, os.path.join(a.targets_dir, tgt, "native.cif"), a.cutoff)
        if tr is None:
            print(f"  ! {tgt}: 정답 결합자리 계산 실패 — 건너뜀"); continue
        true = set(tr[0])

        groups = {}
        for r in sorted(os.listdir(base)):
            if not os.path.isdir(os.path.join(base, r)):
                continue
            comp, _ = split_run(r)
            ps = []
            for c in sorted(glob.glob(os.path.join(base, r, "results", "**", "*.cif"),
                                      recursive=True)):
                ep, _ = EC.pred_epitope(cj, c, a.cutoff)
                if ep:
                    ps.append(ep)
            if not ps:
                continue
            e = ps[0] if len(ps) == 1 else EC.consensus(ps, a.pose_frac)
            if not e:
                e = set().union(*ps)
            groups.setdefault(comp, []).append((e, len(e & true) / len(true) >= a.recall_ok))

        comps = {k: v for k, v in groups.items() if k != FULL}
        full = groups.get(FULL, [])
        if len(comps) < 3:
            print(f"  ! {tgt}: 조성이 {len(comps)}개뿐 — 건너뜀"); continue
        fix = a.fix_comp
        if fix and fix not in comps:
            print(f"  ! {tgt}: 고정할 조성 {fix} 가 없다 (있는 것 {sorted(comps)[:5]}) — 건너뜀")
            continue
        within = {fix: comps[fix]} if fix else comps      # 비우면 조성 전체에 대해 평균
        nmax = min(len(comps), max(len(v) for v in comps.values()))
        if nmax < 2:
            print(f"  ! {tgt}: 비교할 실행 수가 부족 — 건너뜀"); continue

        ac = arm_curve(comps, "across", EC.jac, a.link, nmax, a.nboot)
        wi = arm_curve(within, "within", EC.jac, a.link, nmax, a.nboot)
        fu = arm_curve({FULL: full}, "within", EC.jac, a.link, nmax, a.nboot) if len(full) >= 2 else None

        print("=" * 92)
        print(f"  {tgt} ({depth})   조성 {len(comps)}개 · 조성당 반복 "
              f"{sorted({len(v) for v in comps.values()})} · 원래MSA {len(full)}회 · "
              f"시드군={'조성 ' + fix + ' 고정' if fix else '조성 전체 평균'}")
        print("=" * 92)
        print(f"  {'실행수':>5} | {'조성군':^18} | {'시드군(조성고정)':^18} | {'원래MSA군':^18}")
        print(f"  {'':>5} | {'자리':>7}{'정답도달':>11} | {'자리':>7}{'정답도달':>11} | {'자리':>7}{'정답도달':>11}")
        for i in range(nmax):
            n, m1, r1, _ = ac[i]
            _, m2, r2, _ = wi[i]
            f3 = fu[i] if fu else (0, float("nan"), float("nan"), 0)
            print(f"  {n:>5} | {m1:>7.2f}{r1:>11.2f} | {m2:>7.2f}{r2:>11.2f} | "
                  f"{f3[1]:>7.2f}{f3[2]:>11.2f}")
            rows.append(dict(target=tgt, n_run=n,
                             comp_site=round(m1, 3), comp_reach=round(r1, 3),
                             seed_site=round(m2, 3), seed_reach=round(r2, 3),
                             full_site=round(f3[1], 3), full_reach=round(f3[2], 3),
                             n_comp=len(comps), fix_comp=(fix or "(전체평균)")))
        d_site = ac[-1][1] - wi[-1][1]
        d_reach = ac[-1][2] - wi[-1][2]
        print(f"\n  n={nmax} 에서 조성군 − 시드군:  자리 {d_site:+.2f} · 정답도달 {d_reach:+.2f}")
        summ.append((tgt, nmax, d_site, d_reach))

    if not rows:
        raise SystemExit("!! 아무 타깃도 못 읽었다 — 경로를 확인할 것")

    print("\n" + "=" * 78)
    print("  타깃별 요약 — 같은 실행 수에서 조성군이 시드군보다 얼마나 나은가")
    print("=" * 78)
    print(f"  {'타깃':<12}{'실행수':>6}{'자리 차이':>10}{'정답도달 차이':>14}")
    for t, n, ds, dr in sorted(summ, key=lambda x: -x[2]):
        print(f"  {t:<12}{n:>6}{ds:>+10.2f}{dr:>+14.2f}")
    ds = [x[2] for x in summ]
    dr = [x[3] for x in summ if x[3] == x[3]]
    win = sum(1 for v in ds if v > 0)
    print(f"\n  조성군이 자리를 더 많이 찾은 타깃 {win}/{len(ds)}종 "
          f"· 자리 차이 중앙값 {st.median(ds):+.2f}")
    if dr:
        w2 = sum(1 for v in dr if v > 0)
        print(f"  정답 도달이 더 높은 타깃 {w2}/{len(dr)}종 · 중앙값 {st.median(dr):+.2f}")
    print("\n  ※ 반복이 4회뿐이면 n=4 까지만 볼 수 있다. 더 늘리려면 GPU 로 실행을 추가해야 한다")

    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)
    print(f"→ {a.out}  ({len(summ)}종)")


if __name__ == "__main__":
    main()
