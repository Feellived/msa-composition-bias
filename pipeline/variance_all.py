#!/usr/bin/env python3
"""[전수 분산] "조성이 결합 자리를 정한다"를 중심 사례 5종이 아니라 **타깃 전부**에서 확인한다.

지금까지 4.2절의 논증은 이랬다 — 조성을 바꿨을 때의 차이가 확산 모델의 무작위성 때문이라면,
같은 조성을 반복했을 때의 흔들림과 서로 다른 조성 사이의 흔들림이 비슷해야 한다. 그런데
실제로는 조성 간이 21.7배 컸다. 문제는 이 21.7배가 **5종에서만** 계산됐다는 것이다.

이 스크립트는 같은 계산을 compreps_<타깃>.csv (자세 단위 채점) 전부에 대해 돌린다.
GPU 를 쓰지 않는다. 이미 채점된 CSV 만 읽는다.

흔들림을 재는 기준을 둘 둔다. 둘 다 조성 간과 비교한다.

  (가) 조성 내 반복   같은 조성을 여러 번 돌렸을 때의 실행 간 차이
  (나) 원래 MSA 반복  원래 MSA 를 시드만 바꿔 반복했을 때의 차이
                      ← MSA 가 완전히 같으므로 **순수한 확산 모델의 무작위성**이다.
                        (가)보다 깨끗한 기준이라 이쪽이 더 강한 논증이 된다.

실행 하나의 대표값은 그 실행의 자세 5개 중 최고값이다(--stat 로 중앙값 선택 가능).
한 실행의 자세들은 서로 독립이 아니므로 자세를 표본으로 세지 않는다.

사용:
  python -u variance_all.py --results results
  python -u variance_all.py --results results --metric recall
"""
import argparse, csv, glob, math, os, random, re
import statistics as st

RUN_RE = re.compile(r"^(?P<comp>.+?)_r(?P<rep>\d+)$")
HEAD = ["target", "model", "depth", "seed", "pose", "dockq", "recall", "overrep", "n_contact"]


def load(path):
    """자세 단위 CSV → {조성: {반복: [자세값...]}} · 원래 MSA 는 조성 이름 'seedfull'."""
    rows = list(csv.DictReader(open(path)))
    if not rows or list(rows[0]) != HEAD:
        return None, "스키마가 다르다"
    dm = {(r["model"], r["depth"]) for r in rows}
    return rows, (f"모델·깊이가 여럿이다 {sorted(dm)}" if len(dm) > 1 else None)


def run_values(rows, metric, stat):
    """(조성, 반복) → 실행 대표값."""
    buf = {}
    for r in rows:
        try:
            v = float(r[metric])
        except (ValueError, KeyError):
            continue
        if math.isnan(v):
            continue
        m = RUN_RE.match(r["seed"])
        comp, rep = (m.group("comp"), m.group("rep")) if m else (r["seed"], "0")
        buf.setdefault(comp, {}).setdefault(rep, []).append(v)
    agg = max if stat == "max" else st.median
    return {c: {k: agg(v) for k, v in d.items()} for c, d in buf.items()}


def pooled_sd(groups):
    """여러 집단의 집단 내 표준편차를 자유도로 묶는다. 반복이 2회 이상인 집단만 쓴다."""
    ss, df = 0.0, 0
    for v in groups:
        if len(v) < 2:
            continue
        m = st.mean(v)
        ss += sum((x - m) ** 2 for x in v); df += len(v) - 1
    return (math.sqrt(ss / df) if df else None), df


def ratio_of(comp_runs):
    """조성 간 표준편차 / 조성 내 표준편차."""
    groups = [list(d.values()) for d in comp_runs.values()]
    within, df = pooled_sd(groups)
    means = [st.mean(g) for g in groups if g]
    between = st.stdev(means) if len(means) > 1 else None
    return within, between, df


def perm_p(comp_runs, nperm, seed=0):
    """조성 딱지를 섞어도 이만큼 갈리는가. 반복 구조(집단 크기)는 그대로 둔다."""
    within, between, _ = ratio_of(comp_runs)
    if not within or not between:
        return None, None
    obs = between / within
    flat, sizes = [], []
    for d in comp_runs.values():
        v = list(d.values()); flat += v; sizes.append(len(v))
    rng = random.Random(seed); hit = 0
    for _ in range(nperm):
        rng.shuffle(flat)
        i, groups = 0, []
        for n in sizes:
            groups.append(flat[i:i + n]); i += n
        w, _df = pooled_sd(groups)
        ms = [st.mean(g) for g in groups if g]
        b = st.stdev(ms) if len(ms) > 1 else 0.0
        if w and b / w >= obs:
            hit += 1
    return obs, (hit + 1) / (nperm + 1)


def binom_tail(k, n, p):
    """n 번 중 k 번 이상 나올 확률."""
    return sum(math.comb(n, i) * p ** i * (1 - p) ** (n - i) for i in range(k, n + 1))


def fisher(ps):
    """여러 타깃의 p 값을 하나로 묶는다 (Fisher 방법, 자유도 짝수라 닫힌 형태)."""
    X = -2 * sum(math.log(max(p, 1e-12)) for p in ps)
    k, x = len(ps), X / 2
    return math.exp(-x) * sum(x ** i / math.factorial(i) for i in range(k))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default="results", help="compreps_*.csv 가 있는 폴더")
    ap.add_argument("--metric", default="dockq", choices=["dockq", "recall"])
    ap.add_argument("--stat", default="max", choices=["max", "median"],
                    help="실행 하나의 대표값을 자세 5개 중 무엇으로 잡을지")
    ap.add_argument("--nperm", type=int, default=5000)
    ap.add_argument("--full-name", default="seedfull", help="원래 MSA 의 조성 이름")
    ap.add_argument("--out", default="results/variance_all.csv")
    a = ap.parse_args()

    files = sorted(f for f in glob.glob(os.path.join(a.results, "compreps_*.csv"))
                   if not f.endswith("compreps_summary.csv"))
    if not files:
        raise SystemExit(f"!! {a.results} 에 compreps_<타깃>.csv 가 없다")
    print(f"파일 {len(files)}개 · 지표 {a.metric} · 실행 대표값 = 자세 중 {a.stat}\n")

    rows, skipped = [], []
    for f in files:
        t = os.path.basename(f)[len("compreps_"):-len(".csv")]
        recs, warn = load(f)
        if recs is None:
            skipped.append((t, warn)); continue
        if warn:
            print(f"  ! {t}: {warn}")
        rv = run_values(recs, a.metric, a.stat)
        comps = {c: d for c, d in rv.items() if c != a.full_name}
        full = rv.get(a.full_name, {})
        if len(comps) < 2:
            skipped.append((t, f"조성이 {len(comps)}개뿐")); continue

        within, between, df = ratio_of(comps)
        obs, p = perm_p(comps, a.nperm)
        fullvals = list(full.values())
        sd_full = st.stdev(fullvals) if len(fullvals) > 1 else None
        rows.append(dict(
            target=t, n_comp=len(comps), n_run=sum(len(d) for d in comps.values()),
            n_full=len(fullvals),
            within=(round(within, 4) if within else None),
            between=(round(between, 4) if between else None),
            ratio=(round(between / within, 2) if within and between else None),
            perm_p=(round(p, 5) if p is not None else None),
            sd_full=(round(sd_full, 4) if sd_full else None),
            ratio_full=(round(between / sd_full, 2) if sd_full and between else None)))

    W = [("target", "타깃", 12), ("n_comp", "조성수", 7), ("n_run", "실행", 6),
         ("within", "조성내", 8), ("between", "조성간", 8), ("ratio", "배수", 7),
         ("perm_p", "p", 9), ("sd_full", "원래MSA내", 10), ("ratio_full", "배수2", 7)]
    print("=" * 86)
    print(f"  조성 간 흔들림 대 조성 내 흔들림 — 배수가 큰 순 ({a.metric})")
    print("=" * 86)
    print("  " + "".join(f"{h:>{w}}" for _, h, w in W))
    for r in sorted(rows, key=lambda r: -(r["ratio"] or 0)):
        print("  " + "".join(f"{('-' if r[c] is None else r[c])!s:>{w}}" for c, _, w in W))

    rt = [r["ratio"] for r in rows if r["ratio"]]
    ps = [r["perm_p"] for r in rows if r["perm_p"] is not None]
    print("\n" + "-" * 74)
    print("  전수 요약")
    print("-" * 74)
    print(f"  타깃 {len(rows)}종 · 실행 합계 {sum(r['n_run'] for r in rows)}회")

    # ⚠️ 배수를 1 과 비교하면 안 된다. 조성당 반복이 4회뿐이라 아무 효과가 없어도
    #    조성 평균들은 σ/√4 만큼 흔들린다 — 즉 귀무가설에서 배수의 기댓값이 1 이 아니라
    #    0.5 근처다. 옳은 기준선은 딱지를 섞어 만든 순열 분포이고, 그것이 perm_p 다.
    sig = [r for r in rows if r["perm_p"] is not None and r["perm_p"] < 0.05]
    if ps:
        print(f"\n  ⭐ 타깃별 순열검정 — 조성 딱지를 섞어도 이만큼 갈리는가")
        print(f"    유의한 타깃 {len(sig)}/{len(ps)}종   (우연이라면 {0.05*len(ps):.1f}종 기대)")
        print(f"    이항검정 p = {binom_tail(len(sig), len(ps), 0.05):.3g}")
        print(f"    Fisher 로 30종을 하나로 묶은 p = {fisher(ps):.3g}")
        if sig:
            print(f"    {', '.join(r['target'] for r in sorted(sig, key=lambda r: r['perm_p']))}")
    if rt:
        print(f"\n  참고용 배수(조성 간 / 조성 내) 중앙값 {st.median(rt):.2f} "
              f"· 최소 {min(rt):.2f} · 최대 {max(rt):.2f}")
        print("    이 배수는 타깃끼리 비교하는 용도이지 1 을 기준으로 읽는 값이 아니다.")
    if skipped:
        print(f"\n  ! 제외 {len(skipped)}종")
        for t, why in skipped:
            print(f"    {t:<12} {why}")

    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    print(f"\n→ {a.out}  ({len(rows)}종)")


if __name__ == "__main__":
    main()
