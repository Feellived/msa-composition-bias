#!/usr/bin/env python3
"""[핵심 통계] 조성을 바꾸면 '믿을 만한 다른 결합자리 후보'가 생기는가.

주장은 두 부분이고, 둘 다 이 스크립트가 재는 값이다.

  ① **믿을 만하다(재현된다)** — 같은 조성으로 다시 돌리면 같은 자리가 나온다
       → 조성 **안** 실행끼리의 결합자리 겹침이 높아야 한다.
  ② **다른 후보다** — 조성을 바꾸면 다른 자리가 나온다
       → 조성 **간** 겹침이 낮아야 한다.

  즉 핵심 숫자는 (조성 내 겹침) ÷ (조성 간 겹침). 1에 가까우면 조성은 아무 상관이 없고
  자리는 그냥 실행마다 흔들리는 것이며, 크면 **조성이 자리를 정한다**는 뜻이다.

  ⚠️ 이 비율만으로는 부족하다. 실행 수가 적으면 우연히 커질 수 있으므로 **뒤섞기 검정**을
  같이 한다: 실행에 붙은 조성 딱지를 무작위로 다시 붙여 같은 값을 여러 번 계산하고,
  실제 값이 그 분포의 어디쯤인지 본다(p = 뒤섞기가 실제만큼 커진 비율).

그다음 "**후보가 몇 개 만들어졌나**"를 센다. 조성마다 합의 자리(그 조성의 실행 과반에
나온 잔기)를 구하고, 서로 많이 겹치는 조성끼리 묶어 **서로 구별되는 자리 후보의 개수**와
각 후보가 진짜 결합자리를 얼마나 덮는지를 보고한다. 에피토프 비닝 관점에서 쓸 숫자가 이것이다.

단위 = 실행 1회(자세 5개 중 DockQ 최고를 그 실행의 대표로) — epitope_cluster.py와 동일.

사용(DockQ env):
  python dump_seedrep_full.py --data $DATA/compreps --only 8ulr_HL --csv-out results/compreps_8ulr_HL.csv
  python site_reproducibility.py --csv results/compreps_8ulr_HL.csv --data $DATA/compreps
  python site_reproducibility.py --csv ... --link 0.6 --nperm 5000
"""
import argparse, csv, glob, json, os, random, re
import statistics as st
from collections import Counter, defaultdict
import pose_features as PF
import epitope_cluster as EC          # pred_epitope · jac · consensus 재사용

SUCC = 0.49
RUN_RE = re.compile(r"^(?P<comp>.+?)_r(?P<rep>\d+)$")     # 예: seed3_r2 · full_r7


def split_run(name):
    """실행 폴더 이름 → (조성, 반복). 규칙에 안 맞으면 조성=이름 전체."""
    m = RUN_RE.match(name)
    return (m.group("comp"), m.group("rep")) if m else (name, "0")


def pairs_within(groups):
    """같은 조성 안의 모든 쌍."""
    out = []
    for g in groups.values():
        for i in range(len(g)):
            for j in range(i + 1, len(g)):
                out.append((g[i], g[j]))
    return out


def pairs_between(groups):
    """서로 다른 조성 사이의 모든 쌍."""
    ks = list(groups)
    out = []
    for a in range(len(ks)):
        for b in range(a + 1, len(ks)):
            for x in groups[ks[a]]:
                for y in groups[ks[b]]:
                    out.append((x, y))
    return out


def mean_jac(pairs):
    v = [EC.jac(x, y) for x, y in pairs]
    v = [z for z in v if z == z]
    return (st.mean(v), len(v)) if v else (float("nan"), 0)


def link_clusters(items, thr):
    """items = [(이름, 잔기집합)]. 겹침 ≥ thr 이면 같은 무리로 잇는다(단일연결)."""
    n = len(items)
    par = list(range(n))
    def find(x):
        while par[x] != x:
            par[x] = par[par[x]]; x = par[x]
        return x
    for i in range(n):
        for j in range(i + 1, n):
            if EC.jac(items[i][1], items[j][1]) >= thr:
                a, b = find(i), find(j)
                if a != b: par[a] = b
    cl = defaultdict(list)
    for i in range(n):
        cl[find(i)].append(i)
    return list(cl.values())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", required=True)
    ap.add_argument("--data", default=os.environ.get("DATA", "/mnt/data/admuser/msadepth") + "/compreps")
    ap.add_argument("--targets-dir", default="targets")
    ap.add_argument("--cutoff", type=float, default=5.0)
    ap.add_argument("--link", type=float, default=0.5, help="합의 자리끼리 이 값 이상 겹치면 같은 후보로 묶음")
    ap.add_argument("--nperm", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", default="")
    a = ap.parse_args()

    rows = list(csv.DictReader(open(a.csv)))
    if not rows:
        raise SystemExit("!! CSV가 비었음")
    tgt, model = rows[0]["target"], rows[0]["model"]
    depth = rows[0].get("depth", "")
    cj = json.load(open(os.path.join(a.targets_dir, tgt, "chains.json")))
    tr = PF.native_true(cj, os.path.join(a.targets_dir, tgt, "native.cif"), a.cutoff)
    if tr is None:
        raise SystemExit("!! native 결합자리 계산 실패")
    true = set(tr[0])

    # 실행별 대표 자세(DockQ 최고) — epitope_cluster.py와 같은 규칙
    best = {}
    for r in rows:
        try:
            q = float(r["dockq"])
        except Exception:
            continue
        s = r["seed"]
        if s not in best or q > best[s][0]:
            best[s] = (q, r["pose"])

    base = os.path.join(a.data, "seedrep_cand", model, tgt, depth)
    recs = []
    for s, (q, pose) in sorted(best.items()):
        hits = glob.glob(os.path.join(base, s, "results", "**", pose), recursive=True)
        if not hits:
            print(f"  ! {s}: 자세 파일 못 찾음({pose})"); continue
        ep, _ = EC.pred_epitope(cj, hits[0], a.cutoff)
        if not ep:
            print(f"  ! {s}: 결합자리 계산 실패"); continue
        comp, rep = split_run(s)
        recs.append(dict(run=s, comp=comp, rep=rep, dockq=q, ep=ep))
    if not recs:
        raise SystemExit("!! 계산된 실행이 없음")

    groups = defaultdict(list)
    for r in recs:
        groups[r["comp"]].append(r["ep"])
    multi = {k: v for k, v in groups.items() if len(v) >= 2}

    print(f"■ {tgt} · {model} · {depth}   실행 {len(recs)}개 · 조성 {len(groups)}가지"
          f" (반복 2회 이상인 조성 {len(multi)}가지)")
    print(f"  진짜 결합자리 잔기 {len(true)}개\n")
    if not multi:
        print("  ⚠️ 조성당 반복이 1회뿐이라 '조성 내 재현성'을 잴 수 없다. 반복을 늘려 다시 돌릴 것.")
        return

    win, nw = mean_jac(pairs_within(multi))
    btw, nb = mean_jac(pairs_between(groups))
    print("[조성이 자리를 정하는가 — 결합자리 겹침(자카드 0~1)]")
    print(f"  같은 조성 안  {win:.3f}  (쌍 {nw}개)   ← 높을수록 '같은 조성이면 같은 자리' = 믿을 만함")
    print(f"  다른 조성 사이 {btw:.3f}  (쌍 {nb}개)   ← 낮을수록 '조성을 바꾸면 다른 자리'")
    ratio = win / btw if btw and btw == btw and btw > 0 else float("nan")
    print(f"  비율 {ratio:.2f}배" if ratio == ratio else "  비율 계산 불가")

    # 뒤섞기 검정 — 조성 딱지를 무작위로 재배치해도 이만큼 벌어지나
    rng = random.Random(a.seed)
    sizes = [len(v) for v in groups.values()]
    allep = [r["ep"] for r in recs]
    hit = 0
    obs = win - (btw if btw == btw else 0)
    for _ in range(a.nperm):
        idx = list(range(len(allep))); rng.shuffle(idx)
        gg, k = {}, 0
        for gi, n in enumerate(sizes):
            gg[gi] = [allep[i] for i in idx[k:k + n]]; k += n
        mm = {q: v for q, v in gg.items() if len(v) >= 2}
        if not mm:
            continue
        w2, _ = mean_jac(pairs_within(mm))
        b2, _ = mean_jac(pairs_between(gg))
        if w2 == w2 and (w2 - (b2 if b2 == b2 else 0)) >= obs:
            hit += 1
    p = (hit + 1) / (a.nperm + 1)
    verdict = ("✅ 조성이 자리를 정한다" if p < 0.05 else
               "△ 경계" if p < 0.15 else "✗ 조성과 무관 — 자리는 실행마다 흔들릴 뿐")
    print(f"  뒤섞기 검정 p = {p:.4f}  ({a.nperm}회)   {verdict}\n")

    # 조성별 합의 자리 → 서로 구별되는 후보 몇 개인가
    cons = []
    for k, v in sorted(groups.items()):
        c = EC.consensus(v, 0.5)
        if c:
            cons.append((k, c))
    print(f"[만들어진 자리 후보 — 조성별 합의 자리를 겹침 {a.link} 이상이면 같은 것으로 묶음]")
    cl = link_clusters(cons, a.link)
    cl.sort(key=lambda ix: -len(ix))
    print(f"  조성 {len(cons)}가지 → 서로 구별되는 자리 후보 **{len(cl)}개**\n")
    print(f"  {'후보':5}{'조성수':7}{'잔기':6}{'진짜 자리 덮음':>14}{'예측 중 진짜':>13}   조성")
    print("  " + "-" * 74)
    out = []
    for ci, ix in enumerate(cl, 1):
        u = set()
        for i in ix:
            u |= cons[i][1]
        rec = len(u & true) / len(true) if true else float("nan")
        pre = len(u & true) / len(u) if u else float("nan")
        names = ",".join(cons[i][0] for i in ix)
        star = " ★" if rec >= 0.5 else ""
        print(f"  {ci:<5}{len(ix):<7}{len(u):<6}{rec:>13.2f}{pre:>13.2f}   {names}{star}")
        out.append(dict(target=tgt, model=model, depth=depth, cand=ci, n_comp=len(ix),
                        n_res=len(u), true_covered=round(rec, 4), precision=round(pre, 4),
                        comps=names))
    good = [o for o in out if o["true_covered"] >= 0.5]
    print(f"\n  → 후보 {len(cl)}개 중 진짜 자리를 절반 이상 덮는 것 **{len(good)}개**"
          + (f" (후보 {', '.join(str(o['cand']) for o in good)})" if good else ""))
    print("     ※ 후보를 여러 개 만들어 놓고 고르는 것이 목적이므로, '하나라도 맞으면' 성공이다.")

    path = a.out or f"results/site_repro_{tgt}.csv"
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["target", "model", "depth", "within", "between",
                                           "ratio", "perm_p", "n_cand", "cand", "n_comp",
                                           "n_res", "true_covered", "precision", "comps"])
        w.writeheader()
        for o in out:
            o.update(within=round(win, 4), between=round(btw, 4),
                     ratio=(round(ratio, 3) if ratio == ratio else ""),
                     perm_p=round(p, 5), n_cand=len(cl))
            w.writerow(o)
    print(f"\n→ {path}")


if __name__ == "__main__":
    main()
