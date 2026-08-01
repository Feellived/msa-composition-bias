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

⭐ **정답을 안 보는 판(2026-08-01)**. 단위는 여전히 **실행 1회**지만, 그 실행의 결합자리를
"자세 5개 중 DockQ 최고 하나"가 아니라 **자세 5개의 합의**(절반 이상 자세에 나온 잔기)로 잡는다.
  · 왜 — DockQ는 **정답 구조가 있어야** 계산된다. 그걸로 대표를 고르면 이 파이프라인은 실전에서
    돌 수 없고, "그 자리는 정답을 보고 고른 것 아니냐"는 반박에 답할 수 없다.
  · 표본 단위는 안 바뀐다 — 한 실행의 자세 5개는 서로 상관되어 있어 독립 표본이 아니므로
    **자세를 따로 세지 않고 합쳐서 실행 하나의 자리**로 만든다.
  · 옛 방식이 필요하면 `--legacy-best-pose` (비교용으로만).

후보를 조립할 때도 **합집합이 아니라 투표**를 쓴다(`--merge-frac`). 묶인 조성 중 정해진 비율
이상이 지목한 잔기만 넣는다. 합집합은 한 조성에만 나온 잔기까지 다 넣어 후보가 정답의 2.5배로
부풀었다(62잔기 대 25잔기).

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
    ap.add_argument("--pose-frac", type=float, default=0.5,
                    help="한 실행 안에서 자세 몇 비율에 나와야 그 실행의 결합자리로 볼까 (기본 0.5)")
    ap.add_argument("--merge-frac", type=float, default=0.75,
                    help="묶인 조성 중 몇 비율이 지목해야 후보 잔기로 넣을까. 0 이면 옛 방식(합집합)")
    ap.add_argument("--legacy-best-pose", action="store_true",
                    help="⚠️ 옛 방식 — 실행마다 DockQ(정답) 최고 자세를 대표로. 새 판과 맞대볼 때만.")
    ap.add_argument("--nperm", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", default="")
    ap.add_argument("--dump-sites", nargs="?", const="AUTO", default="",
                    help="후보 자리의 '잔기 목록'을 JSON으로 쓴다(경로 생략 시 results/sites_<타깃>.json). "
                         "CSV에는 잔기 수만 있어 유도 재도킹(guided) 입력을 만들 수 없다.")
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

    # 실행 하나의 결합자리 = 그 실행 자세들의 합의 (정답을 안 본다)
    base = os.path.join(a.data, "seedrep_cand", model, tgt, depth)
    by_run = defaultdict(list)
    for r in rows:
        by_run[r["seed"]].append(r)

    def dq(r):
        try:
            return float(r["dockq"])
        except Exception:
            return float("nan")

    recs, npose = [], []
    for s, rr in sorted(by_run.items()):
        # ⚠️ 옛 방식: DockQ(정답)로 대표 자세를 고른다. 새 판과 맞대볼 때만 쓸 것.
        use = [max(rr, key=lambda r: (dq(r) if dq(r) == dq(r) else -1))] if a.legacy_best_pose else rr
        eps = []
        for r in use:
            hits = glob.glob(os.path.join(base, s, "results", "**", r["pose"]), recursive=True)
            if not hits:
                continue
            ep, _ = EC.pred_epitope(cj, hits[0], a.cutoff)
            if ep:
                eps.append(ep)
        if not eps:
            print(f"  ! {s}: 자세를 못 읽음"); continue
        ep = eps[0] if len(eps) == 1 else EC.consensus(eps, a.pose_frac)
        if not ep:                      # 자세들이 서로 너무 달라 합의가 비면 합집합으로 물러선다
            ep = set().union(*eps)
        comp, rep = split_run(s)
        recs.append(dict(run=s, comp=comp, rep=rep, ep=ep,
                         dockq=max((dq(r) for r in rr if dq(r) == dq(r)), default=float("nan"))))
        npose.append(len(eps))
    if not recs:
        raise SystemExit("!! 계산된 실행이 없음")

    groups = defaultdict(list)
    for r in recs:
        groups[r["comp"]].append(r["ep"])
    multi = {k: v for k, v in groups.items() if len(v) >= 2}

    print(f"■ {tgt} · {model} · {depth}   실행 {len(recs)}개 · 조성 {len(groups)}가지"
          f" (반복 2회 이상인 조성 {len(multi)}가지)")
    print(f"  진짜 결합자리 잔기 {len(true)}개")
    if a.legacy_best_pose:
        print("  ⚠️ 옛 방식 — 실행마다 DockQ(정답) 최고 자세를 대표로 썼다. 실전 파이프라인이 아니다.")
    else:
        print(f"  자리 만드는 법 = 실행당 자세 {st.mean(npose):.1f}개의 합의(비율 {a.pose_frac}) "
              f"· 후보 조립 = 투표(비율 {a.merge_frac})   ← 정답을 보지 않는다")
    print()
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
    sites = []
    for ci, ix in enumerate(cl, 1):
        # 합집합이 아니라 투표 — 묶인 조성 중 merge_frac 이상이 지목한 잔기만 넣는다.
        cnt = Counter()
        for i in ix:
            cnt.update(cons[i][1])
        need = max(1, int(round(a.merge_frac * len(ix))))
        u = {r for r, n in cnt.items() if n >= need}
        if not u:                       # 너무 빡세서 다 걸러지면 최빈 잔기라도 남긴다
            top = max(cnt.values())
            u = {r for r, n in cnt.items() if n == top}
        rec = len(u & true) / len(true) if true else float("nan")
        pre = len(u & true) / len(u) if u else float("nan")
        names = ",".join(cons[i][0] for i in ix)
        star = " ★" if rec >= 0.5 else ""
        print(f"  {ci:<5}{len(ix):<7}{len(u):<6}{rec:>13.2f}{pre:>13.2f}   {names}{star}")
        out.append(dict(target=tgt, model=model, depth=depth, cand=ci, n_comp=len(ix),
                        n_res=len(u), true_covered=round(rec, 4), precision=round(pre, 4),
                        comps=names))
        # 잔기 키 = (항원 사슬 순번, 그 사슬 참조서열에서의 0-based 위치).
        # 유도 재도킹은 여기서 1-based 서열 위치로 바꿔 쓴다(sites_to_pocket.py).
        sites.append(dict(cand=ci, n_comp=len(ix), comps=names.split(","),
                          # ⚠️ 잔기 위치는 posmap(서열정렬)에서 온 numpy 정수다. 그대로 두면
                          #    json 이 "int64 is not JSON serializable" 로 죽는다(2026-07-29).
                          residues=sorted([[int(x) for x in k] for k in u]),
                          # ⚠️ 아래 둘은 정답 구조를 본 값이다. 보고용이며 '고르는 데' 쓰면 안 된다.
                          #    rank_sites.py 는 읽을 때 이 키를 버리고 그 사실을 화면에 알린다.
                          true_covered=round(rec, 4), precision=round(pre, 4),
                          from_full_msa=any(c.startswith("seedfull") for c in names.split(","))))
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

    # ── 유도 재도킹(guided)용 잔기 목록 ───────────────────────────────────────
    if a.dump_sites:
        sp = f"results/sites_{tgt}.json" if a.dump_sites == "AUTO" else a.dump_sites
        os.makedirs(os.path.dirname(sp) or ".", exist_ok=True)
        json.dump(dict(target=tgt, model=model, depth=depth, cutoff=a.cutoff,
                       n_true_res=len(true), perm_p=round(p, 5),
                       within=round(win, 4), between=round(btw, 4),
                       candidates=sites),
                  open(sp, "w"), indent=1,
                  # numpy 스칼라가 또 새어 들어와도 죽지 않도록(위 캐스팅이 1차 방어).
                  default=lambda o: int(o) if hasattr(o, "__int__") else float(o))
        nfull = sum(1 for s in sites if s["from_full_msa"])
        print(f"→ {sp}  (후보 {len(sites)}개 · 잔기 목록 포함 · 원래 MSA가 속한 후보 {nfull}개)")


if __name__ == "__main__":
    main()
