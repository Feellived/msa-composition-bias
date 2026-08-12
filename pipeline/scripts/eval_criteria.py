#!/usr/bin/env python3
"""[최종 판정] 6.4 절에서 결과를 보기 전에 확정한 세 기준으로 본 검정을 채점한다 → 6.7 기입용.

  ① 효과    조성 간 이질성이 유의한 복합체 수 + 우연 확률 + 복합체별 p 의 Fisher 결합
  ② 재현성  뒤섞기(순열) 검정 p < 0.05 인 복합체 수
  ③ 후보생성 정답 자리를 절반 이상 덮는 후보가 하나라도 있는 복합체 수

주 지표는 결합자리 겹침(recall, 성공선 0.4). 자세 정확도(DockQ)는 따로 찍는다.

읽는 파일 (run_analyze_target.sh 가 만든 것):
  results/summary_<타깃>_recall.csv   p_heterogeneity
  results/summary_<타깃>_dockq.csv    성공 수(자세 정확도 축이 살아 있는지 확인용)
  results/site_repro_<타깃>.csv       perm_p · true_covered · comps
  maintest.csv                        명단·군(group)·층(stratum)

⚠️ 이 스크립트는 사전 확정 기준을 적용할 뿐 기준을 만들지 않는다.
   분모가 커지면 우연 확률과 Fisher 결합값이 함께 바뀐다 — 그것도 같이 찍는다.

사용:
  python eval_criteria.py
  python eval_criteria.py --set3 "8k3k_D 8k46_I ..."      # 세트 3 명단(빈도 계산용)
  python eval_criteria.py --exclude 8ulr_HL --csv-out results/criteria.csv
"""
import argparse, csv, math, os
from collections import defaultdict

SIG = 0.05          # 유의 기준
COVER = 0.5         # 후보가 정답을 '절반 이상' 덮는다는 기준


def fnum(x, d=None):
    try:
        v = float(x)
        return d if v != v else v
    except (TypeError, ValueError):
        return d


def fisher_combine(ps):
    """복합체별 p 를 하나로 합친다(Fisher). 자유도가 짝수라 정확식으로 계산한다."""
    ps = [min(max(p, 1e-12), 1.0) for p in ps if p is not None]
    if not ps:
        return None
    x = -2 * sum(math.log(p) for p in ps)
    k = len(ps)
    tail = sum((x / 2) ** i / math.factorial(i) for i in range(k))
    return math.exp(-x / 2) * tail


def binom_tail(k, n, p=SIG):
    """'유의한 것이 k 개 이상 나올' 우연 확률. 각 복합체가 우연히 유의할 확률을 p 로 본다."""
    if n == 0:
        return None
    return sum(math.comb(n, i) * p ** i * (1 - p) ** (n - i) for i in range(k, n + 1))


def bh_survivors(ps, q=SIG):
    """Benjamini-Hochberg 보정으로 살아남는 개수."""
    v = sorted(p for p in ps if p is not None)
    n = len(v)
    m = 0
    for i, p in enumerate(v, 1):
        if p <= q * i / n:
            m = i
    return m


def read_target(t):
    """한 복합체의 판정 재료를 모은다."""
    out = dict(het_p=None, perm_p=None, best_cover=None, full_cover=None,
               dq_succ=None, dq_het=None, n_cand=None)
    p = "results/summary_%s_recall.csv" % t
    if os.path.exists(p):
        for r in csv.DictReader(open(p)):
            out["het_p"] = fnum(r.get("p_heterogeneity")); break
    p = "results/summary_%s_dockq.csv" % t
    if os.path.exists(p):
        for r in csv.DictReader(open(p)):
            out["dq_het"] = fnum(r.get("p_heterogeneity"))
            out["dq_succ"] = (fnum(r.get("succ49_full"), 0) or 0) + (fnum(r.get("succ49_red"), 0) or 0)
            break
    p = "results/site_repro_%s.csv" % t
    if os.path.exists(p):
        cs = list(csv.DictReader(open(p)))
        if cs:
            out["perm_p"] = fnum(cs[0].get("perm_p"))
            out["n_cand"] = len(cs)
            covers = [(fnum(c.get("true_covered"), 0) or 0) for c in cs]
            out["best_cover"] = max(covers)
            for c in cs:
                if "seedfull" in (c.get("comps") or ""):
                    out["full_cover"] = fnum(c.get("true_covered"))
                    break
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--maintest", default="maintest.csv")
    ap.add_argument("--set4", default="sweep_targets.csv", help="세트 4 명단(49종)")
    ap.add_argument("--set3", default="set3_targets.csv",
                    help="세트 3 명단(10종). 파일 경로 또는 공백 구분 문자열. "
                         "PDB ID 만 있어도 되고(복합체 식별자의 앞부분으로 대조), "
                         "세트 4 분모에서는 자동으로 뺀다")
    ap.add_argument("--exclude", default="8ulr_HL",
                    help="설계가 달라 같은 표에 못 넣는 복합체(공백 구분)")
    ap.add_argument("--csv-out", default="results/criteria.csv")
    a = ap.parse_args()
    skip = set(a.exclude.split())

    meta, design = {}, {}
    for r in csv.DictReader(open(a.maintest)):
        if r.get("status") == "run":
            meta[r["target"]] = dict(group=r.get("group", ""), stratum=r.get("stratum", ""))
            # 설계값 = 조성 수 × 반복 수 + 원래 MSA 횟수. DockQ 성공률의 분모로 쓴다.
            design[r["target"]] = (int(r.get("n_comp") or 6) * int(r.get("n_reps") or 4)
                                   + int(r.get("n_full") or 8))
    set4 = set()
    if os.path.exists(a.set4):
        set4 = {r["target"] for r in csv.DictReader(open(a.set4))}
    # 세트 3 — PDB ID 로 적혀 있으므로 복합체 식별자(<pdb>_<사슬>)의 앞부분과 대조한다.
    set3_pdb = set()
    if a.set3:
        if os.path.exists(a.set3):
            for r in csv.DictReader(open(a.set3)):
                if r.get("pdb"):
                    set3_pdb.add(r["pdb"].strip().lower())
        else:
            set3_pdb = {x.strip().lower().split("_")[0] for x in a.set3.split()}
    def in_set3(t):
        return t.split("_")[0].lower() in set3_pdb
    # ⚠️ sweep_targets.csv 는 세트 4(49종)에 세트 3(10종)이 더해진 59행이다.
    #    세트 4 분모는 49 여야 하므로 세트 3 멤버를 뺀다(사전 확정 명단 크기를 지킨다).
    set4 = {t for t in set4 if not in_set3(t)}
    set3 = {t for t in meta if in_set3(t)}

    rows = []
    for t in sorted(meta):
        if t in skip:
            continue
        d = read_target(t)
        d.update(target=t, group=meta[t]["group"], stratum=meta[t]["stratum"])
        d["scored"] = d["het_p"] is not None and d["perm_p"] is not None
        rows.append(d)

    done = [r for r in rows if r["scored"]]
    miss = [r["target"] for r in rows if not r["scored"]]
    n = len(done)
    print("■ 본 검정 판정 — 채점된 복합체 %d개%s"
          % (n, ("  (제외: %s)" % ", ".join(sorted(skip)) if skip else "")))
    if miss:
        print("  ⚠️ 아직 채점 안 된 것 %d개: %s" % (len(miss), ", ".join(miss)))
    if n == 0:
        raise SystemExit("!! 채점된 복합체가 없다 — run_analyze_target.sh 를 먼저 돌릴 것")

    # ── 복합체별 표 ────────────────────────────────────────────────────────────
    print()
    print("  %-12s %-4s %-6s %9s %8s %7s %9s %6s"
          % ("복합체", "군", "층", "이질성p", "순열p", "후보수", "최고덮음", "원래"))
    print("  " + "-" * 68)
    for r in sorted(done, key=lambda x: (x["het_p"] is None, x["het_p"])):
        st = {"rich": "다양O", "poor": "다양X"}.get(r["stratum"], "-")
        het = "%.4f" % r["het_p"] if r["het_p"] is not None else "  -  "
        prm = "%.4f" % r["perm_p"] if r["perm_p"] is not None else "  -  "
        bc = "%.2f" % r["best_cover"] if r["best_cover"] is not None else " - "
        fc = "%.2f" % r["full_cover"] if r["full_cover"] is not None else " - "
        star = " ★" if (r["het_p"] is not None and r["het_p"] < SIG) else ""
        print("  %-12s %-4s %-6s %9s %8s %7s %9s %6s%s"
              % (r["target"], r["group"], st, het, prm, r["n_cand"] or "-", bc, fc, star))

    # ── 세 기준 ───────────────────────────────────────────────────────────────
    sig = [r for r in done if r["het_p"] is not None and r["het_p"] < SIG]
    rep = [r for r in done if r["perm_p"] is not None and r["perm_p"] < SIG]
    gen = [r for r in done if (r["best_cover"] or 0) >= COVER]
    fp = fisher_combine([r["het_p"] for r in done])
    chance = binom_tail(len(sig), n)
    nbh = bh_survivors([r["het_p"] for r in done])

    print()
    print("■ 6.4 사전 기준 대조  (n = %d)" % n)
    print("  ① 효과 — 이질성 유의        %3d / %d    우연 확률 %.4f   %s"
          % (len(sig), n, chance, "통과" if chance < SIG else "미달"))
    print("     유의: %s" % (", ".join(r["target"] for r in sig) if sig else "없음"))
    print("     BH 보정(q=0.05) 생존     %3d / %d" % (nbh, n))
    print("  ① 효과 — Fisher 결합        p = %.6f   %s"
          % (fp, "통과" if fp is not None and fp < SIG else "미달"))
    fail_rep = [r["target"] for r in done if r not in rep]
    print("  ② 재현성 — 순열 p<0.05      %3d / %d %s"
          % (len(rep), n, ("  (실패: %s)" % ", ".join(fail_rep)) if fail_rep else "  실패 없음"))
    print("  ③ 후보생성 — 정답 절반 덮음 %3d / %d" % (len(gen), n))

    # ── 원래 MSA 보다 나은 후보가 생긴 복합체 ─────────────────────────────────
    # ⚠️ 6.6.2 표는 '원래 MSA 보다 더 덮는 후보가 있다'(차이 > 0)로 뽑았고, 그중 일부는
    #    원래 MSA 도 이미 절반 이상 덮고 있었다(9mqr_DE 0.50 · 8siq_HL 0.59 · 8t4d_OQ 0.64).
    #    그래서 두 가지를 함께 찍는다 — 차이가 있는 전부(6.6.2 와 같은 기준)와,
    #    그중 '원래는 절반을 못 덮었는데 후보는 넘긴' 것(더 강한 주장).
    gains = [r for r in done
             if r["full_cover"] is not None and r["best_cover"] is not None
             and r["best_cover"] > r["full_cover"]]
    crossed = [r for r in gains if r["full_cover"] < COVER <= r["best_cover"]]
    print()
    print("■ 원래 MSA 보다 더 덮는 후보가 생긴 복합체 — %d개 (6.6.2 와 같은 기준)" % len(gains))
    if gains:
        print("  %-12s %10s %10s %8s %7s %9s"
              % ("복합체", "원래 덮음", "최고 덮음", "차이", "후보수", "순열p"))
        print("  " + "-" * 62)
        for r in sorted(gains, key=lambda x: -(x["best_cover"] - x["full_cover"])):
            mark = " ★" if r in crossed else ""
            print("  %-12s %10.2f %10.2f %8.2f %7d %9.4f%s"
                  % (r["target"], r["full_cover"], r["best_cover"],
                     r["best_cover"] - r["full_cover"], r["n_cand"] or 0, r["perm_p"], mark))
    print("  ★ = 원래 MSA 는 절반을 못 덮었는데 후보가 넘긴 경우 — %d개. 주장에 쓸 것은 이쪽이다."
          % len(crossed))
    print("  ⚠️ 기준 ③(정답 절반 덮음)은 관대하다 — 후보가 하나뿐인 복합체는 '조성이 새 후보를")
    print("     만든 것이 아니라 전부 같은 자리로 갔다'는 뜻이다.")

    # ── 자세 정확도 축 ────────────────────────────────────────────────────────
    dq_tot = sum(r["dq_succ"] or 0 for r in done)
    dq_any = [r["target"] for r in done if (r["dq_succ"] or 0) > 0]
    dq_sig = [r["target"] for r in done if r["dq_het"] is not None and r["dq_het"] < SIG]
    n_runs = sum(design.get(r["target"], 0) for r in done)
    print()
    print("■ 자세 정확도(DockQ) 축")
    print("  DockQ ≥ 0.49 인 실행 %d회 / 전체 %d회  (복합체 %d개: %s)"
          % (dq_tot, n_runs, len(dq_any), ", ".join(dq_any) if dq_any else "없음"))
    print("  DockQ 이질성이 유의한 복합체 %d개%s"
          % (len(dq_sig), ("  (%s)" % ", ".join(dq_sig)) if dq_sig else ""))
    print("  → 성공이 거의 없으면 이 축의 검정은 정보를 내지 못한다.")
    print("     정확한 서술은 '자리는 찾되 자세는 맞히지 못한다'.")

    # ── 명단·층·군별 빈도 ─────────────────────────────────────────────────────
    print()
    print("■ 빈도 — 이질성이 유의한 복합체 수를 각 명단의 분모로 나눈 값")
    for name, members, denom in (("세트 3", set3, len(set3_pdb)), ("세트 4", set4, len(set4))):
        if not members:
            print("  %s: 명단이 없어 계산하지 않음 (--set3 로 넘길 것)" % name)
            continue
        inset = [r for r in done if r["target"] in members]
        s = sum(1 for r in inset if r["het_p"] is not None and r["het_p"] < SIG)
        print("  %s(분모 %d): 채점 %d개 중 유의 %d  →  %d/%d = %.3f"
              % (name, denom, len(inset), s, s, denom, s / denom))
        if len(inset) < denom:
            print("     ※ 분모는 사전 확정 명단 전체다. 본 검정에 안 들어간 %d개는 '유의하지 않음'으로"
                  " 세는 것과 같다(안전한 방향)." % (denom - len(inset)))
    for key, lab, denom in (("stratum", "층", None), ("group", "군", None)):
        g = defaultdict(lambda: [0, 0])
        for r in done:
            v = r[key] or "-"
            g[v][1] += 1
            if r["het_p"] is not None and r["het_p"] < SIG:
                g[v][0] += 1
        print("  %s별: %s" % (lab, " · ".join("%s %d/%d" % (k, v[0], v[1]) for k, v in sorted(g.items()))))
    print("  ⚠️ 층은 성적이 아니라 입력(Neff80)으로 나눈 것이다. 다만 이 자료에서는 층이 항원 계열과")
    print("     거의 겹쳐 '조성 다양성 때문'과 '항원 계열 때문'을 분리할 수 없다(6.2.2 교란).")

    os.makedirs(os.path.dirname(a.csv_out) or ".", exist_ok=True)
    cols = ["target", "group", "stratum", "het_p", "perm_p", "n_cand",
            "best_cover", "full_cover", "dq_succ", "dq_het", "scored"]
    with open(a.csv_out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    print()
    print("→ %s" % a.csv_out)


if __name__ == "__main__":
    main()
