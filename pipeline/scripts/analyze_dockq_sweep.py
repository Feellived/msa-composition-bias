#!/usr/bin/env python3
"""depth-sweep DockQ 다각도 분석 (boltz + protenix 함께). 순수 stdlib(csv/math/random)만 사용 → 어느 env서도 실행.
입력 = results/dockq_sweep.csv (컬럼: target,group,ab,model,rung,neff80,best_dockq,n_pose).
 ⚠️ best_dockq = rung당 pose 5개 중 최고(best-of-5). rung 11~12개 걸친 max = best-of-~60 → best-of-N 착시 큼 → 순열 null로 통제.
 ⚠️ 이 CSV엔 per-pose·ipTM 없음 → pass@k·ipTM-regret은 별도(confidence 파서) 필요, 여기선 제외.

분석(전부 boltz·protenix 각각):
 ①성공률 3-tier(0.23/0.49/0.80) — full-depth(rung0) vs oracle-over-depth(best rung) vs group별.
 ②depth-민감도 — 타깃별 range(max-min)·std·Spearman(neff,dockq). 모델 간 비교(Boltz 평평 vs Protenix 반응?)=핵심 신규분석.
 ③family×ab 요약 — Δ(저깊이-full)·Spearman.
 ④순열 null — best-of-N 통제(rung 라벨 셔플, max-over-rung 재계산). rescue/gain이 우연 초과인가.
 ⑤평균-기반(비 best-of-N) — 저깊이 평균 vs full. depth 감소가 '평균적으로' 돕나(부호검정).
 ⑥shape 분류 — flat/deep-better/shallow-better/mid-peak(sweet-spot)/spiky.
 ⑦모델 head-to-head + ⑧진짜 rescue 후보(둘다 full<0.23인데 어떤 감소 rung이 ≥0.49) = 가설살리기 case 후보.

사용: python analyze_dockq_sweep.py [--csv results/dockq_sweep.csv] [--outdir results/analysis]
출력: results/analysis/*.csv + summary.txt (+ 그림, matplotlib 있으면).
"""
import argparse, csv, math, os, random
from collections import defaultdict

random.seed(0)

def mean(xs): return sum(xs)/len(xs) if xs else float("nan")
def std(xs):
    if len(xs) < 2: return 0.0
    m = mean(xs); return math.sqrt(sum((x-m)**2 for x in xs)/(len(xs)-1))
def median(xs):
    if not xs: return float("nan")
    s = sorted(xs); n = len(s)
    return s[n//2] if n%2 else (s[n//2-1]+s[n//2])/2

def spearman(x, y):
    n = len(x)
    if n < 3: return float("nan")
    def rank(v):
        order = sorted(range(len(v)), key=lambda i: v[i]); r = [0.0]*len(v); i = 0
        while i < len(v):
            j = i
            while j+1 < len(v) and v[order[j+1]] == v[order[i]]: j += 1
            avg = (i+j)/2.0
            for k in range(i, j+1): r[order[k]] = avg
            i = j+1
        return r
    rx, ry = rank(x), rank(y); mx = mean(rx); my = mean(ry)
    num = sum((a-mx)*(b-my) for a, b in zip(rx, ry))
    dx = math.sqrt(sum((a-mx)**2 for a in rx)); dy = math.sqrt(sum((b-my)**2 for b in ry))
    return num/(dx*dy) if dx*dy else float("nan")

def classify_shape(dq, amp_flat=0.15):
    """dq = rung0(full)..rungN(shallow) 순 DockQ. 반환 (shape, amp, argmax_pos_norm)."""
    if len(dq) < 4: return ("n<4", 0.0, float("nan"))
    amp = max(dq) - min(dq)
    if amp < amp_flat: return ("flat", amp, float("nan"))
    im = max(range(len(dq)), key=lambda i: dq[i])
    pos = im/(len(dq)-1)   # 0=full-depth, 1=shallowest
    thr = min(dq) + 0.5*amp
    peaks = sum(1 for v in dq if v >= thr)
    if peaks >= 4: shape = "spiky"       # 여러 rung이 고르게 높음(=noisy/운)
    elif pos <= 0.25: shape = "deep-better"
    elif pos >= 0.75: shape = "shallow-better"
    else: shape = "mid-peak(sweet-spot)"
    return (shape, amp, pos)

def perm_null(T, key_full=0, fail=0.23, succ=0.49, perms=20000):
    """T = 타깃별 dq 리스트(rung0 first). rung 라벨 셔플로 'full이 특별히 나쁜가' 검정(best-of-N 통제).
       observed: gain=max-full, rescue=(full<fail & max>=succ). null: full 자리에 랜덤 rung."""
    T = [d for d in T if len(d) >= 4]
    if not T: return None
    r0 = [d[key_full] for d in T]; rmax = [max(d) for d in T]
    obs_gain = mean([rmax[i]-r0[i] for i in range(len(T))])
    obs_res = sum(1 for i in range(len(T)) if r0[i] < fail and rmax[i] >= succ)
    ng = [0.0]*perms; nr = [0]*perms
    for p in range(perms):
        r0p = [d[random.randrange(len(d))] for d in T]
        ng[p] = mean([rmax[i]-r0p[i] for i in range(len(T))])
        nr[p] = sum(1 for i in range(len(T)) if r0p[i] < fail and rmax[i] >= succ)
    p_gain = sum(1 for v in ng if v >= obs_gain)/perms
    p_res = sum(1 for v in nr if v >= obs_res)/perms
    return dict(n=len(T), obs_gain=obs_gain, null_gain=mean(ng), p_gain=p_gain,
                obs_res=obs_res, null_res=mean(nr), p_res=p_res)

def sign_test(deltas):
    """저깊이평균-full의 부호검정(중앙값=0 귀무). + 개수 vs - 개수, 이항 근사 p(양측)."""
    pos = sum(1 for d in deltas if d > 1e-9); neg = sum(1 for d in deltas if d < -1e-9)
    n = pos+neg
    if n == 0: return (pos, neg, float("nan"))
    k = max(pos, neg)
    # 이항 양측 p = 2*P(X>=k | n,0.5)
    from math import comb
    tail = sum(comb(n, i) for i in range(k, n+1))/(2**n)
    return (pos, neg, min(1.0, 2*tail))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", default="results/dockq_sweep.csv")
    ap.add_argument("--outdir", default="results/analysis")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    # ── 로드 → (model, target) 별 rung 정렬 ──
    rows = list(csv.DictReader(open(a.csv)))
    key = defaultdict(list)   # (model,target) -> list of (rung, neff, dockq, group, ab)
    for r in rows:
        try:
            key[(r["model"], r["target"])].append(
                (int(r["rung"]), float(r["neff80"]), float(r["best_dockq"]), r["group"], r["ab"]))
        except ValueError:
            continue
    models = sorted({r["model"] for r in rows})
    per = {}   # (model,target) -> dict(group,ab,neff[],dq[](rung순),full,best,amp,shape,spr)
    for (m, t), lst in key.items():
        lst.sort(key=lambda z: z[0])
        neff = [z[1] for z in lst]; dq = [z[2] for z in lst]
        grp, ab = lst[0][3], lst[0][4]
        shape, amp, pos = classify_shape(dq)
        spr = spearman(neff, dq)   # +면 깊을수록 DockQ↑
        per[(m, t)] = dict(model=m, target=t, group=grp, ab=ab, neff=neff, dq=dq,
                           full=dq[0], best=max(dq), best_rung=dq.index(max(dq)),
                           best_neff=neff[dq.index(max(dq))], amp=amp, shape=shape, pos=pos, spr=spr,
                           mean_reduced=mean(dq[1:]) if len(dq) > 1 else float("nan"))

    out = []   # summary.txt 라인
    def P(*s): line = " ".join(str(x) for x in s); print(line); out.append(line)

    # ═══ ① 성공률 3-tier ═══
    P("="*90); P("① 성공률 (타깃 수 / 전체) — full-depth(rung0) vs oracle-over-depth(최고 rung)")
    P(f"{'model':9}{'scope':16}{'n':>4} | full≥.23  full≥.49  full≥.80 | orc≥.23  orc≥.49  orc≥.80")
    succ_rows = []
    for m in models:
        for scope in ["ALL", "A", "B", "C"]:
            ts = [v for (mm, t), v in per.items() if mm == m and (scope == "ALL" or v["ab"] == scope)]
            if not ts: continue
            n = len(ts)
            f23 = sum(v["full"] >= 0.23 for v in ts); f49 = sum(v["full"] >= 0.49 for v in ts); f80 = sum(v["full"] >= 0.80 for v in ts)
            o23 = sum(v["best"] >= 0.23 for v in ts); o49 = sum(v["best"] >= 0.49 for v in ts); o80 = sum(v["best"] >= 0.80 for v in ts)
            P(f"{m:9}{scope:16}{n:>4} | {f23:>7} {f49:>9} {f80:>9} | {o23:>7} {o49:>8} {o80:>8}")
            succ_rows.append(dict(model=m, scope=scope, n=n, full_23=f23, full_49=f49, full_80=f80,
                                  orc_23=o23, orc_49=o49, orc_80=o80))
    _write(os.path.join(a.outdir, "success_rates.csv"), succ_rows)
    P("  ⚠️ orc(oracle-over-depth)=rung 12개 중 최고 = best-of-N 심함 → ④ 순열null로 통제해야 의미.")

    # ═══ ② depth-민감도 (모델 간 비교) ═══
    P("\n" + "="*90); P("② depth-민감도 — 타깃별 range(max-min)·std·|Spearman|. Boltz 평평 vs Protenix 반응 검정")
    sens = defaultdict(list)
    for (m, t), v in per.items():
        sens[m].append((v["amp"], std(v["dq"]), abs(v["spr"]) if v["spr"] == v["spr"] else 0.0))
    P(f"{'model':9}{'n':>4} | mean_range  median_range  mean_std  mean|Spearman|")
    for m in models:
        a_ = [x[0] for x in sens[m]]; s_ = [x[1] for x in sens[m]]; r_ = [x[2] for x in sens[m]]
        P(f"{m:9}{len(a_):>4} | {mean(a_):>10.3f} {median(a_):>13.3f} {mean(s_):>9.3f} {mean(r_):>13.3f}")
    if len(models) == 2:
        m0, m1 = models
        common = [t for (mm, t) in per if mm == m0 and (m1, t) in per]
        d_amp = [per[(m0, t)]["amp"] - per[(m1, t)]["amp"] for t in common]
        pos, neg, p = sign_test(d_amp)
        P(f"  짝지은 range 차({m0}-{m1}) 부호검정: {m0}>{m1} {pos}개 vs 반대 {neg}개, p={p:.4f} "
          f"(mean Δrange={mean(d_amp):+.3f}) → {'유의' if p<0.05 else '비유의'}")
        P(f"  해석: Δrange 부호가 한쪽으로 쏠리면 그 모델이 depth에 더 민감. (가설: protenix가 더 민감 = boltz-protenix<0)")

    # ═══ ③ family×ab 요약 ═══
    P("\n" + "="*90); P("③ family×ab 요약 — Δ(저깊이-full, low-full) & Spearman(neff,dockq) [+개수/n]")
    P(f"{'model':9}{'fam/ab':10}{'n':>4} | mean_full  mean_best | Δ(low-full)  Spearman[+/n]")
    fam_rows = []
    for m in models:
        groups = defaultdict(list)
        for (mm, t), v in per.items():
            if mm == m: groups[(v["group"], v["ab"])].append(v)
        for gkey in sorted(groups):
            vs = groups[gkey]; n = len(vs)
            d_lowfull = [v["dq"][-1] - v["dq"][0] for v in vs]     # 마지막 rung(저깊이) - rung0(full)
            sprs = [v["spr"] for v in vs if v["spr"] == v["spr"]]
            pos = sum(1 for s in sprs if s > 0)
            P(f"{m:9}{gkey[0]+'/'+gkey[1]:10}{n:>4} | {mean([v['full'] for v in vs]):>9.3f} {mean([v['best'] for v in vs]):>9.3f} | "
              f"{mean(d_lowfull):>+10.3f} {mean(sprs) if sprs else float('nan'):>+8.2f}[{pos}/{len(sprs)}]")
            fam_rows.append(dict(model=m, family=gkey[0], ab=gkey[1], n=n,
                                 mean_full=round(mean([v['full'] for v in vs]), 3),
                                 mean_best=round(mean([v['best'] for v in vs]), 3),
                                 delta_lowfull=round(mean(d_lowfull), 3),
                                 mean_spearman=round(mean(sprs), 3) if sprs else "", spearman_pos=f"{pos}/{len(sprs)}"))
    _write(os.path.join(a.outdir, "family_summary.csv"), fam_rows)

    # ═══ ④ 순열 null (best-of-N 통제) ═══
    P("\n" + "="*90); P("④ 순열 null — rung 라벨 셔플, max-over-rung 재계산. 'full이 특별히 나빠 감소가 돕나' 검정")
    P("   (gain=max-full, rescue=full<0.23 & max≥0.49. p<0.05면 우연(best-of-N) 초과=진짜 방향)")
    for m in models:
        T = [v["dq"] for (mm, t), v in per.items() if mm == m]
        res = perm_null(T)
        if res:
            P(f"  [{m}] n={res['n']}  gain 관측={res['obs_gain']:.3f} null={res['null_gain']:.3f} p={res['p_gain']:.3f}  |  "
              f"rescue 관측={res['obs_res']} null평균={res['null_res']:.1f} p={res['p_res']:.3f}")
    # A/B/C 별 순열도
    P("   (family/ab 별)")
    for m in models:
        for ab in ["A", "B", "C"]:
            T = [v["dq"] for (mm, t), v in per.items() if mm == m and v["ab"] == ab]
            res = perm_null(T)
            if res and res["n"] >= 4:
                P(f"    [{m}/{ab}] n={res['n']}  gain p={res['p_gain']:.3f}  rescue {res['obs_res']}/{res['null_res']:.1f} p={res['p_res']:.3f}")

    # ═══ ⑤ 평균-기반(비 best-of-N) ═══
    P("\n" + "="*90); P("⑤ 평균-기반(best-of-N 아님) — 저깊이 rung 평균 vs full. '평균적으로' depth감소가 돕나(부호검정)")
    for m in models:
        vs = [v for (mm, t), v in per.items() if mm == m and len(v["dq"]) > 1]
        deltas = [v["mean_reduced"] - v["full"] for v in vs]
        pos, neg, p = sign_test(deltas)
        P(f"  [{m}] n={len(vs)}  mean(저깊이평균-full)={mean(deltas):+.3f}  개선 {pos} vs 악화 {neg}  p={p:.4f} "
          f"→ {'감소가 평균적으로 도움(유의)' if (p<0.05 and pos>neg) else '감소가 평균적으로 손해(유의)' if (p<0.05 and neg>pos) else '평균으론 효과 없음'}")

    # ═══ ⑥ shape 분류 ═══
    P("\n" + "="*90); P("⑥ depth-response shape 분류 (타깃 수)")
    P(f"{'model':9} | flat  deep-better  shallow-better  mid-peak(sweet)  spiky")
    shape_rows = []
    for m in models:
        cnt = defaultdict(int)
        for (mm, t), v in per.items():
            if mm == m: cnt[v["shape"]] += 1
        P(f"{m:9} | {cnt['flat']:>4} {cnt['deep-better']:>11} {cnt['shallow-better']:>15} {cnt['mid-peak(sweet-spot)']:>16} {cnt['spiky']:>7}")
        for sh, c in cnt.items(): shape_rows.append(dict(model=m, shape=sh, count=c))
    _write(os.path.join(a.outdir, "shape_counts.csv"), shape_rows)

    # ═══ ⑦ 모델 head-to-head ═══
    P("\n" + "="*90); P("⑦ 모델 head-to-head (두 모델 다 있는 타깃) — full & oracle 각각 누가 이기나")
    if len(models) == 2:
        m0, m1 = models
        common = sorted(t for (mm, t) in per if mm == m0 and (m1, t) in per)
        fw = [0, 0, 0]; ow = [0, 0, 0]   # m0승 m1승 무
        for t in common:
            a0, a1 = per[(m0, t)]["full"], per[(m1, t)]["full"]
            fw[0 if a0 > a1+0.02 else 1 if a1 > a0+0.02 else 2] += 1
            b0, b1 = per[(m0, t)]["best"], per[(m1, t)]["best"]
            ow[0 if b0 > b1+0.02 else 1 if b1 > b0+0.02 else 2] += 1
        P(f"  full-depth: {m0} 승 {fw[0]} / {m1} 승 {fw[1]} / 무 {fw[2]}  (n={len(common)})")
        P(f"  oracle    : {m0} 승 {ow[0]} / {m1} 승 {ow[1]} / 무 {ow[2]}")
        P(f"  consensus 상한(둘 중 최고, oracle-over-depth): "
          f"≥.49 {sum(max(per[(m0,t)]['best'],per[(m1,t)]['best'])>=0.49 for t in common)}/{len(common)}  "
          f"≥.80 {sum(max(per[(m0,t)]['best'],per[(m1,t)]['best'])>=0.80 for t in common)}/{len(common)}")

    # ═══ ⑧ 진짜 rescue 후보 (가설살리기 case) ═══
    P("\n" + "="*90); P("⑧ ⭐ rescue 후보 = full-depth 실패인데 감소 rung이 성공 (가설살리기 case-study 후보)")
    P("   (A) 모델별: full<0.23 인데 어떤 rung ≥0.49")
    resc_rows = []
    for m in models:
        for (mm, t), v in per.items():
            if mm != m: continue
            if v["full"] < 0.23 and v["best"] >= 0.49:
                P(f"    [{m}] {t:10} {v['group']}/{v['ab']} full={v['full']:.2f} → best={v['best']:.2f} "
                  f"@rung{v['best_rung']}(neff{v['best_neff']:.0f}) shape={v['shape']}")
                resc_rows.append(dict(model=m, target=t, group=v["group"], ab=v["ab"], full=round(v["full"],3),
                                      best=round(v["best"],3), best_rung=v["best_rung"], best_neff=round(v["best_neff"],1),
                                      shape=v["shape"], kind="model_rescue"))
    P("   (B) ⭐⭐ 강한 후보: 두 모델 다 full<0.23 인데, 어떤 모델의 어떤 rung이 ≥0.49 (아무도 full로 못 푸는데 감소로 건짐)")
    if len(models) == 2:
        m0, m1 = models
        common = sorted(t for (mm, t) in per if mm == m0 and (m1, t) in per)
        for t in common:
            f0, f1 = per[(m0, t)]["full"], per[(m1, t)]["full"]
            b0, b1 = per[(m0, t)]["best"], per[(m1, t)]["best"]
            if f0 < 0.23 and f1 < 0.23 and max(b0, b1) >= 0.49:
                win = m0 if b0 >= b1 else m1
                P(f"    ⭐ {t:10} {per[(m0,t)]['group']}/{per[(m0,t)]['ab']}: 둘다 full<0.23 (b={m0}{f0:.2f}/{m1}{f1:.2f}) "
                  f"→ {win} best={max(b0,b1):.2f} 로 rescue")
                resc_rows.append(dict(model=win, target=t, group=per[(m0,t)]["group"], ab=per[(m0,t)]["ab"],
                                      full=round(min(f0,f1),3), best=round(max(b0,b1),3), best_rung="", best_neff="",
                                      shape="both-full-fail", kind="strong_rescue"))
    _write(os.path.join(a.outdir, "rescue_candidates.csv"), resc_rows)

    # per-target 전체 요약(가설살리기 선별용)
    pt_rows = []
    for (m, t), v in sorted(per.items()):
        pt_rows.append(dict(model=m, target=t, family=v["group"], ab=v["ab"], n_rung=len(v["dq"]),
                            full=round(v["full"],3), best=round(v["best"],3), best_rung=v["best_rung"],
                            best_neff=round(v["best_neff"],1), amp=round(v["amp"],3),
                            spearman=round(v["spr"],3) if v["spr"]==v["spr"] else "", shape=v["shape"]))
    _write(os.path.join(a.outdir, "per_target_summary.csv"), pt_rows)

    with open(os.path.join(a.outdir, "summary.txt"), "w") as f:
        f.write("\n".join(out))
    P("\n" + "="*90)
    P(f"→ 저장: {a.outdir}/ (success_rates·family_summary·shape_counts·rescue_candidates·per_target_summary·summary.txt)")

    # ═══ 그림 (matplotlib 있으면) ═══
    try:
        import matplotlib; matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        # (a) 모델별 depth-민감도(range) 히스토그램
        fig, ax = plt.subplots(figsize=(7, 4.5), dpi=130)
        for m, c in zip(models, ["#0072B2", "#D55E00"]):
            amps = [v["amp"] for (mm, t), v in per.items() if mm == m]
            ax.hist(amps, bins=[i*0.05 for i in range(19)], alpha=0.55, label=f"{m} (median {median(amps):.2f})", color=c)
        ax.set_xlabel("depth-range (max-min DockQ over rungs)"); ax.set_ylabel("타깃 수")
        ax.set_title("depth-민감도: Boltz(평평) vs Protenix"); ax.legend()
        fig.tight_layout(); fig.savefig(os.path.join(a.outdir, "fig_depth_sensitivity.png")); plt.close(fig)
        # (b) rescue 후보 궤적 small-multiples
        cand = sorted({t for r in resc_rows for t in [r["target"]]})
        if cand:
            ncol = 3; nrow = (len(cand)+ncol-1)//ncol
            fig, axs = plt.subplots(nrow, ncol, figsize=(4*ncol, 3*nrow), dpi=120, squeeze=False)
            for i, t in enumerate(cand):
                ax = axs[i//ncol][i%ncol]
                for m, c in zip(models, ["#0072B2", "#D55E00"]):
                    if (m, t) in per:
                        v = per[(m, t)]
                        ax.plot(v["neff"], v["dq"], "o-", ms=4, color=c, label=m)
                ax.set_xscale("log"); ax.axhline(0.49, ls="--", c="gray", lw=0.8); ax.axhline(0.23, ls=":", c="gray", lw=0.8)
                ax.set_title(t, fontsize=9); ax.set_ylim(-0.03, 1.0)
                if i == 0: ax.legend(fontsize=7)
            for j in range(len(cand), nrow*ncol): axs[j//ncol][j%ncol].axis("off")
            fig.suptitle("rescue 후보 궤적 (x=Neff log, --=0.49 --=0.23)")
            fig.tight_layout(); fig.savefig(os.path.join(a.outdir, "fig_rescue_candidates.png")); plt.close(fig)
        P("  그림 저장: fig_depth_sensitivity.png · fig_rescue_candidates.png")
    except Exception as e:
        P(f"  (matplotlib 없거나 그림 실패: {type(e).__name__} — CSV/txt만 저장됨)")

def _write(path, rows):
    if not rows: return
    cols = list(rows[0].keys())
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols); w.writeheader(); w.writerows(rows)

if __name__ == "__main__":
    main()
