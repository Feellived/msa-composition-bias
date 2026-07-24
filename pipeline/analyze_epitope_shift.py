#!/usr/bin/env python3
"""에피토프 위치 편향 다각도 분석 통합 (boltz + protenix 등 여러 모델 비교). 순수 stdlib만 사용.
입력 = epitope_shift.py 출력 CSV 1개 이상(모델별로 나눠 실행했으면 여러 파일 전달, model 컬럼으로 구분).
 예: python analyze_epitope_shift.py --csv results/epitope_shift_boltz.csv results/epitope_shift_protenix.csv

Boltz에서 했던 것과 정확히 같은 5가지 분석 + 모델 2개 이상일 때 자연히 추가되는 교차확인 1가지:
 ①family×ab 요약 — Δ(low-full)·Spearman(neff,·) for over-rep/recall/pop_rank/dcc_pop, 모델별 나란히.
 ②관대한 지표(threshold-free) — true_rank·dcc_true 같은 요약(recall이 너무 엄격한 거 아니냐는 검토용).
 ③이탈-정답 개별연동 — B그룹 안에서 '이탈 정도'(over-rep Δ)와 '정답행 정도'(true_rank Δ)가 타깃별로
   같이 움직이는지 corr. (이탈≠정답행 임을 재확인하거나 뒤집는 용도)
 ④'진짜 편향' 필터 — B라벨인데 full-depth mean_overrep이 native_overrep보다 뚜렷이(excess≥0.3) 높은
   타깃만 골라냄(B라벨 전부를 편향으로 취급하는 오류 방지).
 ⑤flat-궤적 확인 — ④에서 걸러진 '진짜 편향' 타깃들이 깊이를 줄여도(full vs 최저깊이) 그대로인지.
 ⑥⭐모델 간 교차확인(신규, 모델 2개 이상일 때) — 같은 타깃이 두 모델에서 똑같이 편향/flat인지,
   아니면 모델마다 다른지(다르면 = 계열 자체 문제 아니라 모델-특이적 = 재랭커 근거 강화).

사용: python analyze_epitope_shift.py --csv results/epitope_shift_boltz.csv results/epitope_shift_protenix.csv
      [--out-dir results/analysis] [--excess 0.3]
"""
import argparse, csv, math, os
from collections import defaultdict

def mean(xs): return sum(xs) / len(xs) if xs else float("nan")

def spearman(x, y):
    n = len(x)
    if n < 3: return float("nan")
    def rank(v):
        order = sorted(range(len(v)), key=lambda i: v[i]); r = [0.0] * len(v); i = 0
        while i < len(v):
            j = i
            while j + 1 < len(v) and v[order[j + 1]] == v[order[i]]: j += 1
            avg = (i + j) / 2.0
            for k in range(i, j + 1): r[order[k]] = avg
            i = j + 1
        return r
    rx, ry = rank(x), rank(y); mx = mean(rx); my = mean(ry)
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    dx = math.sqrt(sum((a - mx) ** 2 for a in rx)); dy = math.sqrt(sum((b - my) ** 2 for b in ry))
    return num / (dx * dy) if dx * dy else float("nan")

def load_rows(paths):
    rows = []
    for p in paths:
        for r in csv.DictReader(open(p)):
            if "model" not in r or not r["model"]:
                r["model"] = os.path.basename(p).replace("epitope_shift_", "").replace(".csv", "")  # 구식 파일(model 없음) fallback
            rows.append(r)
    return rows

def f(r, k):
    v = r.get(k, "")
    return float(v) if v not in ("", None) else float("nan")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", nargs="+", required=True)
    ap.add_argument("--out-dir", default="results/analysis")
    ap.add_argument("--excess", type=float, default=0.3, help="'진짜 편향' 판정 문턱(full_overrep - native_overrep)")
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    rows = load_rows(a.csv)
    models = sorted({r["model"] for r in rows})
    print(f"로드: {len(rows)}행, 모델={models}")

    # (model,target) -> rung 정렬된 레코드 리스트
    byt = defaultdict(list)
    for r in rows:
        try:
            byt[(r["model"], r["target"])].append(r)
        except KeyError:
            continue
    for k in byt: byt[k].sort(key=lambda z: int(z["rung"]))

    out = []
    def P(*s): line = " ".join(str(x) for x in s); print(line); out.append(line)

    METRICS = ["mean_overrep", "mean_recall", "mean_true_rank", "mean_dcc_true", "mean_pop_rank", "mean_dcc_pop"]

    # 타깃별 요약(Δ,rho) 계산 — 전 분석에서 공용으로 씀
    per = {}   # (model,target) -> dict(family,ab,{metric:(delta_lowfull,rho)}, native_overrep, full_overrep, low_overrep)
    for (m, t), xs in byt.items():
        valid = [x for x in xs if x.get("neff80", "") != ""]
        if len(set(round(float(z["neff80"]), 1) for z in valid)) < 4: continue
        # neff80 오름차순 정렬(=낮은 깊이→높은 깊이) 후에야 low-full이 맞게 계산됨(rung 순서 ≠ neff 순서일 수 있음)
        valid_sorted = sorted(valid, key=lambda z: float(z["neff80"]))
        neff_s = [float(z["neff80"]) for z in valid_sorted]
        d2 = {}
        for key in METRICS:
            has = all(z.get(key, "") != "" for z in valid_sorted)
            if not has: d2[key] = (float("nan"), float("nan")); continue
            vals = [float(z[key]) for z in valid_sorted]
            d2[key] = (vals[0] - vals[-1], spearman(neff_s, vals))
        native_or = f(valid_sorted[0], "native_overrep")
        full_or = f(valid_sorted[-1], "mean_overrep")   # 가장 깊은(neff 최대) = full
        low_or = f(valid_sorted[0], "mean_overrep")     # 가장 얕은 = low
        per[(m, t)] = dict(family=xs[0]["family"], ab=xs[0]["ab"], d=d2,
                           native_overrep=native_or, full_overrep=full_or, low_overrep=low_or,
                           n_pop=int(float(xs[0].get("n_pop", 0) or 0)))

    # ═══ ① family×ab 요약 (모델별 나란히) ═══
    P("=" * 110); P("① family×ab 요약 — Δ(low-full)·Spearman(neff,·) [+개수/n] — 모델별")
    fam_rows = []
    for m in models:
        P(f"\n[{m}]")
        P(f"{'fam/ab':10}{'n':>3} | over-rep Δ,ρ         | pop_rank Δ,ρ         | dcc_pop Δ,ρ           | recall Δ,ρ")
        groups = defaultdict(list)
        for (mm, t), v in per.items():
            if mm == m: groups[(v["family"], v["ab"])].append(v)
        for gkey in sorted(groups):
            vs = groups[gkey]; n = len(vs)
            def col(key):
                ds = [v["d"][key][0] for v in vs if v["d"][key][0] == v["d"][key][0]]
                rs = [v["d"][key][1] for v in vs if v["d"][key][1] == v["d"][key][1]]
                if not ds: return "NA"
                pos = sum(1 for x in rs if x > 0)
                return f"{mean(ds):+.3f},{mean(rs):+.2f}[{pos}/{len(rs)}]" if rs else f"{mean(ds):+.3f},NA"
            P(f"{gkey[0]+'/'+gkey[1]:10}{n:>3} | {col('mean_overrep'):21} | {col('mean_pop_rank'):21} | "
              f"{col('mean_dcc_pop'):21} | {col('mean_recall')}")
            for key in METRICS:
                ds = [v["d"][key][0] for v in vs if v["d"][key][0] == v["d"][key][0]]
                rs = [v["d"][key][1] for v in vs if v["d"][key][1] == v["d"][key][1]]
                fam_rows.append(dict(model=m, family=gkey[0], ab=gkey[1], metric=key, n=n,
                                     delta_lowfull=round(mean(ds), 3) if ds else "",
                                     mean_rho=round(mean(rs), 3) if rs else "",
                                     rho_pos=f"{sum(1 for x in rs if x > 0)}/{len(rs)}" if rs else ""))
    _write(os.path.join(a.out_dir, "epitope_family_summary.csv"), fam_rows)

    # ═══ ② 관대한 지표(threshold-free) 요약 ═══
    P("\n" + "=" * 110); P("② 관대한 지표(true_rank·dcc_true) — recall이 너무 엄격한가 재검토, 모델별")
    for m in models:
        P(f"\n[{m}]")
        groups = defaultdict(list)
        for (mm, t), v in per.items():
            if mm == m: groups[(v["family"], v["ab"])].append(v)
        for gkey in sorted(groups):
            vs = groups[gkey]
            def col(key):
                ds = [v["d"][key][0] for v in vs if v["d"][key][0] == v["d"][key][0]]
                rs = [v["d"][key][1] for v in vs if v["d"][key][1] == v["d"][key][1]]
                pos = sum(1 for x in rs if x > 0)
                return f"{mean(ds):+.3f}({mean(rs):+.2f}[{pos}/{len(rs)}])" if ds else "NA"
            P(f"  {gkey[0]+'/'+gkey[1]:10} true_rank(관대,기대+)={col('mean_true_rank'):22}  "
              f"dcc_true(관대,기대-,Å)={col('mean_dcc_true')}")

    # ═══ ③ 이탈-정답 개별연동 (B그룹만) ═══
    P("\n" + "=" * 110); P("③ B그룹 안에서 '이탈량(over-rep Δ)' vs '정답행(true_rank Δ)' 타깃별 연동 — 모델별")
    for m in models:
        for fam in ["RBD", "HA", "Env"]:
            vs = [v for (mm, t), v in per.items() if mm == m and v["family"] == fam and v["ab"] == "B"]
            xo = [v["d"]["mean_overrep"][0] for v in vs if v["d"]["mean_overrep"][0] == v["d"]["mean_overrep"][0]
                  and v["d"]["mean_true_rank"][0] == v["d"]["mean_true_rank"][0]]
            yt = [v["d"]["mean_true_rank"][0] for v in vs if v["d"]["mean_overrep"][0] == v["d"]["mean_overrep"][0]
                  and v["d"]["mean_true_rank"][0] == v["d"]["mean_true_rank"][0]]
            if len(xo) < 3: continue
            rho = spearman(xo, yt)
            tag = "(이탈↔정답 연동)" if rho < -0.3 else "(반대 연동!)" if rho > 0.3 else "(연동 약함/없음)"
            P(f"  [{m}] {fam}/B (n={len(xo)}): corr(over-rep Δ, true_rank Δ) = {rho:+.2f} {tag}")

    # ═══ ④ '진짜 편향' 필터 ═══
    P("\n" + "=" * 110); P(f"④ '진짜 편향' 필터 — B라벨인데 full-depth overrep이 native보다 뚜렷이(excess≥{a.excess}) 높은 타깃")
    biased_rows = []
    biased = defaultdict(set)   # model -> set(target)
    for m in models:
        P(f"\n[{m}]")
        for (mm, t), v in per.items():
            if mm != m or v["ab"] != "B" or v["family"] not in ("RBD", "HA", "Env"): continue
            if v["native_overrep"] != v["native_overrep"] or v["full_overrep"] != v["full_overrep"]: continue
            excess = v["full_overrep"] - v["native_overrep"]
            tag = "★편향" if excess >= a.excess else ("약함" if excess >= 0.1 else "편향아님")
            P(f"  {t:12} {v['family']:4} native={v['native_overrep']:.3f} full={v['full_overrep']:.3f} "
              f"excess={excess:+.3f}  {tag}")
            biased_rows.append(dict(model=m, target=t, family=v["family"], native_overrep=round(v["native_overrep"], 3),
                                    full_overrep=round(v["full_overrep"], 3), excess=round(excess, 3), tag=tag))
            if tag == "★편향": biased[m].add(t)
        P(f"  → [{m}] 진짜편향 {len(biased[m])}개: {sorted(biased[m])}")
    _write(os.path.join(a.out_dir, "epitope_biased_filter.csv"), biased_rows)

    # ═══ ⑤ flat-궤적 확인 (진짜편향 타깃만) ═══
    P("\n" + "=" * 110); P("⑤ '진짜 편향' 타깃의 깊이별 궤적 — 깊이 줄여도(full→최저) over-rep이 그대로인가")
    flat_rows = []
    for m in models:
        for t in sorted(biased[m]):
            v = per[(m, t)]
            moved = v["full_overrep"] - v["low_overrep"]
            shape = "flat(안 움직임)" if abs(moved) < 0.05 else "이탈(감소)" if moved > 0 else "반대(증가)"
            P(f"  [{m}] {t:12} full={v['full_overrep']:.3f} → low={v['low_overrep']:.3f}  Δ={-moved:+.3f}  {shape}")
            flat_rows.append(dict(model=m, target=t, full_overrep=round(v["full_overrep"], 3),
                                  low_overrep=round(v["low_overrep"], 3), delta=round(-moved, 3), shape=shape))
    _write(os.path.join(a.out_dir, "epitope_biased_trajectory.csv"), flat_rows)

    # ═══ ⑥ 모델 간 교차확인 (모델 2개 이상) ═══
    if len(models) >= 2:
        P("\n" + "=" * 110); P("⑥ ⭐ 모델 간 교차확인 — 같은 타깃이 두 모델에서 똑같이 편향/flat인가")
        all_targets = sorted(set().union(*[biased[m] for m in models]))
        cross_rows = []
        for t in all_targets:
            flags = {m: (t in biased[m]) for m in models}
            both = all(flags.values())
            tag = "양쪽 다 편향(계열/항원 자체 문제 가능성)" if both else "한쪽만 편향(모델-특이적, 재랭커로 해결 가능성↑)"
            P(f"  {t:12} " + " / ".join(f"{m}={'편향' if flags[m] else '-'}" for m in models) + f"  → {tag}")
            cross_rows.append(dict(target=t, **{f"biased_{m}": flags[m] for m in models}, both_biased=both))
        _write(os.path.join(a.out_dir, "epitope_cross_model.csv"), cross_rows)
        n_both = sum(1 for t in all_targets if all(t in biased[m] for m in models))
        P(f"\n  → {len(all_targets)}개 편향타깃 중 양쪽 다 편향={n_both}, 한쪽만={len(all_targets)-n_both}")
        P("  해석: '한쪽만 편향'이 많으면 편향이 모델-특이적(학습데이터·아키텍처 차이) = 재랭커/consensus로 해결 여지.")
        P("        '양쪽 다 편향'이 많으면 그 항원·에피토프 자체가 구조적으로 어려운 것(모델 불문) = 근본적 난제.")

    with open(os.path.join(a.out_dir, "epitope_summary.txt"), "w") as fh:
        fh.write("\n".join(out))
    P(f"\n→ 저장: {a.out_dir}/epitope_*.csv, epitope_summary.txt")

def _write(path, rows):
    if not rows: return
    cols = list(rows[0].keys())
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols); w.writeheader(); w.writerows(rows)

if __name__ == "__main__":
    main()
