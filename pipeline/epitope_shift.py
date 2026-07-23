#!/usr/bin/env python3
"""에피토프 '인기 자리 이동' 측정 — MSA 깊이 ↓ 일 때 예측 에피토프가 과대표집(인기) 자리에서 벗어나나.
왜 recall만으론 부족한가: recall = '진짜 에피토프에 가까워졌나'. 하지만 '인기 자리에서 벗어남'과 '진짜 자리로 감'은
다른 사건 → 주 지표 = over-rep overlap(예측 접촉잔기 중 인기영역 비율), 보조 = recall(이탈이 옳은 방향인지 판별).

인기영역(=dominant/A 영역, epitope_defs 재사용):
 - RBD  = RBM(437-508 + K417), Wuhan spike(P0DTC2) 넘버링. 항원 ref서열을 P0DTC2에 정렬해 위치 매핑.
 - HA   = head globular(P03437 참조 56-306).  stem/HA2 = 인기 아님.
 - Env  = 구조적 supersite(CD4bs ∪ N332 ∪ V2apex, HXB2/P04578).  FP/MPER/gp120gp41 = 인기 아님.
 - group C(대조 항원: CD38·B7-H3·PD-1 등) = 인기영역 정의 없음 → over-rep = NA(재랭커 음성대조로만).

읽는 법(도장):
 - 깊이(Neff80)와 over-rep overlap 이 양(+)의 상관 = 깊을수록 인기 자리로 쏠림 = 'MSA 깊이 = 편향 통로'.
 - B 라벨(진짜=인기 바깥): 깊이↓ 시 over-rep↓(이탈) & recall↑(정답행) = 진짜 편향 통로.
 - A 라벨(진짜=인기): 깊이↓ 시 over-rep↓ & recall↓(역대조) = 이탈이 손해 → 이탈이 '깊이 탓'임을 반증.
 - ⚠️ over-rep/recall 은 rung당 pose 5개의 '평균'을 씀(best-of-N max 착시 회피). 추세는 per-target Spearman.

경로: pose=$DATA/<model>/<t>/rung<r>/results/**.cif · native=targets/<t>/native.cif · chains.json=targets/<t>/chains.json
사용(biopython+scipy env; classify_epitope 가 UniProt fetch):
  python epitope_shift.py [--models boltz] [--rungs 12] [--cutoff 5.0] [--out results/epitope_shift.csv]
"""
import argparse, csv, glob, math, os
import numpy as np
import epitope_defs as E
from epitope_recall import (load, best, scored_epitope,
                            antigen_refs, antibody_refs, native_true, neff_of)
from classify_epitope import fetch_ref, map_to_ref

# ── 인기영역(과대표집) 정의: family -> (UniProt ref acc, popular position set in ref numbering) ──
HA_HEAD = set(range(E.HA_HEAD_REFRANGE[0], E.HA_HEAD_REFRANGE[1] + 1))          # P03437 56-306
ENV_POP = E.ENV_CLASSES["CD4bs"] | E.ENV_CLASSES["N332"] | E.ENV_CLASSES["V2apex"]
FAMILY_POP = {
    "RBD": ("P0DTC2", E.RBD_RBM),    # 全spike 넘버링 = RBD_RBM(437-508+417) 그대로
    "HA":  ("P03437", HA_HEAD),
    "Env": ("P04578", ENV_POP),
}
_REFCACHE = {}

def popular_refset(cj, family):
    """항원 사슬별 chains.json ref서열 index 중 '인기영역'에 매핑되는 (i, idx) 집합. C(미지원)면 None.
    좌표계 = native_true/scored_epitope 와 동일(항원사슬 index i, ref서열 0-based index)."""
    if family not in FAMILY_POP:
        return None
    acc, popset = FAMILY_POP[family]
    if acc not in _REFCACHE:
        _REFCACHE[acc] = fetch_ref(acc)
    ref = _REFCACHE[acc]
    out = set()
    for i, agseq in enumerate(antigen_refs(cj)):
        # map_to_ref(qseq, qids, ref) -> {qid: ref_pos(1-based)}. qids=0..len → 키=ref서열 index.
        m = map_to_ref(agseq, list(range(len(agseq))), ref)
        for idx, fpos in m.items():
            if fpos in popset:
                out.add((i, idx))
    return out

def pose_pred(cj, pose_path, cutoff):
    """pose에서 예측 에피토프 접촉잔기 집합(i, ref_idx). 실패 None. (epitope_recall.pose_metrics 와 동일 배정)"""
    m = load(pose_path); used = set(); ag = []
    for i, ref in enumerate(antigen_refs(cj)):
        cid, _, rr = best(m, ref, exclude=used)
        if cid is None: continue
        used.add(cid); ag.append((i, rr, ref))
    ab = []
    for ref in antibody_refs(cj):
        cid, _, rr = best(m, ref, exclude=used)
        if cid: used.add(cid); ab.extend(rr)
    if not ag or not ab: return None
    pred, _ = scored_epitope(ag, ab, cutoff)
    return pred

def frac(pred, region):
    return len(pred & region) / len(pred) if pred else float("nan")

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
    rx, ry = rank(x), rank(y); mx = sum(rx) / n; my = sum(ry) / n
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    dx = math.sqrt(sum((a - mx) ** 2 for a in rx)); dy = math.sqrt(sum((b - my) ** 2 for b in ry))
    return num / (dx * dy) if dx * dy else float("nan")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--list", default="sweep_targets.csv")
    ap.add_argument("--targets-dir", default="targets")
    ap.add_argument("--data", default=os.environ.get("DATA", "/mnt/data/admuser/msadepth"))
    ap.add_argument("--models", nargs="+", default=["boltz"])
    ap.add_argument("--rungs", type=int, default=12)
    ap.add_argument("--cutoff", type=float, default=5.0)
    ap.add_argument("--out", default="results/epitope_shift.csv")
    a = ap.parse_args()
    rows = []
    for r in csv.DictReader(open(a.list)):
        tgt = r["target"]; fam = r["group"]; ab = r.get("ab", "")
        cjp = os.path.join(a.targets_dir, tgt, "chains.json")
        native = os.path.join(a.targets_dir, tgt, "native.cif")
        if not os.path.exists(cjp): continue
        import json; cj = json.load(open(cjp))
        tr = native_true(cj, native, a.cutoff)
        if tr is None: print(f"{tgt}: native epitope 실패"); continue
        true, _ = tr
        try:
            pop = popular_refset(cj, fam)
        except Exception as e:
            print(f"{tgt}: popular_refset 실패({type(e).__name__}) → over-rep NA"); pop = None
        native_or = frac(true, pop) if pop is not None else float("nan")   # 진짜 에피토프의 인기겹침(A≈높음/B≈낮음)
        n_pop = len(pop) if pop is not None else 0
        nmap = neff_of(tgt, os.path.join(a.data, "ladders"))
        for model in a.models:
            for rung in range(a.rungs):
                poses = glob.glob(os.path.join(a.data, model, tgt, f"rung{rung}", "results", "**", "*.cif"), recursive=True)
                if not poses: continue
                ors, recs, brec, bor = [], [], -1.0, float("nan")
                for pose in poses:
                    try:
                        pred = pose_pred(cj, pose, a.cutoff)
                    except Exception:
                        continue
                    if not pred: continue
                    rc = frac(pred, true); recs.append(rc)
                    if pop is not None: ors.append(frac(pred, pop))
                    if rc > brec: brec, bor = rc, (frac(pred, pop) if pop is not None else float("nan"))
                if not recs: continue
                mean_or = (sum(o for o in ors if not math.isnan(o)) / len([o for o in ors if not math.isnan(o)])) if ors else float("nan")
                rows.append(dict(target=tgt, family=fam, ab=ab, label=r.get("label", ""), rung=rung,
                                 neff80=nmap.get(rung, ""), n_pose=len(recs), n_pop=n_pop,
                                 native_overrep=round(native_or, 3) if not math.isnan(native_or) else "",
                                 mean_overrep=round(mean_or, 3) if not math.isnan(mean_or) else "",
                                 best_overrep=round(bor, 3) if not math.isnan(bor) else "",
                                 mean_recall=round(sum(recs) / len(recs), 3), best_recall=round(max(recs), 3)))
                print(f"  {tgt:12} {fam:4}{ab:2} r{rung:<2} Neff80={nmap.get(rung,'?'):>7} "
                      f"over-rep(mean)={mean_or:.2f} recall(mean)={sum(recs)/len(recs):.2f} (n_pop={n_pop})")
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    cols = ["target","family","ab","label","rung","neff80","n_pose","n_pop",
            "native_overrep","mean_overrep","best_overrep","mean_recall","best_recall"]
    with open(a.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols); w.writeheader(); w.writerows(rows)

    # ── 요약(도장): family×ab 별 full→low 변화 + per-target Spearman(neff, over-rep) ──
    from collections import defaultdict
    byt = defaultdict(list)
    for x in rows: byt[x["target"]].append(x)
    summ = defaultdict(list)   # (family, ab) -> [(d_or, d_rec, rho_or)]
    for tgt, xs in byt.items():
        xs = [x for x in xs if x["neff80"] != "" and x["mean_overrep"] != ""]
        if len(xs) < 4: continue
        xs.sort(key=lambda z: float(z["neff80"]))                # 낮은 깊이 → 높은 깊이
        neff = [float(z["neff80"]) for z in xs]
        orr = [float(z["mean_overrep"]) for z in xs]
        rec = [float(z["mean_recall"]) for z in xs]
        d_or = orr[0] - orr[-1]      # low - full : 음수면 깊이↓ 시 인기겹침 감소(이탈)
        d_rec = rec[0] - rec[-1]     # low - full : 양수면 깊이↓ 시 recall 증가
        rho = spearman(neff, orr)    # +면 깊을수록 인기겹침↑ = 편향통로
        summ[(xs[0]["family"], xs[0]["ab"])].append((d_or, d_rec, rho))
    print("\n" + "=" * 92)
    print("요약 — 깊이(Neff80) full→최저 변화(low−full) & per-target Spearman(neff, over-rep)")
    print(f"{'family/ab':10}{'n':>3} | over-rep Δ(low−full) | recall Δ(low−full) | rho(neff,overrep) [+개수]")
    for key in sorted(summ):
        v = summ[key]; n = len(v)
        dor = np.mean([z[0] for z in v]); drec = np.mean([z[1] for z in v])
        rhos = [z[2] for z in v if not math.isnan(z[2])]
        mrho = np.mean(rhos) if rhos else float("nan"); pos = sum(1 for z in rhos if z > 0)
        print(f"{key[0]+'/'+key[1]:10}{n:>3} |   {dor:+.3f}            |   {drec:+.3f}          |   {mrho:+.2f}  [{pos}/{len(rhos)}]")
    print("\n도장 판정:")
    print("  B/* 행: over-rep Δ < 0 (깊이↓ 시 인기자리 이탈) 이고 recall Δ > 0 (진짜 자리로) = 편향 통로 확증.")
    print("  A/* 행: over-rep Δ < 0 이지만 recall Δ < 0 (역대조) = 이탈이 '깊이 탓'임을 뒷받침.")
    print("  rho(neff,overrep) 다수 + = 깊을수록 인기자리 쏠림. C/* = 인기영역 없음(over-rep NA).")
    print(f"\n→ {a.out} ({len(rows)}행). 다음: mean_overrep vs log(Neff) 그림 + 순열/짝-Wilcoxon 로 추세 유의성.")

if __name__ == "__main__":
    main()
