#!/usr/bin/env python3
"""에피토프 '인기 자리 이동' 측정 — MSA 깊이 ↓ 일 때 예측 에피토프가 과대표집(인기) 자리에서 벗어나나.
왜 recall(접촉잔기 겹침 %) 하나로는 부족한가 — 문헌이 실제로 쓰는 두 축을 추가로 본다:

  ① 순위점수(rank score, DiscoTope-3.0 방식) — threshold(5Å 컷오프) 없이 보는 버전.
     각 항원 잔기를 '항체까지 최소거리'로 가까운 순~먼 순으로 줄세운 뒤, 관심 잔기집합(진짜 에피토프 or
     인기영역)이 그 줄에서 평균 몇 퍼센타일(가까운 쪽=1.0)에 있나를 봄. recall처럼 '컷오프 넘었나/안 넘었나'
     이분법이 아니라 연속값이라, 컷오프 값 하나에 결론이 흔들리는 걸 피함.
  ② 중심점 거리(DCC, distance between centers; 결합주머니 예측 평가에서 쓰는 방식) — Å 단위 실측 거리.
     '예측 접촉잔기 중심'과 '진짜 에피토프 중심'/'인기영역 중심' 사이의 3차원 거리. ⚠️ 이 거리는 같은
     pose(같은 구조) 안에서만 의미가 있음(서로 다른 pose·모델끼리 좌표를 정렬 없이 비교하면 안 됨) —
     그래서 항상 '예측접촉 ↔ X' 형태로, 한 pose 내부 비교로만 씀. dcc_true↓=진짜자리에 가까움(좋음),
     dcc_pop↑=인기자리에서 멀어짐(이탈).

인기영역(=dominant/A 영역, epitope_defs 재사용):
 - RBD  = RBM(437-508 + K417), Wuhan spike(P0DTC2) 넘버링.
 - HA   = head globular(P03437 참조 56-306).  stem/HA2 = 인기 아님.
 - Env  = 구조적 supersite(CD4bs ∪ N332 ∪ V2apex, HXB2/P04578).  FP/MPER/gp120gp41 = 인기 아님.
 - group C(대조 항원) = 인기영역 정의 없음 → over-rep/pop_rank/dcc_pop = NA(재랭커 음성대조로만).

oracle 정의(프로젝트 관행 재사용 — Phase 0의 'oracle vs ipTM-pick'과 동일 방식):
  rung 안 pose 5개 중 recall이 가장 좋은 pose **하나**를 기준점으로 고정하고, 그 pose에서 나머지 지표
  (overrep·rank·dcc)도 같이 읽음(지표마다 제각각 최댓값 pose를 따로 고르면 통계적으로 앞뒤가 안 맞음
  — '어떤 pose를 뽑았는지'가 지표마다 달라지는 cherry-pick을 막기 위함). mean_* = pose 5개 평균
  (best-of-N 착시 회피용 — 도장 판정은 이 mean 기준 추세로).

읽는 법(도장):
 - 깊이(Neff80)와 over-rep/pop_rank 가 양(+)의 상관, dcc_pop 이 음(-)의 상관 = 깊을수록 인기 자리로 쏠림.
 - B 라벨(진짜=인기 바깥): 깊이↓ 시 over-rep↓·pop_rank↓·dcc_pop↑(이탈) & recall↑(정답행) = 편향 통로 확증.
 - A 라벨(진짜=인기, 역대조): 깊이↓ 시 이탈은 하는데 recall도 같이 ↓ = 이탈이 '깊이 탓'임을 뒷받침.
 - ⚠️ 전부 rung당 pose 5개의 '평균'을 주 지표로 씀(best-of-N max 착시 회피). 추세는 per-target Spearman.

경로: pose=$DATA/<model>/<t>/rung<r>/results/**.cif · native=targets/<t>/native.cif · chains.json=targets/<t>/chains.json
사용(biopython+scipy env; classify_epitope 가 UniProt fetch):
  python eval_epitope_shift.py [--models boltz] [--rungs 12] [--cutoff 5.0] [--out results/epitope_shift.csv]
  ⚠️ 모델별로 --out 파일을 따로 두고 한 번에 한 모델씩 실행 권장(예: epitope_shift_boltz.csv·epitope_shift_protenix.csv).
     각 행에 model 컬럼이 있어 나중에 analyze_epitope_shift.py로 합쳐서 모델 간 비교 가능.
"""
import argparse, csv, glob, math, os
import numpy as np
from scipy.spatial import cKDTree
import lib_epitope_defs as E
from lib_epitope_recall import (load, best, scored_epitope, posmap, T3,
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
        m = map_to_ref(agseq, list(range(len(agseq))), ref)
        for idx, fpos in m.items():
            if fpos in popset:
                out.add((i, idx))
    return out

def scored_epitope_full(ag_chains, ab_residues, cutoff):
    """scored_epitope(epitope_recall) 확장판 — 접촉집합·거리 외에 잔기 3D 중심좌표(coord)도 반환.
    좌표는 이 pose 자신의 프레임 그대로(다른 pose·모델과 정렬 불필요 — 이 pose 내부 거리비교에만 씀)."""
    ab_coords = np.array([a.coord for r in ab_residues for a in r if a.element != "H"], dtype=float)
    if len(ab_coords) == 0: return set(), {}, {}
    tree = cKDTree(ab_coords)
    contacts, dist, coord = set(), {}, {}
    for ckey, residues, ref in ag_chains:
        seq = "".join(T3.get(r.get_resname().upper(), "X") for r in residues)
        pm = posmap(seq, ref)
        for i, res in enumerate(residues):
            if i not in pm: continue
            coords = np.array([a.coord for a in res if a.element != "H"], dtype=float)
            if len(coords) == 0: continue
            md = float(tree.query(coords, k=1)[0].min())
            key = (ckey, pm[i])
            dist[key] = md
            coord[key] = coords.mean(axis=0)
            if md <= cutoff: contacts.add(key)
    return contacts, dist, coord

def rank_score(dist, region):
    """DiscoTope식 '에피토프 순위점수' — 항원 잔기를 항체까지 거리(가까운→먼) 순으로 줄세워, region이
    그 줄에서 평균 몇 퍼센타일(1.0=제일 가까움/접촉스러움, 0.0=제일 멂)에 있나. 컷오프(threshold) 불필요."""
    keys = list(dist.keys())
    if not keys or not region: return float("nan")
    order = sorted(keys, key=lambda k: dist[k])
    n = len(order)
    pct_top = {}
    i = 0
    while i < n:
        j = i
        while j + 1 < n and dist[order[j + 1]] == dist[order[i]]: j += 1
        avg_rank = (i + j) / 2.0
        for k in range(i, j + 1): pct_top[order[k]] = (1.0 - avg_rank / (n - 1)) if n > 1 else 1.0
        i = j + 1
    vals = [pct_top[k] for k in region if k in pct_top]
    return sum(vals) / len(vals) if vals else float("nan")

def centroid_dist(coord, region_a, region_b):
    """DCC(중심점간 거리)식 — region_a 중심좌표 ↔ region_b 중심좌표 사이 거리(Å). 같은 pose 프레임 안에서만
    비교(다른 pose·모델 간 절대 좌표 비교는 무의미 — 항상 이 pose의 '예측접촉 ↔ X' 형태로만 씀)."""
    ca = [coord[k] for k in region_a if k in coord]
    cb = [coord[k] for k in region_b if k in coord]
    if not ca or not cb: return float("nan")
    a = np.mean(ca, axis=0); b = np.mean(cb, axis=0)
    return float(np.linalg.norm(a - b))

def frac(pred, region):
    return len(pred & region) / len(pred) if pred else float("nan")

def nanmean(xs):
    xs = [x for x in xs if not (isinstance(x, float) and math.isnan(x))]
    return sum(xs) / len(xs) if xs else float("nan")

def fmt(x):
    return "" if (x is None or x != x) else round(x, 3)     # x!=x ⇔ NaN

def pose_metrics_full(cj, pose_path, cutoff, true, pop):
    """pose 하나 채점 — recall/over-rep(기존, 컷오프 기반) + 순위점수(threshold-free) + 중심거리(DCC, Å)."""
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
    pred, dist, coord = scored_epitope_full(ag, ab, cutoff)
    if not pred: return None
    return dict(
        recall=frac(pred, true),
        overrep=frac(pred, pop) if pop is not None else float("nan"),
        true_rank=rank_score(dist, true),
        pop_rank=rank_score(dist, pop) if pop is not None else float("nan"),
        dcc_true=centroid_dist(coord, pred, true),
        dcc_pop=centroid_dist(coord, pred, pop) if pop is not None else float("nan"),
    )

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
        native_or = frac(true, pop) if pop is not None else float("nan")
        n_pop = len(pop) if pop is not None else 0
        nmap = neff_of(tgt, os.path.join(a.data, "ladders"))
        for model in a.models:
            for rung in range(a.rungs):
                poses = glob.glob(os.path.join(a.data, model, tgt, f"rung{rung}", "results", "**", "*.cif"), recursive=True)
                if not poses: continue
                ms = []
                for pose in poses:
                    try:
                        mm = pose_metrics_full(cj, pose, a.cutoff, true, pop)
                    except Exception:
                        continue
                    if mm: ms.append(mm)
                if not ms: continue
                oracle = max(ms, key=lambda z: z["recall"])       # recall 기준 고정 pose (cherry-pick 방지)
                row = dict(target=tgt, model=model, family=fam, ab=ab, label=r.get("label", ""), rung=rung,
                          neff80=nmap.get(rung, ""), n_pose=len(ms), n_pop=n_pop,
                          native_overrep=fmt(native_or),
                          mean_overrep=fmt(nanmean([m["overrep"] for m in ms])),
                          oracle_overrep=fmt(oracle["overrep"]),
                          mean_recall=fmt(nanmean([m["recall"] for m in ms])),
                          oracle_recall=fmt(oracle["recall"]),
                          mean_true_rank=fmt(nanmean([m["true_rank"] for m in ms])),
                          oracle_true_rank=fmt(oracle["true_rank"]),
                          mean_pop_rank=fmt(nanmean([m["pop_rank"] for m in ms])),
                          oracle_pop_rank=fmt(oracle["pop_rank"]),
                          mean_dcc_true=fmt(nanmean([m["dcc_true"] for m in ms])),
                          oracle_dcc_true=fmt(oracle["dcc_true"]),
                          mean_dcc_pop=fmt(nanmean([m["dcc_pop"] for m in ms])),
                          oracle_dcc_pop=fmt(oracle["dcc_pop"]))
                rows.append(row)
                print(f"  {tgt:12} {fam:4}{ab:2} r{rung:<2} Neff80={nmap.get(rung,'?'):>7} "
                      f"overrep(mean)={row['mean_overrep']} recall(mean)={row['mean_recall']} "
                      f"pop_rank(mean)={row['mean_pop_rank']} dcc_pop(mean)={row['mean_dcc_pop']} (n_pop={n_pop})")
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    cols = ["target","model","family","ab","label","rung","neff80","n_pose","n_pop",
            "native_overrep","mean_overrep","oracle_overrep","mean_recall","oracle_recall",
            "mean_true_rank","oracle_true_rank","mean_pop_rank","oracle_pop_rank",
            "mean_dcc_true","oracle_dcc_true","mean_dcc_pop","oracle_dcc_pop"]
    with open(a.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols); w.writeheader(); w.writerows(rows)

    # ── 요약(도장): family×ab 별 full→low 변화 + per-target Spearman(neff, ·) ──
    from collections import defaultdict
    byt = defaultdict(list)
    for x in rows: byt[x["target"]].append(x)

    def series(xs, key):
        return [float(z[key]) for z in xs if z[key] != ""]

    summ = defaultdict(list)   # (family, ab) -> list of per-target dicts
    for tgt, xs in byt.items():
        xs = [x for x in xs if x["neff80"] != "" and x["mean_overrep"] != ""]
        if len(xs) < 4: continue
        xs.sort(key=lambda z: float(z["neff80"]))                # 낮은 깊이 → 높은 깊이
        neff = [float(z["neff80"]) for z in xs]
        d = {}
        for key in ("mean_overrep", "mean_recall", "mean_pop_rank", "mean_dcc_pop"):
            vals = series(xs, key)
            if len(vals) != len(xs): d[key] = (float("nan"), float("nan")); continue
            d[key] = (vals[0] - vals[-1], spearman(neff, vals))   # (low-full, rho)
        summ[(xs[0]["family"], xs[0]["ab"])].append(d)

    print("\n" + "=" * 100)
    print("요약 — 깊이(Neff80) full→최저 변화(low−full) & per-target Spearman(neff, ·)")
    print("  (over-rep/pop_rank: ↓ 저깊이=이탈. dcc_pop: ↑ 저깊이=이탈(거리라 부호 반대). recall: ↑ 저깊이=정답행.)")
    print(f"{'family/ab':10}{'n':>3} | over-rep Δ,ρ        | pop_rank Δ,ρ        | dcc_pop Δ,ρ          | recall Δ,ρ")
    for key in sorted(summ):
        v = summ[key]; n = len(v)
        def col(k):
            ds = [z[k][0] for z in v if not math.isnan(z[k][0])]
            rs = [z[k][1] for z in v if not math.isnan(z[k][1])]
            if not ds: return "NA"
            pos = sum(1 for x in rs if x > 0)
            return f"{np.mean(ds):+.3f},{np.mean(rs):+.2f}[{pos}/{len(rs)}]" if rs else f"{np.mean(ds):+.3f},NA"
        print(f"{key[0]+'/'+key[1]:10}{n:>3} | {col('mean_overrep'):20} | {col('mean_pop_rank'):20} | "
              f"{col('mean_dcc_pop'):20} | {col('mean_recall')}")
    print("\n도장 판정:")
    print("  B/* 행: over-rep↓·pop_rank↓·dcc_pop↑ (깊이↓ 시 인기자리 이탈, 3개 지표 다 같은 방향 가리킴)")
    print("         + recall↑(진짜 자리로) = 편향 통로 확증 (threshold 하나에 안 걸림 — 3중 확인).")
    print("  A/* 행: 이탈은 하는데 recall도 같이 ↓(역대조) = 이탈이 '깊이 탓'임을 뒷받침.")
    print("  rho(neff,·) 다수 해당 부호 = 깊을수록 인기자리 쏠림. C/* = 인기영역 없음(NA).")
    print(f"\n→ {a.out} ({len(rows)}행). true_rank/dcc_true = recall의 threshold-free 짝(QC용, 요약 미포함).")
    print("  다음: mean_pop_rank·mean_dcc_pop vs log(Neff) 그림 + 순열/짝-Wilcoxon 로 추세 유의성.")

if __name__ == "__main__":
    main()
