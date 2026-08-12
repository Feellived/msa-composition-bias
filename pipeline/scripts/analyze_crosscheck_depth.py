#!/usr/bin/env python3
"""[교차확인] DockQ 개선 깊이(rung)와 epitope(recall/over-rep) 개선 깊이가 겹치나 — MSA깊이 가설 케이스스터디.

가설의 핵심 예측: "MSA 깊이를 줄이면 예측이 흔한(over-represented) 자리에서 벗어나 진짜 에피토프로 가고,
그래서 DockQ가 오른다." → 그렇다면 **DockQ가 최고인 깊이(rung)와 recall이 최고인 깊이가 같아야** 한다.
반대로 DockQ만 오르고 recall은 딴 깊이에서 오르면 = DockQ 개선이 에피토프 교정 때문이 아님(fold/우연).

⚠️ 이건 재랭커(pose 선택)가 아니라 **깊이-반응 케이스 스터디**다. pose_features.csv 하나로 됨(두 CSV join 불필요):
   같은 파일에 dockq·recall·overrep·dcc_true가 pose 단위로 다 있고(Boltz+Protenix), recall/overrep은
   eval_epitope_shift.py와 동일 계산(pose_features가 그걸 import). rung0 = full MSA(가장 깊음), rung↑ = 얕음.

각 (target, model)마다:
  · rung마다 best DockQ pose를 잡고, 그 pose의 recall(=rec@dqpk)·전 pose max recall·min over-rep 집계
  · full(rung0) 대비 peak: DockQ가 어느 깊이서 최고인지, recall이 어느 깊이서 최고인지
  · coincide = DockQ-peak 깊이와 recall-peak 깊이가 ±1 이내로 겹치나(가설 지지 신호)
  · responsive = dq_gain(=peak−full) 크고 peak가 얕은 깊이(rung>0)에서 발생 = "깊이 줄여 좋아진" 후보

판정(관대하되 정직):
  · dq_gain 크고 + peak가 얕은 깊이 + recall도 같은 깊이서 개선(coincide) = 유력 후보(→ seed-복제로 확인).
  · dq만 오르고 recall 안 겹침 = 에피토프 교정 아님(fold/배치 운일 수). ⚠️ 어느 쪽도 이걸로 '증명'은 아님 — seed-복제·개별순열 통제가 남음.

사용(stdlib only, 어디서나):
  python analyze_crosscheck_depth.py                         # results/pose_features.csv, 전 타깃
  python analyze_crosscheck_depth.py --focus                 # 후보만(8txu_HL·9y0a_AB·9azr_HL·RBD앵커)
  python analyze_crosscheck_depth.py --min-gain 0.15         # responsive 문턱
  python analyze_crosscheck_depth.py --targets "8txu_HL 9y0a_AB"
→ results/crosscheck_depth.csv 도 씀(seed-복제 우선순위표).
"""
import argparse, csv, math
from collections import defaultdict

CAND = ["8txu_HL", "9y0a_AB", "9azr_HL", "8wpy_AB", "8k3k_D", "8k46_I", "8y6a_CD", "8ulr_HL"]
OUTCOLS = ["target", "model", "group", "ab", "n_rung",
           "dq_full", "dq_peak", "dq_peak_rung", "dq_peak_neff", "dq_gain",
           "rec_full", "rec_peak", "rec_peak_rung", "rec_gain", "rec_at_dqpk",
           "over_full", "over_min", "over_min_rung", "over_drop",
           "coincide", "responsive"]


def f(x):
    try:
        v = float(x)
        return None if math.isnan(v) else v
    except Exception:
        return None


def load(path):
    rows = []
    for r in csv.DictReader(open(path)):
        rows.append(dict(target=r["target"], group=r.get("group", ""), ab=r.get("ab", ""),
                         model=r["model"], rung=int(float(r["rung"])), neff=f(r.get("neff80", "")),
                         dockq=f(r.get("dockq", "")), recall=f(r.get("recall", "")),
                         overrep=f(r.get("overrep", "")), dcc_true=f(r.get("dcc_true", ""))))
    return rows


def rung_stats(poses):
    by = defaultdict(list)
    for p in poses:
        by[p["rung"]].append(p)
    out = {}
    for rg, ps in by.items():
        dqs = [(p["dockq"], p) for p in ps if p["dockq"] is not None]
        if not dqs:
            continue
        best_dq, best_pose = max(dqs, key=lambda x: x[0])
        recs = [p["recall"] for p in ps if p["recall"] is not None]
        overs = [p["overrep"] for p in ps if p["overrep"] is not None]
        out[rg] = dict(neff=ps[0]["neff"], best_dq=best_dq, rec_at_dq=best_pose["recall"],
                       max_rec=max(recs) if recs else None, min_over=min(overs) if overs else None)
    return out


def peak(d, key, want_max=True):
    items = [(rg, s[key]) for rg, s in d.items() if s.get(key) is not None]
    if not items:
        return None, None, None
    rg, v = (max if want_max else min)(items, key=lambda x: x[1])
    return rg, v, d[rg]["neff"]


def summarize(poses, min_gain):
    st = rung_stats(poses)
    if 0 not in st or not st:
        return None
    dq_pk_rg, dq_pk, dq_pk_neff = peak(st, "best_dq", True)
    rec_pk_rg, rec_pk, _ = peak(st, "max_rec", True)
    ov_min_rg, ov_min, _ = peak(st, "min_over", False)
    dq_full = st[0]["best_dq"]
    rec_full = st[0]["max_rec"]
    ov_full = st[0]["min_over"]
    dq_gain = dq_pk - dq_full if dq_full is not None else None
    rec_gain = (rec_pk - rec_full) if (rec_pk is not None and rec_full is not None) else None
    ov_drop = (ov_full - ov_min) if (ov_full is not None and ov_min is not None) else None
    coincide = (rec_pk_rg is not None and dq_pk_rg is not None and abs(dq_pk_rg - rec_pk_rg) <= 1)
    responsive = (dq_gain is not None and dq_gain >= min_gain and dq_pk_rg is not None and dq_pk_rg > 0)
    return dict(n_rung=len(st), dq_full=dq_full, dq_peak=dq_pk, dq_peak_rung=dq_pk_rg, dq_peak_neff=dq_pk_neff,
                dq_gain=dq_gain, rec_full=rec_full, rec_peak=rec_pk, rec_peak_rung=rec_pk_rg, rec_gain=rec_gain,
                rec_at_dqpk=st[dq_pk_rg]["rec_at_dq"] if dq_pk_rg in st else None,
                over_full=ov_full, over_min=ov_min, over_min_rung=ov_min_rg, over_drop=ov_drop,
                coincide=coincide, responsive=responsive)


def g(v, n=3):
    return f"{v:.{n}f}" if isinstance(v, float) else ("-" if v is None else str(v))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", default="results/pose_features.csv")
    ap.add_argument("--targets", default="")
    ap.add_argument("--focus", action="store_true", help="후보만(8txu_HL·9y0a_AB·9azr_HL·RBD앵커)")
    ap.add_argument("--min-gain", type=float, default=0.15)
    ap.add_argument("--out", default="results/crosscheck_depth.csv")
    a = ap.parse_args()
    rows = load(a.csv)
    want = set(a.targets.replace(",", " ").split()) if a.targets else (set(CAND) if a.focus else None)

    grp = defaultdict(list)
    meta = {}
    for r in rows:
        grp[(r["target"], r["model"])].append(r)
        meta[(r["target"], r["model"])] = (r["group"], r["ab"])
    recs = []
    for (t, m), ps in grp.items():
        if want and t not in want:
            continue
        s = summarize(ps, a.min_gain)
        if s is None:
            continue
        gr, ab = meta[(t, m)]
        recs.append(dict(target=t, model=m, group=gr, ab=ab, **s))

    recs.sort(key=lambda r: -(r["dq_gain"] or -9))
    # 표
    print(f"== DockQ×epitope 교차확인 | {a.csv} | responsive 문턱 dq_gain≥{a.min_gain} ==")
    print("  깊이: rung0=full MSA(깊음), rung↑=얕음. coincide=DockQ-peak와 recall-peak 깊이가 ±1 이내.")
    print(f"\n  {'target':10}{'mdl':9}{'DockQ full→peak@rung(neff)':30}{'gain':>7}"
          f"   {'recall f→pk@r':16}{'r@dqpk':>7}{'coin':>6}{'resp':>6}")
    print("  " + "-" * 96)
    for r in recs:
        dqs = f"{g(r['dq_full'],2)}→{g(r['dq_peak'],2)}@r{r['dq_peak_rung']}(n{g(r['dq_peak_neff'],0)})"
        rcs = f"{g(r['rec_full'],2)}→{g(r['rec_peak'],2)}@r{r['rec_peak_rung']}"
        print(f"  {r['target']:10}{r['model']:9}{dqs:30}{g(r['dq_gain'],3):>7}   {rcs:16}"
              f"{g(r['rec_at_dqpk'],2):>7}{('YES' if r['coincide'] else '·'):>6}{('★' if r['responsive'] else '·'):>6}")

    # 요약
    resp = [r for r in recs if r["responsive"]]
    coin = [r for r in resp if r["coincide"]]
    print(f"\n[요약] responsive(★, 얕은 깊이서 DockQ +{a.min_gain}↑) = {len(resp)}/{len(recs)}")
    print(f"       그중 recall도 같은 깊이서 개선(coincide) = {len(coin)}/{len(resp)}"
          f"  ← 가설('DockQ 개선=에피토프 교정') 직접 지지 후보")
    if coin:
        print("       coincide 후보:", ", ".join(f"{r['target']}/{r['model']}" for r in coin))
    print("\n⚠️ coincide여도 '증명'은 아님 — best-of-5 운일 수 있음. 다음 = seed-복제(8txu_HL 최우선)·개별 순열.")

    import os
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=OUTCOLS)
        w.writeheader()
        for r in sorted(recs, key=lambda r: -(r["dq_gain"] or -9)):
            w.writerow({k: (round(r[k], 3) if isinstance(r[k], float) else r[k]) for k in OUTCOLS})
    print(f"\n→ {a.out} (seed-복제 우선순위표, {len(recs)}행)")


if __name__ == "__main__":
    main()
