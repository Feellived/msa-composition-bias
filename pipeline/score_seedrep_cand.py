#!/usr/bin/env python3
"""[채점] seed-복제 케이스스터디 — depth-rescue가 재현되나 vs best-of-5 단발 운.

각 후보 (model, target)에 대해 peak 깊이 × N seed 예측을 DockQ+recall로 채점하고,
seed마다 best(5샘플 중)를 잡아 **N seed 분포**를 관측 peak(cross-check의 obs_dq)과 대조한다.

판정(관대하되 정직):
  · dq_mean이 obs 근처 + sd 작음  = 그 깊이서 rescue가 **재현** = 깊이(개수)에 강건 = 진짜 신호.
  · sd 큼(seed마다 딴판)          = **조성(어느 서열이 뽑혔나)** 이 원인 → 여전히 실재하나 조성 의존(Exp2 nested).
  · 전부 낮음(mean≪obs)           = 관측 peak가 best-of-5 **단발 운** = 그 후보는 기각.
(seed_replicate.py 설계 그대로 — 개수 고정·조성만 흔들어 '개수 vs 조성' 판별.)

DockQ env(biopython + DockQ + scipy)에서:
  python score_seedrep_cand.py
"""
import argparse, csv, glob, json, os, tempfile
import statistics as st
import pose_features as PF   # 채점 로직 재사용(native_merged·pose_merged·dockq·pose_all_metrics·native_true·ES)


def load_group(list_csv):
    g = {}
    for r in csv.DictReader(open(list_csv)):
        g[r["target"]] = r.get("group", "")
    return g


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cand", default="seedrep_cand.csv")
    ap.add_argument("--list", default="sweep_targets.csv")
    ap.add_argument("--data", default=os.environ.get("DATA", "/mnt/data/admuser/msadepth"))
    ap.add_argument("--targets-dir", default="targets")
    ap.add_argument("--cutoff", type=float, default=5.0)
    ap.add_argument("--out", default="results/seedrep_cand_scored.csv")
    a = ap.parse_args()
    grp = load_group(a.list)
    outroot = os.path.join(a.data, "seedrep_cand")

    print(f"{'target':10}{'model':9}{'depth':>7}{'nseed':>6}  "
          f"{'DockQ mean±sd (min~max)':26}{'≥.49':>6}  {'recall m':>9}   obs_dq")
    print("-" * 92)
    rows = []
    for r in csv.DictReader(open(a.cand)):
        t, model = r["target"], r["model"]
        obs = r.get("obs_dq", "")
        cjp = os.path.join(a.targets_dir, t, "chains.json")
        native = os.path.join(a.targets_dir, t, "native.cif")
        if not os.path.exists(cjp):
            print(f"{t:10}{model:9}  chains.json 없음 skip"); continue
        cj = json.load(open(cjp))
        tr = PF.native_true(cj, native, a.cutoff)
        if tr is None:
            print(f"{t:10}{model:9}  native epitope 실패 skip"); continue
        true, _ = tr
        try:
            pop = PF.ES.popular_refset(cj, grp.get(t, ""))
        except Exception:
            pop = None
        with tempfile.TemporaryDirectory() as td:
            natm = PF.native_merged(cj, native, td)
            if natm is None:
                print(f"{t:10}{model:9}  native merge 실패 skip"); continue
            base = os.path.join(outroot, model, t)
            best_dq, best_rec = [], []
            depth_lbl = ""
            for ddir in sorted(glob.glob(os.path.join(base, "d*"))):
                depth_lbl = os.path.basename(ddir)
                for sdir in sorted(glob.glob(os.path.join(ddir, "seed*"))):
                    poses = glob.glob(os.path.join(sdir, "results", "**", "*.cif"), recursive=True)
                    bdq = brec = None
                    for p in poses:
                        try:
                            met = PF.pose_all_metrics(cj, p, a.cutoff, true, pop)
                            pm = PF.pose_merged(cj, p, td)
                            q = PF.dockq(pm, natm) if pm else None
                        except Exception:
                            continue
                        if q is not None and (bdq is None or q > bdq):
                            bdq = q
                        if met and (brec is None or met["recall"] > brec):
                            brec = met["recall"]
                    if bdq is not None:
                        best_dq.append(bdq)
                    if brec is not None:
                        best_rec.append(brec)
            n = len(best_dq)
            if n == 0:
                print(f"{t:10}{model:9}{depth_lbl:>7}   예측 없음(미실행?) skip"); continue
            m = sum(best_dq) / n
            sd = st.pstdev(best_dq) if n > 1 else 0.0
            mn, mx = min(best_dq), max(best_dq)
            ge = sum(1 for x in best_dq if x >= 0.49)
            rm = (sum(best_rec) / len(best_rec)) if best_rec else float("nan")
            dqs = f"{m:.2f}±{sd:.2f} ({mn:.2f}~{mx:.2f})"
            print(f"{t:10}{model:9}{depth_lbl:>7}{n:>6}  {dqs:26}{ge}/{n:<4}  {rm:>9.2f}   {obs}")
            rows.append(dict(target=t, model=model, depth=depth_lbl, n_seed=n,
                             dq_mean=round(m, 3), dq_sd=round(sd, 3), dq_min=round(mn, 3), dq_max=round(mx, 3),
                             n_ge49=ge, rec_mean=round(rm, 3), obs_dq=obs))

    if rows:
        os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
        with open(a.out, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            w.writeheader(); w.writerows(rows)
        print(f"\n→ {a.out}")
    print("\n판정: dq_mean≈obs·sd작음=재현(깊이 강건, 진짜) / sd큼=조성 의존(Exp2) / dq_mean≪obs=peak가 운(기각).")


if __name__ == "__main__":
    main()
