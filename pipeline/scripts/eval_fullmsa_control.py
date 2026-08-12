#!/usr/bin/env python3
"""[통제 채점] full MSA에서 자세를 많이 뽑았을 때도 정답이 나오나 — 얕은 깊이와 자세 단위로 맞대결.

판정 질문: 얕은 깊이에서 나온 성공이 '깊이가 연 것'인가, '자세를 많이 뽑아서 걸린 것'인가.
  · full MSA는 조성이 하나뿐이므로, 여기서 자세를 40개 뽑아도 성공이 없으면 → 깊이가 연 것(지지).
  · full MSA에서도 비슷한 비율로 성공이 나오면 → 표본 수 문제(기각).

두 쪽 다 '자세 단위 성공 개수 / 전체 자세 수'로 비교하고, 개수가 달라도 되도록 Fisher 정확검정
(단측)으로 우연 확률을 같이 낸다. 표본이 작으니 p는 참고용 — 방향과 크기를 먼저 본다.

  얕은 깊이 자료 = results/seedrep_poses.csv (eval_dump_seedrep.py 출력)
  통제 자료      = $DATA/fullmsa_ctl/<model>/<target>/full_n*/results/**/*.cif

사용(DockQ env):
  python eval_fullmsa_control.py
  python eval_fullmsa_control.py --only 8ulr_HL
"""
import argparse, csv, glob, json, math, os, tempfile
from collections import defaultdict
import lib_pose_features as PF

SUCC = 0.49


def fisher_one_sided(a, b, c, d):
    """[[a,b],[c,d]] 에서 1행의 성공비율이 더 높을 단측 p (초기하 꼬리합)."""
    n1, n2, k = a + b, c + d, a + c
    if n1 == 0 or n2 == 0:
        return float("nan")
    tot = math.comb(n1 + n2, k)
    if tot == 0:
        return float("nan")
    p = 0
    for x in range(a, min(n1, k) + 1):
        p += math.comb(n1, x) * math.comb(n2, k - x)
    return p / tot


def load_seedrep(path):
    """target,model -> (성공 자세 수, 전체 자세 수)"""
    out = defaultdict(lambda: [0, 0])
    if not os.path.exists(path):
        return out
    for r in csv.DictReader(open(path)):
        try:
            q = float(r["dockq"])
        except Exception:
            continue
        k = (r["target"], r["model"])
        out[k][1] += 1
        if q >= SUCC:
            out[k][0] += 1
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cand", default="seedrep_cand.csv")
    ap.add_argument("--list", default="sweep_targets.csv")
    ap.add_argument("--pf", default="results/pose_features.csv")
    ap.add_argument("--seedrep", default="results/seedrep_poses.csv")
    ap.add_argument("--data", default=os.environ.get("DATA", "/mnt/data/admuser/msadepth"))
    ap.add_argument("--targets-dir", default="targets")
    ap.add_argument("--cutoff", type=float, default=5.0)
    ap.add_argument("--only", default="")
    ap.add_argument("--out", default="results/fullmsa_ctl_scored.csv")
    a = ap.parse_args()

    grp = {r["target"]: r.get("group", "") for r in csv.DictReader(open(a.list))}
    sr = load_seedrep(a.seedrep)
    # full(rung0) best-of-5 참고값
    full5 = {}
    if os.path.exists(a.pf):
        tmp = defaultdict(list)
        for r in csv.DictReader(open(a.pf)):
            if int(float(r["rung"])) == 0:
                try:
                    tmp[(r["target"], r["model"])].append(float(r["dockq"]))
                except Exception:
                    pass
        full5 = {k: max(v) for k, v in tmp.items() if v}

    outroot = os.path.join(a.data, "fullmsa_ctl")
    rows = []
    print(f"성공 기준 DockQ ≥ {SUCC}\n")
    print(f"{'target':10}{'model':9}{'[통제] full MSA 자세':>22}{'[비교] 얕은 깊이':>20}"
          f"{'full best5':>12}{'단측 p':>9}  판정")
    print("-" * 104)

    for r in csv.DictReader(open(a.cand)):
        t, model = r["target"], r["model"]
        if a.only and t != a.only:
            continue
        base = os.path.join(outroot, model, t)
        poses = sorted(glob.glob(os.path.join(base, "full_n*", "results", "**", "*.cif"), recursive=True))
        if not poses:
            print(f"{t:10}{model:9}  통제 예측 없음(미실행) — run_fullmsa_control.sh 필요")
            continue
        cjp = os.path.join(a.targets_dir, t, "chains.json")
        native = os.path.join(a.targets_dir, t, "native.cif")
        if not os.path.exists(cjp):
            print(f"{t:10}{model:9}  chains.json 없음 skip")
            continue
        cj = json.load(open(cjp))
        tr = PF.native_true(cj, native, a.cutoff)
        if tr is None:
            print(f"{t:10}{model:9}  native epitope 실패 skip")
            continue
        true, _ = tr
        try:
            pop = PF.ES.popular_refset(cj, grp.get(t, ""))
        except Exception:
            pop = None

        dqs, recs, hi_rec = [], [], []
        with tempfile.TemporaryDirectory() as td:
            natm = PF.native_merged(cj, native, td)
            if natm is None:
                print(f"{t:10}{model:9}  native merge 실패 skip")
                continue
            for p in poses:
                try:
                    met = PF.pose_all_metrics(cj, p, a.cutoff, true, pop)
                    pm = PF.pose_merged(cj, p, td)
                    q = PF.dockq(pm, natm) if pm else None
                except Exception:
                    continue
                if q is None:
                    continue
                dqs.append(q)
                if met:
                    recs.append(met["recall"])
                    if q >= SUCC:
                        hi_rec.append(met["recall"])

        if not dqs:
            print(f"{t:10}{model:9}  채점 실패(자세 {len(poses)}개)")
            continue
        c_ok, c_n = sum(1 for x in dqs if x >= SUCC), len(dqs)
        s_ok, s_n = sr.get((t, model), [0, 0])
        p = fisher_one_sided(s_ok, s_n - s_ok, c_ok, c_n - c_ok) if s_n else float("nan")
        f5 = full5.get((t, model))

        if s_n == 0:
            verdict = "얕은 깊이 자료 없음(dump 먼저)"
        elif c_ok == 0 and s_ok > 0:
            verdict = "깊이가 연 것 — 지지"
        elif c_ok / c_n >= s_ok / s_n:
            verdict = "표본 수 문제 — 기각"
        else:
            verdict = "통제서도 일부 나옴 — 약화"
        ctl = f"{c_ok}/{c_n} (최고 {max(dqs):.2f})"
        cmp_ = f"{s_ok}/{s_n}" if s_n else "-"
        print(f"{t:10}{model:9}{ctl:>22}{cmp_:>20}"
              f"{(f'{f5:.3f}' if f5 is not None else '-'):>12}{p:>9.3f}  {verdict}")
        rows.append(dict(target=t, model=model, ctl_ok=c_ok, ctl_n=c_n,
                         ctl_best=round(max(dqs), 3),
                         ctl_rec_at_succ=(round(sum(hi_rec) / len(hi_rec), 3) if hi_rec else ""),
                         seedrep_ok=s_ok, seedrep_n=s_n,
                         full_best5=(round(f5, 3) if f5 is not None else ""),
                         p_one_sided=(round(p, 4) if p == p else ""), verdict=verdict))

    if rows:
        os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
        with open(a.out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            w.writeheader(); w.writerows(rows)
        print(f"\n→ {a.out}")
    print("\n※ 자세 수가 적으면 p는 쉽게 커진다(예: 3/40 vs 0/40 이면 p≈0.24). 방향·크기를 먼저 보고 p는 참고로.")


if __name__ == "__main__":
    main()
