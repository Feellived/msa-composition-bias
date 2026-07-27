#!/usr/bin/env python3
"""[전체 덤프] seed-복제 결과를 요약이 아니라 자세 하나하나까지 전부 출력.

score_seedrep_cand.py는 평균±표준편차만 보여줘서 '어떤 seed가 성공하고 어떤 seed가
바닥인지'(분포 모양)를 가린다. 여기서는 세 층을 나란히 찍는다:

  ① full(rung0)     — 원래 MSA 전부. 조성이 하나뿐이라 비교의 기준선.
  ② ladder(peak)    — 우리가 처음 관측한 그 깊이. ⚠️ 이것도 무작위 추첨 1회다
                       (neff_ladder.py:82 = seed_replicate.py:39, 같은 추첨식).
  ③ seed 복제 ×N     — 같은 깊이에서 조성만 다시 뽑은 것. 즉 ②와 같은 모집단의 표본들.

따라서 ②는 ③의 '0번째 표본'으로 봐야 하며, "②가 ③ 분포 안에 있나"가 아니라
**"③ 분포 자체가 ①보다 위에 있나"**가 판정 질문이다.

표본예산 맞춤(matched-N): seed 하나 = 자세 5개 중 최고(best-of-5). full도 자세 5개 중 최고.
그래서 seed별 best-of-5 값들(n=N)을 full의 best-of-5 값(n=1)과 견주면 예산이 맞는다.
자세 단위 성공률도 같이 찍어 best-of-N 부풀림을 확인한다.

사용(DockQ env):
  python dump_seedrep_full.py                      # 전체
  python dump_seedrep_full.py --only 8ulr_HL       # 한 후보만
  python dump_seedrep_full.py --csv-out results/seedrep_poses.csv
"""
import argparse, csv, glob, json, math, os, tempfile
import statistics as st
import pose_features as PF

SUCC = 0.49


def f(x):
    try:
        v = float(x)
        return None if math.isnan(v) else v
    except Exception:
        return None


def g(v, n=3):
    return f"{v:.{n}f}" if isinstance(v, float) else "  -  "


def load_pf(path):
    """pose_features.csv → (target, model, rung) -> [ {pose,dockq,recall}, ... ]"""
    out = {}
    if not os.path.exists(path):
        return out
    for r in csv.DictReader(open(path)):
        key = (r["target"], r["model"], int(float(r["rung"])))
        out.setdefault(key, []).append(
            dict(pose=r.get("pose", ""), dockq=f(r.get("dockq")), recall=f(r.get("recall")),
                 neff=f(r.get("neff80"))))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cand", default="seedrep_cand.csv")
    ap.add_argument("--list", default="sweep_targets.csv")
    ap.add_argument("--pf", default="results/pose_features.csv")
    ap.add_argument("--data", default=os.environ.get("DATA", "/mnt/data/admuser/msadepth"))
    ap.add_argument("--targets-dir", default="targets")
    ap.add_argument("--cutoff", type=float, default=5.0)
    ap.add_argument("--only", default="")
    ap.add_argument("--csv-out", default="results/seedrep_poses.csv")
    a = ap.parse_args()

    grp = {r["target"]: r.get("group", "") for r in csv.DictReader(open(a.list))}
    pf = load_pf(a.pf)
    outroot = os.path.join(a.data, "seedrep_cand")
    allrows = []

    for r in csv.DictReader(open(a.cand)):
        t, model, peak = r["target"], r["model"], int(r["peak_rung"])
        if a.only and t != a.only:
            continue
        obs = r.get("obs_dq", "")
        print("=" * 100)
        print(f"■ {t} · {model}   (관측 peak = 깊이단계 {peak}, obs_dq={obs})")
        print("=" * 100)

        # ── ① full(rung0) ─────────────────────────────────────────────
        full = pf.get((t, model, 0), [])
        print(f"\n[① full MSA · 깊이단계 0]  조성 1가지(전부) · 자세 {len(full)}개")
        if full:
            print("     " + "  ".join(f"{g(p['dockq'],2)}" for p in full)
                  + f"   → best {g(max((p['dockq'] for p in full if p['dockq'] is not None), default=None),3)}")
            print("     recall: " + "  ".join(f"{g(p['recall'],2)}" for p in full))
        else:
            print("     (pose_features.csv에 없음)")

        # ── ② ladder(peak rung) = 무작위 추첨 1회 ──────────────────────
        lad = pf.get((t, model, peak), [])
        print(f"\n[② 사다리 깊이단계 {peak}]  ⚠️ 이것도 무작위 추첨 1회(= seed 복제의 0번째 표본) · 자세 {len(lad)}개")
        if lad:
            print("     " + "  ".join(f"{g(p['dockq'],2)}" for p in lad)
                  + f"   → best {g(max((p['dockq'] for p in lad if p['dockq'] is not None), default=None),3)}")
            print("     recall: " + "  ".join(f"{g(p['recall'],2)}" for p in lad))

        # ── ③ seed 복제 전체 ──────────────────────────────────────────
        cjp = os.path.join(a.targets_dir, t, "chains.json")
        native = os.path.join(a.targets_dir, t, "native.cif")
        if not os.path.exists(cjp):
            print(f"\n[③ seed 복제] chains.json 없음 — 건너뜀\n")
            continue
        cj = json.load(open(cjp))
        tr = PF.native_true(cj, native, a.cutoff)
        if tr is None:
            print(f"\n[③ seed 복제] native epitope 실패 — 건너뜀\n")
            continue
        true, _ = tr
        try:
            pop = PF.ES.popular_refset(cj, grp.get(t, ""))
        except Exception:
            pop = None

        with tempfile.TemporaryDirectory() as td:
            natm = PF.native_merged(cj, native, td)
            if natm is None:
                print(f"\n[③ seed 복제] native merge 실패 — 건너뜀\n")
                continue
            base = os.path.join(outroot, model, t)
            ddirs = sorted(glob.glob(os.path.join(base, "d*")))
            if not ddirs:
                print(f"\n[③ seed 복제] 예측 폴더 없음: {base}")
                print("     → 미실행. run_seedrep_cand.sh 로그에서 이 후보의 메시지를 확인할 것.\n")
                continue
            for ddir in ddirs:
                depth = os.path.basename(ddir)
                sdirs = sorted(glob.glob(os.path.join(ddir, "seed*")))
                print(f"\n[③ seed 복제 · {depth}]  seed {len(sdirs)}개 × 자세 5개 = 전부 나열")
                print(f"     {'seed':7}{'자세별 DockQ (5개)':46}{'best':>7}{'recall@best':>12}{'max rec':>9}")
                print("     " + "-" * 82)
                bests, brecs = [], []
                for sdir in sdirs:
                    s = os.path.basename(sdir)
                    poses = sorted(glob.glob(os.path.join(sdir, "results", "**", "*.cif"), recursive=True))
                    vals = []
                    for p in poses:
                        try:
                            met = PF.pose_all_metrics(cj, p, a.cutoff, true, pop)
                            pm = PF.pose_merged(cj, p, td)
                            q = PF.dockq(pm, natm) if pm else None
                        except Exception:
                            continue
                        vals.append((q, (met or {}).get("recall"), os.path.basename(p)))
                        allrows.append(dict(target=t, model=model, depth=depth, seed=s,
                                            pose=os.path.basename(p),
                                            dockq=(round(q, 4) if q is not None else ""),
                                            recall=(round(met["recall"], 4) if met else "")))
                    if not vals:
                        print(f"     {s:7}(자세 없음)")
                        continue
                    dqs = [v[0] for v in vals if v[0] is not None]
                    if not dqs:
                        print(f"     {s:7}(채점 실패)")
                        continue
                    bi = max(range(len(vals)), key=lambda i: (vals[i][0] is not None, vals[i][0] or -1))
                    bq, brec = vals[bi][0], vals[bi][1]
                    recs = [v[1] for v in vals if v[1] is not None]
                    mark = " ★" if (bq is not None and bq >= SUCC) else "  "
                    print(f"     {s:7}" + "  ".join(f"{g(v[0],2)}" for v in vals).ljust(46)
                          + f"{g(bq,3):>7}{mark}{g(brec,2):>10}{g(max(recs) if recs else None,2):>9}")
                    bests.append(bq)
                    if brec is not None:
                        brecs.append((bq, brec))
                if bests:
                    n = len(bests)
                    fullbest = max((p["dockq"] for p in full if p["dockq"] is not None), default=None)
                    ge = sum(1 for x in bests if x >= SUCC)
                    print("     " + "-" * 82)
                    print(f"     분포: 중앙값 {st.median(bests):.3f} · 최소 {min(bests):.3f} · 최대 {max(bests):.3f}"
                          f" · ≥{SUCC} {ge}/{n}")
                    print(f"     [표본예산 맞춤] seed별 best-of-5 (n={n})  vs  full의 best-of-5 = {g(fullbest,3)}"
                          f"  → full 초과 {sum(1 for x in bests if fullbest is not None and x > fullbest)}/{n}")
                    if brecs:
                        hi = [rc for dq, rc in brecs if dq >= SUCC]
                        lo = [rc for dq, rc in brecs if dq < SUCC]
                        print(f"     [짝지어 보기] DockQ 성공한 seed의 결합자리 회복 평균 "
                              f"{(sum(hi)/len(hi)) if hi else float('nan'):.2f} (n={len(hi)})  vs  "
                              f"실패한 seed {(sum(lo)/len(lo)) if lo else float('nan'):.2f} (n={len(lo)})")
        print()

    if allrows:
        os.makedirs(os.path.dirname(a.csv_out) or ".", exist_ok=True)
        with open(a.csv_out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(allrows[0].keys()))
            w.writeheader(); w.writerows(allrows)
        print(f"→ 자세 단위 원자료 {a.csv_out} ({len(allrows)}행)")
    print("\n판정 질문: ②가 ③ 안에 있나(당연 — 같은 추첨)가 아니라, ③ 분포가 ①보다 위에 있나.")


if __name__ == "__main__":
    main()
