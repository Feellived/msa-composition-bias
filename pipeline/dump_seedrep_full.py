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
    ap.add_argument("--maintest", default="maintest.csv",
                    help="본 검정 명단. --only 는 이 명단(status=run)을 먼저 본다")
    ap.add_argument("--cand-first", action="store_true",
                    help="옛 후보 명단(--cand)을 우선한다. 옛 arm 결과를 다시 볼 때만 쓸 것")
    ap.add_argument("--depth", default="",
                    help="깊이 폴더를 하나로 못박는다(예: d90 또는 90). 폴더가 여럿이면 필수")
    ap.add_argument("--all-depths", action="store_true",
                    help="깊이 폴더 여러 개를 한 표에 합친다. 설계가 다른 실행이 섞이므로 권장하지 않음")
    ap.add_argument("--csv-out", default="results/seedrep_poses.csv")
    a = ap.parse_args()

    # --only 를 쓰면 그 후보 행만 남아 원자료를 통째로 덮어쓴다(다른 후보 채점 결과 소실).
    # → 기본 파일명일 때는 후보별 파일로 자동 분리한다.
    if a.only and a.csv_out == ap.get_default("csv_out"):
        a.csv_out = f"results/seedrep_poses_{a.only}.csv"
        print(f"[안내] --only 사용 → 원자료를 {a.csv_out} 로 따로 씁니다"
              f" (results/seedrep_poses.csv 보존). 전체 갱신은 --only 없이 실행.\n")

    # 채점 대상 만들기.
    #   ⚠️ 명단이 둘이다 — seedrep_cand.csv(옛 후보 5개, 일부는 폐기된 boltz arm)와
    #      maintest.csv(본 검정). 두 명단에 같은 타깃이 있으면 어느 쪽을 쓰느냐로
    #      모델·깊이가 달라진다. 예전에는 옛 후보를 우선해서, 본 검정 타깃인데도
    #      없어진 boltz 폴더를 찾다가 "자료 없음"으로 끝났다(2026-07-29 9y0a_AB).
    #   → --only 는 본 검정 명단을 먼저 본다. 옛 arm을 보려면 --cand-first.
    cand_rows = list(csv.DictReader(open(a.cand))) if os.path.exists(a.cand) else []
    known = {r["target"] for r in cand_rows}

    def from_maintest(name):
        if not os.path.exists(a.maintest):
            return None
        for r in csv.DictReader(open(a.maintest)):
            if r.get("target") == name and r.get("status") == "run":
                return dict(target=r["target"], model=r.get("model") or "protenix",
                            peak_rung=r.get("rung") or "0",
                            replicas=r.get("n_comp") or "6", obs_dq="", obs_rec="",
                            # 본 검정이 쓴 깊이(서열 수). 같은 타깃 밑에 옛 arm 깊이 폴더가
                            # 같이 있을 때 어느 쪽이 본 검정인지 가리는 데 쓴다.
                            depth_hint=(r.get("n_rows") or ""))
        return None

    if a.only:
        got = None if a.cand_first else from_maintest(a.only)
        if got is not None:
            if a.only in known:
                print(f"[안내] {a.only} 는 두 명단에 다 있다 → 본 검정({a.maintest}) 행을 쓴다."
                      f" 옛 후보 행으로 보려면 --cand-first.")
                cand_rows = [r for r in cand_rows if r["target"] != a.only]
            cand_rows.append(got)
        elif a.only not in known:
            print(f"!! {a.only} 는 {a.cand} 에도 {a.maintest}(status=run) 에도 없다 — 채점할 대상이 없다.")
            print(f"   본 검정 타깃이면 pick_maintest_depth.py 로 {a.maintest} 를 먼저 만들 것.")
            raise SystemExit(2)

    grp = {r["target"]: r.get("group", "") for r in csv.DictReader(open(a.list))}
    pf = load_pf(a.pf)
    outroot = os.path.join(a.data, "seedrep_cand")
    allrows = []

    for r in cand_rows:
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
            # ⚠️ 깊이 폴더가 둘 이상이면 설계가 다른 실행(옛 arm + 본 검정)이 한 타깃 밑에
            #    같이 있는 것이다. 아래 채점기들은 깊이를 보지 않고 폴더 이름(seedfull/c#_r#)으로만
            #    묶으므로, 그대로 두면 두 설계가 한 표에 섞인 채 이질성 검정에 들어간다
            #    (2026-07-29 8k5g_HL d58+d90 · 8q7s_C d35+d86). 조용히 섞지 말고 여기서 멈춘다.
            if len(ddirs) > 1 and not a.all_depths:
                want = (a.depth or r.get("depth_hint") or "").strip()
                if want and not want.startswith("d"):
                    want = "d" + want
                keep = [d for d in ddirs if os.path.basename(d) == want]
                if keep:
                    src = "--depth" if a.depth else f"{a.maintest} 의 n_rows"
                    print(f"\n[③ seed 복제] 깊이 폴더 {len(ddirs)}개 발견 "
                          f"({', '.join(os.path.basename(d) for d in ddirs)})"
                          f" → {want} 만 쓴다 ({src}).")
                    ddirs = keep
                else:
                    print(f"\n[③ seed 복제] !! 깊이 폴더가 여러 개다: "
                          f"{', '.join(os.path.basename(d) for d in ddirs)}")
                    print("     설계가 다른 실행이 섞여 있을 수 있어 그대로 채점하지 않는다.")
                    print(f"     어느 것이 본 검정인지 정해서 다시 실행할 것 — 예: --depth {os.path.basename(ddirs[-1])}")
                    print("     (정말 합쳐서 보려면 --all-depths)")
                    raise SystemExit(4)
            for ddir in ddirs:
                depth = os.path.basename(ddir)
                sdirs = sorted(glob.glob(os.path.join(ddir, "seed*")))
                print(f"\n[③ seed 복제 · {depth}]  seed {len(sdirs)}개 × 자세 5개 = 전부 나열")
                print(f"     {'실행':13}{'자세별 DockQ (5개)':46}{'best':>7}{'recall@best':>12}{'max rec':>9}")
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
                                            recall=(round(met["recall"], 4) if met else ""),
                                            overrep=(round(met["overrep"], 4) if met and met["overrep"] == met["overrep"] else ""),
                                            n_contact=(met["n_contact"] if met else "")))
                    if not vals:
                        print(f"     {s:13}(자세 없음)")
                        continue
                    dqs = [v[0] for v in vals if v[0] is not None]
                    if not dqs:
                        print(f"     {s:13}(채점 실패)")
                        continue
                    bi = max(range(len(vals)), key=lambda i: (vals[i][0] is not None, vals[i][0] or -1))
                    bq, brec = vals[bi][0], vals[bi][1]
                    recs = [v[1] for v in vals if v[1] is not None]
                    mark = " ★" if (bq is not None and bq >= SUCC) else "  "
                    print(f"     {s:13}" + "  ".join(f"{g(v[0],2)}" for v in vals).ljust(46)
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

    # ⚠️ DockQ 도구가 없으면(= DockQ conda 환경 밖) dockq 열만 통째로 비고, recall·overrep 은
    #    자체 계산이라 멀쩡히 나온다. 그대로 CSV를 쓰면 뒤 단계가 조용히 무너진다 —
    #    score_compreps 는 빈 리스트에 max()로 죽고, epitope_cluster·site_reproducibility 는
    #    대표 자세를 DockQ 최고값으로 고르므로 "계산된 실행이 없음"이 된다(2026-07-29 실제 발생).
    #    반쪽짜리 파일을 남기면 다음 실행이 재계산을 건너뛰므로, 아예 쓰지 않고 멈춘다.
    if allrows and not any(r["dockq"] != "" for r in allrows):
        print(f"!! 자세 {len(allrows)}개를 훑었는데 DockQ 값이 한 개도 없다 (recall 은 계산됨).")
        print("   → DockQ 도구를 못 찾은 것이다. conda activate DockQ 후 다시 실행할 것.")
        print("   반쪽짜리 원자료를 남기지 않으려고 파일을 쓰지 않았다.")
        raise SystemExit(5)

    if allrows:
        os.makedirs(os.path.dirname(a.csv_out) or ".", exist_ok=True)
        with open(a.csv_out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(allrows[0].keys()))
            w.writeheader(); w.writerows(allrows)
        print(f"→ 자세 단위 원자료 {a.csv_out} ({len(allrows)}행)")
    else:
        print("!! 채점된 자세가 한 개도 없다 — 결과 파일을 쓰지 않았다.")
        print("   흔한 원인: ① 예측 폴더가 비었다 ② targets/<타깃>/native.cif 가 없다")
        print("             ③ chains.json 의 사슬 배정이 예측 산출물과 안 맞는다")
        raise SystemExit(3)
    print("\n판정 질문: ②가 ③ 안에 있나(당연 — 같은 추첨)가 아니라, ③ 분포가 ①보다 위에 있나.")


if __name__ == "__main__":
    main()
