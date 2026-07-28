#!/usr/bin/env python3
"""[판독] 사다리 결과를 읽어 본 검정에 넣을 타깃과 깊이를 규칙대로 고른다.

⚠️ 규칙은 **결과를 보기 전(2026-07-28)에 확정**한 것이다. 여기서 바꾸지 말 것.

  각 칸마다 구조 5개의 '진짜 결합자리 덮음'(recall)을 보고 0.4 이상인 개수를 센다.
    · 12칸 60개 구조에서 0.4 이상이 **하나도 없음**  → 본 검정 안 함. 비율의 분모에는
      남기고 **실패로 계상**(안 돌린 것을 성공으로 세지 않으므로 안전한 방향).
    · **모든 칸에서 5개 전부** 0.4 이상          → 본 검정 안 함. '구제할 것이 없던 경우'로 별도 표시.
    · 그 밖                                    → **본 검정 대상.** 깊이는 0.4 이상이
      **1~4개인 칸**(되기도 하고 안 되기도 하는 중간 지대) 중 **서열이 가장 많은 칸**.
      그런 칸이 없으면 서열 수가 **1746개**(8ulr에서 효과가 확인된 수)에 가장 가까운 칸.
    · ⚠️ rung0은 후보에서 제외한다(2026-07-28 추가). rung0 = 원래 MSA 전체이고
      comp_x_reps.sh가 그것을 그대로 대조군(seedfull.a3m)으로 쓴다. 거기서 조성을 재추첨하면
      전체에서 전체를 뽑는 것이라 조성이 전부 같아지고 조성군=대조군이 되어 이질성 검정이
      정의되지 않는다. 이는 결과를 보고 기준을 옮기는 것이 아니라 **실행 불가능한 선택지를
      막는 것**이며, 검정 방식·성공 기준·설계값(6조성×4반복+원래8)은 그대로다.

왜 결과를 보고 깊이를 골라도 되나: 조성 간 이질성 검정은 **그 실험의 총 성공 횟수를
조건으로 삼기 때문에**, 그 깊이에서 성공률이 얼마인지는 검정에서 상쇄된다. 게다가 본 검정은
조성과 반복을 전부 새로 뽑으므로, 깊이를 고를 때 쓴 사다리 자료가 새 자료에 들어오지 않는다.
→ 선별이 영향을 주는 것은 '몇 개 중 몇 개'라는 비율뿐이며, 검정 자체는 타당하다.
   (그래서 절대 성공률이나 '얕은 쪽이 원래보다 낫다'는 비교는 주장하지 않는다.)

사용(stdlib only):
  python pick_maintest_depth.py                        # 전체
  python pick_maintest_depth.py --group RBD            # 코로나 세트만
  python pick_maintest_depth.py --only 8sis_HL 9zdu_HL
  python pick_maintest_depth.py --out maintest.csv     # 실행기가 읽을 명단

--only 로 지정한 이름이 표에 안 나오면 이유를 구분해 알린다(조용히 빈 표를 내지 않는다):
  · sweep_targets.csv 에는 있으나 pose_features.csv 에 없음 → 아직 채점 안 됨(먼저 pose_features.py)
  · sweep_targets.csv 에 아예 없음                        → 이름 오타
  · --group 에 걸러짐
대상이 하나도 안 남으면 종료 코드 1로 멈춘다. (--only 에 아무것도 안 넘기면 필터가 통째로
무시되어 전체가 나오므로, 셸에서 인자가 비어 전달되는 경우를 위 오타 경고가 잡아준다.)
"""
import argparse, csv, os
from collections import defaultdict

SUCC_RECALL = 0.40      # 진짜 결합자리를 이만큼 덮으면 성공 (사전 확정 = 판정 기준)
SUCC_DOCKQ = 0.23       # 참고 표시용 문턱(자세 품질). 판정에는 쓰지 않는다.
SUCC_DOCKQ2 = 0.49      # 참고 표시용 문턱(CAPRI Medium)
ANCHOR_ROWS = 1746      # 8ulr에서 효과가 확인된 서열 수
N_COMP, N_REPS, N_FULL = 6, 4, 8       # 본 검정 설계 (2026-07-28 확정)


def f(x):
    try:
        v = float(x)
        return None if v != v else v
    except Exception:
        return None


def ladder_rows(data, target):
    """{rung: 서열수} — $DATA/ladders/<타깃>/<사슬>/neff.tsv 에서. 사슬이 여럿이면 첫 사슬."""
    base = os.path.join(data, "ladders", target)
    if not os.path.isdir(base):
        return {}
    for ch in sorted(os.listdir(base)):
        p = os.path.join(base, ch, "neff.tsv")
        if not os.path.exists(p):
            continue
        out = {}
        for i, ln in enumerate(open(p)):
            if i == 0:
                continue
            q = ln.split()
            if len(q) >= 2:
                out[int(q[0])] = int(q[1])
        if out:
            return out
    return {}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pf", default="results/pose_features.csv")
    ap.add_argument("--list", default="sweep_targets.csv")
    ap.add_argument("--data", default=os.environ.get("DATA", "/mnt/data/admuser/msadepth"))
    ap.add_argument("--model", default="protenix")
    ap.add_argument("--rungs", type=int, default=12)
    ap.add_argument("--group", default="")
    ap.add_argument("--only", nargs="*", default=[])
    ap.add_argument("--out", default="maintest.csv")
    ap.add_argument("--metric", default="recall", choices=["recall", "dockq"],
                    help="판정에 쓸 지표. 기본 recall(사전 확정). dockq는 민감도 확인 전용")
    ap.add_argument("--thr", type=float, default=None,
                    help="성공 문턱. 생략 시 recall=0.40 / dockq=0.49")
    a = ap.parse_args()

    meta = {r["target"]: r for r in csv.DictReader(open(a.list))}
    per = defaultdict(lambda: defaultdict(list))       # target -> rung -> [recall]
    for r in csv.DictReader(open(a.pf)):
        if r["model"] != a.model:
            continue
        per[r["target"]][int(float(r["rung"]))].append(
            {k: f(r.get(k)) for k in ("recall", "dockq", "overrep")})

    tgts = [t for t in meta if t in per]
    if a.group:
        tgts = [t for t in tgts if meta[t].get("group") == a.group]
    if a.only:
        tgts = [t for t in tgts if t in a.only]
        # 조용히 빈 표를 내지 않도록, 요청한 이름이 왜 빠졌는지를 구분해 알린다.
        no_score = [t for t in a.only if t in meta and t not in per]
        if no_score:
            print(f"⚠️ 채점 결과가 없는 타깃 {len(no_score)}개 — 먼저 "
                  f"pose_features.py --models {a.model} 를 돌릴 것:")
            print("   " + " ".join(no_score) + "\n")
        unknown = [t for t in a.only if t not in meta]
        if unknown:
            print(f"⚠️ {a.list} 에 없는 이름 {len(unknown)}개: {' '.join(unknown)}\n")
        dropped = [t for t in a.only if t in per and t in meta and t not in tgts]
        if dropped and a.group:
            print(f"⚠️ --group {a.group} 에 걸러진 타깃 {len(dropped)}개: {' '.join(dropped)}\n")
    tgts.sort()
    if not tgts:
        raise SystemExit(
            "!! 판독할 타깃이 없다. 위 경고를 먼저 해결할 것."
            if a.only else
            f"!! 판독할 타깃이 없다. {a.pf} 에 model={a.model} 인 행이 없거나 "
            f"--group {a.group!r} 에 걸리는 타깃이 없다.\n"
            "   (참고: --only 에 아무것도 안 넘기면 필터가 통째로 무시되어 전체가 나온다. "
            "이름을 넘겼는데 이 메시지가 보이면 셸에서 인자가 빈 채로 전달된 것이다.)")

    MET = a.metric
    THR = a.thr if a.thr is not None else (SUCC_RECALL if MET == "recall" else SUCC_DOCKQ2)
    label = {"recall": "결합자리 덮음(recall)", "dockq": "자세 품질(DockQ)"}[MET]
    print(f"판독 규칙: 칸마다 구조 5개 중 {label} ≥{THR} 인 개수를 센다.")
    print(f"           1~4개인 칸(중간 지대) 중 서열이 가장 많은 칸을 고른다. rung0(원래 MSA)은 제외.")
    if MET != "recall" or a.thr is not None:
        print("⚠️ 사전 확정 기준(recall ≥ 0.40)이 아니다 — 민감도 확인용 실행이다. "
              "이 결과로 본 검정 명단을 바꾸지 말 것.")
    print()
    print(f"{'타깃':13}{'군':4}{'칸':4}{'칸별 성공 구조수 (0~11칸)':30}{'고른 칸':>8}{'서열':>8}   판정")
    print("-" * 104)

    rows, n_run, n_skip_none, n_skip_all, n_incomplete = [], 0, 0, 0, 0
    ref, picked = {}, {}
    for t in tgts:
        nr = ladder_rows(a.data, t)
        hits, missing = {}, []
        d23, d49, rec_all, dq_all, ov = {}, {}, [], [], {}
        for k in range(a.rungs):
            rows_k = per[t].get(k, [])
            v = [x[MET] for x in rows_k if x.get(MET) is not None]
            if not v:
                missing.append(k); continue
            hits[k] = sum(1 for x in v if x >= THR)
            rr = [x["recall"] for x in rows_k if x.get("recall") is not None]
            dd = [x["dockq"] for x in rows_k if x.get("dockq") is not None]
            oo = [x["overrep"] for x in rows_k if x.get("overrep") is not None]
            rec_all += rr; dq_all += dd
            d23[k] = sum(1 for x in dd if x >= SUCC_DOCKQ)
            d49[k] = sum(1 for x in dd if x >= SUCC_DOCKQ2)
            if oo:
                ov[k] = sum(oo) / len(oo)
        ref[t] = dict(rec_max=(max(rec_all) if rec_all else float("nan")),
                      dq_max=(max(dq_all) if dq_all else float("nan")),
                      d23=d23, d49=d49, ov=ov)
        grp = meta[t].get("group", "")
        strip = "".join(str(hits.get(k, "·")) for k in range(a.rungs))
        if missing:
            n_incomplete += 1
            print(f"{t:13}{grp:4}{len(hits):<4}{strip:30}{'-':>8}{'-':>8}   ⚠️ 미완성 — 빠진 칸 {missing}")
            continue

        tot = sum(hits.values())
        if tot == 0:
            n_skip_none += 1
            print(f"{t:13}{grp:4}{len(hits):<4}{strip:30}{'-':>8}{'-':>8}   건너뜀 · 어느 깊이에서도 안 됨(실패로 계상)")
            R = ref.get(t, {})
            rows.append(dict(target=t, group=grp, model=a.model, rung="", n_rows="",
                             n_comp="", n_reps="", n_full="", status="no_response",
                             metric=MET, thr=THR,
                             recall_max=(f"{R.get('rec_max', float('nan')):.2f}"),
                             dockq_max=(f"{R.get('dq_max', float('nan')):.2f}"),
                             dq23_pick="", dq49_pick="",
                             overrep_full=(f"{R['ov'][0]:.2f}" if R.get("ov", {}).get(0) is not None else ""),
                             overrep_pick=""))
            continue
        if all(h == 5 for h in hits.values()):
            n_skip_all += 1
            print(f"{t:13}{grp:4}{len(hits):<4}{strip:30}{'-':>8}{'-':>8}   건너뜀 · 모든 깊이에서 성공(구제 불필요)")
            R = ref.get(t, {})
            rows.append(dict(target=t, group=grp, model=a.model, rung="", n_rows="",
                             n_comp="", n_reps="", n_full="", status="always_ok",
                             metric=MET, thr=THR,
                             recall_max=(f"{R.get('rec_max', float('nan')):.2f}"),
                             dockq_max=(f"{R.get('dq_max', float('nan')):.2f}"),
                             dq23_pick="", dq49_pick="",
                             overrep_full=(f"{R['ov'][0]:.2f}" if R.get("ov", {}).get(0) is not None else ""),
                             overrep_pick=""))
            continue

        # rung0 = 원래 MSA 전체(comp_x_reps.sh가 rung0.a3m을 그대로 seedfull.a3m으로 복사한다).
        # 거기서 조성을 재추첨하면 전체에서 전체를 뽑는 것이라 여섯 조성이 같은 집합이 되고,
        # 조성군과 대조군도 같아져 이질성 검정이 정의되지 않는다. 그래서 후보에서 제외한다.
        # (결과를 보고 기준을 옮기는 것이 아니라, 규칙이 실행 불가능한 선택지를 허용하던 것을
        #  막는 것이다. 검정 방식·성공 기준·설계값은 그대로다.)
        cand = [k for k in hits if k >= 1]
        mid = [k for k in cand if 1 <= hits[k] <= 4]
        if mid:
            pick = max(mid, key=lambda k: nr.get(k, 0))
            why = "중간 지대"
            if 1 <= hits.get(0, 0) <= 4 and nr.get(0, 0) > nr.get(pick, 0):
                why += " (rung0은 원래 MSA라 제외)"
        elif cand:
            pick = min(cand, key=lambda k: abs(nr.get(k, 10**9) - ANCHOR_ROWS))
            why = f"중간 지대 없음 → 서열 {ANCHOR_ROWS}에 가장 가까운 칸"
            if hits.get(0, 0) > 0:
                why += " ⚠️ 성공이 rung0(원래 MSA)에만 있음"
        else:
            print(f"{t:13}{grp:4}{len(hits):<4}{strip:30}{'-':>8}{'-':>8}   "
                  f"건너뜀 · rung1 이상이 없음(사다리 미완성)")
            n_incomplete += 1
            continue
        n_run += 1
        print(f"{t:13}{grp:4}{len(hits):<4}{strip:30}{('rung'+str(pick)):>8}{nr.get(pick,0):>8}   본 검정 · {why}")
        picked[t] = pick
        R = ref.get(t, {})
        rows.append(dict(target=t, group=grp, model=a.model, rung=pick, n_rows=nr.get(pick, ""),
                         n_comp=N_COMP, n_reps=N_REPS, n_full=N_FULL, status="run",
                         metric=MET, thr=THR,
                         recall_max=(f"{R.get('rec_max', float('nan')):.2f}"),
                         dockq_max=(f"{R.get('dq_max', float('nan')):.2f}"),
                         dq23_pick=R.get("d23", {}).get(pick, ""),
                         dq49_pick=R.get("d49", {}).get(pick, ""),
                         overrep_full=(f"{R['ov'][0]:.2f}" if R.get("ov", {}).get(0) is not None else ""),
                         overrep_pick=(f"{R['ov'][pick]:.2f}" if R.get("ov", {}).get(pick) is not None else "")))

    # ── 참고 지표: 판정에 쓰지 않지만 선별을 눈으로 확인하기 위한 전 후보 일람 ──────
    if ref:
        print("\n" + "-" * 116)
        print("[참고 지표 — 판정에는 쓰지 않는다] 후보 전체. 두 지표가 어긋나면 그 자체가 보고 대상이다.")
        print(f"{'타깃':13}{'군':4}{'recall최고':>9}{'DockQ최고':>10}   "
              f"{'DockQ≥0.23 성공수 (0~11칸)':32}{'인기자리 겹침 rung0→고른칸':>26}   상태")
        print("-" * 116)
        for t in tgts:
            R = ref.get(t)
            if not R:
                continue
            grp = meta[t].get("group", "")
            d23s = "".join(str(R["d23"].get(k, "·")) for k in range(a.rungs))
            pk = picked.get(t)
            o0 = R["ov"].get(0)
            op = R["ov"].get(pk) if pk is not None else None
            ovs = (f"{o0:.2f} → {op:.2f}" if (o0 is not None and op is not None)
                   else (f"{o0:.2f} → -" if o0 is not None else "-"))
            st = ("본 검정 rung%s" % pk) if pk is not None else "제외"
            print(f"{t:13}{grp:4}{R['rec_max']:>9.2f}{R['dq_max']:>10.2f}   "
                  f"{d23s:32}{ovs:>26}   {st}")
        print("\n  · recall최고 / DockQ최고 = 그 타깃의 전 칸·전 구조 통틀어 최고값. "
              "판정 문턱을 못 넘은 이유가 '거의 근접'인지 '전혀 못 찾음'인지를 가른다.")
        print("  · DockQ≥0.23 성공수 = 자세 품질 축의 칸별 성공 개수. recall 쪽 띠와 비교해 축이 어긋나는지 본다.")
        print("  · 인기자리 겹침 = 예측 접촉 중 과대표집 부위 비율의 칸 평균. 내려가면 인기 자리에서 벗어난 것"
              "(C군은 정의가 없어 빈칸).")

    print("\n" + "=" * 104)
    print(f"본 검정 대상 {n_run}개 · 어느 깊이에서도 안 됨 {n_skip_none}개 · 모든 깊이 성공 {n_skip_all}개"
          + (f" · 미완성 {n_incomplete}개" if n_incomplete else ""))
    if n_run:
        print(f"예상 실험 횟수: {n_run} × ({N_COMP}조성 × {N_REPS}반복 + 원래 {N_FULL}회)"
              f" = {n_run * (N_COMP * N_REPS + N_FULL)}회")
    if n_incomplete:
        print("⚠️ 미완성 타깃이 있다 — sweep이 아직 안 끝났거나 실패한 칸이 있다. 먼저 확인할 것.")
    if rows:
        with open(a.out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            w.writeheader(); w.writerows(rows)
        print(f"\n→ {a.out} ({len(rows)}행) — 본 검정 실행기가 읽을 명단")
    print("\n※ 비율을 낼 때 '어느 깊이에서도 안 됨'은 분모에 남기고 실패로 센다(안전한 방향).")
    print("   '모든 깊이 성공'은 구제할 것이 없던 경우이므로 비율에서 따로 표시한다.")


if __name__ == "__main__":
    main()
