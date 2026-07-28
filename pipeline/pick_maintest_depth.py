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

SUCC_RECALL = 0.40      # 진짜 결합자리를 이만큼 덮으면 성공
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
    a = ap.parse_args()

    meta = {r["target"]: r for r in csv.DictReader(open(a.list))}
    per = defaultdict(lambda: defaultdict(list))       # target -> rung -> [recall]
    for r in csv.DictReader(open(a.pf)):
        if r["model"] != a.model:
            continue
        v = f(r.get("recall"))
        per[r["target"]][int(float(r["rung"]))].append(v)

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

    print(f"판독 규칙: 칸마다 구조 5개 중 결합자리 덮음 ≥{SUCC_RECALL} 인 개수를 센다.")
    print(f"           1~4개인 칸(중간 지대) 중 서열이 가장 많은 칸을 고른다.\n")
    print(f"{'타깃':13}{'군':4}{'칸':4}{'칸별 성공 구조수 (0~11칸)':30}{'고른 칸':>8}{'서열':>8}   판정")
    print("-" * 104)

    rows, n_run, n_skip_none, n_skip_all, n_incomplete = [], 0, 0, 0, 0
    for t in tgts:
        nr = ladder_rows(a.data, t)
        hits, missing = {}, []
        for k in range(a.rungs):
            v = [x for x in per[t].get(k, []) if x is not None]
            if not v:
                missing.append(k); continue
            hits[k] = sum(1 for x in v if x >= SUCC_RECALL)
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
            rows.append(dict(target=t, group=grp, model=a.model, rung="", n_rows="",
                             n_comp="", n_reps="", n_full="", status="no_response"))
            continue
        if all(h == 5 for h in hits.values()):
            n_skip_all += 1
            print(f"{t:13}{grp:4}{len(hits):<4}{strip:30}{'-':>8}{'-':>8}   건너뜀 · 모든 깊이에서 성공(구제 불필요)")
            rows.append(dict(target=t, group=grp, model=a.model, rung="", n_rows="",
                             n_comp="", n_reps="", n_full="", status="always_ok"))
            continue

        mid = [k for k, h in hits.items() if 1 <= h <= 4]
        if mid:
            pick = max(mid, key=lambda k: nr.get(k, 0))
            why = "중간 지대"
        else:
            pick = min(hits, key=lambda k: abs(nr.get(k, 10**9) - ANCHOR_ROWS))
            why = f"중간 지대 없음 → 서열 {ANCHOR_ROWS}에 가장 가까운 칸"
        n_run += 1
        print(f"{t:13}{grp:4}{len(hits):<4}{strip:30}{('rung'+str(pick)):>8}{nr.get(pick,0):>8}   본 검정 · {why}")
        rows.append(dict(target=t, group=grp, model=a.model, rung=pick, n_rows=nr.get(pick, ""),
                         n_comp=N_COMP, n_reps=N_REPS, n_full=N_FULL, status="run"))

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
