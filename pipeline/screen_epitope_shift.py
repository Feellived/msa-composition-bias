#!/usr/bin/env python3
"""[후보 탐색 2] '흔한 자리 → 진짜 자리'로 예측이 옮겨간 복합체를 찾는다.

기존 선별(screen_candidates.py)은 DockQ(자세 점수)로만 골랐다. 그런데 결합 자리를 맞혀도
방향이 틀리면 DockQ는 낮게 나온다. '편향에서 벗어나 정답 자리를 찾았다'는 이야기를 하려면
**자리 자체의 이동**으로 골라야 한다.

pose_features.csv의 두 값을 쓴다(둘 다 '예측한 접촉면 중 몇 %가 그 영역인가'):
  recall   = 진짜 결합자리와 겹치는 비율      ← 오를수록 정답 쪽
  overrep  = 그 항원에서 흔한 자리와 겹치는 비율 ← 내릴수록 흔한 자리에서 벗어남

찾는 것: MSA를 줄인 어떤 칸에서 **recall은 오르고 overrep은 내려간** 복합체.
        (= 흔한 자리에 붙어 있다가 진짜 자리로 옮겨감)

⚠️ A군은 진짜 자리가 곧 흔한 자리라 이 패턴이 원리상 안 나온다(8ulr이 A군). **B군이 본선.**

사용(stdlib only):
  python screen_epitope_shift.py                      # protenix
  python screen_epitope_shift.py --model boltz        # boltz 재실행 후
  python screen_epitope_shift.py --min-shift 0.15
"""
import argparse, csv, json, math, os
from collections import defaultdict


def f(x):
    try:
        v = float(x)
        return None if math.isnan(v) else v
    except Exception:
        return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", default="results/pose_features.csv")
    ap.add_argument("--model", default="protenix")
    ap.add_argument("--targets-dir", default="targets")
    ap.add_argument("--min-shift", type=float, default=0.10,
                    help="recall 상승분 + overrep 하락분의 합 문턱")
    ap.add_argument("--out", default="results/screen_epitope_shift.csv")
    a = ap.parse_args()

    # (target, rung) -> 그 칸에서 '진짜 자리와 가장 많이 겹친' pose의 (recall, overrep, dockq)
    best = defaultdict(dict)
    for r in csv.DictReader(open(a.csv)):
        if r["model"] != a.model:
            continue
        rc, ov = f(r.get("recall")), f(r.get("overrep"))
        if rc is None or ov is None:
            continue                      # 흔한 자리 정보가 없는 항원은 제외
        t, rg = r["target"], int(float(r["rung"]))
        cur = best[t].get(rg)
        if cur is None or rc > cur[0]:
            best[t][rg] = (rc, ov, f(r.get("dockq")) or 0.0)

    def meta(t):
        try:
            cj = json.load(open(os.path.join(a.targets_dir, t, "chains.json")))
            return str(cj.get("AB", "?")), str(cj.get("label", "")), str(cj.get("antigen_grp", ""))
        except Exception:
            return "?", "", ""

    rows = []
    for t, bk in best.items():
        if 0 not in bk or len(bk) < 2:
            continue
        r0, o0, d0 = bk[0]
        cand = []
        for rg, (rc, ov, dq) in bk.items():
            if rg == 0:
                continue
            shift = (rc - r0) + (o0 - ov)          # 진짜 쪽으로 + 흔한 자리에서 멀어짐
            if rc > r0 and ov < o0:                # 두 방향 모두 만족해야 함
                cand.append((shift, rg, rc, ov, dq))
        if not cand:
            continue
        shift, rg, rc, ov, dq = max(cand)
        g, lab, agrp = meta(t)
        rows.append(dict(target=t, grp=g, site=lab, antigen=agrp,
                         full_recall=round(r0, 3), full_overrep=round(o0, 3), full_dockq=round(d0, 3),
                         best_rung=rg, recall=round(rc, 3), overrep=round(ov, 3), dockq=round(dq, 3),
                         d_recall=round(rc - r0, 3), d_overrep=round(ov - o0, 3),
                         shift=round(shift, 3), n_rung_ok=len(cand)))

    rows.sort(key=lambda r: ((0 if r["grp"] == "B" else 1), -r["shift"]))
    print(f"모델 {a.model} · '흔한 자리 → 진짜 자리' 이동 탐색 (문턱 {a.min_shift})")
    print("  recall=예측 중 진짜 자리 비율 ↑ / overrep=예측 중 흔한 자리 비율 ↓ 이면 이동\n")
    print(f"{'target':11}{'군':>3}{'결합자리':>10}{'recall(full→최적)':>20}"
          f"{'overrep(full→최적)':>21}{'이동':>7}{'@rung':>7}{'DockQ':>7}{'칸수':>5}")
    print("-" * 92)
    hits = []
    for r in rows:
        if r["shift"] < a.min_shift:
            continue
        star = "  ★B군" if r["grp"] == "B" else ""
        rc_s = f"{r['full_recall']:.2f} → {r['recall']:.2f}"
        ov_s = f"{r['full_overrep']:.2f} → {r['overrep']:.2f}"
        print(f"{r['target']:11}{r['grp']:>3}{r['site'][:9]:>10}{rc_s:>20}{ov_s:>21}"
              f"{r['shift']:>7.3f}{r['best_rung']:>7}{r['dockq']:>7.2f}{r['n_rung_ok']:>5}{star}")
        hits.append(r)
    nb = [r for r in hits if r["grp"] == "B"]
    print(f"\n문턱 넘은 복합체 {len(hits)}개 · 그중 B군(드문 자리) {len(nb)}개"
          + (f" → {', '.join(r['target'] for r in nb)}" if nb else " (없음)"))
    print("⭐ B군이 본선: 진짜 자리가 흔한 자리와 달라 '흔한 자리로 쏠렸다가 벗어나 정답을 찾았다'를 직접 보일 수 있다.")
    print("   A군은 진짜 자리 = 흔한 자리라 이 서사가 원리상 성립하지 않는다.")
    print("⚠️ 여기 뜬 것은 '한 칸에서 그렇게 보였다'까지다. 반드시 comp_x_reps로 반복 실행해 확인할 것"
          "\n   (한 칸은 추첨 1회 — 9azr가 그렇게 떨어졌다).")

    if rows:
        os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
        with open(a.out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            w.writeheader(); w.writerows(rows)
        print(f"→ {a.out}")


if __name__ == "__main__":
    main()
