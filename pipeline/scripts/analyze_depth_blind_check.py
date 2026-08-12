#!/usr/bin/env python3
"""[깊이 검증] 사다리 깊이를 **답을 보지 않고** 골라도 되는가 — 이미 있는 채점 결과만으로 확인한다.

■ 왜 이걸 묻나
  본 검정에 쓴 깊이는 `prep_pick_depth.py` 가 **정답 대비 recall 을 보고** 골랐다(칸마다
  자세 5개 중 0.4 를 넘는 것이 1~4개인 '되기도 안 되기도 하는' 칸). 실전에서는 정답이 없으므로
  이 규칙을 못 쓴다. 그래서 묻는다:

      "위 몇 칸을 그냥 다 돌린다"는 **눈먼 규칙**으로도 같은 것을 건지는가?

  건진다면 깊이는 계산으로 살 수 있는 문제고, 남는 문제는 선택뿐이다.

■ 무엇을 재나
  '건졌다' = 그 전략이 뽑은 자세들 중 **정답 결합자리를 0.4 이상 덮는 자세가 하나라도 있다**.
  0.4 는 `prep_pick_depth.py` 의 SUCC_RECALL 과 같은 값이며 **결과를 보기 전에 확정**됐다.
  ⚠️ 이것은 '자세가 맞았다'가 아니라 **'후보 안에 정답 자리가 들어왔다'**(= 생성 성공)이다.
  고르는 것은 별개 문제이고 이 스크립트는 거기에 답하지 않는다.

■ ⭐ 자세 예산을 맞춘다 (이 스크립트의 핵심)
  칸 4개를 쓰면 자세가 20개, 답을 보고 고른 칸 1개는 5개다. 그냥 비교하면 **"많이 뽑아서 이긴 것"**과
  구별이 안 된다(이 프로젝트에서 여러 번 데인 함정). 그래서 두 가지를 함께 낸다.

    · 그대로     — 전략이 뽑은 자세 전부를 쓴다(생성 예산이 다름을 인정하고 보는 값)
    · 예산 맞춤  — 그 자세 더미에서 **무작위로 5개만** 뽑았을 때 성공할 확률.
                   부트스트랩이 아니라 정확식:  P = 1 − C(n−k, 5) / C(n, 5)
                   (n = 자세 수, k = 성공 자세 수).  n < 5 면 그대로 쓴다.
      → '예산 맞춤' 합계가 답을 본 칸과 비슷하면, **깊이는 눈감고 골라도 손해가 없다**는 뜻이다.

■ 전략
    원래MSA    rung 0 만            (아무것도 안 한 기준선)
    답을본칸   maintest.csv 의 rung  (지금 쓴 규칙 — 정답을 봄)
    눈먼 1~2   rung 1,2
    눈먼 1~4   rung 1,2,3,4
    사다리전체 rung ≥ 1

■ ⚠️ 이 검증의 한계 (반드시 함께 말할 것)
  · 사다리는 칸마다 **1회(자세 5개)**뿐이다. 본 검정처럼 조성을 재추첨한 자료가 아니다.
    따라서 이 표는 '깊이를 눈감고 골라도 되나'의 **1차 확인**이지 본 검정의 재현이 아니다.
  · boltz 사다리 자료는 2026-07-27 a3m 사고로 **무효**다. 기본값이 protenix 인 이유다.
  · '후보 목록이 몇 개로 줄어드나'는 자리 군집이 필요해 여기서 못 잰다(구조 파일이 있어야 한다).

사용 (stdlib only · CPU 몇 초):
    cd ~/projects/msa-composition-bias/pipeline
    python analyze_depth_blind_check.py
    python analyze_depth_blind_check.py --thr 0.40 --model protenix --out results/depth_blind_check.csv
"""
import argparse
import csv
import math
import os
from collections import defaultdict

SIG = {"8k3k_D", "8k46_I", "8k5h_HL", "8siq_HL", "9mqr_DE", "9zdu_HL"}          # 본 검정 유의 6종
STAR = {"8k3k_D", "8sit_HL", "8ume_HL", "8k46_I"}                               # ★ 4종


def fnum(x):
    try:
        v = float(x)
        return None if v != v else v
    except Exception:
        return None


def p_budget(n, k, budget):
    """자세 n개 중 성공이 k개일 때, 무작위 budget개를 뽑아 하나라도 성공할 확률(정확식)."""
    if n <= 0:
        return 0.0
    if n <= budget:
        return 1.0 if k > 0 else 0.0
    if k == 0:
        return 0.0
    return 1.0 - math.comb(n - k, budget) / math.comb(n, budget)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pf", default="results/pose_features.csv")
    ap.add_argument("--maintest", default="maintest.csv")
    ap.add_argument("--model", default="protenix",
                    help="boltz 사다리는 a3m 사고로 무효 — 기본은 protenix")
    ap.add_argument("--thr", type=float, default=0.40,
                    help="정답 결합자리를 이만큼 덮으면 성공 (사전 확정값 0.40)")
    ap.add_argument("--budget", type=int, default=5, help="예산 맞춤 비교에 쓸 자세 수")
    ap.add_argument("--out", default="results/depth_blind_check.csv")
    a = ap.parse_args()

    if not os.path.exists(a.pf):
        raise SystemExit(f"!! {a.pf} 가 없다. pipeline 폴더에서 실행할 것.")
    if not os.path.exists(a.maintest):
        raise SystemExit(f"!! {a.maintest} 가 없다.")

    picked, meta = {}, {}
    for r in csv.DictReader(open(a.maintest)):
        if r.get("status") != "run":
            continue
        picked[r["target"]] = int(float(r["rung"]))
        meta[r["target"]] = (r.get("group", ""), r.get("n_rows", ""))
    if not picked:
        raise SystemExit("!! maintest.csv 에 status=run 인 행이 없다.")

    # target → rung → [recall, ...]
    rec = defaultdict(lambda: defaultdict(list))
    skipped_model = 0
    for r in csv.DictReader(open(a.pf)):
        if r.get("model") != a.model:
            skipped_model += 1
            continue
        t = r["target"]
        if t not in picked:
            continue
        v, rg = fnum(r.get("recall")), fnum(r.get("rung"))
        if v is None or rg is None:
            continue
        rec[t][int(rg)].append(v)

    missing = [t for t in picked if t not in rec]
    if missing:
        print(f"⚠️ 사다리 자료가 없는 복합체 {len(missing)}종 — 분모에서 뺀다: {' '.join(sorted(missing))}\n")

    STRATS = [
        ("원래MSA",    lambda t: [0]),
        ("답을본칸",   lambda t: [picked[t]]),
        ("눈먼 1~2",   lambda t: [1, 2]),
        ("눈먼 1~4",   lambda t: [1, 2, 3, 4]),
        ("사다리전체", lambda t: sorted(g for g in rec[t] if g >= 1)),
    ]

    rows, agg = [], {n: dict(hit=0, exp=0.0, npose=0, ncov=0) for n, _ in STRATS}
    tgts = sorted(t for t in picked if t in rec)

    print(f"모델 {a.model} · 성공선 recall ≥ {a.thr} · 예산 맞춤 = 자세 {a.budget}개\n")
    print("'건졌다' = 그 전략의 자세 중 정답 결합자리를 0.4 이상 덮는 것이 하나라도 있음(= 생성 성공).")
    print("괄호 안 = 그 더미에서 자세 5개만 무작위로 뽑았을 때 성공할 확률.\n")
    hdr = f"{'복합체':<11}{'칸':>3}{'서열':>6}  " + "".join(f"{n:>15}" for n, _ in STRATS)
    print(hdr)
    print("-" * 96)

    for t in tgts:
        row = dict(target=t, group=meta[t][0], n_rows=meta[t][1], picked_rung=picked[t])
        cells = []
        for name, sel in STRATS:
            rg = [g for g in sel(t) if g in rec[t]]
            poses = [v for g in rg for v in rec[t][g]]
            n = len(poses)
            k = sum(1 for v in poses if v >= a.thr)
            hit = 1 if k > 0 else 0
            pb = p_budget(n, k, a.budget)
            if n:
                agg[name]["ncov"] += 1
                agg[name]["hit"] += hit
                agg[name]["exp"] += pb
                agg[name]["npose"] += n
            row[f"{name}_hit"] = hit
            row[f"{name}_p5"] = round(pb, 3)
            row[f"{name}_npose"] = n
            cells.append(f"{'○' if hit else '·'} ({pb:.2f})" if n else "  -   ")
        mark = "★" if t in STAR else ("✓" if t in SIG else " ")
        print(f"{mark}{t:<10}{picked[t]:>3}{meta[t][1]:>6}  " + "".join(f"{c:>15}" for c in cells))
        rows.append(row)

    N = len(tgts)
    print("\n" + "=" * 96)
    print(f"{'전략':<12}{'건진 수':>10}{'예산맞춤 기대':>16}{'자세/복합체':>13}   해석")
    print("-" * 96)
    base = None
    for name, _ in STRATS:
        d = agg[name]
        if not d["ncov"]:
            continue
        pp = d["npose"] / d["ncov"]
        note = ""
        if name == "답을본칸":
            base = d["exp"]
            note = "← 지금 쓴 규칙(정답을 봄)"
        elif base is not None:
            note = f"답을본칸 대비 {d['exp'] - base:+.1f}"
        print(f"{name:<12}{d['hit']:>7}/{N:<3}{d['exp']:>13.1f}/{N:<3}{pp:>13.1f}   {note}")

    print("\n[읽는 법]")
    print("  · '건진 수'만 보면 칸을 늘린 쪽이 당연히 유리하다(자세를 더 뽑았으니까).")
    print("  · ⭐ 판정은 **'예산맞춤 기대'**로 한다. 자세 5개로 맞춰 놓고도 눈먼 규칙이")
    print("    답을 본 칸과 비슷하면 → 깊이는 눈감고 골라도 손해가 없다.")
    print("    뚜렷이 낮으면 → 깊이 선택에 정답 정보가 실제로 들어가 있었다는 뜻이다.")
    print("\n⚠️ 사다리는 칸마다 1회(자세 5개)뿐이라 본 검정의 재현이 아니라 1차 확인이다.")
    print("⚠️ 이 표는 '후보 안에 정답 자리가 있나'까지만 답한다. 그중 하나를 고르는 문제는 별개다.")

    if rows:
        os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
        with open(a.out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)
        print(f"\n→ {a.out}")


if __name__ == "__main__":
    main()
