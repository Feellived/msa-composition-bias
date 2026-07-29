#!/usr/bin/env python3
"""[결과 모으기] 타깃마다 흩어진 채점 결과를 한 파일로 합친다.

analyze_target.sh 를 돌리면 타깃당 파일이 대여섯 개씩 생겨 29개 타깃이면 백 개가 넘는다.
그대로는 훑어보기도, 옮기기도 어렵다. 여기서 두 개로 줄인다.

  results/maintest_summary.csv   타깃 x 지표 한 줄씩 — 이질성 p, 성공 수, 검정 p, 판정
  results/maintest_poses.csv     자세 단위 원자료 전부(타깃 열 추가)

지표 세 가지의 뜻:
  dockq   자세가 정답 구조와 얼마나 맞나(높을수록 좋음, 성공선 0.49)
  recall  진짜 결합자리를 얼마나 덮나(높을수록 좋음)
  overrep 흔한 자리에 얼마나 붙나(낮을수록 편향에서 벗어난 것)

화면에는 요약표를 찍고, 조성 간 이질성이 유의한 타깃을 따로 모아 보여준다.
⚠️ 이건 모아 보기용이지 판정이 아니다. 판정은 인수인계서 Ⅱ 6.4절 기준으로만 한다.

사용:
  python collect_results.py
  python collect_results.py --metric recall     # 화면 표만 이 지표로
"""
import argparse, csv, glob, os, re


def num(x, d=None):
    try:
        return float(x)
    except (TypeError, ValueError):
        return d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default="results")
    ap.add_argument("--metric", default="recall", choices=["dockq", "recall", "overrep"],
                    help="화면 표에 쓸 지표(파일에는 세 가지 모두 들어간다)")
    a = ap.parse_args()

    # ── ① 지표별 요약 합치기 ──────────────────────────────────────────────
    rows = []
    pat = re.compile(r"summary_(.+)_(dockq|recall|overrep)\.csv$")
    for f in sorted(glob.glob(os.path.join(a.dir, "summary_*_*.csv"))):
        m = pat.search(os.path.basename(f))
        if not m:
            continue
        tgt, metric = m.group(1), m.group(2)
        for r in csv.DictReader(open(f)):
            r = dict(r)
            r["metric"] = metric
            r.setdefault("target", tgt)
            rows.append(r)

    if not rows:
        raise SystemExit(f"!! {a.dir} 에 summary_*_*.csv 가 없다.\n"
                         "   analyze_target.sh 를 먼저 돌릴 것"
                         "(옛 버전으로 돌렸다면 요약이 한 파일에 덮어써졌을 수 있다 → 다시 돌릴 것).")

    cols = ["target", "metric"] + [c for c in rows[0] if c not in ("target", "metric")]
    out1 = os.path.join(a.dir, "maintest_summary.csv")
    with open(out1, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(sorted(rows, key=lambda r: (r["target"], r["metric"])))
    print(f"→ {out1}  ({len(rows)}행, 타깃 {len({r['target'] for r in rows})}개)\n")

    # ── ② 자세 단위 원자료 합치기 ─────────────────────────────────────────
    poses, npose = [], 0
    for f in sorted(glob.glob(os.path.join(a.dir, "compreps_*.csv"))):
        base = os.path.basename(f)
        if base.startswith("compreps_summary"):
            continue
        for r in csv.DictReader(open(f)):
            poses.append(r)
            npose += 1
    if poses:
        out2 = os.path.join(a.dir, "maintest_poses.csv")
        with open(out2, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(poses[0].keys()), extrasaction="ignore")
            w.writeheader()
            w.writerows(poses)
        print(f"→ {out2}  ({npose}행)\n")

    # ── ③ 화면 표 ────────────────────────────────────────────────────────
    sel = [r for r in rows if r["metric"] == a.metric]
    print(f"지표 = {a.metric}   (단위 = 실행 1회, 자세 중 최고)")
    print(f"{'타깃':13}{'원래':>6}{'조성':>6}{'중앙값(원래→조성)':>20}"
          f"{'이질성 p':>10}{'순위검정 p':>11}  판정")
    print("-" * 88)
    het = []
    for r in sorted(sel, key=lambda x: x["target"]):
        ph = num(r.get("p_heterogeneity"))
        mark = " ⭐" if (ph is not None and ph < 0.05) else ""
        if mark:
            het.append((r["target"], ph))
        med = "{} → {}".format(r.get("med_full", ""), r.get("med_red", ""))
        pstr = f"{ph:.4f}" if ph is not None else "-"
        print(f"{r['target']:13}{r.get('n_full',''):>6}{r.get('n_red',''):>6}"
              f"{med:>20}{pstr:>10}"
              f"{str(r.get('p_ranktest','-')):>11}  {r.get('verdict','')}{mark}")

    print("-" * 88)
    if het:
        print(f"조성 간 이질성이 유의한 타깃 {len(het)}개 (p < 0.05):")
        for t, p in sorted(het, key=lambda x: x[1]):
            print(f"   {t}  p = {p:.4f}")
    else:
        print("조성 간 이질성이 유의한 타깃 없음.")
    print("\n※ 이 표는 모아 보기용이다. 판정은 결과를 보기 전에 정해둔 기준(인수인계서 Ⅱ 6.4절)으로만 한다.")
    print("※ 여러 타깃의 p 를 하나로 합치지 말 것 — 반응 없는 타깃이 섞이면 희석된다.")


if __name__ == "__main__":
    main()
