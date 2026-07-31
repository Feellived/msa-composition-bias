#!/usr/bin/env python3
"""[그림] 복합체 × 조성 성공률 히트맵 — 발표 S02(한 장 요약)용.

예전 F17_per_composition_heatmap.png 를 만든 스크립트가 남아 있지 않아 자료에서 다시 그린다.

═══ 예전 그림의 문제 (이 스크립트가 고치는 것) ═══
조성을 8가지 돌린 것은 8ulr_HL 하나뿐인데 열을 8개로 그려서, 나머지 29행의 7·8열이 비었다.
그 빈 칸이 컬러맵의 0.0(전부 실패)과 똑같이 흰색으로 보여 "성공률 0"과 구별되지 않았다.
캡션에 "빈 칸은 미실행이지 0이 아니다"라고 적어 글로 때우고 있었다.
→ 여기서는 **조성 열을 6개로 고정**한다. 8ulr_HL 의 조성 7·8번은 별도 슬라이드에서 다루므로
  정보 손실이 없고, 남는 흰 칸은 전부 진짜 0.0 이 된다.

═══ 값의 정의 ═══
  행   = 복합체 (results/compreps_<타깃>.csv 하나)
  열   = 조성 0~5, 세로 구분선 오른쪽에 원래 MSA 한 칸
  칸   = 그 조건의 **실행 단위 성공률**
         실행 하나 = seed 폴더 하나(seed<조성>_r<반복> · seedfull_r<n>)
         실행의 점수 = 그 실행의 자세들 중 recall 최댓값
         성공 = 점수 ≥ 0.4 (6.4 사전 확정 기준)
         칸 값 = 성공한 실행 수 ÷ 그 조건의 실행 수

⚠️ 자료가 없는 칸이 생기면 회색 빗금으로 그리고 개수를 화면에 알린다 — 0.0 과 섞이지 않게.

사용:
  python report/fig_composition_heatmap.py                      # 기본 경로
  python report/fig_composition_heatmap.py --results pipeline/results --outdir report/figures
  python report/fig_composition_heatmap.py --thr 0.4 --max-comp 6
출력: composition_heatmap.png(라벨 있음) · composition_heatmap_nolabel.png(발표용) · .csv(원자료)
"""
import argparse, csv, glob, os, re
from collections import defaultdict

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

RUN_RE = re.compile(r"^seed(?P<comp>full|\d+)_r(?P<rep>\d+)$")


def load_target(path, thr):
    """compreps CSV 하나 → {조성키: 성공률}. 조성키는 int(0~) 또는 'full'."""
    runs = defaultdict(list)                    # (조성, 반복) → [recall...]
    tgt = ""
    for r in csv.DictReader(open(path)):
        tgt = r.get("target") or tgt
        m = RUN_RE.match(r.get("seed", ""))
        if not m:
            continue
        try:
            v = float(r["recall"])
        except (TypeError, ValueError, KeyError):
            continue
        if v != v:
            continue
        comp = m.group("comp")
        comp = comp if comp == "full" else int(comp)
        runs[(comp, m.group("rep"))].append(v)
    by = defaultdict(list)                      # 조성 → [실행 점수...]
    for (comp, _), vals in runs.items():
        by[comp].append(max(vals))              # 실행의 점수 = 자세 중 최댓값
    return tgt, {c: (sum(1 for x in v if x >= thr) / len(v)) for c, v in by.items() if v}, \
           {c: len(v) for c, v in by.items()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default="pipeline/results")
    ap.add_argument("--outdir", default="report/figures")
    ap.add_argument("--thr", type=float, default=0.4, help="성공 기준(epitope recall)")
    ap.add_argument("--max-comp", type=int, default=6, help="그릴 조성 열 수")
    ap.add_argument("--order", default="diff", choices=["diff", "comp", "name"],
                    help="행 정렬 — diff(조성−원래, 기본) · comp(조성 평균) · name")
    a = ap.parse_args()

    rows = []
    for p in sorted(glob.glob(os.path.join(a.results, "compreps_*.csv"))):
        base = os.path.basename(p)[len("compreps_"):-4]
        if base.endswith("_summary") or base == "summary":
            continue
        tgt, rate, n = load_target(p, a.thr)
        if not rate:
            print(f"  ! {base}: 실행을 못 읽었다 — 건너뜀")
            continue
        rows.append((tgt or base, rate, n))
    if not rows:
        raise SystemExit(f"!! {a.results} 에서 compreps_*.csv 를 못 찾았다")

    # 조성 수가 기본값을 넘는 복합체를 알린다(예전 그림이 빈 열을 만든 원인)
    extra = [(t, max(c for c in r if c != "full") + 1) for t, r, _ in rows
             if any(isinstance(c, int) and c >= a.max_comp for c in r)]
    if extra:
        print("조성이 %d개를 넘는 복합체(초과분은 그리지 않는다):" % a.max_comp)
        for t, k in extra:
            print(f"  {t}: 조성 {k}가지 → 0~{a.max_comp-1} 만 사용")

    def key(x):
        t, r, _ = x
        comp = [r[c] for c in range(a.max_comp) if c in r]
        cm = sum(comp) / len(comp) if comp else 0.0
        fu = r.get("full", 0.0)
        return {"diff": -(cm - fu), "comp": -cm, "name": t}[a.order]
    rows.sort(key=key)

    M = np.full((len(rows), a.max_comp + 1), np.nan)     # 마지막 열 = 원래 MSA
    for i, (_, r, _) in enumerate(rows):
        for c in range(a.max_comp):
            if c in r:
                M[i, c] = r[c]
        if "full" in r:
            M[i, -1] = r["full"]
    nmiss = int(np.isnan(M).sum())

    # ── 그리기 ────────────────────────────────────────────────────────────────
    os.makedirs(a.outdir, exist_ok=True)
    for labeled in (True, False):
        h = max(4.0, 0.22 * len(rows) + 1.2)
        fig, ax = plt.subplots(figsize=(5.6, h))
        cmap = plt.get_cmap("Blues").copy()
        cmap.set_bad("#e8e8e8")                          # 자료 없음 = 회색(0.0과 구별)
        im = ax.imshow(np.ma.masked_invalid(M), cmap=cmap, vmin=0, vmax=1,
                       aspect="auto", interpolation="nearest")
        # 조성과 원래 MSA 사이 구분선
        ax.axvline(a.max_comp - 0.5, color="#333333", lw=1.4)
        ax.set_xticks(list(range(a.max_comp + 1)))
        if labeled:
            ax.set_xticklabels([f"C{c}" for c in range(a.max_comp)] + ["original\nMSA"], fontsize=8)
            ax.set_yticks(range(len(rows)))
            ax.set_yticklabels([t for t, _, _ in rows], fontsize=6)
            ax.set_xlabel("MSA composition", fontsize=9)
        else:
            ax.set_xticklabels([])
            ax.set_yticks([])
        for s in ("top", "right", "left", "bottom"):
            ax.spines[s].set_visible(False)
        ax.tick_params(length=2)
        cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cb.set_label("Run-level success rate", fontsize=9)
        cb.ax.tick_params(labelsize=8)
        fig.tight_layout()
        name = "composition_heatmap" + ("" if labeled else "_nolabel") + ".png"
        fig.savefig(os.path.join(a.outdir, name), dpi=300, bbox_inches="tight")
        plt.close(fig)
        print("→", os.path.join(a.outdir, name))

    # 원자료도 남긴다 — 그림의 숫자를 표로 검산할 수 있게
    cp = os.path.join(a.outdir, "composition_heatmap.csv")
    with open(cp, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["target"] + [f"C{c}" for c in range(a.max_comp)] + ["original_MSA"]
                   + [f"n_C{c}" for c in range(a.max_comp)] + ["n_original"])
        for (t, r, n), mrow in zip(rows, M):
            w.writerow([t] + [("" if x != x else round(x, 4)) for x in mrow]
                       + [n.get(c, 0) for c in range(a.max_comp)] + [n.get("full", 0)])
    print("→", cp)
    print(f"\n복합체 {len(rows)}종 · 조성 열 {a.max_comp} + 원래 MSA 1"
          f" · 성공 기준 recall ≥ {a.thr}")
    if nmiss:
        print(f"⚠️ 자료 없는 칸 {nmiss}개 — 회색으로 그렸다(0.0 과 구별). 위 CSV 의 빈칸과 같다.")
    else:
        print("자료 없는 칸 없음 — 흰색은 전부 성공률 0.0 이다.")


if __name__ == "__main__":
    main()
