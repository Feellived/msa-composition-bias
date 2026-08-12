#!/usr/bin/env python3
"""[Phase 1 확인] 무학습 재랭커 make-or-break — 배포가능 피처가 pose 선택기로 ipTM을 이기나.

lib_pose_features.py docstring이 지정한 Phase 1: "복합체 내부에서 각 피처로 pose 순위 → ipTM·무작위·
best-of-N null 대비 비교". 학습 없이(무학습) 각 후보 피처로 타깃별 pose를 하나 골라, 그 pose의
라벨(DockQ/recall)을 ipTM 선택·무작위·oracle과 비교한다.

배포가능 선택기(정답 안 씀, 실배포에서 계산 가능):
  iptm·ptm·plddt(신뢰도) · n_contact(접촉수) · dcc_pop(popular 자리서 멀리=off-hotspot)
  · overrep_lo(popular 적게) · pop_rank_lo
라벨(정답, 평가용): dockq / recall.  oracle=상한, random=무작위 단일선택 기대값(best-of-1 null).

지표(단위=타깃):
  mean   = 선택 pose 라벨 평균
  hit@th = 선택 pose가 성공(≥th)인 타깃 비율
  regret = oracle − 선택 (0=완벽)
  win    = 선택 > ipTM-선택 인 타깃 비율(>0.5면 ipTM을 이김)
진단 = 모델별 Spearman(iptm, label): 낮거나 모델마다 다르면 ipTM은 cross-model 비교 불가.

판정:
  배포가능(✅) 중 iptm보다 mean↑·regret↓·win>0.5 인 피처 있으면 → 무학습 재랭커 성립(GO).
  전부 iptm 근처거나 못 이기면 → 무학습 한계 → 학습 재랭커(Phase 2)로.

사용(stdlib only):
  python analyze_phase1_rerank.py                                   # dockq, scope=rung0(배포)
  python analyze_phase1_rerank.py --label recall                    # 에피토프 위치 라벨
  python analyze_phase1_rerank.py --scope all                       # 전 pose(깊이 섭동 포함=다양성↑)
  python analyze_phase1_rerank.py --csv 다른.csv --succ 0.23
"""
import argparse, csv, math
from collections import defaultdict

FEAT = ("dockq", "recall", "n_contact", "overrep", "pop_rank", "dcc_pop",
        "dcc_true", "true_rank", "iptm", "ptm", "plddt")

# (이름, 피처키, maximize, 배포가능)
SELECTORS = [
    ("oracle",     "__oracle__", None, False),
    ("random",     "__random__", None, False),
    ("iptm",       "iptm",       True, True),
    ("ptm",        "ptm",        True, True),
    ("plddt",      "plddt",      True, True),
    ("n_contact",  "n_contact",  True, True),
    ("dcc_pop",    "dcc_pop",    True, True),   # popular서 멀수록(off-hotspot) native일 것이라는 가설
    ("overrep_lo", "overrep",    False, True),  # popular 잔기 적게
    ("pop_rank_lo","pop_rank",   False, True),
]


def f(x):
    try:
        v = float(x)
        return None if math.isnan(v) else v
    except Exception:
        return None


def load(path):
    rows = []
    for r in csv.DictReader(open(path)):
        d = dict(target=r["target"], model=r["model"], rung=str(r.get("rung", "")), group=r.get("group", ""))
        for k in FEAT:
            d[k] = f(r.get(k, ""))
        rows.append(d)
    return rows


def spearman(pairs):
    xy = [(a, b) for a, b in pairs if a is not None and b is not None]
    n = len(xy)
    if n < 3:
        return None
    def ranks(vals):
        order = sorted(range(n), key=lambda i: vals[i])
        rk = [0.0] * n
        i = 0
        while i < n:
            j = i
            while j + 1 < n and vals[order[j + 1]] == vals[order[i]]:
                j += 1
            avg = (i + j) / 2.0 + 1
            for k in range(i, j + 1):
                rk[order[k]] = avg
            i = j + 1
        return rk
    xs = [a for a, _ in xy]; ys = [b for _, b in xy]
    rx = ranks(xs); ry = ranks(ys)
    mx = sum(rx) / n; my = sum(ry) / n
    num = sum((rx[i] - mx) * (ry[i] - my) for i in range(n))
    dx = math.sqrt(sum((rx[i] - mx) ** 2 for i in range(n)))
    dy = math.sqrt(sum((ry[i] - my) ** 2 for i in range(n)))
    return num / (dx * dy) if dx and dy else None


def pick(poses, feat, label, maximize):
    """argmax(또는 argmin) 피처 pose들의 라벨 평균(동점=무작위 tiebreak 기대값)."""
    cand = [(p[feat], p[label]) for p in poses if p.get(feat) is not None and p[label] is not None]
    if not cand:
        return None
    ext = max(c[0] for c in cand) if maximize else min(c[0] for c in cand)
    tied = [c[1] for c in cand if c[0] == ext]
    return sum(tied) / len(tied)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", default="results/pose_features.csv")
    ap.add_argument("--label", default="dockq", choices=["dockq", "recall"])
    ap.add_argument("--scope", default="rung0", choices=["rung0", "all"],
                    help="rung0=배포(full MSA)만 / all=전 pose(깊이 섭동 포함)")
    ap.add_argument("--targets", default="",
                    help="쉼표/공백 구분. 지정 시 그 복합체만(=케이스 스터디 렌즈). "
                         "blanket 49타깃은 rescue 없는 타깃에 희석됨 — 앵커만 보려면 여기 지정.")
    ap.add_argument("--detail", action="store_true", help="타깃별 선택값 표(특정 복합체 깊이 보기)")
    ap.add_argument("--succ", type=float, default=0.49)
    a = ap.parse_args()
    rows = load(a.csv)
    if a.scope == "rung0":
        rows = [r for r in rows if r["rung"] in ("0", "0.0")]
    if a.targets:
        want = set(a.targets.replace(",", " ").split())
        rows = [r for r in rows if r["target"] in want]
    if not rows:
        raise SystemExit("!! 해당 scope/targets에 pose 없음")
    lab = a.label

    print(f"== Phase 1 make-or-break | label={lab} | scope={a.scope} | 성공≥{a.succ} ==")
    print(f"[진단] 모델별 Spearman(iptm, {lab}) — 낮거나 모델마다 다르면 ipTM은 못 믿을 선택기:")
    for m in sorted({r["model"] for r in rows}):
        rs = spearman([(r["iptm"], r[lab]) for r in rows if r["model"] == m])
        print(f"    {m:9} rho = {rs:.3f}" if rs is not None else f"    {m:9} rho = NA")

    byt = defaultdict(list)
    for r in rows:
        byt[r["target"]].append(r)
    targets = sorted(byt)
    orc = {t: max((r[lab] for r in byt[t] if r[lab] is not None), default=None) for t in targets}
    ipk = {t: pick(byt[t], "iptm", lab, True) for t in targets}

    print(f"\n[선택기 비교] 단위=타깃 · 값=선택 pose의 {lab}")
    print(f"  {'selector':13}{'n':>4}{'mean':>8}{'hit@'+format(a.succ,'g'):>9}"
          f"{'hit@0.23':>9}{'regret':>8}{'win_vs_iptm':>12}  배포?")
    print("  " + "-" * 70)
    allpicks = {}
    for name, feat, mx, dep in SELECTORS:
        picks = {}
        for t in targets:
            if name == "oracle":
                picks[t] = orc[t]
            elif name == "random":
                vs = [r[lab] for r in byt[t] if r[lab] is not None]
                picks[t] = sum(vs) / len(vs) if vs else None
            else:
                picks[t] = pick(byt[t], feat, lab, mx)
        allpicks[name] = picks
        ok = [t for t in targets if picks[t] is not None]
        if not ok:
            continue
        vals = [picks[t] for t in ok]
        mean = sum(vals) / len(vals)
        h49 = sum(1 for v in vals if v >= a.succ) / len(vals)
        h23 = sum(1 for v in vals if v >= 0.23) / len(vals)
        reg = [orc[t] - picks[t] for t in ok if orc[t] is not None]
        mreg = sum(reg) / len(reg) if reg else float("nan")
        wt = [t for t in ok if ipk[t] is not None]
        win = (sum(1 for t in wt if picks[t] > ipk[t] + 1e-9) / len(wt)) if wt else float("nan")
        print(f"  {name:13}{len(ok):>4}{mean:8.3f}{h49:9.2f}{h23:9.2f}{mreg:8.3f}{win:12.2f}  {'✅' if dep else '—'}")

    if a.detail:
        cols = [c for c in ("oracle", "iptm", "ptm", "plddt", "n_contact",
                            "dcc_pop", "overrep_lo", "pop_rank_lo") if c in allpicks]
        print(f"\n[타깃별 상세] 값 = 각 선택기가 고른 pose의 {lab}  (oracle=천장, iptm보다 높은 배포피처가 그 복합체의 rescue)")
        print("  " + "target".ljust(11) + "".join(c[:8].rjust(9) for c in cols))
        for t in targets:
            def g(nm):
                v = allpicks[nm].get(t)
                return f"{v:.3f}" if v is not None else "-"
            print("  " + t.ljust(11) + "".join(g(c).rjust(9) for c in cols))

    print("\n판정 가이드:")
    print("  · 배포가능(✅) 피처 중 iptm행보다 mean↑ · regret↓ · win>0.50 인 것 → 무학습 재랭커 성립(GO).")
    print("  · 전부 iptm 근처거나 못 이김 → 무학습 한계 → 학습 재랭커(Phase 2, DeepRank-GNN-esm 등).")
    print("  ⚠️ dcc_pop/overrep/pop_rank는 popular set 있는 타깃(RBD/HA/Env)만 평가됨(n 열 참고).")


if __name__ == "__main__":
    main()
