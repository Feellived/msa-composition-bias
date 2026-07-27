#!/usr/bin/env python3
"""[기제 분석] 실패한 예측들은 한곳으로 모이는가 — '모델이 선호하는 잘못된 자리'가 있는지 본다.

왜 이 형태인가: "PDB에서 인기 있는 자리로 쏠린다"는 주장은 우리 데이터가 반박한다
(8ulr의 진짜 에피토프 CD4bs가 곧 인기 자리라, 성공할수록 인기 자리와 더 겹친다).
그래서 외부 인기도를 끌어오지 않고 **이 실험 안에서만** 묻는다:

  · 실패한 실행들의 예측 결합자리가 서로 겹치나?  → 겹치면 '선호하는 잘못된 자리'가 존재
  · 성공 실행들의 자리와 실패 실행들의 자리는 다른가?
  · 실패 자리는 진짜 에피토프에서 얼마나 떨어져 있나?

겹침은 자카드 지수(두 집합의 교집합 ÷ 합집합, 0~1)로 잰다.

단위 = 실행 1회(자세 5개 중 DockQ 최고인 것을 그 실행의 대표로).

사용(DockQ env):
  python dump_seedrep_full.py --data $DATA/compreps --only 8ulr_HL --csv-out results/compreps_8ulr_HL.csv
  python epitope_cluster.py --csv results/compreps_8ulr_HL.csv --data $DATA/compreps
"""
import argparse, csv, glob, json, os
import statistics as st
from collections import Counter, defaultdict
import pose_features as PF

SUCC = 0.49


def pred_epitope(cj, pose_path, cutoff):
    """pose_all_metrics와 같은 방식으로 '예측된 항원 접촉 잔기 집합'을 꺼낸다."""
    m = PF.er_load(pose_path)
    used, ag = set(), []
    for i, ref in enumerate(PF.antigen_refs(cj)):
        cid, _, rr = PF.er_best(m, ref, exclude=used)
        if cid is None:
            continue
        used.add(cid); ag.append((i, rr, ref))
    ab = []
    for ref in PF.antibody_refs(cj):
        cid, _, rr = PF.er_best(m, ref, exclude=used)
        if cid:
            used.add(cid); ab.extend(rr)
    if not ag or not ab:
        return None, None
    pred, dist, coord = PF.ES.scored_epitope_full(ag, ab, cutoff)
    return (set(pred) if pred else None), coord


def jac(a, b):
    if not a or not b:
        return float("nan")
    return len(a & b) / len(a | b)


def mean_pairwise(sets):
    vals = [jac(sets[i], sets[j]) for i in range(len(sets)) for j in range(i + 1, len(sets))]
    vals = [v for v in vals if v == v]
    return (st.mean(vals), len(vals)) if vals else (float("nan"), 0)


def consensus(sets, frac=0.5):
    if not sets:
        return set()
    c = Counter()
    for s in sets:
        c.update(s)
    need = frac * len(sets)
    return {r for r, n in c.items() if n >= need}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", required=True)
    ap.add_argument("--data", default=os.environ.get("DATA", "/mnt/data/admuser/msadepth") + "/compreps")
    ap.add_argument("--targets-dir", default="targets")
    ap.add_argument("--cutoff", type=float, default=5.0)
    ap.add_argument("--out", default="")
    a = ap.parse_args()

    rows = list(csv.DictReader(open(a.csv)))
    if not rows:
        raise SystemExit("!! CSV가 비었음")
    tgt, model = rows[0]["target"], rows[0]["model"]
    depth = rows[0].get("depth", "")
    cj = json.load(open(os.path.join(a.targets_dir, tgt, "chains.json")))
    native = os.path.join(a.targets_dir, tgt, "native.cif")
    tr = PF.native_true(cj, native, a.cutoff)
    if tr is None:
        raise SystemExit("!! native epitope 계산 실패")
    true = set(tr[0])

    # 실행별 대표 자세(= DockQ 최고)
    best = {}
    for r in rows:
        try:
            q = float(r["dockq"])
        except Exception:
            continue
        s = r["seed"]
        if s not in best or q > best[s][0]:
            best[s] = (q, r["pose"], float(r.get("recall") or "nan"))

    base = os.path.join(a.data, "seedrep_cand", model, tgt, depth)
    recs = []
    for s, (q, pose, rc) in sorted(best.items()):
        hits = glob.glob(os.path.join(base, s, "results", "**", pose), recursive=True)
        if not hits:
            print(f"  ! {s}: 자세 파일 못 찾음({pose})"); continue
        ep, _ = pred_epitope(cj, hits[0], a.cutoff)
        if not ep:
            print(f"  ! {s}: 결합자리 계산 실패"); continue
        recs.append(dict(run=s, dockq=q, recall=rc, ep=ep, ok=q >= SUCC))
    if not recs:
        raise SystemExit("!! 계산된 실행이 없음")

    ok = [r for r in recs if r["ok"]]
    ng = [r for r in recs if not r["ok"]]
    print(f"■ {tgt} · {model} · {depth}   실행 {len(recs)}개 (성공 {len(ok)} · 실패 {len(ng)})")
    print(f"  진짜 결합자리 잔기 {len(true)}개\n")

    mo, no_ = mean_pairwise([r["ep"] for r in ok])
    mn, nn_ = mean_pairwise([r["ep"] for r in ng])
    cross = [jac(x["ep"], y["ep"]) for x in ok for y in ng]
    cross = [v for v in cross if v == v]
    # 우연 수준: 같은 크기 접촉면을 항원 표면에서 무작위로 잡았을 때의 기대 겹침
    nag = sum(len(c.get("seq", "")) for c in cj.get("chains", []) if c.get("role") == "antigen") or 0
    def chance(sets):
        if not sets or not nag:
            return float("nan")
        m = st.mean(len(x) for x in sets)
        inter = m * m / nag
        return inter / (2 * m - inter) if (2 * m - inter) > 0 else float("nan")
    ch_ok, ch_ng = chance([r["ep"] for r in ok]), chance([r["ep"] for r in ng])

    print(f"[예측 자리끼리 얼마나 겹치나 — 자카드 0~1]   (항원 {nag}잔기 기준 우연 수준과 비교)")
    def x(v, c):
        return f"  (우연 {c:.3f}의 {v/c:.1f}배)" if c == c and c > 0 else ""
    print(f"  성공 실행끼리   {mo:.3f}  (쌍 {no_}개){x(mo, ch_ok)}")
    print(f"  실패 실행끼리   {mn:.3f}  (쌍 {nn_}개){x(mn, ch_ng)}   ← 우연보다 높으면 '선호하는 잘못된 자리' 존재")
    print(f"  성공 vs 실패    {(st.mean(cross) if cross else float('nan')):.3f}  ← 낮으면 서로 다른 자리\n")

    # 흔한 자리(그 항원에서 항체가 자주 붙는 부위)와의 비교 — B군에서 '편향 이탈'을 직접 본다
    try:
        pop = set(PF.ES.popular_refset(cj, cj.get("antigen_grp", "")) or [])
    except Exception:
        pop = set()

    cs, cn = consensus([r["ep"] for r in ok]), consensus([r["ep"] for r in ng])
    print(f"  접촉면 평균 크기  성공 {st.mean(len(r['ep']) for r in ok) if ok else float('nan'):.0f}잔기"
          f" · 실패 {st.mean(len(r['ep']) for r in ng) if ng else float('nan'):.0f}잔기"
          f"   ← 실패가 훨씬 넓으면 '펼쳐 붙음'\n")
    print("[합의 자리 — 그 무리의 절반 이상 실행에 등장한 잔기]")
    print(f"  성공 합의자리 {len(cs):3d}개 · 진짜와 겹침 {jac(cs,true):.3f} · 진짜 잔기 포함률 {(len(cs&true)/len(true) if true else 0):.3f}")
    print(f"  실패 합의자리 {len(cn):3d}개 · 진짜와 겹침 {jac(cn,true):.3f} · 진짜 잔기 포함률 {(len(cn&true)/len(true) if true else 0):.3f}")
    print(f"  성공 합의자리 vs 실패 합의자리 겹침 {jac(cs,cn):.3f}\n")

    if pop:
        ab = str(cj.get("AB", "?"))
        print(f"[흔한 자리와의 비교]  이 복합체는 {ab}군"
              + ("(진짜 자리 = 흔한 자리 → '편향 이탈' 서사 부적합)" if ab == "A"
                 else "(진짜 자리가 흔한 자리와 다름 → '편향 이탈'을 직접 볼 수 있음)" if ab == "B" else ""))
        print(f"  흔한 자리 {len(pop)}잔기 · 진짜 자리와 겹침 {jac(pop, true):.3f}")
        mo_p = st.mean([jac(r["ep"], pop) for r in ok]) if ok else float("nan")
        mn_p = st.mean([jac(r["ep"], pop) for r in ng]) if ng else float("nan")
        print(f"  성공 실행의 예측 ↔ 흔한 자리  {mo_p:.3f}")
        print(f"  실패 실행의 예측 ↔ 흔한 자리  {mn_p:.3f}"
              + ("   ← 실패가 더 높으면 '흔한 자리로 쏠렸다'" if mn_p > mo_p else "   (실패가 더 낮음 = 흔한 자리 쏠림 아님)"))
        print()

    verdict = []
    if mn >= 0.5 and mn > mo * 0.8:
        verdict.append("실패가 한곳으로 모임 = 모델이 선호하는 잘못된 자리가 있다")
    elif mn >= 0.3:
        verdict.append("실패가 어느 정도 모임(중간)")
    else:
        verdict.append("실패가 흩어짐 = 특정 자리로의 쏠림이라 보기 어렵다")
    if jac(cs, cn) < 0.2:
        verdict.append("성공/실패가 서로 다른 자리")
    print("판정: " + " · ".join(verdict))

    out = a.out or f"results/epitope_cluster_{tgt}.csv"
    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["run", "dockq", "recall", "n_pred", "jac_to_true",
                    "jac_to_fail_consensus", "jac_to_succ_consensus", "jac_to_popular", "success"])
        for r in sorted(recs, key=lambda x: -x["dockq"]):
            w.writerow([r["run"], round(r["dockq"], 3), round(r["recall"], 3) if r["recall"] == r["recall"] else "",
                        len(r["ep"]), round(jac(r["ep"], true), 3),
                        round(jac(r["ep"], cn), 3), round(jac(r["ep"], cs), 3),
                        (round(jac(r["ep"], pop), 3) if pop else ""), int(r["ok"])])
    print(f"\n→ {out}  (그림용: '실패 자리와의 겹침' vs DockQ 로 이동을 그릴 수 있음)")


if __name__ == "__main__":
    main()
