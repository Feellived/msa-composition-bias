#!/usr/bin/env python3
"""[진단] 우리 자리를 제약으로 준 것이 해로웠던 타깃의 원인을 셋으로 가른다.

  (i)  선택 실패 — 후보 목록에 정답을 잘 덮는 후보가 **있었는데** 다른 것을 골랐다
  (ii) 생성 실패 — 후보 목록에 정답을 덮는 후보가 **아예 없었다**
  (iii) 제약 방해 — 제대로 골랐는데도 DockQ 가 내려갔다

셋은 처방이 완전히 다르다. (i) 은 선택기 문제, (ii) 는 생성 폭 문제,
(iii) 은 **애초에 개입하면 안 되는 구간**이라는 뜻이다.

⚠️ 판정 문턱(`--cov-ok`, `--sel-gap`, `--eps`)은 **사람이 정한 값**이다. 결과가 문턱에
   민감하면 그 판정은 보고서에 쓸 수 없다. 그래서 `--sensitivity` 로 문턱을 격자로
   흔들어 **종수가 뒤집히는지** 먼저 확인한다.

⚠️ 후보의 '덮음'만 보면 안 된다. 덮음 1.0 인데 잔기 107개짜리 후보가 실제로 있었다
   (8sis_HL). 정답을 다 덮지만 그 네 배 넓이를 덮는 것이고, 그런 제약은 모델을
   넓은 영역 아무 데나 자신 있게 붙여놓는다. 그래서 **정밀도를 항상 같이 본다.**

사용 (CSV/JSON 만 읽는다):
  python -u diagnose_harm.py --sensitivity     # 문턱을 흔들어 본다 (먼저 이것부터)
  python -u diagnose_harm.py
"""
import argparse, csv, glob, json, os
import statistics as st

AB = os.path.expanduser("~/projects/bk21-antibody-ml/consensus_docking")
CAUSES = ("이득", "무변화", "제약방해", "선택실패", "생성실패")


def load_pick(path):
    """선택기가 고른 후보 번호. 열 이름이 판마다 달라 유연하게 찾는다."""
    rows = list(csv.DictReader(open(path)))
    if not rows:
        raise SystemExit(f"!! 비어 있음: {path}")
    cols = list(rows[0])
    tcol = next((c for c in cols if c.lower() == "target"), None)
    ccol = next((c for c in cols if c.lower() in
                 ("cand", "candidate", "pick", "cand_id", "chosen", "sel")), None)
    if not tcol or not ccol:
        raise SystemExit(f"!! pick CSV 에서 target/cand 열을 못 찾음: {cols}")
    out = {}
    for r in rows:
        if r.get(ccol) in (None, ""):
            continue
        try:
            out[r[tcol]] = int(float(r[ccol]))
        except ValueError:
            print(f"  ! {r[tcol]}: 후보 번호를 숫자로 못 읽었다({r[ccol]!r}) — 건너뜀")
    return out


def load_dockq(path):
    """(타깃, 팔) → 행. 같은 칸이 두 번 나오면 알리고 뒤엣것을 쓴다."""
    rows = list(csv.DictReader(open(path)))
    if not rows:
        raise SystemExit(f"!! 비어 있음: {path}")
    need = {"target", "arm", "dockq_max"}
    miss = need - set(rows[0])
    if miss:
        raise SystemExit(f"!! {path} 에 열이 없다: {sorted(miss)}\n   있는 열 = {list(rows[0])}")
    out, dup = {}, []
    for r in rows:
        k = (r["target"], r["arm"])
        if k in out:
            dup.append(k)
        out[k] = r
    if dup:
        print(f"  ! 중복 행 {len(dup)}개 — 예: {dup[:3]}  (모델이 섞였을 수 있다)")
    return out


def f(r, col):
    try:
        return float(r[col])
    except (TypeError, ValueError, KeyError):
        return float("nan")


def band(q):
    if q != q:
        return "?"
    if q < 0.13:
        return "매우낮음"
    if q < 0.23:
        return "애매"
    return "이미성공"


def classify(rec, cov_ok, sel_gap, eps):
    """원인 판정. 순서가 중요하다 — 후보가 없었으면 선택기 탓을 할 수 없다."""
    if rec["best_cov"] < cov_ok:
        return "생성실패"
    if rec["cov"] < rec["best_cov"] - sel_gap:
        return "선택실패"
    if rec["delta"] < -eps:
        return "제약방해"
    if rec["delta"] > eps:
        return "이득"
    return "무변화"


def collect(a):
    """타깃마다 판정에 필요한 값을 모은다(문턱과 무관한 부분만)."""
    pick = load_pick(a.pick or os.path.join(a.sites, "pick_abepi_max.csv"))
    dq = load_dockq(a.dockq)
    print(f"선택 기록 {len(pick)}종 · DockQ 칸 {len(dq)}개")

    recs = []
    for jf in sorted(glob.glob(os.path.join(a.sites, "sites_*.json"))):
        d = json.load(open(jf))
        t = d["target"]
        nc, ou = dq.get((t, "noconstraint")), dq.get((t, "ours"))
        if not nc or not ou:
            print(f"  ! {t}: 무제약/우리자리 행이 없어 건너뜀"); continue
        cands = {c["cand"]: c for c in d["candidates"]}
        ci = pick.get(t)
        if ci not in cands:
            print(f"  ! {t}: 고른 후보 {ci} 가 후보 목록에 없다 — 건너뜀"); continue
        ch = cands[ci]
        best = max(cands.values(), key=lambda c: c["true_covered"])
        q0, q1 = f(nc, "dockq_max"), f(ou, "dockq_max")
        recs.append(dict(
            target=t, band=band(q0), dq_no=round(q0, 3), dq_ours=round(q1, 3),
            delta=round(q1 - q0, 3), cand=ci, n_cand=len(cands),
            n_res=len(ch["residues"]), n_true=d.get("n_true_res", ""),
            n_comp=ch["n_comp"], cov=round(ch["true_covered"], 3),
            prec=round(ch["precision"], 3), best_cand=best["cand"],
            best_cov=round(best["true_covered"], 3),
            from_full_msa=ch.get("from_full_msa", ""),
            rec_no=round(f(nc, "recall_max"), 3), rec_ours=round(f(ou, "recall_max"), 3),
            perm_p=d.get("perm_p", "")))
    if not recs:
        raise SystemExit("!! 판정된 타깃이 없다 — 경로를 확인할 것")
    return recs


def sensitivity(recs, a):
    """문턱을 격자로 흔들어 종수와 핵심 주장이 버티는지 본다."""
    print("\n" + "=" * 96)
    print("  문턱 민감도 — 사람이 정한 값이 결과를 만들고 있지 않은지 확인한다")
    print("=" * 96)
    big = sorted(recs, key=lambda r: r["delta"])[:6]
    print(f"  '크게 망가진 6종' = {', '.join(r['target'] for r in big)}"
          f"   (변화 {big[0]['delta']} ~ {big[-1]['delta']})\n")
    print(f"  {'덮음문턱':>8}{'차이문턱':>8}{'변화문턱':>8}  "
          + "".join(f"{c:>7}" for c in CAUSES) + "   큰손해6종이 전부 선택실패인가")
    base_ok = True
    for cov_ok in (0.3, 0.4, 0.5, 0.6, 0.7):
        for sel_gap in (0.10, 0.15, 0.20, 0.25):
            for eps in (0.02,):
                cnt = {c: 0 for c in CAUSES}
                for r in recs:
                    cnt[classify(r, cov_ok, sel_gap, eps)] += 1
                allsel = all(classify(r, cov_ok, sel_gap, eps) == "선택실패" for r in big)
                nsel = sum(1 for r in big if classify(r, cov_ok, sel_gap, eps) == "선택실패")
                base_ok &= allsel
                print(f"  {cov_ok:>8}{sel_gap:>8}{eps:>8}  "
                      + "".join(f"{cnt[c]:>7}" for c in CAUSES)
                      + f"   {'예' if allsel else f'아니오 ({nsel}/6)'}")
    print()
    if base_ok:
        print("  ⭐ 모든 문턱 조합에서 '큰 손해 6종은 전부 선택 실패'가 유지된다 — 문턱에 안 기댄다")
    else:
        print("  ⚠️ 일부 문턱에서 뒤집힌다 — 이 주장은 문턱에 기대고 있다. 보고서에 조건을 명시할 것")
    print("  ※ 종수 자체는 문턱에 따라 움직이는 것이 정상이다. 봐야 할 것은 주장이 버티는가이다")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sites", default="results/honest",
                    help="정직한 자리 JSON 과 pick CSV 가 있는 폴더")
    ap.add_argument("--pick", default="", help="기본 = <sites>/pick_abepi_max.csv")
    ap.add_argument("--dockq", default=os.path.join(AB, "results/dockq_honest_merged.csv"))
    ap.add_argument("--cov-ok", type=float, default=0.5,
                    help="후보가 정답을 이만큼 덮으면 '쓸 만한 후보'로 본다")
    ap.add_argument("--sel-gap", type=float, default=0.15,
                    help="고른 후보가 최선 후보보다 이만큼 덜 덮으면 선택 실패로 본다")
    ap.add_argument("--eps", type=float, default=0.02, help="DockQ 변화가 이 미만이면 무변화")
    ap.add_argument("--sensitivity", action="store_true", help="문턱을 격자로 흔들어 본다")
    ap.add_argument("--out", default="results/harm_diag.csv")
    a = ap.parse_args()

    recs = collect(a)
    if a.sensitivity:
        sensitivity(recs, a)
        return

    for r in recs:
        r["why"] = classify(r, a.cov_ok, a.sel_gap, a.eps)
    recs.sort(key=lambda r: r["delta"])

    W = ("target", "band", "why", "dq_no", "dq_ours", "delta",
         "cov", "prec", "n_res", "n_true", "best_cov", "n_cand", "rec_no", "rec_ours")
    hdr = {"target": "타깃", "band": "구간", "why": "원인", "dq_no": "무제약", "dq_ours": "우리",
           "delta": "변화", "cov": "덮음", "prec": "정밀도", "n_res": "후보잔기",
           "n_true": "정답잔기", "best_cov": "최선덮음", "n_cand": "후보수",
           "rec_no": "무제약자리", "rec_ours": "우리자리"}
    print("\n" + "=" * 118)
    print("  타깃별 판정 — 변화량 오름차순 (덮음·정밀도는 **고른 후보**의 값이다)")
    print("=" * 118)
    print("  ".join(f"{hdr[c]:>8}" for c in W))
    for r in recs:
        print("  ".join(f"{r[c]!s:>8}" for c in W))

    print("\n" + "-" * 60)
    print("  원인별 종수")
    print("-" * 60)
    for w in CAUSES:
        s = [r for r in recs if r["why"] == w]
        if s:
            print(f"  {w:<6} {len(s):>2}종   {', '.join(r['target'] for r in s)}")

    unk = [r for r in recs if r["band"] == "?"]
    if unk:
        print(f"\n  ! 무제약 DockQ 를 못 읽어 구간을 못 정한 {len(unk)}종: "
              f"{', '.join(r['target'] for r in unk)} — 아래 구간 표에서 빠진다")

    print("\n" + "-" * 78)
    print(f"  무제약 DockQ 구간 × 방향  (변화 > +{a.eps} = 오름, < -{a.eps} = 내림)")
    print("-" * 78)
    print(f"  {'구간':<10}{'종수':>5}{'오름':>6}{'내림':>6}{'무변화':>8}   {'중앙 변화량':>10}")
    for b in ("매우낮음", "애매", "이미성공"):
        s = [r for r in recs if r["band"] == b]
        if not s:
            continue
        up = sum(1 for r in s if r["delta"] > a.eps)
        dn = sum(1 for r in s if r["delta"] < -a.eps)
        print(f"  {b:<10}{len(s):>5}{up:>6}{dn:>6}{len(s)-up-dn:>8}"
              f"   {st.median(r['delta'] for r in s):>10.3f}")

    print("\n" + "-" * 78)
    print("  ⭐ 고른 후보의 덮음 × DockQ 변화")
    print("-" * 78)
    print(f"  {'덮음 구간':<14}{'종수':>5}{'오름':>6}{'내림':>6}   {'중앙 변화량':>10}   {'중앙 정밀도':>10}")
    for lo, hi, nm in ((0.0, 0.3, "낮음 <0.3"), (0.3, 0.7, "중간 0.3~0.7"), (0.7, 1.01, "높음 ≥0.7")):
        s = [r for r in recs if lo <= r["cov"] < hi]
        if not s:
            continue
        up = sum(1 for r in s if r["delta"] > a.eps)
        dn = sum(1 for r in s if r["delta"] < -a.eps)
        print(f"  {nm:<14}{len(s):>5}{up:>6}{dn:>6}"
              f"   {st.median(r['delta'] for r in s):>10.3f}"
              f"   {st.median(r['prec'] for r in s):>10.3f}")

    print("\n  ※ 판정 문턱이 결과를 만들고 있지 않은지 --sensitivity 로 확인할 것")

    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(recs[0]))
        w.writeheader(); w.writerows(recs)
    print(f"\n→ {a.out}  ({len(recs)}종)")


if __name__ == "__main__":
    main()
