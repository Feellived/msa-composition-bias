#!/usr/bin/env python3
"""[진단] 우리 자리를 제약으로 준 것이 해로웠던 타깃의 원인을 셋으로 가른다.

  (i)  선택 실패 — 후보 목록에 정답을 잘 덮는 후보가 **있었는데** 다른 것을 골랐다
  (ii) 생성 실패 — 후보 목록에 정답을 덮는 후보가 **아예 없었다**
  (iii) 제약 방해 — 제대로 골랐는데도 DockQ 가 내려갔다
        (무제약 예측이 이미 정답 근처였는데, 우리 제약이 정답보다 거칠어 되레 밀어냈다)

셋은 처방이 완전히 다르다. (i) 은 선택기 문제, (ii) 는 생성 폭 문제,
(iii) 은 **애초에 개입하면 안 되는 구간**이라는 뜻이다.

같이 뽑는 그림 데이터 (results/harm_diag.csv):
   가로 = 고른 후보가 정답을 덮는 비율 · 세로 = DockQ 변화량
   점 크기 = 후보 잔기 수 · 색 = 무제약 DockQ 구간

⚠️ 후보의 '덮음'만 보면 안 된다. 덮음 1.0 인데 잔기 107개짜리 후보가 실제로 있었다
   (8sis_HL). 정답을 다 덮지만 그 네 배 넓이를 덮는 것이고, 그런 제약은 모델을
   넓은 영역 아무 데나 자신 있게 붙여놓는다. 그래서 **정밀도를 항상 같이 본다.**

사용 (DockQ env 필요 없음 — CSV/JSON 만 읽는다):
  python -u diagnose_harm.py
  python -u diagnose_harm.py --cov-ok 0.4 --out results/harm_diag.csv
"""
import argparse, csv, glob, json, os

AB = os.path.expanduser("~/projects/bk21-antibody-ml/consensus_docking")


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
    return {r[tcol]: int(float(r[ccol])) for r in rows if r.get(ccol) not in (None, "")}


def load_dockq(path):
    """(타깃, 팔) → 행. 같은 칸이 두 번 나오면 알리고 뒤엣것을 쓴다."""
    out, dup = {}, []
    for r in csv.DictReader(open(path)):
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
    ap.add_argument("--out", default="results/harm_diag.csv")
    a = ap.parse_args()

    pick = load_pick(a.pick or os.path.join(a.sites, "pick_abepi_max.csv"))
    dq = load_dockq(a.dockq)
    print(f"선택 기록 {len(pick)}종 · DockQ 칸 {len(dq)}개\n")

    rows = []
    for jf in sorted(glob.glob(os.path.join(a.sites, "sites_*.json"))):
        d = json.load(open(jf))
        t = d["target"]
        nc, ou = dq.get((t, "noconstraint")), dq.get((t, "ours"))
        if not nc or not ou:
            print(f"  ! {t}: 무제약/우리자리 행이 없어 건너뜀")
            continue
        cands = {c["cand"]: c for c in d["candidates"]}
        ci = pick.get(t)
        if ci not in cands:
            print(f"  ! {t}: 고른 후보 {ci} 가 후보 목록에 없다 — 건너뜀")
            continue
        ch = cands[ci]
        best = max(cands.values(), key=lambda c: c["true_covered"])

        q0, q1 = f(nc, "dockq_max"), f(ou, "dockq_max")
        dlt = q1 - q0
        cov, prec = ch["true_covered"], ch["precision"]

        # 원인 판정 — 순서가 중요하다. 후보가 없었으면 선택기 탓을 할 수 없다.
        if best["true_covered"] < a.cov_ok:
            why = "생성실패"
        elif cov < best["true_covered"] - a.sel_gap:
            why = "선택실패"
        elif dlt < -a.eps:
            why = "제약방해"
        elif dlt > a.eps:
            why = "이득"
        else:
            why = "무변화"

        rows.append(dict(
            target=t, band=band(q0), why=why,
            dq_no=round(q0, 3), dq_ours=round(q1, 3), delta=round(dlt, 3),
            cand=ci, n_cand=len(cands), n_res=len(ch["residues"]),
            n_true=d.get("n_true_res", ""), n_comp=ch["n_comp"],
            cov=round(cov, 3), prec=round(prec, 3),
            best_cand=best["cand"], best_cov=round(best["true_covered"], 3),
            from_full_msa=ch.get("from_full_msa", ""),
            rec_no=round(f(nc, "recall_max"), 3), rec_ours=round(f(ou, "recall_max"), 3),
            perm_p=d.get("perm_p", "")))

    if not rows:
        raise SystemExit("!! 판정된 타깃이 없다 — 경로를 확인할 것")

    rows.sort(key=lambda r: r["delta"])
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
    for r in rows:
        print("  ".join(f"{r[c]!s:>8}" for c in W))

    print("\n" + "-" * 60)
    print("  원인별 종수")
    print("-" * 60)
    for w in ("이득", "무변화", "제약방해", "선택실패", "생성실패"):
        s = [r for r in rows if r["why"] == w]
        if s:
            print(f"  {w:<6} {len(s):>2}종   {', '.join(r['target'] for r in s)}")

    print("\n" + "-" * 78)
    print("  무제약 DockQ 구간 × 방향  (변화 > +0.02 = 오름, < -0.02 = 내림)")
    print("-" * 78)
    print(f"  {'구간':<10}{'종수':>5}{'오름':>6}{'내림':>6}{'무변화':>8}   {'중앙 변화량':>10}")
    for b in ("매우낮음", "애매", "이미성공"):
        s = [r for r in rows if r["band"] == b]
        if not s:
            continue
        up = sum(1 for r in s if r["delta"] > a.eps)
        dn = sum(1 for r in s if r["delta"] < -a.eps)
        md = sorted(r["delta"] for r in s)[len(s) // 2]
        print(f"  {b:<10}{len(s):>5}{up:>6}{dn:>6}{len(s)-up-dn:>8}   {md:>10.3f}")

    # 핵심 확인 — 덮음이 높을수록 이득인가. 덮음 구간별 변화량 중앙값.
    print("\n" + "-" * 78)
    print("  ⭐ 고른 후보의 덮음 × DockQ 변화  — 이 표가 (iii) 가설의 핵심이다")
    print("-" * 78)
    print(f"  {'덮음 구간':<14}{'종수':>5}{'오름':>6}{'내림':>6}   {'중앙 변화량':>10}   {'중앙 정밀도':>10}")
    for lo, hi, nm in ((0.0, 0.3, "낮음 <0.3"), (0.3, 0.7, "중간 0.3~0.7"), (0.7, 1.01, "높음 ≥0.7")):
        s = [r for r in rows if lo <= r["cov"] < hi]
        if not s:
            continue
        up = sum(1 for r in s if r["delta"] > a.eps)
        dn = sum(1 for r in s if r["delta"] < -a.eps)
        md = sorted(r["delta"] for r in s)[len(s) // 2]
        mp = sorted(r["prec"] for r in s)[len(s) // 2]
        print(f"  {nm:<14}{len(s):>5}{up:>6}{dn:>6}   {md:>10.3f}   {mp:>10.3f}")

    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"\n→ {a.out}  ({len(rows)}종)")


if __name__ == "__main__":
    main()
