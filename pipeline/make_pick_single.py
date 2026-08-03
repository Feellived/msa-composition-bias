#!/usr/bin/env python3
"""[①] 후보가 1개뿐이라 선택 목록에서 빠진 타깃을 되살린다.

eval_selectors.py 는 후보가 2개 이상인 타깃만 다룬다 — "고를 것이 없다"는 이유였다.
그런데 **후보가 1개라는 것은 조성들이 완전히 합의했다는 뜻**이고, 후보 개수 규칙에서는
오히려 가장 확신 있는 경우다. 안 돌린 것이 착오였다.

  · 고를 필요가 없으므로 선택기도 필요 없다. 그 하나를 그대로 쓴다.
  · 정답을 쓰지 않는다 — 덮음·정밀도·F1 칸은 **비워 둔다**(어차피 실행에 안 쓰인다).
    JSON 안에는 그 값이 들어 있지만 여기 옮기면 "정답을 본 목록"이 되어 버린다.

⚠️ 원본 pick_abepi_max.csv 는 **건드리지 않는다.** 지금까지 나온 23종 결과의 근거 파일이라
   덮으면 재현이 깨진다. 원본 + 새 줄을 합친 **새 파일**을 만든다.

사용 (pipeline/ 에서):
  python -u make_pick_single.py
  → results/honest/pick_abepi_max_plus.csv

그다음:
  PICKCSV=$PWD/results/honest/pick_abepi_max_plus.csv \
  ONLY="..." ARMS="noconstraint ours" bash scripts/run_honest_guided.sh --apply
"""
import argparse, csv, glob, json, os

COLS = ["target", "selector", "cand", "n_comp", "n_res",
        "true_covered", "precision", "f1", "from_full_msa"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sites", default="results/honest")
    ap.add_argument("--pick", default="results/honest/pick_abepi_max.csv")
    ap.add_argument("--out", default="results/honest/pick_abepi_max_plus.csv")
    ap.add_argument("--max-cand", type=int, default=1,
                    help="후보가 이 개수 이하인 타깃만 되살린다")
    a = ap.parse_args()

    if os.path.abspath(a.out) == os.path.abspath(a.pick):
        raise SystemExit("!! 원본을 덮으려 한다. --out 을 다른 이름으로 줄 것")
    if not os.path.exists(a.pick):
        raise SystemExit(f"!! 원본 목록이 없다: {a.pick}")

    orig = list(csv.DictReader(open(a.pick)))
    have = {r["target"] for r in orig}
    print(f"원본 {len(orig)}종: {a.pick}")

    added, skipped = [], []
    for jf in sorted(glob.glob(os.path.join(a.sites, "sites_*.json"))):
        d = json.load(open(jf))
        t = d["target"]
        if t in have:
            continue
        cands = d.get("candidates", [])
        if not cands:
            skipped.append(f"{t}(후보 0)"); continue
        if len(cands) > a.max_cand:
            skipped.append(f"{t}(후보 {len(cands)})"); continue
        c = cands[0]
        added.append({"target": t, "selector": "single", "cand": c["cand"],
                      "n_comp": c["n_comp"], "n_res": len(c["residues"]),
                      "true_covered": "", "precision": "", "f1": "",
                      "from_full_msa": c.get("from_full_msa", "")})

    if skipped:
        print(f"  ! 되살리지 않음 {len(skipped)}종 — {', '.join(skipped)}")
    if not added:
        raise SystemExit("!! 되살릴 타깃이 없다 — 이미 다 들어 있거나 후보 수 조건에 안 맞는다")

    print(f"\n되살린 {len(added)}종 (후보가 {a.max_cand}개 이하 · 조성이 완전히 합의한 경우)")
    print(f"  {'타깃':<12}{'후보':>4}{'조성수':>6}{'잔기':>6}{'원래MSA포함':>12}")
    for r in added:
        print(f"  {r['target']:<12}{r['cand']:>4}{r['n_comp']:>6}{r['n_res']:>6}"
              f"{str(r['from_full_msa']):>12}")

    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS)
        w.writeheader()
        for r in orig:
            w.writerow({k: r.get(k, "") for k in COLS})
        w.writerows(added)
    print(f"\n→ {a.out}   원본 {len(orig)} + 추가 {len(added)} = {len(orig)+len(added)}종")
    print("   ⚠️ 원본은 그대로 두었다.")


if __name__ == "__main__":
    main()
