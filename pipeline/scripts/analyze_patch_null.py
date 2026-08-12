#!/usr/bin/env python3
"""[기준선 근거] "자리 겹침 40%"가 우연히 넘기 어려운 선인지 실제로 계산한다.

문제:
  결합자리 겹침(recall) 기준을 0.4로 두었는데, 그 근거가 "항원 표면은 넓고 결합자리는
  좁으니 어렵다"는 말뿐이었다. 몇 퍼센트인지 재지 않으면 기준이 임의로 보인다.

방법 (귀무분포를 손으로 만들지 않고 구조에서 직접):
  ① native 에서 항원만 떼어 표면 잔기를 구한다(Shrake-Rupley SASA ≥ --sasa Å²).
  ② 진짜 결합자리 = native 항체와 --cutoff Å 이내 항원 잔기.
  ③ **표면 잔기 하나하나를 중심으로**, 결합자리와 같은 개수의 이웃 표면 잔기를 모아
     '항체 크기의 패치'를 만든다 = 임의 위치에 항체를 놓아 본 것.
  ④ 그 패치가 진짜 결합자리를 몇 % 덮는지 계산 → 표면 전체에 대한 분포.
  → "40%를 넘는 위치가 표면의 몇 %인가" = 우연히 넘을 확률.

같은 크기 패치를 쓰므로 '크게 잡아서 덮었다'가 배제된다(데모의 sizematch 팔과 같은 취지).
접촉 잔기가 실제로 이어져 있다는 점(연속성)까지 반영되므로, 잔기를 무작위로 뽑는
초기하분포보다 현실에 가깝다. 비교용으로 초기하 값도 함께 낸다.

사용 (biopython+scipy env):
  python analyze_patch_null.py                       # targets/ 전체
  python analyze_patch_null.py --targets "8k3k_D 8ulr_HL" --thresh 0.4
"""
import argparse, csv, glob, json, os
import numpy as np
from Bio.PDB import MMCIFParser, PDBParser
from Bio.PDB.SASA import ShrakeRupley

try:
    from scipy.stats import hypergeom
except Exception:
    hypergeom = None


def load(p):
    par = MMCIFParser(QUIET=True) if p.endswith(".cif") else PDBParser(QUIET=True)
    return par.get_structure("x", p)[0]


def heavy(res):
    return np.array([a.coord for a in res if a.element != "H"], dtype=float)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=".", help="targets/ 가 있는 폴더")
    ap.add_argument("--targets", default="", help="비우면 targets/* 전부")
    ap.add_argument("--cutoff", type=float, default=5.0, help="결합자리 접촉 기준 Å")
    ap.add_argument("--sasa", type=float, default=10.0, help="표면 잔기 판정 SASA Å²")
    ap.add_argument("--thresh", default="0.4 0.2 0.6", help="검사할 겹침 기준들")
    ap.add_argument("--out", default="results/patch_null.csv")
    a = ap.parse_args()

    ths = [float(x) for x in a.thresh.replace(",", " ").split()]
    ts = a.targets.replace(",", " ").split() or sorted(
        os.path.basename(d) for d in glob.glob(os.path.join(a.root, "targets", "*"))
        if os.path.isdir(d))
    sr = ShrakeRupley()
    rows = []

    hdr = f"{'복합체':<11}{'표면':>6}{'결합자리':>8}{'비율':>8}" + \
          "".join(f"{f'≥{t:.0%} 위치':>12}" for t in ths)
    print(hdr); print("-" * 78)

    for t in ts:
        nat = os.path.join(a.root, "targets", t, "native.cif")
        cjp = os.path.join(a.root, "targets", t, "chains.json")
        if not (os.path.exists(nat) and os.path.exists(cjp)):
            continue
        cj = json.load(open(cjp))
        agc = cj["antigen"]; agc = [agc] if isinstance(agc, str) else list(agc)
        abc = cj["antibody"]; abc = [abc] if isinstance(abc, str) else list(abc)
        try:
            m = load(nat)
            ag = [r for c in agc if c in m for r in m[c] if r.id[0] == " " and "CA" in r]
            ab = [r for c in abc if c in m for r in m[c] if r.id[0] == " "]
            if len(ag) < 20 or not ab:
                print(f"{t:<11} 사슬 확인 필요 — 건너뜀"); continue

            # ② 진짜 결합자리
            abx = np.vstack([heavy(r) for r in ab])
            true = set()
            for i, r in enumerate(ag):
                v = heavy(r)
                if v.size and np.min(np.linalg.norm(v[:, None, :] - abx[None, :, :], axis=2)) <= a.cutoff:
                    true.add(i)
            if len(true) < 5:
                print(f"{t:<11} 결합자리 {len(true)}잔기 — 건너뜀"); continue

            # ① 표면 잔기 (항체를 뗀 항원 단독 상태에서)
            import copy
            solo = copy.deepcopy(m)
            for c in list(solo):
                if c.id not in agc:
                    solo.detach_child(c.id)
            sr.compute(solo, level="R")
            acc = {}
            for c in solo:
                for r in c:
                    if r.id[0] == " ":
                        acc[(c.id, r.id)] = r.sasa
            surf = [i for i, r in enumerate(ag)
                    if acc.get((r.get_parent().id, r.id), 0.0) >= a.sasa]
            if len(surf) < 30:
                print(f"{t:<11} 표면 {len(surf)}잔기 — 건너뜀"); continue

            # ③ 표면 잔기마다 '결합자리와 같은 개수'의 이웃을 모아 패치
            E = len(true)
            ca = np.array([ag[i]["CA"].coord for i in surf], dtype=float)
            d = np.linalg.norm(ca[:, None, :] - ca[None, :, :], axis=2)
            order = np.argsort(d, axis=1)[:, :E]           # 자기 자신 포함 E개
            tset = np.array([1 if surf[j] in true else 0 for j in range(len(surf))])
            rec = tset[order].sum(axis=1) / E              # 위치마다 겹침 비율

            fr = [float((rec >= th).mean()) for th in ths]
            print(f"{t:<11}{len(surf):>6}{E:>8}{E/len(surf):>7.1%}" +
                  "".join(f"{f:>11.2%}" for f in fr))
            rows.append(dict(target=t, n_surface=len(surf), n_epitope=E,
                             epitope_frac=round(E / len(surf), 4),
                             mean_recall=round(float(rec.mean()), 4),
                             **{f"frac_ge_{th}": round(f, 5) for th, f in zip(ths, fr)}))
        except Exception as e:
            print(f"{t:<11} 실패: {type(e).__name__} {e}")

    if not rows:
        raise SystemExit("!! 계산된 복합체가 없다. --root 를 확인할 것.")
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

    S = np.array([r["n_surface"] for r in rows], float)
    E = np.array([r["n_epitope"] for r in rows], float)
    print(f"\n■ 전체 {len(rows)}종")
    print(f"  표면 잔기 중앙값 {np.median(S):.0f} · 결합자리 중앙값 {np.median(E):.0f} "
          f"({np.median(E/S):.1%})")
    for th, k in zip(ths, [f"frac_ge_{t}" for t in ths]):
        v = np.array([r[k] for r in rows], float)
        print(f"  겹침 {th:.0%} 이상인 위치의 비율 — 중앙값 {np.median(v):.2%} · "
              f"가장 관대한 복합체 {v.max():.2%} · 30종 평균 {v.mean():.2%}")
    if hypergeom is not None:
        s, e = np.median(S), np.median(E)
        p = hypergeom.sf(np.ceil(ths[0] * e) - 1, int(s), int(e), int(e))
        print(f"\n  (참고) 잔기를 이어붙이지 않고 무작위로 {int(e)}개 뽑았을 때 "
              f"{ths[0]:.0%}를 넘길 확률 = {p:.2%}")
        print("   위 실측이 이보다 크면, 패치가 이어져 있어 결합자리 근처를 통째로 덮는 위치가")
        print("   있기 때문이다. 기준의 근거로는 실측값을 쓸 것.")
    print(f"\n→ {a.out}")


if __name__ == "__main__":
    main()
