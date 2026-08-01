#!/usr/bin/env python3
"""[후보 정제 sweep] 자리 후보를 좁히는 손잡이 넷을 훑어 F1 천장이 얼마나 오르는지 본다.

■ 왜 (2026-08-01 발견)
  지금 후보는 **정답보다 평균 2.3배 넓다** (정답 25잔기 · 최고 후보 58잔기).
  그래서 F1 천장이 0.535 에 갇혀 있고, 선택기를 아무리 잘 만들어도 그 위로 못 간다.
      선택기 개선 여지   0.366 → 0.405 → 0.535   +0.17
      후보 정제 여지     0.535 → 0.856(낙관 상한)  +0.32   ← 2배
  넓어지는 직접 원인은 site_reproducibility.py 의 후보 조립이
  **묶인 조성들의 합의 자리를 전부 union 하기 때문**이다 (`u |= cons[i][1]`).

■ 손잡이 넷
  --cons-frac   조성 안에서 몇 %의 실행에 나와야 그 조성의 합의 자리에 넣나 (현재 0.5)
  --merge-frac  묶인 조성 중 몇 %에 나와야 후보에 넣나. **1.0 이면 교집합, 0 이면 지금의 union**
  --max-res     후보를 잔기 수 상한으로 자른다(등장 빈도 높은 순)
  --main-patch  흩어진 잔기를 버리고 가장 큰 연결 덩어리만 남긴다  (구조 필요)

■ 무엇을 재나
  후보마다 덮음·정밀도·F1 을 계산하고, **타깃별 최고 F1 = 그 설정의 천장**을 낸다.
  ⚠️ 이건 천장이지 선택기 성능이 아니다. 실제로 그 후보를 고를 수 있는지는 별개다.

사용 (msa-depth pipeline 폴더에서, DockQ env 불필요·CPU):
  python refine_sites.py --targets "8k3k_D 8ulr_HL 8sis_HL"        # 기본 sweep
  python refine_sites.py --all --out results/refine_sweep.csv
  python refine_sites.py --all --dump-sites results/sites_refined --use 0.7,0.5,30
      → 그 설정으로 sites_*.json 재생성. 이어서
        python eval_selectors.py --sites results/sites_refined --abepi <abepiscore_all.csv>
      로 "정제하면 선택기가 실제로 좋아지나"를 잰다(지금까지의 F1 은 천장일 뿐이다).
"""
import argparse
import csv
import glob
import json
import os
import sys
from collections import Counter, defaultdict

import epitope_cluster as EC            # pred_epitope · jac · consensus
import pose_features as PF              # native_true
from site_reproducibility import split_run, link_clusters


def load_runs(tgt, csv_path, data, targets_dir, cutoff):
    """실행별 대표 자세(DockQ 최고)의 예측 결합자리를 조성별로 모은다."""
    rows = list(csv.DictReader(open(csv_path)))
    if not rows:
        return None
    # ⚠️ results/ 에는 compreps_summary.csv 처럼 형식이 다른 파일도 있다. 컬럼으로 걸러낸다.
    need = {"model", "seed", "pose", "dockq"}
    miss = need - set(rows[0])
    if miss:
        print(f"  ! {tgt}: 자세 단위 CSV 가 아님(없는 열 {sorted(miss)}) — 건너뜀")
        return None
    model, depth = rows[0]["model"], rows[0].get("depth", "")
    cj = json.load(open(os.path.join(targets_dir, tgt, "chains.json")))
    tr = PF.native_true(cj, os.path.join(targets_dir, tgt, "native.cif"), cutoff)
    if tr is None:
        return None
    true = set(tr[0])

    best = {}
    for r in rows:
        try:
            q = float(r["dockq"])
        except Exception:
            continue
        s = r["seed"]
        if s not in best or q > best[s][0]:
            best[s] = (q, r["pose"])

    base = os.path.join(data, "seedrep_cand", model, tgt, depth)
    groups = defaultdict(list)
    for s, (q, pose) in sorted(best.items()):
        hits = glob.glob(os.path.join(base, s, "results", "**", pose), recursive=True)
        if not hits:
            continue
        ep, _ = EC.pred_epitope(cj, hits[0], cutoff)
        if ep:
            groups[split_run(s)[0]].append(ep)
    return dict(target=tgt, model=model, depth=depth, true=true, groups=dict(groups))


def build(groups, cons_frac, merge_frac, max_res, link):
    """조성별 합의 → 묶기 → 후보 잔기 집합. merge_frac 이 union 을 대체하는 핵심."""
    cons = []
    for k, v in sorted(groups.items()):
        c = EC.consensus(v, cons_frac)
        if c:
            cons.append((k, c))
    if not cons:
        return [], 0
    out = []
    for ix in link_clusters(cons, link):
        cnt = Counter()
        for i in ix:
            cnt.update(cons[i][1])
        need = max(1, int(round(merge_frac * len(ix))))
        res = {r for r, n in cnt.items() if n >= need}
        if not res:                                  # 너무 빡세면 최빈 잔기라도 남긴다
            top = max(cnt.values())
            res = {r for r, n in cnt.items() if n == top}
        if max_res and len(res) > max_res:           # 등장 빈도 높은 순으로 자른다
            res = set(sorted(res, key=lambda r: -cnt[r])[:max_res])
        out.append(dict(comps=[cons[i][0] for i in ix], res=res))
    out.sort(key=lambda c: -len(c["comps"]))
    return out, len(cons)


def score(cands, true, n_cons=0, n_group=0):
    """⚠️ 천장(F1)만 보면 안 된다. 후보가 1개로 줄면 고를 것이 없어져 선택기가 죽는다."""
    best = dict(f1=0.0, rec=0.0, pre=0.0, n=0,
                n_cand=len(cands), n_cons=n_cons, n_group=n_group)
    for c in cands:
        u = c["res"]
        rec = len(u & true) / len(true) if true else 0.0
        pre = len(u & true) / len(u) if u else 0.0
        f1 = 0.0 if pre + rec == 0 else 2 * pre * rec / (pre + rec)
        if f1 > best["f1"]:
            best.update(f1=f1, rec=rec, pre=pre, n=len(u))
    return best


def loo(per, base_key):
    """한 타깃을 빼고 나머지에서 최고 조합을 고른 뒤, 뺀 타깃에서 채점한다.

    ⚠️ 36조합을 같은 30종에서 고르면 +0.104 는 낙관적이다. 여기서 나오는 값이 정직한 상승폭.
    계산이 아니라 위에서 이미 만든 점수판을 다시 세는 것이라 비용이 없다.
    """
    tab = defaultdict(dict)                       # 조합 → {타깃: (f1, n_cand)}
    for r in per:
        tab[(r["cons_frac"], r["merge_frac"], r["max_res"])][r["target"]] = (r["f1"], r["n_cand"])
    tgts = sorted({r["target"] for r in per})
    held, chosen = [], Counter()
    for t in tgts:
        pick, bestv = None, -1.0
        for k, d in tab.items():
            tr = [v for tt, v in d.items() if tt != t]
            if sum(c for _, c in tr) / len(tr) < 1.8:      # 후보가 말라 선택기가 죽는 설정은 배제
                continue
            m = sum(f for f, _ in tr) / len(tr)
            if m > bestv:
                pick, bestv = k, m
        if pick is None:
            continue
        held.append(tab[pick][t][0])
        chosen[pick] += 1
    if not held:
        return
    cur = [tab[base_key][t][0] for t in tgts] if base_key in tab else []
    print(f"\n■ 한 타깃 빼기(LOO) — 조합을 나머지 {len(tgts)-1}종에서 고르고 뺀 1종에서 채점")
    print(f"  정직한 F1 천장 {sum(held)/len(held):.3f}", end="")
    if cur:
        print(f"   (현재 설정 {sum(cur)/len(cur):.3f} 대비 {sum(held)/len(held)-sum(cur)/len(cur):+.3f})")
    else:
        print()
    top = chosen.most_common(3)
    print("  고른 조합 " + " · ".join(f"cons{k[0]} merge{k[1]} max{k[2] or '-'} ×{n}" for k, n in top))
    if len(chosen) == 1:
        print(f"  → {len(held)}번 다 같은 조합을 골랐다. 표면이 고원이라는 뜻이고 전이 위험이 낮다.")
    else:
        print(f"  ⚠️ 조합이 {len(chosen)}가지로 갈렸다 — 최적점이 타깃에 따라 흔들린다는 신호.")


def dump_sites(loaded, cf, mf, mr, link, outdir):
    """정제한 후보를 site_reproducibility.py --dump-sites 와 **같은 형식**으로 쓴다.

    이걸 eval_selectors.py 에 넘겨야 "정제하면 선택기가 실제로 좋아지나"를 잴 수 있다.
    지금까지 잰 F1 은 천장(후보 중 최고)이지 고를 수 있는 값이 아니다.
    """
    os.makedirs(outdir, exist_ok=True)
    n_ok = 0
    for d in loaded:
        cands, _ = build(d["groups"], cf, mf, mr, link)
        true, sites = d["true"], []
        for ci, c in enumerate(cands, 1):
            u = c["res"]
            rec = len(u & true) / len(true) if true else float("nan")
            pre = len(u & true) / len(u) if u else float("nan")
            sites.append(dict(
                cand=ci, n_comp=len(c["comps"]), comps=list(c["comps"]),
                # ⚠️ 잔기 키는 posmap 에서 온 numpy 정수 — 캐스팅 없으면 json 이 죽는다(2026-07-29).
                residues=sorted([[int(x) for x in k] for k in u]),
                # ⚠️ 아래 둘은 정답 구조를 본 값이다. 보고용이고 '고르는 데' 쓰면 안 된다.
                true_covered=round(rec, 4), precision=round(pre, 4),
                from_full_msa=any(str(x).startswith("seedfull") for x in c["comps"])))
        if len(sites) >= 1:
            n_ok += 1
        json.dump(dict(target=d["target"], model=d["model"], depth=d["depth"],
                       n_true_res=len(true),
                       refined=dict(cons_frac=cf, merge_frac=mf, max_res=mr, link=link),
                       candidates=sites),
                  open(os.path.join(outdir, f"sites_{d['target']}.json"), "w"), indent=1,
                  default=lambda o: int(o) if hasattr(o, "__int__") else float(o))
    print(f"\n→ {outdir}/sites_*.json  ({n_ok}종 · cons={cf} merge={mf} max_res={mr or '-'})")
    print("  다음: python eval_selectors.py --sites "
          f"{outdir} --abepi <abepiscore_all.csv>  ← 정제 전 결과와 나란히 놓고 본다")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--targets", default="")
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--dir", default="results")
    ap.add_argument("--data", default=os.environ.get("DATA", "/mnt/data/admuser/msadepth") + "/compreps")
    ap.add_argument("--targets-dir", default="targets")
    ap.add_argument("--cutoff", type=float, default=5.0)
    ap.add_argument("--link", type=float, default=0.5)
    ap.add_argument("--cons-fracs", default="0.5 0.7 0.8")
    ap.add_argument("--merge-fracs", default="0 0.5 0.75 1.0")
    ap.add_argument("--max-res-list", default="0 40 30")
    ap.add_argument("--out", default="results/refine_sweep.csv")
    ap.add_argument("--per-out", default="results/refine_per_target.csv",
                    help="조합×타깃 점수판. 이게 있으면 한-타깃-빼기(LOO)가 계산 없이 된다")
    ap.add_argument("--dump-sites", default="",
                    help="정제한 후보를 sites_<타깃>.json 으로 이 폴더에 쓴다 "
                         "(site_reproducibility.py --dump-sites 와 같은 형식 → eval_selectors.py 로 바로 넘어간다)")
    ap.add_argument("--use", default="",
                    help="떨굴 설정 'cons,merge,maxres' (예: 0.7,0.5,30). 생략하면 후보고갈을 통과한 최고점")
    a = ap.parse_args()

    tg = a.targets.replace(",", " ").split()
    if a.all or not tg:
        SKIP = {"summary", "all"}
        tg = sorted(t for t in (os.path.basename(p)[9:-4]
                                for p in glob.glob(os.path.join(a.dir, "compreps_*.csv")))
                    if t not in SKIP)
    loaded = []
    for t in tg:
        p = os.path.join(a.dir, f"compreps_{t}.csv")
        if not os.path.exists(p):
            print(f"  ! {t}: {p} 없음 — 건너뜀"); continue
        d = load_runs(t, p, a.data, a.targets_dir, a.cutoff)
        if d and d["groups"]:
            loaded.append(d)
            print(f"  {t:<11} 조성 {len(d['groups']):>2}가지 · 정답 {len(d['true']):>3}잔기")
        else:
            print(f"  ! {t}: 실행/자세를 못 읽음 — 건너뜀")
    if not loaded:
        sys.exit("!! 읽은 타깃이 없다")

    CF = [float(x) for x in a.cons_fracs.split()]
    MF = [float(x) for x in a.merge_fracs.split()]
    MR = [int(x) for x in a.max_res_list.split()]
    rows, per, base_key = [], [], (0.5, 0.0, 0)
    print(f"\n{'cons':>5}{'merge':>7}{'maxres':>8}{'F1천장':>9}{'덮음':>8}{'정밀도':>8}"
          f"{'후보크기':>9}{'후보수':>7}{'살아남은조성':>13}")
    print("-" * 78)
    for cf in CF:
        for mf in MF:
            for mr in MR:
                sc = []
                for d in loaded:
                    cands, ncons = build(d["groups"], cf, mf, mr, a.link)
                    s = score(cands, d["true"], ncons, len(d["groups"]))
                    sc.append(s)
                    per.append(dict(cons_frac=cf, merge_frac=mf, max_res=mr,
                                    target=d["target"], model=d["model"], depth=d["depth"],
                                    f1=round(s["f1"], 4), recall=round(s["rec"], 4),
                                    precision=round(s["pre"], 4), n_res=s["n"],
                                    n_cand=s["n_cand"], n_cons=s["n_cons"],
                                    n_group=s["n_group"], n_true=len(d["true"])))
                m = lambda k: sum(s[k] for s in sc) / len(sc)
                tag = "  ← 현재" if (cf, mf, mr) == base_key else ""
                warn = "  ⚠️ 후보 부족" if m("n_cand") < 1.8 else ""
                print(f"{cf:>5.1f}{mf:>7.2f}{mr if mr else '-':>8}"
                      f"{m('f1'):>9.3f}{m('rec'):>8.3f}{m('pre'):>8.3f}{m('n'):>9.1f}"
                      f"{m('n_cand'):>7.1f}{m('n_cons'):>7.1f}/{m('n_group'):<5.1f}{tag}{warn}")
                rows.append(dict(cons_frac=cf, merge_frac=mf, max_res=mr,
                                 f1=round(m("f1"), 4), recall=round(m("rec"), 4),
                                 precision=round(m("pre"), 4), n_res=round(m("n"), 1),
                                 n_cand=round(m("n_cand"), 2), n_cons=round(m("n_cons"), 2),
                                 n_group=round(m("n_group"), 2), n_target=len(loaded)))
    rows.sort(key=lambda r: -r["f1"])
    ok = [r for r in rows if r["n_cand"] >= 1.8]      # 후보가 평균 2개 미만이면 고를 것이 없다
    b = ok[0] if ok else rows[0]
    if ok and ok[0] is not rows[0]:
        t = rows[0]
        print(f"\n⚠️ F1 최고는 cons={t['cons_frac']} merge={t['merge_frac']} "
              f"max_res={t['max_res'] or '-'} ({t['f1']:.3f}) 인데 후보가 평균 {t['n_cand']:.1f}개뿐이라 제외했다.")
    print(f"\n■ 최고 설정  cons={b['cons_frac']} merge={b['merge_frac']} "
          f"max_res={b['max_res'] or '-'}  →  F1 천장 {b['f1']:.3f}")
    cur = next(r for r in rows if (r["cons_frac"], r["merge_frac"], r["max_res"]) == base_key)
    print(f"  현재 설정 대비 {b['f1'] - cur['f1']:+.3f}  ({cur['f1']:.3f} → {b['f1']:.3f})")
    print("\n⚠️ 이건 천장이다. 실제로 그 후보를 고를 수 있는지는 eval_selectors.py 가 따로 답한다.")

    loo(per, base_key)

    if a.dump_sites:
        if a.use:
            cf, mf, mr = [float(x) for x in a.use.replace(",", " ").split()]
            mr = int(mr)
        else:
            cf, mf, mr = b["cons_frac"], b["merge_frac"], b["max_res"]
        dump_sites(loaded, cf, mf, mr, a.link, a.dump_sites)

    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    print(f"→ {a.out}")
    if a.per_out:
        os.makedirs(os.path.dirname(a.per_out) or ".", exist_ok=True)
        with open(a.per_out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(per[0])); w.writeheader(); w.writerows(per)
        print(f"→ {a.per_out}  (조합×타깃 점수판 — 재실행 없이 LOO·부트스트랩 가능)")


if __name__ == "__main__":
    main()
