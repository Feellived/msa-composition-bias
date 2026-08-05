#!/usr/bin/env python3
"""[A] 같은 항원에 붙는 서로 다른 항체를 우리 파이프라인이 가려내는가 — 비닝 예행연습.

■ 왜 이게 필요한가 (2026-08-05)
  지금까지 모든 시험은 **복합체 하나당 항체 하나**였다. 그래서 한 번도 확인하지 않은 것이 있다.

      같은 항원에 다른 항체를 넣으면 우리 답이 달라지는가?

  안 달라진다면 우리는 항체-항원 결합을 예측하는 게 아니라 **항원 표면의 인기 있는 곳**을
  찍고 있는 것이다. 그건 DiscoTope 과 같은 물건이고(그래서 선택 단계에 쓰면 순환이라고
  적어 두었다), 파이프라인의 전제가 무너진다.

■ 새 예측이 필요 없다
  이미 끝낸 30종 안에 같은 항원군 항체가 여럿 있다 — RBD 8 · HA 8 · Env 3 → 쌍 59개.
  특히 8q7s_C · 8q7s_H · 8q7s_O 는 **같은 PDB 의 같은 항원**이라 정렬조차 필요 없다.

■ 무엇을 재나 — wet epitope binning 과 같은 판으로 만든다
  wet binning 은 구조를 보는 실험이 아니라 경쟁 결합 실험이다. 항체 A 를 먼저 붙이고 B 가
  더 붙는지 봐서 **쌍마다 차단/비차단**을 적는다. 그 판을 그대로 흉내 낸다.

    ① 정답 차단 행렬 : 결정 구조의 진짜 에피토프끼리 겹치면(자카드 ≥ --block) '차단'
    ② 예측 차단 행렬 : 우리가 고른 후보 자리끼리의 겹침
    ③ 판별력(AUC)   : 예측 겹침이 정답 차단쌍을 위로 올리는가. 0.5 = 무작위
    ④ 순열검정      : 항체 딱지를 무작위로 다시 붙여도(Mantel 식) 이만큼 나오나

■ 기준선을 반드시 같이 본다
  **항체 서열 동일도만으로 다 맞혀지면 이 파이프라인은 필요 없다.** 그래서 같이 잰다.
  우리 쪽 차별점은 "서열이 달라도 같은 자리에 붙는 항체를 묶는다"이므로,
  서열 기준선이 낮은데 우리가 높아야 의미가 있다.

■ 좌표 맞추기
  잔기 키는 (항원 사슬 순번, 그 타깃 참조서열에서의 0-based 위치)라 타깃끼리 그대로 비교할
  수 없다. 쌍마다 항원 서열을 국소 정렬해 **상대 타깃의 좌표계로 옮긴 뒤** 비교한다.
  정렬 동일도가 --min-id 미만인 쌍은 빼고, 몇 쌍을 뺐는지 끝에 보고한다(조용히 비지 않게).

■ 판정 기준은 plan/PLAN_202608_final8days.md §6 에 **결과를 보기 전에** 적어 두었다.
  성공 = AUC ≥ 0.65 이고 순열 p < 0.05 이고 서열 기준선보다 높을 것.

사용 (conda activate boltz · pipeline/ 에서):
  python -u cross_antibody.py --group RBD
  python -u cross_antibody.py --group all --out results/cross_antibody.csv
  python -u cross_antibody.py --targets "8q7s_C 8q7s_H 8q7s_O"     # 같은 PDB, 가장 깨끗
"""
import argparse
import csv
import glob
import json
import os
import random
import statistics as st
from collections import defaultdict

from Bio.Align import PairwiseAligner

import pose_features as PF

_aln = PairwiseAligner()
_aln.mode = "local"
_aln.match_score = 2
_aln.mismatch_score = -1
_aln.open_gap_score = -10
_aln.extend_gap_score = -0.5


# ── 서열 정렬로 좌표 옮기기 ──────────────────────────────────────────────────
def pos_map(seq_a, seq_b):
    """A 의 위치 → B 의 위치. 국소 정렬로 대응되는 구간만. (map, 동일도, 대응길이)."""
    if not seq_a or not seq_b:
        return {}, 0.0, 0
    try:
        a = _aln.align(seq_a, seq_b)[0]
    except Exception:
        return {}, 0.0, 0
    m, same, n = {}, 0, 0
    for (s0, e0), (s1, e1) in zip(a.aligned[0], a.aligned[1]):
        for k in range(e0 - s0):
            m[s0 + k] = s1 + k
            n += 1
            if seq_a[s0 + k] == seq_b[s1 + k]:
                same += 1
    return m, (same / n if n else 0.0), n


def chain_pairing(seqs_a, seqs_b):
    """A 의 항원 사슬마다 대응하는 B 의 사슬을 고른다(대응 잔기 많은 순 탐욕 배정).

    반환 {A사슬순번: (B사슬순번, 위치맵, 동일도)}. 항원이 한 사슬인 타깃이 대부분이라
    보통은 자명하지만, HA1+HA2 처럼 다사슬 항원이 있어 일반적으로 처리한다.
    """
    cand = []
    for i, sa in enumerate(seqs_a):
        for j, sb in enumerate(seqs_b):
            m, ident, n = pos_map(sa, sb)
            cand.append((ident * n, i, j, m, ident))
    cand.sort(key=lambda t: -t[0])
    out, ui, uj = {}, set(), set()
    for _, i, j, m, ident in cand:
        if i in ui or j in uj or not m:
            continue
        out[i] = (j, m, ident)
        ui.add(i); uj.add(j)
    return out


def project(res, pairing):
    """A 좌표의 잔기집합 → B 좌표. 옮기지 못한 잔기 수도 함께 반환."""
    out, drop = set(), 0
    for ci, p in res:
        pr = pairing.get(ci)
        if pr is None or p not in pr[1]:
            drop += 1
            continue
        out.add((pr[0], pr[1][p]))
    return out, drop


def jac(a, b):
    if not a or not b:
        return float("nan")
    return len(a & b) / len(a | b)


# ── 통계 ────────────────────────────────────────────────────────────────────
def auc(scores, labels):
    """Mann-Whitney 판별력. 양성(차단)이 음성보다 위에 올 확률. 0.5 = 무작위."""
    pos = [s for s, l in zip(scores, labels) if l]
    neg = [s for s, l in zip(scores, labels) if not l]
    if not pos or not neg:
        return float("nan")
    tot = 0.0
    for p in pos:
        for q in neg:
            tot += 1.0 if p > q else 0.5 if p == q else 0.0
    return tot / (len(pos) * len(neg))


def _ranks(v):
    order = sorted(range(len(v)), key=lambda i: v[i])
    r = [0.0] * len(v)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
            j += 1
        avg = (i + j) / 2.0 + 1
        for k in range(i, j + 1):
            r[order[k]] = avg
        i = j + 1
    return r


def spearman(a, b):
    """순위 상관. 예측 겹침이 정답 겹침을 따라가는가(이진화 없이)."""
    if len(a) < 3:
        return float("nan")
    ra, rb = _ranks(a), _ranks(b)
    ma, mb = st.mean(ra), st.mean(rb)
    num = sum((x - ma) * (y - mb) for x, y in zip(ra, rb))
    da = sum((x - ma) ** 2 for x in ra) ** 0.5
    db = sum((y - mb) ** 2 for y in rb) ** 0.5
    return num / (da * db) if da and db else float("nan")


# ── 타깃 읽기 ───────────────────────────────────────────────────────────────
def pick_site(d, how):
    """sites_<타깃>.json 에서 '이 항체가 붙는다고 우리가 말한 자리' 하나를 정한다."""
    cs = d.get("candidates") or []
    if not cs:
        return None
    if how == "union":
        u = set()
        for c in cs:
            u |= {tuple(r) for r in c["residues"]}
        return u
    if how == "first":
        c = cs[0]
    elif how == "oracle":                      # ⚠️ 정답 사용 — 천장 확인용
        def f1(c):
            p, r = c.get("precision") or 0.0, c.get("true_covered") or 0.0
            return 0.0 if p + r == 0 else 2 * p * r / (p + r)
        c = max(cs, key=f1)
    else:                                      # ncomp — 보고서의 선택기와 같은 규칙
        c = max(cs, key=lambda c: (c["n_comp"], -len(c["residues"])))
    return {tuple(r) for r in c["residues"]}


def load_target(t, sites_dir, targets_dir, cutoff, pick):
    sp = os.path.join(sites_dir, f"sites_{t}.json")
    cp = os.path.join(targets_dir, t, "chains.json")
    np_ = os.path.join(targets_dir, t, "native.cif")
    for p in (sp, cp, np_):
        if not os.path.exists(p):
            print(f"  ! {t}: {p} 없음 — 건너뜀")
            return None
    cj = json.load(open(cp))
    tr = PF.native_true(cj, np_, cutoff)
    if tr is None:
        print(f"  ! {t}: 정답 결합자리 계산 실패 — 건너뜀")
        return None
    pred = pick_site(json.load(open(sp)), pick)
    if not pred:
        print(f"  ! {t}: 후보가 없음 — 건너뜀")
        return None
    ag = [c["seq"] for c in cj["chains"] if c["role"] == "antigen"]
    ab = "".join(c["seq"] for c in cj["chains"] if c["role"] in ("heavy", "light"))
    return dict(target=t, ag_seqs=ag, ab_seq=ab,
                true={tuple(x) for x in tr[0]}, pred=pred)


# ── 한 항원군 분석 ──────────────────────────────────────────────────────────
def analyse(group, TS, block, min_id, nperm, seed, verbose=True):
    n = len(TS)
    pairs, skipped = [], 0
    for i in range(n):
        for j in range(i + 1, n):
            A, B = TS[i], TS[j]
            pairing = chain_pairing(A["ag_seqs"], B["ag_seqs"])
            if not pairing:
                skipped += 1
                continue
            ident = min(v[2] for v in pairing.values())
            if ident < min_id:
                skipped += 1
                continue
            tA, dt = project(A["true"], pairing)
            pA, dp = project(A["pred"], pairing)
            tj, pj = jac(tA, B["true"]), jac(pA, B["pred"])
            if tj != tj or pj != pj:
                skipped += 1
                continue
            _, sid, _ = pos_map(A["ab_seq"], B["ab_seq"])
            pairs.append(dict(group=group, a=A["target"], b=B["target"], i=i, j=j,
                              ag_ident=round(ident, 3), true_jac=round(tj, 4),
                              pred_jac=round(pj, 4), ab_seqid=round(sid, 4),
                              blocked=int(tj >= block),
                              drop_true=dt, drop_pred=dp))
    if len(pairs) < 3:
        print(f"  ! {group}: 쓸 수 있는 쌍이 {len(pairs)}개뿐 — 건너뜀 (뺀 쌍 {skipped})")
        return pairs, None

    tv = [p["true_jac"] for p in pairs]
    pv = [p["pred_jac"] for p in pairs]
    sv = [p["ab_seqid"] for p in pairs]
    lb = [p["blocked"] for p in pairs]
    nb = sum(lb)

    res = dict(group=group, n_target=n, n_pair=len(pairs), n_skip=skipped,
               n_blocked=nb, pred_sd=round(st.pstdev(pv), 4),
               pred_mean=round(st.mean(pv), 4), true_mean=round(st.mean(tv), 4),
               rho=round(spearman(pv, tv), 4))

    if 0 < nb < len(pairs):
        obs = auc(pv, lb)
        res["auc"] = round(obs, 4)
        res["auc_seqid"] = round(auc(sv, lb), 4)
        # 순열 — 항체 딱지를 무작위로 다시 붙인다(Mantel 식 노드 재배치)
        P = defaultdict(float)
        for p in pairs:
            P[(p["i"], p["j"])] = p["pred_jac"]
            P[(p["j"], p["i"])] = p["pred_jac"]
        rng, hit, done = random.Random(seed), 0, 0
        idx = list(range(n))
        for _ in range(nperm):
            rng.shuffle(idx)
            sc = [P.get((idx[p["i"]], idx[p["j"]])) for p in pairs]
            if any(s is None for s in sc):
                continue
            done += 1
            if auc(sc, lb) >= obs:
                hit += 1
        res["perm_p"] = round((hit + 1) / (done + 1), 5) if done else float("nan")
    else:
        res["auc"] = res["auc_seqid"] = res["perm_p"] = float("nan")
        if verbose:
            print(f"  ! {group}: 차단쌍이 {nb}/{len(pairs)} 이라 판별력을 못 잰다 "
                  f"(--block 을 조정하거나 항체를 더 넣을 것)")
    return pairs, res


def report(res, block):
    g = res["group"]
    print(f"\n■ {g}  항체 {res['n_target']}개 · 쌍 {res['n_pair']}개"
          + (f" · 정렬 불량으로 뺀 쌍 {res['n_skip']}개" if res["n_skip"] else ""))
    print(f"   정답 기준 차단쌍 {res['n_blocked']}/{res['n_pair']}  (자카드 ≥ {block})")
    print(f"   진짜 에피토프 겹침 평균 {res['true_mean']:.3f}"
          f"   ·  우리 예측 겹침 평균 {res['pred_mean']:.3f}")

    sd = res["pred_sd"]
    if sd <= 0.05:
        print(f"   ⚠️ 예측 겹침의 표준편차 {sd:.3f} — 거의 상수다. "
              f"**항체를 바꿔도 같은 자리를 내놓고 있다**는 뜻이고, 이 경우 아래 판별력은 의미가 없다.")
    else:
        print(f"   예측 겹침 표준편차 {sd:.3f}  → 항체에 따라 답이 갈리기는 한다")

    a, s, p = res["auc"], res["auc_seqid"], res["perm_p"]
    if a == a:
        print(f"   판별력 AUC  우리 {a:.3f}   |  항체 서열 동일도 기준선 {s:.3f}   |  무작위 0.500")
        print(f"   순열검정 p = {p:.4f}")
        ok = a >= 0.65 and p < 0.05 and a > s and sd > 0.05
        edge = 0.55 <= a < 0.65
        print("   판정: " + ("✅ 성공 — 정답 구조 없이 항체끼리를 가른다" if ok else
                            "△ 경계" if edge else
                            "✗ 실패 — 우리 답은 항체와 무관하다"))
        if a >= 0.65 and a <= s:
            print("   ⚠️ 다만 항체 서열만 봐도 같거나 더 잘 맞힌다 — 이 파이프라인이 보태는 게 없다.")
    print(f"   순위 상관(이진화 없이) rho = {res['rho']:.3f}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--group", default="all", help="RBD · HA · Env · all")
    ap.add_argument("--targets", default="", help="항원군 대신 타깃을 직접 지정")
    ap.add_argument("--sites", default="results")
    ap.add_argument("--targets-dir", default="targets")
    ap.add_argument("--manifest", default="targets_manifest.csv")
    ap.add_argument("--cutoff", type=float, default=5.0)
    ap.add_argument("--pick", default="ncomp", choices=["ncomp", "first", "union", "oracle"],
                    help="어느 후보를 '우리 답'으로 볼까. oracle 은 정답을 보므로 천장 확인용")
    ap.add_argument("--block", type=float, default=0.10,
                    help="진짜 에피토프끼리 이 이상 겹치면 '차단'으로 본다(사전 고정 0.10)")
    ap.add_argument("--min-id", type=float, default=0.5,
                    help="항원 서열 동일도가 이보다 낮은 쌍은 좌표를 못 믿으므로 뺀다")
    ap.add_argument("--nperm", type=int, default=20000)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", default="results/cross_antibody.csv")
    ap.add_argument("--summary-out", default="results/cross_antibody_summary.csv")
    a = ap.parse_args()

    have = sorted(os.path.basename(p)[6:-5] for p in glob.glob(os.path.join(a.sites, "sites_*.json")))
    if a.targets:
        groups = {"직접지정": [t for t in a.targets.replace(",", " ").split() if t in have]}
    else:
        man = {r["target"]: r for r in csv.DictReader(open(a.manifest))}
        g = defaultdict(list)
        for t in have:
            if t in man:
                g[man[t]["antigen_grp"]].append(t)
        groups = ({k: v for k, v in g.items() if len(v) >= 3} if a.group == "all"
                  else {a.group: g.get(a.group, [])})
    if not groups:
        raise SystemExit("!! 분석할 항원군이 없다 — --group 또는 --targets 확인")

    print(f"후보 선택 = {a.pick}"
          + ("   ⚠️ 정답을 본다(천장 확인용)" if a.pick == "oracle" else "   (정답을 보지 않는다)"))
    print(f"차단 판정 문턱 = 진짜 에피토프 자카드 ≥ {a.block}   (사전 고정)\n")

    allpairs, summ = [], []
    for gname, ts in sorted(groups.items()):
        print(f"── {gname}: {' '.join(ts)}")
        TS = [x for x in (load_target(t, a.sites, a.targets_dir, a.cutoff, a.pick) for t in ts) if x]
        if len(TS) < 3:
            print(f"  ! {gname}: 읽힌 타깃이 {len(TS)}종뿐 — 건너뜀")
            continue
        pairs, res = analyse(gname, TS, a.block, a.min_id, a.nperm, a.seed)
        allpairs += pairs
        if res:
            summ.append(res)
            report(res, a.block)

    # 항원군을 합친 전체 판 — 표본이 가장 크다
    if len(summ) > 1 and allpairs:
        pv = [p["pred_jac"] for p in allpairs]
        sv = [p["ab_seqid"] for p in allpairs]
        tv = [p["true_jac"] for p in allpairs]
        lb = [p["blocked"] for p in allpairs]
        if 0 < sum(lb) < len(lb):
            pooled = dict(group="전체합침", n_target=sum(s["n_target"] for s in summ),
                          n_pair=len(allpairs), n_skip=sum(s["n_skip"] for s in summ),
                          n_blocked=sum(lb), pred_sd=round(st.pstdev(pv), 4),
                          pred_mean=round(st.mean(pv), 4), true_mean=round(st.mean(tv), 4),
                          rho=round(spearman(pv, tv), 4),
                          auc=round(auc(pv, lb), 4), auc_seqid=round(auc(sv, lb), 4),
                          perm_p=float("nan"))
            summ.append(pooled)
            report(pooled, a.block)
            print("   ⚠️ 합친 판의 순열검정은 항원군마다 좌표계가 달라 따로 내지 않는다. "
                  "군별 p 를 볼 것.")

    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    if allpairs:
        with open(a.out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(allpairs[0])); w.writeheader(); w.writerows(allpairs)
        print(f"\n→ {a.out}  (쌍 {len(allpairs)}개)")
    if summ:
        keys = sorted({k for s in summ for k in s})
        with open(a.summary_out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=keys); w.writeheader(); w.writerows(summ)
        print(f"→ {a.summary_out}")


if __name__ == "__main__":
    main()
