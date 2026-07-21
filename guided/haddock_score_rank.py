#!/usr/bin/env python3
"""HADDOCK '생성 O · 선택 X' 정량화: raw pool best-DockQ pose가 HADDOCK 자기 점수로 몇 등인가.
   ── co-folder의 'ipTM은 best pose 못 고른다'에 대응하는, 물리 모델(HADDOCK) 버전 진단.

설계(적대검증 반영 — 함정 회피):
  · 스테이지 하나 고정(--stage 4_emref). score와 DockQ를 '같은 emref survivor pool'에서 계산
    → DockQ 히트맵 pool ≠ score pool / 스테이지 혼합(가중치 다름) / 파일명index=rank 오해 원천 차단.
  · score는 clustfcc TSV의 `score` 컬럼으로 직접 정렬(더 음수=좋음). 파일명 정수 안 씀.
  · rank 두 종을 별도 컬럼으로: (주지표) per-model 전역 score rank / (HADDOCK 실보고) cluster rank.
  · 'HADDOCK이 실제 내놓는 값' 3렌즈로 명시 분리:
      (A) clust_pick = rank-1 클러스터의 best-SCORE pose DockQ  ← HADDOCK 기본(cluster) 모드 실배포값
      (B) score1     = 전역 per-model best-score pose DockQ      ← 무클러스터(single-model) 렌즈
      (C) topclust_oracle = rank-1 클러스터 멤버 중 MAX DockQ     ← 클러스터 내 상한(과대평가, 참고용)
    헤드라인 regret = oracle_best − clust_pick(A).  best-DockQ pose가 unclustered면 = 클러스터 배제=미배포.
  · rank는 백분위(k/N) + hit@1/5/10 + regret(ΔDockQ)로 pool-size 정규화.
  · DockQ = merged-chain(A=항원, B=항체 H+L 병합), oracle=MAX 규약(Phase-0 min과 혼용 금지).

⚠️ 범위·정직성:
  · 이 pool은 emref 'survivor'(seletop 점수컷 통과분, ~200). seletop서 '검열된' near-native는 여기 없음
    → 그건 별도 질문(1_rigidbody oracle = haddock_pool_oracle.py). 여기 best는 그 상한이 아님.
  · Spearman ρ는 survivor(범위축소) 부분집합 조건부 → |ρ|가 range restriction으로 과소평가될 수 있음.
  · HADDOCK score는 Eair(제약 위반 페널티) 포함. 이 런은 ab-initio(CDR/표면 제약, 참값-계면 아님)
    → 점수가 정답 위치를 직접 인코딩하는 순환 leak 없음. 그래도 Eair 포함임을 명시.

사용(boltz env, repo 루트 consensus_docking, DockQ 설치):
  python scripts/haddock_score_rank.py --targets "8P5M 8SDF 8SIQ 8SIS 8SIT 8XSI 9ML8 9ML9 9SBB 9ZDU"
  (resumable: 포즈별 DockQ를 results/haddock_score_rank_cache.csv 에 캐시)"""
import argparse, csv, glob, os, shutil, sys, tempfile
from collections import defaultdict
import build_msafree_summary as B
import haddock_dockq as HD

ACC, MED = 0.23, 0.49   # CAPRI Acceptable / Medium


def qualifying_tsvs(rundir):
    """clustfcc TSV(model_name/score/cluster_id) 후보 전부(중복 제거, clustfcc 우선)."""
    cands = glob.glob(os.path.join(rundir, "*clustfcc*", "*.tsv")) + \
            glob.glob(os.path.join(rundir, "**", "*.tsv"), recursive=True)
    seen, out = set(), []
    for f in cands:
        if f in seen: continue
        seen.add(f)
        try:
            with open(f) as fh:
                hdr = next(csv.reader(fh, delimiter="\t"))
            if {"model_name", "score", "cluster_id"} <= set(hdr): out.append(f)
        except Exception: pass
    return out


def _is_float(s):
    try: float(s); return True
    except (TypeError, ValueError): return False


def spearman(x, y):
    """scipy 없이 Spearman ρ. score↑=나쁨, dockq↑=좋음 → 이상적 선택기면 ρ<0(음), 0 근처면 선택 고장."""
    n = len(x)
    if n < 3: return float("nan")
    def ranks(v):
        order = sorted(range(n), key=lambda i: v[i]); r = [0.0] * n; i = 0
        while i < n:
            j = i
            while j + 1 < n and v[order[j + 1]] == v[order[i]]: j += 1
            avg = (i + j) / 2.0 + 1
            for k in range(i, j + 1): r[order[k]] = avg
            i = j + 1
        return r
    rx, ry = ranks(x), ranks(y)
    mx, my = sum(rx) / n, sum(ry) / n
    num = sum((rx[i] - mx) * (ry[i] - my) for i in range(n))
    dx = sum((rx[i] - mx) ** 2 for i in range(n)) ** 0.5
    dy = sum((ry[i] - my) ** 2 for i in range(n)) ** 0.5
    return num / (dx * dy) if dx > 0 and dy > 0 else float("nan")


def load_cache(path):
    c = defaultdict(dict)   # (target) -> {model_name: dockq_or_None}
    if not os.path.exists(path): return c
    with open(path) as fh:
        for r in csv.DictReader(fh):
            v = (r.get("dockq") or "").strip()          # 손상/절단 라인에 강건(resume 중단 안 함)
            c[r["target"]][r["model_name"]] = float(v) if _is_float(v) else None
    return c


def dockq_pose(pose_path, ni, natm, td):
    """merged-chain DockQ (float) or None. 역할판별은 HD.id_and_recall 재사용."""
    try:
        pm = HD.load_gz(pose_path)
        _, pid, abids = HD.id_and_recall(pm, ni)
        if not abids: return None
        mp = os.path.join(td, "m.pdb"); HD.write_merged(pm, pid, abids, mp)
        return HD.dockq(mp, natm)
    except Exception:
        return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--targets", default=" ".join(B.TARGETS))
    ap.add_argument("--stage", default="4_emref", help="점수·DockQ를 계산할 단일 스테이지(고정)")
    ap.add_argument("--haddock", default="haddock")
    ap.add_argument("--out", default="results/haddock_score_rank.csv")
    ap.add_argument("--cache", default="results/haddock_score_rank_cache.csv")
    ap.add_argument("--limit", type=int, default=0, help=">0이면 타깃당 포즈 수 제한(스모크용, best-score 우선)")
    a = ap.parse_args()
    if shutil.which("DockQ") is None:   # silent 전부-None 방지: DockQ 없으면 즉시 중단
        sys.exit("✗ DockQ가 PATH에 없음 — `conda activate boltz` 후 `which DockQ` 확인하고 실행하세요.")
    os.makedirs("results", exist_ok=True)
    cache = load_cache(a.cache)
    cf = open(a.cache, "a", newline=""); cw = csv.writer(cf)
    if os.stat(a.cache).st_size == 0: cw.writerow(["target", "model_name", "cluster_id", "score", "dockq"])

    hdr = ["target", "class", "stage", "n_tsv", "n_scored", "coverage",
           "best_dockq", "best_pose", "best_pose_score", "best_pose_unclustered",
           "model_rank", "model_pct", "cluster_rank", "n_clusters",
           "clustpick_pose", "clustpick_dockq", "regret_clustmode",
           "score1_pose", "score1_dockq", "regret_permodel", "topclust_oracle_dockq",
           "spearman_score_dockq", "hit1_023", "hit5_023", "hit10_023", "hit1_049", "hit5_049", "hit10_049"]
    rows_out = []
    print(f"[stage={a.stage}]  '생성 O·선택 X' 진단 — best-DockQ pose가 HADDOCK score로 몇 등인가")
    print("  HADpick = rank-1 클러스터의 best-SCORE pose DockQ (HADDOCK cluster 모드 실배포값)")
    print(f"{'tgt':6}{'class':10}{'N/tsv':>8}{'best':>7}{'model_rk':>10}{'clust':>7}{'HADpick':>9}{'regret':>8}{'ρ':>7}  hit@1/5/10(.49)")
    print("-" * 100)

    for t in a.targets.split():
        t = t.upper(); ni = B.native_info(t)
        if ni is None: print(f"{t:6} native 없음"); continue
        rundir = os.path.join(a.haddock, t.lower(), "run")
        tsvs = qualifying_tsvs(rundir)
        if not tsvs: print(f"{t:6} clustfcc TSV 없음(HADDOCK 미완료?) skip"); continue
        if len(tsvs) > 1: print(f"   ⚠️ {t}: clustfcc TSV {len(tsvs)}개 발견 → 첫 번째 사용({os.path.relpath(tsvs[0], rundir)}); 스테이지 확인 권장")
        tsv = tsvs[0]
        emdir = os.path.join(rundir, a.stage)
        recs = [(r["model_name"], float(r["score"]), (r.get("cluster_id") or "-").strip())
                for r in csv.DictReader(open(tsv), delimiter="\t")]
        if a.limit > 0: recs = sorted(recs, key=lambda x: x[1])[:a.limit]
        if not recs: print(f"{t:6} TSV 행 없음 skip"); continue
        # 가드: TSV model_name이 emdir 실제 파일로 해석되나(score-pool↔pose-pool 스테이지 일치)
        nres = sum(1 for mn, _, _ in recs
                   if os.path.exists(os.path.join(emdir, mn + ".gz")) or os.path.exists(os.path.join(emdir, mn)))
        if nres == 0:
            print(f"{t:6} ✗ TSV model_name이 {a.stage}/ 파일로 하나도 안 풀림 → 스테이지 mismatch skip"); continue
        if nres < 0.5 * len(recs):
            print(f"   ⚠️ {t}: TSV {len(recs)}행 중 {nres}개만 {a.stage}/ 파일 해석 → 부분 mismatch 의심")

        with tempfile.TemporaryDirectory() as td:
            natm = HD.native_merged(t, td)
            if natm is None: print(f"{t:6} native merge 실패 skip"); continue
            data = []   # (model_name, score, cid, dockq)
            done = cache.get(t, {})
            for i, (mn, sc, cid) in enumerate(recs):
                if mn in done and done[mn] is not None:   # None-캐시(과거 실패)는 자동 재계산
                    q = done[mn]
                else:
                    pose = os.path.join(emdir, mn + ".gz")
                    if not os.path.exists(pose): pose = os.path.join(emdir, mn)
                    q = dockq_pose(pose, ni, natm, td) if os.path.exists(pose) else None
                    cw.writerow([t, mn, cid, f"{sc:.3f}", ("" if q is None else f"{q:.4f}")]); cf.flush()
                data.append((mn, sc, cid, q))
                if (i + 1) % 50 == 0: print(f"   {t}: {i+1}/{len(recs)} DockQ …")

        scored = [(mn, sc, cid, q) for (mn, sc, cid, q) in data if q is not None]
        M = len(data); N = len(scored)
        if not scored:
            print(f"{t:6} DockQ 전부 실패 skip (N/tsv=0/{M})"); continue
        cov = N / M
        if cov < 0.8:
            print(f"   ⚠️ {t}: DockQ 커버리지 {N}/{M}={cov:.0%} (낮음 → rank%/hit@k가 부분 pool 기준임 유의)")
        by_score = sorted(scored, key=lambda x: x[1])          # score 오름차순(1등=최고점)
        srank = {mn: i + 1 for i, (mn, _, _, _) in enumerate(by_score)}
        # oracle best = MAX DockQ; 동점이면 score 최저(순위 앞) pose로 = HADDOCK에 benefit-of-doubt
        maxq = max(q for _, _, _, q in scored)
        best_mn, best_sc, best_cid, best_q = min(((mn, sc, cid, q) for (mn, sc, cid, q) in scored if q == maxq),
                                                 key=lambda x: x[1])
        mrank = srank[best_mn]; mpct = mrank / N
        best_unclust = best_cid in ("", "-")
        s1_mn, s1_sc, s1_cid, s1_q = by_score[0]                # 전역 per-model top1(무클러스터 렌즈)

        # cluster rank: 각 클러스터 top-4 score 평균 오름차순(HADDOCK 규약)
        cl = defaultdict(list)
        for mn, sc, cid, q in scored:
            if cid and cid != "-": cl[cid].append((sc, q, mn))
        clkey = {cid: sum(s for s, _, _ in sorted(v)[:4]) / min(4, len(v)) for cid, v in cl.items()}
        clorder = sorted(clkey, key=lambda c: clkey[c])
        clrank = {c: i + 1 for i, c in enumerate(clorder)}
        best_clrank = clrank.get(best_cid)                     # best-DockQ pose 소속 클러스터 순위
        # (A) HADDOCK cluster 모드 실배포값 = rank-1 클러스터의 best-SCORE pose DockQ
        clust_pick_q = clust_pick_mn = None; topclust_oracle_q = None
        if clorder:
            mem = sorted(cl[clorder[0]])                       # (sc,q,mn) score 오름차순
            clust_pick_mn, _sc0, clust_pick_q = mem[0][2], mem[0][0], mem[0][1]
            topclust_oracle_q = max(q for _, q, _ in mem)      # (C) 클러스터 내 상한

        def hit(k, thr): return int(any(q >= thr for _, _, _, q in by_score[:k]))
        rho = spearman([sc for _, sc, _, _ in scored], [q for _, _, _, q in scored])

        clrank_str = (f"{best_clrank}/{len(clorder)}" if best_clrank else ("미배포" if best_unclust else "-"))
        hp = clust_pick_q if clust_pick_q is not None else float("nan")
        reg_clust = (best_q - clust_pick_q) if clust_pick_q is not None else float("nan")
        print(f"{t:6}{B.EPICLASS.get(t,''):10}{f'{N}/{M}':>8}{best_q:>7.2f}"
              f"{f'{mrank}/{N}':>10}{clrank_str:>7}{hp:>9.2f}{reg_clust:>8.2f}{rho:>7.2f}  "
              f"{hit(1,MED)}/{hit(5,MED)}/{hit(10,MED)}")

        def f(v, nd=4): return "" if v is None else f"{v:.{nd}f}"
        rows_out.append({
            "target": t, "class": B.EPICLASS.get(t, ""), "stage": a.stage, "n_tsv": M, "n_scored": N,
            "coverage": f"{cov:.3f}",
            "best_dockq": f"{best_q:.4f}", "best_pose": best_mn, "best_pose_score": f"{best_sc:.2f}",
            "best_pose_unclustered": int(best_unclust),
            "model_rank": mrank, "model_pct": f"{mpct:.3f}",
            "cluster_rank": (best_clrank if best_clrank else ("unclustered" if best_unclust else "")),
            "n_clusters": len(clorder),
            "clustpick_pose": (clust_pick_mn or ""), "clustpick_dockq": f(clust_pick_q),
            "regret_clustmode": f(reg_clust if clust_pick_q is not None else None),
            "score1_pose": s1_mn, "score1_dockq": f"{s1_q:.4f}", "regret_permodel": f"{best_q - s1_q:.4f}",
            "topclust_oracle_dockq": f(topclust_oracle_q),
            "spearman_score_dockq": f"{rho:.4f}",
            "hit1_023": hit(1, ACC), "hit5_023": hit(5, ACC), "hit10_023": hit(10, ACC),
            "hit1_049": hit(1, MED), "hit5_049": hit(5, MED), "hit10_049": hit(10, MED),
        })
    cf.close()
    with open(a.out, "w", newline="") as fo:
        w = csv.DictWriter(fo, fieldnames=hdr); w.writeheader(); w.writerows(rows_out)

    # 집계
    if rows_out:
        n = len(rows_out)
        orc = sum(1 for r in rows_out if float(r["best_dockq"]) >= MED)
        cpick = sum(1 for r in rows_out if r["clustpick_dockq"] and float(r["clustpick_dockq"]) >= MED)
        s1 = sum(1 for r in rows_out if float(r["score1_dockq"]) >= MED)
        h5 = sum(r["hit5_049"] for r in rows_out)
        print("-" * 100)
        print(f"집계(N타깃={n}, Medium 0.49): oracle(pool내 best)={orc}/{n}  |  "
              f"HADDOCK cluster-pick(실배포)={cpick}/{n}  per-model score-top1={s1}/{n}  score-top5 hit={h5}/{n}")
        print(f"→ {a.out} , 포즈별 캐시 → {a.cache}")
    print("\n[해석] oracle ≫ cluster-pick(=HADDOCK 실배포) → 'near-native는 만들되(생성 O) 점수로 못 고름(선택 X)'.")
    print("       model_rk 뒤(예 #40/200)+regret 큼+ρ≈0(음이 아님) → 물리 채점이 pose 품질과 무관 → ipTM처럼 고장.")
    print("       → consensus/재랭커가 co-folder뿐 아니라 HADDOCK pose에도 필요(단일 물리점수로 안 됨).")
    print("       주의: best는 emref survivor 상한(seletop 검열분 제외). ρ는 survivor 조건부(range restriction).")
    print("             score는 Eair(ab-initio 제약) 포함. topclust_oracle_dockq=클러스터 내 상한(참고, 과대평가).")


if __name__ == "__main__": main()
