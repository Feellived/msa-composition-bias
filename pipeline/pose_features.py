#!/usr/bin/env python3
"""[Phase 0] pose 단위 라벨·피처 추출 — 재랭커(선택기) 학습/평가의 전제 데이터.

════════ 왜 이 스크립트가 필요한가 (나중에 이것만 봐도 알도록) ════════
기존 산출물은 pose를 못 가린다:
  · dockq_sweep.csv   = (target, model, rung) 한 행에 **5 pose 중 최고값**(best_dockq, n_pose=5)
  · epitope_shift.csv = 같은 단위의 mean_/oracle_ **집계값**
그런데 재랭커가 푸는 문제는 "같은 복합체의 후보 pose 중 어느 게 나은가"이므로
**pose 하나하나의 라벨**이 있어야 순위를 배우고/평가할 수 있다. → 이 스크립트가 그걸 만든다.

라벨 2종(둘 다 저장):
  · dockq  = 구조 품질(엄격, 표준 성공기준 0.23/0.49/0.80)
  · recall = 에피토프 위치 정확도(관대) — 편향 서사와 직결되고, n이 작을 때 더 배우기 쉬울 수 있음

피처(pose마다 달라지는 것만):
  ⚠️ 같은 복합체의 pose들은 **서열이 전부 같다** → ESM-2 같은 서열 전용 임베딩은 pose를 구분 못 한다.
     pose를 가르는 건 기하(어느 잔기가 접촉하나·거리·중심) + 모델 자기신뢰도다. 그래서 여기선 기하·신뢰도만 뽑는다.
  n_contact · overrep · true_rank · pop_rank · dcc_true · dcc_pop  (+ 있으면 iptm/ptm/plddt)

⭐ 부수 목적 = 게이트 확인: co-folder가 남긴 confidence(JSON)가 보존됐는지 같이 집계해서 출력한다.
   보존됐으면 ipTM 대비 실험(Phase 1)과 PAE 계열 피처가 가능하고, 없으면 기하 피처만으로 가야 한다.

경로 규약은 dockq_sweep.py / epitope_shift.py와 동일:
  pose = $DATA/<model>/<t>/rung<r>/results/**/*.cif · native·chains.json = targets/<t>/

사용(DockQ + biopython + scipy env):
  python pose_features.py --models boltz protenix --out results/pose_features.csv
  ⭐ 이어달리기: 이미 CSV에 있는 pose는 skip(--rescore로 강제). 타깃마다 flush(중단 안전).
"""
import argparse, csv, glob, json, os, tempfile

import epitope_shift as ES                      # scored_epitope_full, rank_score, centroid_dist, frac, popular_refset, fmt
from epitope_recall import (load as er_load, best as er_best,
                            antigen_refs, antibody_refs, native_true, neff_of)
from dockq_sweep import native_merged, pose_merged, dockq

COLS = ["target", "group", "ab", "label", "model", "rung", "neff80", "pose",
        "dockq", "recall", "n_contact", "overrep", "true_rank", "pop_rank", "dcc_true", "dcc_pop",
        "iptm", "ptm", "plddt"]

CONF_KEYS = {"iptm": "iptm", "ptm": "ptm", "complexplddt": "plddt", "plddt": "plddt"}


def find_confidence(pose_path):
    """pose 옆에 남은 confidence JSON에서 스칼라 신뢰도만 best-effort로 회수.
    모델마다 파일명·키가 달라 '있으면 줍는' 방식. 하나도 없으면 {} → 게이트 판정용으로 집계된다."""
    out = {}
    d = os.path.dirname(pose_path)
    stem = os.path.splitext(os.path.basename(pose_path))[0]
    cands = glob.glob(os.path.join(d, "*.json")) + glob.glob(os.path.join(os.path.dirname(d), "*.json"))
    for j in cands:
        b = os.path.basename(j).lower()
        if "confidence" not in b and stem.lower() not in b:
            continue
        try:
            data = json.load(open(j))
        except Exception:
            continue
        if not isinstance(data, dict):
            continue
        for k, v in data.items():
            kl = k.lower().replace("_", "")
            if kl in CONF_KEYS and isinstance(v, (int, float)):
                out.setdefault(CONF_KEYS[kl], float(v))
    return out


def pose_all_metrics(cj, pose_path, cutoff, true, pop):
    """epitope_shift.pose_metrics_full 과 동일 계산 + 접촉 잔기 수(n_contact) 추가."""
    m = er_load(pose_path)
    used = set(); ag = []
    for i, ref in enumerate(antigen_refs(cj)):
        cid, _, rr = er_best(m, ref, exclude=used)
        if cid is None:
            continue
        used.add(cid); ag.append((i, rr, ref))
    ab = []
    for ref in antibody_refs(cj):
        cid, _, rr = er_best(m, ref, exclude=used)
        if cid:
            used.add(cid); ab.extend(rr)
    if not ag or not ab:
        return None
    pred, dist, coord = ES.scored_epitope_full(ag, ab, cutoff)
    if not pred:
        return None
    nan = float("nan")
    return dict(
        n_contact=len(pred),
        recall=ES.frac(pred, true),
        overrep=ES.frac(pred, pop) if pop is not None else nan,
        true_rank=ES.rank_score(dist, true),
        pop_rank=ES.rank_score(dist, pop) if pop is not None else nan,
        dcc_true=ES.centroid_dist(coord, pred, true),
        dcc_pop=ES.centroid_dist(coord, pred, pop) if pop is not None else nan,
    )


def write_csv(out, rows):
    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLS); w.writeheader()
        for k in sorted(rows):
            w.writerow(rows[k])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--list", default="sweep_targets.csv")
    ap.add_argument("--targets-dir", default="targets")
    ap.add_argument("--data", default=os.environ.get("DATA", "/mnt/data/admuser/msadepth"))
    ap.add_argument("--models", nargs="+", default=["boltz", "protenix"])
    ap.add_argument("--rungs", type=int, default=12)
    ap.add_argument("--cutoff", type=float, default=5.0)
    ap.add_argument("--out", default="results/pose_features.csv")
    ap.add_argument("--rescore", action="store_true")
    a = ap.parse_args()
    import subprocess
    if subprocess.run(["which", "DockQ"], capture_output=True).returncode != 0:
        raise SystemExit("!! DockQ 없음 — DockQ 있는 env에서 실행")

    rows = {}
    if os.path.exists(a.out) and not a.rescore:
        for r in csv.DictReader(open(a.out)):
            rows[(r["target"], r["model"], int(r["rung"]), r["pose"])] = r
        print(f"기존 {len(rows)} pose 로드 → skip (강제 재채점=--rescore)")

    n_conf = n_tot = 0
    for r in csv.DictReader(open(a.list)):
        tgt = r["target"]; fam = r["group"]; abl = r.get("ab", "")
        cjp = os.path.join(a.targets_dir, tgt, "chains.json")
        native = os.path.join(a.targets_dir, tgt, "native.cif")
        if not os.path.exists(cjp):
            continue
        cj = json.load(open(cjp))
        tr = native_true(cj, native, a.cutoff)
        if tr is None:
            print(f"{tgt}: native epitope 실패 skip"); continue
        true, _ = tr
        try:
            pop = ES.popular_refset(cj, fam)
        except Exception:
            pop = None
        nmap = neff_of(tgt, os.path.join(a.data, "ladders"))

        with tempfile.TemporaryDirectory() as td:
            natm = native_merged(cj, native, td)
            if natm is None:
                print(f"{tgt}: native merge 실패 skip"); continue
            for model in a.models:
                for rung in range(a.rungs):
                    poses = sorted(glob.glob(os.path.join(
                        a.data, model, tgt, f"rung{rung}", "results", "**", "*.cif"), recursive=True))
                    for pose in poses:
                        pid = os.path.basename(pose)
                        key = (tgt, model, rung, pid)
                        if not a.rescore and key in rows:
                            continue
                        try:
                            met = pose_all_metrics(cj, pose, a.cutoff, true, pop)
                            if met is None:
                                continue
                            pm = pose_merged(cj, pose, td)
                            q = dockq(pm, natm) if pm else None
                        except Exception:
                            continue
                        conf = find_confidence(pose)
                        n_tot += 1; n_conf += 1 if conf else 0
                        rows[key] = dict(
                            target=tgt, group=fam, ab=abl, label=r.get("label", ""), model=model,
                            rung=rung, neff80=nmap.get(rung, ""), pose=pid,
                            dockq=ES.fmt(q), recall=ES.fmt(met["recall"]), n_contact=met["n_contact"],
                            overrep=ES.fmt(met["overrep"]), true_rank=ES.fmt(met["true_rank"]),
                            pop_rank=ES.fmt(met["pop_rank"]), dcc_true=ES.fmt(met["dcc_true"]),
                            dcc_pop=ES.fmt(met["dcc_pop"]),
                            iptm=ES.fmt(conf.get("iptm")), ptm=ES.fmt(conf.get("ptm")),
                            plddt=ES.fmt(conf.get("plddt")))
                    if poses:
                        print(f"  {tgt:14} {model:8} rung{rung:<2} pose {len(poses)} 처리 (누적 {len(rows)})")
        write_csv(a.out, rows)          # 타깃마다 flush

    write_csv(a.out, rows)
    print(f"\n→ {a.out} ({len(rows)} pose 행)")
    if n_tot:
        print(f"⭐ 게이트 — confidence 회수된 pose: {n_conf}/{n_tot} ({100*n_conf/n_tot:.0f}%)")
        print("   높으면: ipTM 대비 실험(Phase 1) + 신뢰도 피처 가능")
        print("   0%면 : co-folder confidence가 보존 안 됨 → 기하 피처(n_contact·dcc·rank)만으로 진행")
    print("다음(Phase 1): 복합체 내부에서 각 피처로 pose 순위 → ipTM·무작위·best-of-N null 대비 비교")


if __name__ == "__main__":
    main()
