#!/usr/bin/env python3
"""[D] 가짜 쌍 음성 대조 — '붙는 것끼리'만 넣어 온 판에 '안 붙는 쌍'을 넣는다.

■ 왜
  지금까지 모든 시험은 **붙는다는 게 이미 확정된 항체-항원 쌍**만 넣었다. 그런데 랩 실무에서
  정작 필요한 건 "이게 붙긴 붙나?"다. Beacon 스크리닝에서 나온 후보 중 진짜 결합체를
  미리 거를 수 있다면 그것만으로 쓸모가 있다.

■ 방법
  항원은 그대로 두고 **다른 항원군의 항체**를 붙인다(RBD 항원 + HA 항체 등). 정의상 안 붙는다.
  그다음 우리 지표(모델 신뢰도 · 조성 간 합의 정도)가 진짜 쌍과 가짜 쌍을 가르는지 본다.

■ 이 설계가 값싼 이유
  · **정답 구조가 필요 없다** — 가짜 쌍은 정의상 답이 없다. native.cif 를 만들지 않는다.
  · **항원 MSA 를 그대로 재사용한다** — $DATA/ladders/<항원타깃>/ 이 이미 있으므로
    MSA 생성(make_msa.sh)을 다시 돌릴 필요가 없다. 바뀌는 건 항체 사슬뿐이다.
  · 진짜 쌍 쪽 값은 이미 다 있다. 새로 돌리는 건 가짜 쌍뿐이다.

■ 판정 기준 (plan/PLAN_202608_final8days.md §6, 결과 보기 전 고정)
  성공 = 진짜 쌍 대 가짜 쌍의 판별력 AUC ≥ 0.70 · 가짜 쌍 24개 이상

■ ⚠️ 안전
  이 스크립트는 **기본이 미리보기(dry-run)** 다. 실제로 폴더를 만들려면 --write 를 준다.
  만드는 건 targets/decoy_* 뿐이고 기존 타깃 폴더는 건드리지 않는다.
  그다음 예측은 GPU 를 오래 쓴다 — **돌리기 전에 반드시 허가를 받을 것.**

사용 (pipeline/ 에서):
  python -u prep_decoy_pairs.py --n 24                       # 미리보기
  python -u prep_decoy_pairs.py --n 24 --write               # 실제 생성
"""
import argparse
import csv
import json
import os
import random
import shutil
from collections import defaultdict


def chains_of(cj, roles):
    return [c for c in cj["chains"] if c["role"] in roles]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", default="targets_manifest.csv")
    ap.add_argument("--targets-dir", default="targets")
    ap.add_argument("--sites", default="results",
                    help="여기 sites_*.json 이 있는 타깃만 쓴다(예측이 끝난 것)")
    ap.add_argument("--n", type=int, default=24, help="만들 가짜 쌍 개수")
    ap.add_argument("--per-antigen", type=int, default=2, help="항원 하나당 가짜 항체 몇 개")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out-manifest", default="results/decoy_manifest.csv")
    ap.add_argument("--write", action="store_true", help="실제로 폴더를 만든다(기본은 미리보기)")
    a = ap.parse_args()

    man = {r["target"]: r for r in csv.DictReader(open(a.manifest))}
    ok = []
    for t, r in man.items():
        if not os.path.exists(os.path.join(a.targets_dir, t, "chains.json")):
            continue
        if not os.path.exists(os.path.join(a.sites, f"sites_{t}.json")):
            continue
        ok.append(t)
    if not ok:
        raise SystemExit(f"!! 쓸 타깃이 없다 — {a.targets_dir}/<t>/chains.json 과 "
                         f"{a.sites}/sites_<t>.json 이 둘 다 있어야 한다")

    grp = defaultdict(list)
    for t in sorted(ok):
        grp[man[t]["antigen_grp"]].append(t)
    grp = {g: v for g, v in grp.items() if g and g != "??"}
    print("항원군: " + " · ".join(f"{g}({len(v)})" for g, v in sorted(grp.items())))
    if len(grp) < 2:
        raise SystemExit("!! 항원군이 2개 미만이라 다른 군의 항체를 빌려올 수 없다")

    rng = random.Random(a.seed)
    made, seen = [], set()
    ags = [t for g in sorted(grp) for t in grp[g]]
    rng.shuffle(ags)
    for ag in ags:
        if len(made) >= a.n:
            break
        gag = man[ag]["antigen_grp"]
        pool = [t for g, v in grp.items() if g != gag for t in v]
        rng.shuffle(pool)
        k = 0
        for ab in pool:
            if k >= a.per_antigen or len(made) >= a.n:
                break
            did = f"decoy_{ag}_x_{ab}"
            if did in seen:
                continue
            seen.add(did)
            made.append(dict(decoy=did, ag_target=ag, ag_group=gag,
                             ab_target=ab, ab_group=man[ab]["antigen_grp"]))
            k += 1

    if len(made) < a.n:
        print(f"⚠️ 요청한 {a.n}개 중 {len(made)}개만 만들 수 있었다 "
              f"(항원군 조합이 부족). --per-antigen 을 올려볼 것.")

    print(f"\n{'가짜 쌍':<34}{'항원(군)':<22}{'빌려온 항체(군)':<22}")
    print("-" * 78)
    for m in made:
        print(f"{m['decoy']:<34}{m['ag_target']+' ('+m['ag_group']+')':<22}"
              f"{m['ab_target']+' ('+m['ab_group']+')':<22}")

    if not a.write:
        print(f"\n※ 미리보기다. 실제로 만들려면 --write 를 붙일 것.")
    else:
        n_ok = 0
        for m in made:
            src = os.path.join(a.targets_dir, m["ag_target"])
            dst = os.path.join(a.targets_dir, m["decoy"])
            cj_ag = json.load(open(os.path.join(src, "chains.json")))
            cj_ab = json.load(open(os.path.join(a.targets_dir, m["ab_target"], "chains.json")))
            ab = chains_of(cj_ab, ("heavy", "light"))
            if not ab:
                print(f"  ! {m['decoy']}: 빌려올 항체 사슬이 없다 — 건너뜀"); continue
            # 항원 사슬은 순서까지 그대로 둔다 — 잔기 키가 원래 타깃과 같아야 비교가 된다.
            ag_chains = chains_of(cj_ag, ("antigen",))
            used_ids = {c["id"] for c in ag_chains}
            # ⚠️ 사슬 id 는 prep_targets.py 가 'A,B,...'로 위치 기준 배정 — 다사슬 항원(HA1+HA2 등)이면
            #    빌려온 항체 id 가 항원 id 와 겹칠 수 있다. 겹치면 새 id 로 다시 붙인다.
            ab_new, relabel = [], {}
            next_ord = ord("A") + len(ag_chains)
            for c in ab:
                cid = c["id"]
                if cid in used_ids:
                    while chr(next_ord) in used_ids:
                        next_ord += 1
                    relabel[cid] = chr(next_ord)
                    used_ids.add(chr(next_ord)); next_ord += 1
                else:
                    used_ids.add(cid)
                c2 = dict(c)
                c2["id"] = relabel.get(cid, cid)
                ab_new.append(c2)
            new = dict(cj_ag)
            new["chains"] = ag_chains + ab_new
            new["antigen"] = [c["id"] for c in ag_chains]        # ⚠️ make_input.py 가 이 키로 항원 사슬을 찾는다
            new["antibody"] = [c["id"] for c in ab_new]          # ⚠️ 원래 ag_target 의 항체 id 를 그대로 남기면
                                                                  #    make_input.py 가 존재하지 않는 id 를 찾다 KeyError
            new["decoy_of"] = m["ag_target"]
            new["antibody_from"] = m["ab_target"]
            new.pop("src_chains", None)         # 원본 항원의 항체 대응 정보라 가짜 쌍엔 의미가 없다
            # ⚠️ 정답이 없다는 사실을 파일에 박아 둔다. 채점 스크립트가 실수로 native 를 찾지 않게.
            new["no_native"] = True
            if relabel:
                print(f"  · {m['decoy']}: 항체 사슬 id 재배정 {relabel} (항원과 겹쳐서)")
            os.makedirs(dst, exist_ok=True)
            json.dump(new, open(os.path.join(dst, "chains.json"), "w"), indent=1)
            fa = os.path.join(src, "antigen.fasta")
            if os.path.exists(fa):
                shutil.copy(fa, os.path.join(dst, "antigen.fasta"))
            n_ok += 1
        print(f"\n→ {a.targets_dir}/decoy_*  ({n_ok}개 생성 · native.cif 없음 = 정답 없음)")

    os.makedirs(os.path.dirname(a.out_manifest) or ".", exist_ok=True)
    with open(a.out_manifest, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(made[0])); w.writeheader(); w.writerows(made)
    print(f"→ {a.out_manifest}")

    print("\n■ 다음 단계 (⚠️ GPU 를 오래 쓴다 — 돌리기 전에 허가를 받을 것)")
    print("  1) 항원 MSA 사다리는 새로 만들 필요가 없다. 가짜 쌍의 항원은 원래 타깃과 같으므로")
    print("     $DATA/ladders/<ag_target>/ 을 그대로 쓴다 (manifest 의 ag_target 열).")
    print("  2) 조성 3가지 × 반복 2회 정도로 시작한다. 진짜 쌍 쪽 값은 이미 있다.")
    print("  3) 채점은 정답 없이 한다 — 모델 신뢰도(ipTM)와 조성 간 합의만 본다.")
    print("     DockQ·결합자리 덮음은 계산하지 말 것(정답이 없다).")


if __name__ == "__main__":
    main()
