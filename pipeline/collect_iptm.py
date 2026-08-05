#!/usr/bin/env python3
"""[선택기 추가] 예측 자세마다 Protenix 가 남긴 ipTM 을 모아 abepiscore_all.csv 와
같은 형식(path,target,run,pose,score)으로 저장한다 → eval_selectors.py 가 그대로 먹는다.

왜 — 4.6절은 "선택기 넷 중 ipTM 은 우연과 구별되지 않았다"고 쓰지만, 선택기 비교표
(eval_selectors.py)에는 ipTM 항목이 없다. 같은 판에서 재보려면 자세별 ipTM 이 필요하다.

Protenix 는 예측 폴더에 자세(.cif)와 나란히 신뢰도 json 을 남긴다. 파일 이름은 버전마다
다르므로(보통 *summary_confidence_sample_N.json) 이름을 고정하지 않고, 같은 폴더에서
sample 번호가 일치하는 json 을 찾아 iptm 키를 재귀로 뒤진다. 못 찾은 자세는 건너뛰고
몇 개를 건너뛰었는지 끝에 보고한다(조용히 비지 않게).

사용 (conda activate boltz · pipeline/ 에서):
  python -u collect_iptm.py --data $DATA/compreps/seedrep_cand/protenix \
      --out results/iptm_all.csv
"""
import argparse
import csv
import glob
import json
import os
import re

IPTM_KEYS = ("iptm", "interface_ptm", "ptm_interface")


def find_iptm(obj):
    """중첩 dict 어디에 있든 ipTM 값을 찾아낸다. chain 쌍별 행렬이면 평균."""
    if isinstance(obj, dict):
        for k, v in obj.items():
            if k.lower() in IPTM_KEYS:
                if isinstance(v, (int, float)):
                    return float(v)
                if isinstance(v, list):
                    flat = [x for x in _flatten(v) if isinstance(x, (int, float))]
                    if flat:
                        return sum(flat) / len(flat)
            got = find_iptm(v)
            if got is not None:
                return got
    elif isinstance(obj, list):
        for v in obj:
            got = find_iptm(v)
            if got is not None:
                return got
    return None


def _flatten(x):
    if isinstance(x, list):
        for v in x:
            yield from _flatten(v)
    else:
        yield x


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True, help="…/compreps/seedrep_cand/protenix")
    ap.add_argument("--out", default="results/iptm_all.csv")
    a = ap.parse_args()

    if not os.path.isdir(a.data):
        raise SystemExit(f"!! 폴더가 없다: {a.data}")

    rows, miss, seen_json = [], 0, set()
    cifs = sorted(glob.glob(os.path.join(a.data, "*", "*", "*", "results", "**", "*.cif"),
                            recursive=True))
    print(f"자세 파일 {len(cifs)}개를 훑는다\n")
    for cif in cifs:
        # 경로: …/protenix/<타깃>/<깊이>/<실행>/results/…/predictions/<이름>_sample_N.cif
        m = re.search(r"/protenix/([^/]+)/[^/]+/([^/]+)/results/", cif)
        if not m:
            continue
        target, run = m.group(1), m.group(2)
        sm = re.search(r"sample_(\d+)\.cif$", os.path.basename(cif))
        idx = sm.group(1) if sm else None

        d = os.path.dirname(cif)
        cands = [p for p in glob.glob(os.path.join(d, "*.json"))
                 if idx is None or re.search(rf"sample_{idx}\b", os.path.basename(p))]
        cands.sort(key=lambda p: ("summary" not in os.path.basename(p).lower(), len(p)))

        val = None
        for j in cands:
            try:
                val = find_iptm(json.load(open(j)))
            except Exception:
                val = None
            if val is not None:
                seen_json.add(re.sub(r"sample_\d+", "sample_N", os.path.basename(j)))
                break
        if val is None:
            miss += 1
            continue
        rows.append(dict(path=cif, target=target, run=run,
                         pose=os.path.basename(cif), score=round(val, 5)))

    if not rows:
        raise SystemExit("!! ipTM 을 하나도 못 찾았다 — 예측 폴더에 신뢰도 json 이 있는지 확인할 것")

    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["path", "target", "run", "pose", "score"])
        w.writeheader(); w.writerows(rows)

    tg = sorted({r["target"] for r in rows})
    print(f"찾은 json 이름 꼴: {sorted(seen_json)}")
    print(f"ipTM 수집 {len(rows)}자세 · 타깃 {len(tg)}종 · 못 찾은 자세 {miss}개")
    print(f"→ {a.out}")


if __name__ == "__main__":
    main()
