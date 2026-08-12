#!/usr/bin/env python3
"""[진단] 4.8절 "DockQ를 잴 수 없었던 9b7g_QP" 의 원인을 사슬 단위로 들여다본다.

eval_dockq_one.py 가 쓰는 native_merged / pose_merged 는 항원·항체 사슬을
chains.json 의 서열과 가장 잘 맞는 사슬로 골라 붙인다(eval_dockq_sweep.py 의 best()).
문턱이 없어 아무리 안 맞아도 "제일 나은 것"을 고르므로, 실패는 대개 둘 중 하나다.
  ① AF3 사슬 자체가 항원의 일부만 담고 있어 정렬 점수가 크게 낮다(길이·서열이 안 맞음)
  ② 병합/DockQ CLI 단계에서 예외가 나 조용히 None 이 된다

이 스크립트는 항원 사슬마다 (매칭된 사슬 ID, 정렬 점수, native 길이, pose 길이,
pose 가 native 를 덮는 비율)을 출력해 어느 쪽인지 가른다.

사용(DockQ+biopython env · pipeline/ 에서):
  python -u analyze_dockq_fail.py --target 9b7g_QP \
      --pose <AF3 예측 cif 경로> \
      [--pose <자세 2> --pose <자세 3> ...]

--pose 를 안 주면 targets/9b7g_QP/ 밑에서 *.cif 를 찾는다(AF3 결과를 거기 둔 경우).
"""
import argparse
import glob
import json
import os
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dockq_sweep import load, seqca, best, native_merged, pose_merged, dockq  # noqa: E402


def describe(model, cj, td):
    ag_seqs = [(c["seq"], c.get("crop")) for c in cj["chains"] if c["role"] == "antigen"]
    used = set()
    for i, (seq, crop) in enumerate(ag_seqs):
        cid, score, rr = best(model, seq, exclude=used)
        if cid is None:
            print(f"    항원#{i}: 맞는 사슬을 못 찾음 (모델에 5잔기 이상인 사슬이 없음)")
            continue
        used.add(cid)
        cov = len(rr) / max(len(seq), 1)
        print(f"    항원#{i}: 사슬 {cid}  정렬점수 {score:+.3f}  "
              f"정답서열길이 {len(seq)}  매칭잔기수 {len(rr)}  덮음비율 {cov:.2f}"
              + ("  ⚠️ 낮음(부분예측 의심)" if score < 0.5 or cov < 0.5 else ""))
    for c in cj["chains"]:
        if c["role"] in ("heavy", "light"):
            cid, score, rr = best(model, c["seq"], exclude=used)
            tag = c["role"]
            if cid is None:
                print(f"    {tag}: 맞는 사슬을 못 찾음"); continue
            used.add(cid)
            print(f"    {tag}: 사슬 {cid}  정렬점수 {score:+.3f}  매칭잔기수 {len(rr)}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", required=True)
    ap.add_argument("--pose", action="append", default=[])
    ap.add_argument("--targets-dir", default="targets")
    a = ap.parse_args()

    cjp = os.path.join(a.targets_dir, a.target, "chains.json")
    natp = os.path.join(a.targets_dir, a.target, "native.cif")
    if not os.path.exists(cjp):
        raise SystemExit(f"!! 없음: {cjp}")
    cj = json.load(open(cjp))

    poses = a.pose or sorted(glob.glob(os.path.join(a.targets_dir, a.target, "*.cif")))
    if not poses:
        raise SystemExit("!! --pose 를 안 줬고 targets/<타깃>/ 밑에도 .cif 가 없다")

    print(f"=== {a.target} · native ===")
    nm = load(natp)
    with tempfile.TemporaryDirectory() as td:
        describe(nm, cj, td)

        natm = native_merged(cj, natp, td)
        print(f"\nnative_merged: {'성공 → ' + natm if natm else '❌ 실패(None)'}")

        for p in poses:
            print(f"\n=== 자세: {p} ===")
            model = load(p)
            describe(model, cj, td)
            pm = pose_merged(cj, p, td)
            print(f"  pose_merged: {'성공' if pm else '❌ 실패(None)'}")
            if pm and natm:
                try:
                    q = dockq(pm, natm)
                    print(f"  DockQ = {q if q is not None else '❌ CLI 가 값을 못 냄'}")
                except Exception as e:
                    print(f"  ❌ DockQ 실행 중 예외: {e}")


if __name__ == "__main__":
    main()
