#!/usr/bin/env python3
"""[9b7g_QP 전용 시도] AF3 예측의 항체 사슬이 정답 구조보다 길 때 DockQ를 재시도한다.

analyze_dockq_fail.py 로 확인한 원인: 9b7g_QP 의 항원은 잘 맞는데(덮음 0.96),
AF3 에 제출한 항체 서열이 정답 구조가 담은 가변 영역(H 115 · L 109잔기)보다
약 2배 길다(H 242 · L 216잔기, 정렬점수 0.19~0.23). eval_dockq_sweep.py 의 pose_merged
는 항원에만 "정답 범위로 자르기"가 있고 항체에는 없어서, 긴 사슬을 통째로 붙여
DockQ CLI 가 잔기 수 불일치로 값을 못 낸다(native_merged/pose_merged 는 그대로 둔다
— 다른 29종은 이미 정상 작동하므로 건드리지 않는다).

이 스크립트는 항체 사슬만 국소 정렬(local alignment)로 정답과 대응하는 구간을
찾아 잘라낸 뒤 병합해 DockQ를 다시 시도한다. 항원은 기존 방식(전체 사용) 그대로다.

사용(DockQ+biopython env · pipeline/ 에서):
  python -u eval_dockq_af3_crop.py --target 9b7g_QP \
      --pose /home/user/projects/epitope-guided-docking/pipeline/fold_9b7g_qp_model_0.cif
"""
import argparse
import json
import os
import tempfile

from Bio.Align import PairwiseAligner
from dockq_sweep import load, seqca, best, merged_pdb, dockq, native_merged  # noqa: E402

_local = PairwiseAligner()
_local.mode = "local"
_local.match_score = 2
_local.mismatch_score = -1
_local.open_gap_score = -10
_local.extend_gap_score = -0.5


def crop_to_target(rr, seq, tgt_seq):
    """rr·seq = 모델 사슬 전체(잔기 목록·서열). tgt_seq(정답 서열)와 국소 정렬해
    실제로 대응하는 구간의 잔기만 추린다."""
    aln = _local.align(seq, tgt_seq)[0]
    idx = []
    for s0, e0 in aln.aligned[0]:
        idx.extend(range(s0, e0))
    return [rr[i] for i in idx], len(idx) / max(len(tgt_seq), 1)


def pose_merged_ab_crop(cj, pose_path, td):
    model = load(pose_path)
    ag_seqs = [c["seq"] for c in cj["chains"] if c["role"] == "antigen"]
    ab_seqs = [(c["role"], c["seq"]) for c in cj["chains"] if c["role"] in ("heavy", "light")]
    used = set()
    ag = []
    for seq in ag_seqs:
        cid, score, rr = best(model, seq, exclude=used)
        if cid is None:
            print(f"  ! 항원 매칭 실패"); return None
        used.add(cid); ag.append(rr)
        print(f"  항원: 사슬 {cid} 정렬점수 {score:+.3f} (그대로 사용, {len(rr)}잔기)")
    ab = []
    for role, seq in ab_seqs:
        cid, score, rr = best(model, seq, exclude=used)
        if cid is None:
            print(f"  ! {role} 매칭 실패"); continue
        used.add(cid)
        cropped, cov = crop_to_target(rr, seqca(model[cid])[0], seq)
        print(f"  {role}: 사슬 {cid} 전체 {len(rr)}잔기 → 정답과 국소 정렬해 {len(cropped)}잔기로 자름 "
              f"(정답 {len(seq)}잔기 대비 덮음 {cov:.2f})")
        if cropped:
            ab.append(cropped)
    if not ab or not any(ag):
        return None
    out = os.path.join(td, "pose_cropped.pdb")
    merged_pdb(ag, ab, out)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", required=True)
    ap.add_argument("--pose", required=True)
    ap.add_argument("--targets-dir", default="targets")
    a = ap.parse_args()

    cj = json.load(open(os.path.join(a.targets_dir, a.target, "chains.json")))
    native_path = os.path.join(a.targets_dir, a.target, "native.cif")

    with tempfile.TemporaryDirectory() as td:
        natm = native_merged(cj, native_path, td)
        if natm is None:
            raise SystemExit("!! native_merged 실패")
        print(f"native_merged: 성공")

        pm = pose_merged_ab_crop(cj, a.pose, td)
        if pm is None:
            raise SystemExit("!! pose 병합 실패")

        q = dockq(pm, natm)
        print(f"\nDockQ (항체 자르기 적용) = {q if q is not None else '❌ 여전히 값을 못 냄'}")


if __name__ == "__main__":
    main()
