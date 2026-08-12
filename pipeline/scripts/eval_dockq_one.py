#!/usr/bin/env python3
"""[낱개 채점] 자세 파일 몇 개의 DockQ 를 한 번에 잰다.

본 검정 표에 없는 예측(예: AlphaFold Server 결과)을 같은 기준으로 채점할 때 쓴다.
채점은 dockq_sweep 의 병합 규칙을 그대로 쓴다 — 항원을 사슬 A, 항체를 사슬 B 로 합친 뒤
DockQ CLI 에 넘긴다. 자세 쪽은 **서열 매칭**으로 사슬을 찾으므로 모델마다 사슬 기호가
달라도(AlphaFold 는 A/B, Boltz 는 또 다름) 그대로 넣으면 된다.

⚠️ DockQ CLI 가 있는 env 에서 돌린다(= conda activate boltz).

사용:
  python -u eval_dockq_one.py --target 8q7s_H \\
      --poses fold_8q7s_h_model_0.cif runs/.../noconstraint/....cif runs/.../ours/....cif
"""
import argparse, json, os, sys, tempfile


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", required=True, help="targets/<이름>/ 의 이름")
    ap.add_argument("--poses", nargs="+", required=True)
    ap.add_argument("--msad", default=os.path.expanduser("~/projects/bk21-msa-depth-bias/pipeline"))
    ap.add_argument("--labels", nargs="*", default=None, help="출력에 쓸 이름(자세와 같은 순서)")
    a = ap.parse_args()

    sys.path.insert(0, a.msad)
    from dockq_sweep import native_merged, pose_merged, dockq   # noqa: E402

    tdir = os.path.join(a.msad, "targets", a.target)
    cjp, natp = os.path.join(tdir, "chains.json"), os.path.join(tdir, "native.cif")
    for p in (cjp, natp):
        if not os.path.exists(p):
            raise SystemExit(f"!! 없음: {p}")
    cj = json.load(open(cjp))
    labs = a.labels if a.labels and len(a.labels) == len(a.poses) else [
        os.path.basename(p) for p in a.poses]

    with tempfile.TemporaryDirectory() as td:
        nat = native_merged(cj, natp, td)
        if not nat:
            raise SystemExit("!! native 병합 실패 — chains.json 의 src_chains 확인")
        print(f"{'자세':<34}{'DockQ':>8}")
        print("-" * 42)
        for lab, p in zip(labs, a.poses):
            if not os.path.exists(p):
                print(f"{lab:<34}{'파일없음':>8}"); continue
            pm = pose_merged(cj, p, td)
            if not pm:
                print(f"{lab:<34}{'병합실패':>8}"); continue
            q = dockq(pm, nat)
            print(f"{lab:<34}{(f'{q:.3f}' if q is not None else '  -  '):>8}")


if __name__ == "__main__":
    main()
