#!/usr/bin/env python3
"""[복구] 사다리·seed a3m의 질의행에 붙은 메타 주석(`#<길이>\\t<개수>`)을 떼어낸다.

원인: `prep_ladder_neff.py`의 옛 `read_raw`가 ColabFold a3m 첫 줄의 메타 주석을 서열로 취급하고
첫 `>`에서 버퍼를 비우지 않아, 주석이 **질의 서열 앞에 붙어** 저장됨.
  >101
  #440\t1NLWVTVYYGVPVWK...      ← 앞 6글자가 군더더기
결과: boltz는 "MSA does not match input sequence"로 **MSA를 통째로 버림**(더미 대체),
Protenix는 질의행만 밀린 정렬을 그대로 사용.

서열 자체는 멀쩡하므로 **질의행 앞 군더더기만 잘라내면 복구된다**(MSA 재생성 불필요).
같은 파일에 두 번 돌려도 안전(멱등).

⚠️ 기본은 dry-run이다. 실제 수정은 `--apply`. 내 폴더 밖은 건드리지 않는다.

사용:
  python prep_a3m_fix_query.py                        # dry-run(기본): 몇 개가 어떻게 바뀌는지만
  python prep_a3m_fix_query.py --apply                # 실제 수정
  python prep_a3m_fix_query.py --roots /경로1 /경로2   # 대상 폴더 지정
"""
import argparse, os, re, sys

PAT = re.compile(r"^#\d+\t\d+")          # 관측된 형태: #440\t1
PAT_LOOSE = re.compile(r"^#\S*\s*\d*")   # 변형 대비(보고만 하고 자동 수정 안 함)


def first_seq_line_idx(lines):
    """첫 '>' 헤더 다음의 첫 서열 줄 index (없으면 None)."""
    seen_header = False
    for i, ln in enumerate(lines):
        if ln.startswith(">"):
            if seen_header:
                return None          # 서열 없이 헤더가 연달아 나옴
            seen_header = True
            continue
        if seen_header:
            return i
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--roots", nargs="*", default=None,
                    help="검사할 폴더들(기본: $DATA/ladders, $DATA/seedrep, ./seedrep_cand)")
    ap.add_argument("--data", default=os.environ.get("DATA", "/mnt/data/admuser/msadepth"))
    ap.add_argument("--apply", action="store_true", help="실제로 파일을 고친다(기본은 dry-run)")
    a = ap.parse_args()

    roots = a.roots or [os.path.join(a.data, "ladders"),
                        os.path.join(a.data, "seedrep"),
                        "seedrep_cand"]
    roots = [r for r in roots if os.path.isdir(r)]
    if not roots:
        sys.exit("!! 검사할 폴더가 없음. --roots 로 지정할 것.")
    print(("[실제 수정 모드]" if a.apply else "[dry-run — 파일을 고치지 않음]")
          + "  대상: " + ", ".join(roots) + "\n")

    files = []
    for r in roots:
        for dp, _, fns in os.walk(r):
            files += [os.path.join(dp, f) for f in fns if f.endswith(".a3m")]
    files.sort()

    n_ok = n_fix = n_weird = 0
    sample = None
    for p in files:
        try:
            with open(p) as f:
                lines = f.read().split("\n")
        except Exception as e:
            print(f"  !! 읽기 실패 {p}: {e}"); continue
        i = first_seq_line_idx(lines)
        if i is None:
            n_weird += 1; continue
        ln = lines[i]
        if not ln.startswith("#"):
            n_ok += 1; continue
        m = PAT.match(ln)
        if not m:
            n_weird += 1
            print(f"  ?? 예상과 다른 머리말, 손대지 않음: {p}\n     {ln[:50]}")
            continue
        new = ln[m.end():]
        if sample is None:
            sample = (p, ln[:46], new[:40], m.end())
        n_fix += 1
        if a.apply:
            lines[i] = new
            with open(p, "w") as f:
                f.write("\n".join(lines))

    print(f"검사 {len(files)}개 · 이미 정상 {n_ok} · {'고침' if a.apply else '고칠 것'} {n_fix}"
          + (f" · 확인필요 {n_weird}" if n_weird else ""))
    if sample:
        p, before, after, k = sample
        print(f"\n예시: {p}\n  전: {before}...\n  후: {after}...   (앞 {k}글자 제거)")
    if n_fix and not a.apply:
        print("\n→ 실제로 고치려면 같은 명령에 --apply 를 붙여 실행.")
    if a.apply and n_fix:
        print("\n✅ 복구 완료. 검증: python prep_a3m_check_match.py")
        print("⚠️ 예측 결과는 오염된 입력으로 만들어진 것이므로, 쓰려는 실험은 다시 돌려야 함.")


if __name__ == "__main__":
    main()
