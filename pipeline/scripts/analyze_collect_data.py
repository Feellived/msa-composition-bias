#!/usr/bin/env python3
"""[결과 모으기] 흩어진 채점·분석 CSV 를 레포 한 곳에 모아 그림을 그릴 수 있게 한다.

왜 필요한가 — 결과가 두 레포와 서버 데이터 폴더에 흩어져 있고, 정작 그림을 그리려면
숫자 표만 있으면 된다. 구조 파일(.cif)은 크고 레포에 들어갈 물건이 아니다.
이 스크립트는 **작은 표만** 골라 report/data/ 로 복사하고, 무엇을 어디서 가져왔는지
MANIFEST.md 에 적는다. 나중에 그림을 다시 그릴 때 이 폴더만 보면 된다.

  · 기본이 dry-run 이다. 무엇을 가져올지 크기와 함께 보여주기만 한다.
  · 용량 상한(기본 5 MB)을 넘는 파일은 가져오지 않고 목록에만 남긴다.
  · 구조 파일 확장자(.cif .pdb .a3m .npz .pkl)는 아예 후보에서 뺀다.
  · 같은 이름이 겹치면 앞에 출처 이름을 붙여 구분한다(demo_dockq.csv → cd_demo_dockq.csv).

사용:
  python -u analyze_collect_data.py                      # 무엇을 가져올지만 출력
  python -u analyze_collect_data.py --apply
  python -u analyze_collect_data.py --apply --max-mb 20  # 큰 표까지
"""
import argparse, csv, glob, hashlib, os, shutil

HOME = os.path.expanduser("~")
CD = os.path.join(HOME, "projects/bk21-antibody-ml/pipeline")
MS = os.path.join(HOME, "projects/bk21-msa-depth-bias/pipeline")

# (앞에 붙일 이름, 설명, 찾을 곳) — 앞에 붙일 이름은 파일명이 겹칠 때만 쓴다
SOURCES = [
    ("cd", "유도 재도킹 채점 (네 팔 · 30종)", f"{CD}/results/*.csv"),
    ("cd", "후보 전수 재도킹 채점 (후보별)", f"{CD}/results/allcand/*.csv"),
    ("ms", "본 검정·선택기·신호 분석", f"{MS}/results/*.csv"),
    ("ms", "정답을 제거한 판 (본 검정)", f"{MS}/results/honest/*.csv"),
    ("ms", "후보 자리 정의 (정답 제거판)", f"{MS}/results/honest/sites_*.json"),
]
SKIP_EXT = {".cif", ".pdb", ".ent", ".a3m", ".fasta", ".npz", ".pkl", ".pt", ".log"}


def sha8(p):
    h = hashlib.sha1()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()[:8]


def shape(p):
    """CSV 면 (줄 수, 열 이름). 아니면 (None, None)."""
    if not p.endswith(".csv"):
        return None, None
    try:
        with open(p, newline="") as fh:
            r = csv.reader(fh)
            head = next(r, [])
            n = sum(1 for _ in r)
        return n, head
    except Exception:
        return None, None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(HOME, "projects/bk21-msa-depth-bias/report/data"))
    ap.add_argument("--max-mb", type=float, default=5.0)
    ap.add_argument("--apply", action="store_true")
    a = ap.parse_args()
    cap = a.max_mb * (1 << 20)

    seen, take, toobig, dup = set(), [], [], 0
    for pre, desc, pat in SOURCES:
        for p in sorted(glob.glob(pat)):
            if not os.path.isfile(p):
                continue
            if os.path.splitext(p)[1].lower() in SKIP_EXT:
                continue
            sz = os.path.getsize(p)
            base = os.path.basename(p)
            if sz > cap:
                toobig.append((p, sz)); continue
            name = base
            if name in seen:                       # 이름이 겹치면 출처를 붙인다
                name = f"{pre}_{base}"; dup += 1
            if name in seen:                       # 그래도 겹치면 폴더까지
                name = f"{pre}_{os.path.basename(os.path.dirname(p))}_{base}"
            seen.add(name)
            take.append((name, p, sz, desc))

    print(f"모을 파일 {len(take)}개 · 이름 충돌로 접두어 붙인 것 {dup}개 "
          f"· 상한({a.max_mb} MB) 초과로 제외 {len(toobig)}개")
    print(f"보낼 곳 {a.out}\n")
    print(f"  {'파일':<34}{'KB':>8}{'줄':>8}  설명")
    print("  " + "-" * 78)
    tot = 0
    for name, p, sz, desc in take:
        n, _ = shape(p)
        tot += sz
        print(f"  {name:<34}{sz/1024:>8.1f}{(n if n is not None else '-'):>8}  {desc}")
    print(f"\n  합계 {tot/1024:.1f} KB")
    if toobig:
        print(f"\n  ! 상한 초과로 제외 (필요하면 --max-mb 를 키울 것)")
        for p, sz in toobig:
            print(f"    {sz/(1<<20):>8.1f} MB  {p}")
    if not a.apply:
        print("\n[dry-run — 아무것도 복사하지 않음. 가져오려면 --apply]")
        return

    os.makedirs(a.out, exist_ok=True)
    lines = ["# 그림용 데이터", "",
             "`pipeline/analyze_collect_data.py --apply` 가 만든 폴더다. 표만 모여 있고 구조 파일은 없다.",
             "", "| 파일 | 줄 | 열 | 원본 | 해시 |", "|---|---|---|---|---|"]
    for name, p, sz, desc in take:
        shutil.copy2(p, os.path.join(a.out, name))
        n, head = shape(p)
        cols = ", ".join(head[:8]) + (" …" if head and len(head) > 8 else "") if head else "-"
        src = p.replace(HOME, "~")
        lines.append(f"| `{name}` | {n if n is not None else '-'} | {cols} | `{src}` | `{sha8(p)}` |")
    lines += ["", "## 출처별 설명", ""]
    for pre, desc, pat in SOURCES:
        lines.append(f"- {desc} — `{pat.replace(HOME, '~')}`")
    with open(os.path.join(a.out, "MANIFEST.md"), "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"\n→ {a.out}  ({len(take)}개 + MANIFEST.md)")
    print("  다음: git add report/data && git commit && git push")


if __name__ == "__main__":
    main()
