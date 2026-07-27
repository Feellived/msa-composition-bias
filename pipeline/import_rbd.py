#!/usr/bin/env python3
"""[가져오기] Consensus Docking 쪽 runs_rbd(결합부위 밖 항체 세트)를 깊이실험 형식으로 등록.

왜: 이 세트는 "인기 부위를 피해 붙는 항체"만 모아둔 것이라, 우리가 44개에서 겨우 8개
찾아낸 조건(B군)을 11개 전부가 만족한다. 항원 MSA도 이미 만들어져 있어(msa_<t>/A_env),
가장 비싼 단계를 건너뛸 수 있다. 없는 것은 정답 구조뿐이라 RCSB에서 받는다.

하는 일 (전부 '추가'만 — 기존 파일을 지우거나 덮어쓰지 않는다):
  ① runs_rbd/<pdb>/chains.json 읽어 항원/항체 사슬 판별 → 타깃 이름 <pdb>_<항체사슬>
  ② targets/<이름>/chains.json 작성 (AB 표식 포함)
  ③ targets/<이름>/native.cif 를 RCSB에서 내려받음 (이미 있으면 건너뜀)
  ④ 기존 항원 a3m 을 $DATA/ladders/<이름>/<사슬>/rung0.a3m 으로 복사 (머리말 오염 제거)
  ⑤ sweep_targets.csv 에 행 추가

⚠️ 기본은 dry-run. 실제 기록은 --apply.

사용:
  python import_rbd.py                          # 무엇을 할지만 출력
  python import_rbd.py --apply
  python import_rbd.py --src ~/projects/bk21-antibody-ml/consensus_docking/runs_rbd
"""
import argparse, csv, glob, json, os, re, shutil, sys, urllib.request

# 노션 1.3절 표에 기록된 '붙는 자리' 분류 → A/B 표식.
#   B = 진짜 결합자리가 인기 부위 **밖** (편향이 풀리면 좋아져야 함)
#   A = 진짜 결합자리가 인기 부위와 겹침 (편향이 풀리면 오히려 나빠져야 함 = 반대방향 대조군)
SITE = {
    "9zdu": ("코어",              "B"),
    "9ml9": ("코어",              "B"),
    "8siq": ("코어",              "B"),
    "8sit": ("코어",              "B"),
    "8sis": ("코어와 그 밖",       "B"),
    "8xsi": ("숨은면",            "B"),
    "9ml8": ("숨은면",            "B"),
    "8sdf": ("숨은면",            "B"),
    "9sbb": ("그 밖",             "B"),
    "8p5m": ("인기 부위 가장자리", "A"),   # ★ 반대방향 대조군
    "8sdh": ("미기록",            "?"),   # 노션 표에 없던 폴더 — 분류 확인 필요
}

AA = set("ACDEFGHIKLMNPQRSTVWYXBZUO")
PAT = re.compile(r"^#\d+\t\d+")


def norm(s):
    return "".join(c for c in s.upper() if c not in "-.")


def read_a3m(path):
    """a3m 전체를 (헤더, 서열) 목록으로. 질의행 머리말 오염(2026-07-27 버그)도 제거."""
    heads, seqs, cur, h = [], [], [], None
    for ln in open(path):
        ln = ln.rstrip("\n")
        if not ln:
            continue
        if ln[0] == "#":
            rest = PAT.sub("", ln)
            if not rest:
                continue
            ln = rest
        if ln[0] == ">":
            if h is not None:
                heads.append(h); seqs.append("".join(cur))
            cur = []; h = ln[1:]
        else:
            cur.append(ln)
    if h is not None:
        heads.append(h); seqs.append("".join(cur))
    return heads, seqs


def chain_seqs(cj):
    """chains.json에서 {사슬ID: 서열} 추출 (형식이 여러 가지라 방어적으로)."""
    seqs = {}

    def looks(v):
        return isinstance(v, str) and len(v) >= 20 and set(v.upper()) <= AA

    def walk(node):
        if isinstance(node, dict):
            for k, v in node.items():
                if looks(v) and len(str(k)) <= 4:
                    seqs[str(k)] = v.upper()
                elif isinstance(v, dict) and looks(v.get("sequence", "")) and len(str(k)) <= 4:
                    seqs[str(k)] = v["sequence"].upper()
                else:
                    walk(v)
        elif isinstance(node, list):
            for it in node:
                if isinstance(it, dict):
                    cid = it.get("id") or it.get("chain") or it.get("chain_id")
                    sq = it.get("sequence") or it.get("seq")
                    if cid and looks(sq or ""):
                        seqs[str(cid)] = sq.upper()
                    else:
                        walk(it)
                else:
                    walk(it)
    walk(cj)
    return seqs


def pick_roles(cj, seqs):
    """항원/항체 사슬 목록. chains.json에 표기가 있으면 쓰고, 없으면 길이로 추정."""
    ag, ab = [], []
    if isinstance(cj, dict):
        for k, v in cj.items():
            kl = str(k).lower()
            vals = ([str(x) for x in v] if isinstance(v, list)
                    else [str(x) for x in v.keys()] if isinstance(v, dict)
                    else [str(v)])
            vals = [x for x in vals if len(x) <= 4 and x in seqs]
            if not vals:
                continue
            if "antigen" in kl or kl in ("ag", "ag_chains"):
                ag += vals
            elif "antibody" in kl or kl in ("ab", "ab_chains", "fab", "fv"):
                ab += vals
    ag, ab = sorted(set(ag)), sorted(set(ab))
    if not ag or not ab:
        # 표기가 없으면: 항체 가변부는 보통 110~250, 항원은 그 밖 — 가장 긴 것을 항원으로
        rest = sorted(seqs, key=lambda c: -len(seqs[c]))
        ag = ag or rest[:1]
        ab = ab or [c for c in sorted(seqs) if c not in ag]
    return ag, ab


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--src", default=os.path.expanduser(
        "~/projects/bk21-antibody-ml/consensus_docking/runs_rbd"))
    ap.add_argument("--targets-dir", default="targets")
    ap.add_argument("--data", default=os.environ.get("DATA", "/mnt/data/admuser/msadepth"))
    ap.add_argument("--list", default="sweep_targets.csv")
    ap.add_argument("--group", default="", help="sweep_targets.csv의 group 값(비우면 기존 파일에서 추정)")
    ap.add_argument("--apply", action="store_true")
    a = ap.parse_args()

    if not os.path.isdir(a.src):
        sys.exit(f"!! 원본 폴더 없음: {a.src}")

    rows = list(csv.DictReader(open(a.list))) if os.path.exists(a.list) else []
    fields = list(rows[0].keys()) if rows else ["target", "pdb", "group", "ab_chains", "AB", "ag_chains"]
    have = {r["target"] for r in rows}
    groups = sorted({r.get("group", "") for r in rows if r.get("group")})
    print(f"[{'실제 기록' if a.apply else 'dry-run — 아무것도 쓰지 않음'}]  원본 {a.src}")
    print(f"sweep_targets.csv 열: {fields}")
    print(f"기존 group 값들: {groups}\n")
    grp = a.group or (groups[0] if len(groups) == 1 else "")

    new = []
    for pdb in sorted(os.listdir(a.src)):
        d = os.path.join(a.src, pdb)
        cjp = os.path.join(d, "chains.json")
        if not os.path.isdir(d) or not os.path.exists(cjp):
            continue
        try:
            cj = json.load(open(cjp))
        except Exception as e:
            print(f"  !! {pdb}: chains.json 읽기 실패 {e}"); continue
        seqs = chain_seqs(cj)
        if not seqs:
            print(f"  !! {pdb}: 서열을 못 찾음 — chains.json 최상위 키 {list(cj)[:8] if isinstance(cj,dict) else type(cj)}")
            continue
        ag, ab = pick_roles(cj, seqs)
        name = f"{pdb}_{''.join(ab)}"
        site, abflag = SITE.get(pdb, ("?", "?"))

        # 항원 a3m 찾기
        a3ms = {}
        for c in ag:
            cands = (glob.glob(os.path.join(d, f"msa_{pdb}", f"{c}_*", "*.a3m"))
                     + glob.glob(os.path.join(d, f"msa_{pdb}", f"{c}_*", "**", "*.a3m"), recursive=True)
                     + glob.glob(os.path.join(d, f"msa_{pdb}", "**", "*.a3m"), recursive=True))
            best = None
            for p in cands:
                hs, ss = read_a3m(p)
                if ss and norm(ss[0]) == seqs[c]:
                    if best is None or len(ss) > best[1]:
                        best = (p, len(ss))
            a3ms[c] = best

        mark = "이미 등록됨" if name in have else "추가"
        print(f"■ {name:14} {mark:8} 붙는자리={site:14} 표식={abflag}")
        print(f"    항원 {ag} (길이 {[len(seqs[c]) for c in ag]}) · 항체 {ab} (길이 {[len(seqs[c]) for c in ab]})")
        for c in ag:
            if a3ms[c]:
                print(f"    a3m  {c}: {a3ms[c][1]:>6}개 서열  ←  {os.path.relpath(a3ms[c][0], d)}")
            else:
                print(f"    a3m  {c}: ❌ 질의서열과 맞는 a3m을 못 찾음 — 새로 만들어야 함")
        if name in have:
            print()
            continue

        if a.apply:
            td = os.path.join(a.targets_dir, name)
            os.makedirs(td, exist_ok=True)
            out = dict(cj)
            if isinstance(out, dict):
                out.setdefault("AB", abflag)
                out.setdefault("site_class", site)
                out.setdefault("source", "runs_rbd")
            json.dump(out, open(os.path.join(td, "chains.json"), "w"), ensure_ascii=False, indent=1)

            nat = os.path.join(td, "native.cif")
            if not os.path.exists(nat):
                url = f"https://files.rcsb.org/download/{pdb.upper()}.cif"
                try:
                    urllib.request.urlretrieve(url, nat)
                    print(f"    native.cif 받음 ({os.path.getsize(nat)//1024} KB)")
                except Exception as e:
                    print(f"    !! native.cif 실패: {e}")

            for c in ag:
                if not a3ms[c]:
                    continue
                dst = os.path.join(a.data, "ladders", name, c)
                os.makedirs(dst, exist_ok=True)
                r0 = os.path.join(dst, "rung0.a3m")
                if os.path.exists(r0):
                    print(f"    rung0.a3m 이미 있음 — 건드리지 않음")
                    continue
                hs, ss = read_a3m(a3ms[c][0])           # 읽는 순간 머리말 오염이 제거됨
                with open(r0, "w") as fh:
                    for h, s in zip(hs, ss):
                        fh.write(f">{h}\n{s}\n")
                print(f"    rung0.a3m 기록 ({len(ss)}개 서열)")

        row = {k: "" for k in fields}
        for k, v in (("target", name), ("pdb", pdb), ("group", grp),
                     ("ab_chains", "|".join(ab)), ("ag_chains", "|".join(ag)),
                     ("AB", abflag), ("site_class", site)):
            if k in row:
                row[k] = v
        new.append(row)
        print()

    print("=" * 78)
    if not new:
        print("추가할 타깃 없음 (전부 이미 등록되어 있거나 읽기 실패).")
        return
    print(f"새로 추가할 타깃 {len(new)}개: " + ", ".join(r["target"] for r in new))
    if not grp:
        print("⚠️ group 값을 못 정했다 — --group 으로 지정할 것 (인기자리 기준집합을 고르는 데 쓰임).")
    if a.apply:
        with open(a.list, "a", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=fields)
            for r in new:
                w.writerow(r)
        print(f"→ {a.list} 에 {len(new)}행 추가")
        print("\n다음: python neff_ladder.py  (사다리 생성)  →  bash run_sweep.sh")
    else:
        print("\n→ 실제로 등록하려면 --apply (필요하면 --group 도 함께).")


if __name__ == "__main__":
    main()
