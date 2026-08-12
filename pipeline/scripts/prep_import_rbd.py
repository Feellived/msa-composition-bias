#!/usr/bin/env python3
"""[가져오기] Consensus Docking 쪽 runs_rbd(결합부위 밖 항체 세트)를 깊이실험에 편입.

왜: 이 세트는 "인기 부위를 피해 붙는 항체"만 모은 것이라, 44개에서 8개밖에 못 찾던 조건(B군)을
10개가 전부 만족한다. 항원 MSA도 이미 있어(msa_<pdb>/A.a3m) 가장 비싼 단계를 건너뛴다.
8p5m만 A군(인기 부위 가장자리) = 반대방향 대조군. 8sdh는 분류 미기록이라 제외.

⚠️ chains.json을 직접 만들지 않는다. 이 저장소의 정규 경로인 `prep_targets.py`에 태운다
   (사슬을 항원=A·중쇄=B·경쇄=C로 개명하고, native.cif를 structures/로 심링크하며,
    RBD 그룹은 400잔기 초과 사슬만 319-541로 크롭 — 우리 항원은 195잔기라 크롭 없음).
   그래야 채점 코드가 기대하는 형식과 정확히 같아진다.

단계 (각각 기본 dry-run, 실제 기록은 --apply):
  --stage csv   ① runs_rbd에서 사슬 역할을 읽어 prep용 CSV 작성 + structures/<pdb>.cif 확보
     (여기서 사용자가:  python prep_targets.py --csv rbd_offhot.csv --struct structures --outdir targets)
  --stage msa   ② prep이 만든 chains.json의 개명된 항원 사슬에 기존 a3m을 붙여
                   $DATA/ladders/<타깃>/<개명ID>/rung0.a3m 작성 + sweep_targets.csv 행 추가
  --stage undo  ③ 2026-07-27 잘못된 형식으로 만들어졌던 산출물 정리(우리가 만든 것만)

사용:
  python prep_import_rbd.py --stage undo            # 먼저 정리(확인 후 --apply)
  python prep_import_rbd.py --stage csv --apply
  python prep_targets.py --csv rbd_offhot.csv --struct structures --outdir targets
  python prep_import_rbd.py --stage msa --apply
  python prep_a3m_check_match.py --only 8sis_HL ...  # 반드시 게이트로 확인
"""
import argparse, csv, glob, json, os, re, shutil, sys, urllib.request

# 노션 인수인계서 1.3절 "결합부위 밖 항체 10종" 표의 붙는자리 분류 → A/B 표식.
#   B = 진짜 결합자리가 인기 부위 **밖** (편향이 풀리면 좋아져야 함)
#   A = 진짜 결합자리가 인기 부위와 겹침 (편향이 풀리면 오히려 나빠져야 함)
SITE = {
    "9zdu": ("offhot:코어",            "B"),
    "9ml9": ("offhot:코어",            "B"),
    "8siq": ("offhot:코어",            "B"),
    "8sit": ("offhot:코어",            "B"),
    "8sis": ("offhot:코어와그밖",       "B"),
    "8xsi": ("offhot:숨은면",          "B"),
    "9ml8": ("offhot:숨은면",          "B"),
    "8sdf": ("offhot:숨은면",          "B"),
    "9sbb": ("offhot:그밖",            "B"),
    "8p5m": ("offhot:인기부위가장자리", "A"),   # ★ 반대방향 대조군
    # 8sdh = 노션 표에 없고 분류 미기록 → 2026-07-27 사용자 지시로 제외.
}
GROUP = "RBD"
AA = set("ACDEFGHIKLMNPQRSTVWYXBZUO")
PAT = re.compile(r"^#\d+\t\d+")
# 항체 J 영역 모티프 — 중쇄는 WGxG, 경쇄는 FGxG 로 끝난다. 길이보다 신뢰도가 높다.
HEAVY_MOTIF = re.compile(r"WG[QRKAE]G")
LIGHT_MOTIF = re.compile(r"FG[QGSTA]G")


def norm(s):
    return "".join(c for c in s.upper() if c not in "-.")


def read_a3m(path):
    """a3m을 (헤더, 서열) 목록으로. 질의행 머리말 오염(2026-07-27 버그)도 제거."""
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
    """chains.json에서 {사슬ID: 서열} 추출(형식이 여러 가지라 방어적으로)."""
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
    """(항원 사슬들, 중쇄, 경쇄). 중/경쇄는 J 영역 모티프로 가르고, 없으면 길이로."""
    ag, ab = [], []
    if isinstance(cj, dict):
        for k, v in cj.items():
            kl = str(k).lower()
            vals = ([str(x) for x in v] if isinstance(v, list)
                    else [str(x) for x in v.keys()] if isinstance(v, dict) else [str(v)])
            vals = [x for x in vals if len(x) <= 4 and x in seqs]
            if not vals:
                continue
            if "antigen" in kl or kl in ("ag", "ag_chains"):
                ag += vals
            elif "antibody" in kl or kl in ("ab", "ab_chains", "fab", "fv"):
                ab += vals
    ag, ab = sorted(set(ag)), sorted(set(ab))
    if not ag or not ab:
        rest = sorted(seqs, key=lambda c: -len(seqs[c]))
        ag = ag or rest[:1]
        ab = ab or [c for c in sorted(seqs) if c not in ag]

    heavy = light = None
    scored = []
    for c in ab:
        s = seqs[c]
        scored.append((c, bool(HEAVY_MOTIF.search(s)), bool(LIGHT_MOTIF.search(s)), len(s)))
    hs = [c for c, h, l, _ in scored if h and not l]
    ls = [c for c, h, l, _ in scored if l and not h]
    if len(hs) == 1 and len(ls) == 1:
        heavy, light = hs[0], ls[0]
    elif len(ab) == 2:                       # 모티프가 안 갈리면 긴 쪽을 중쇄로
        a, b = sorted(ab, key=lambda c: -len(seqs[c]))
        heavy, light = a, b
    elif len(ab) == 1:
        heavy, light = ab[0], ""
    return ag, (heavy or ""), (light or ""), scored


def src_targets(src):
    """runs_rbd 훑어 (pdb, 타깃이름, 항원들, H, L, 서열, 판정근거) 목록."""
    out = []
    for pdb in sorted(os.listdir(src)):
        d = os.path.join(src, pdb)
        cjp = os.path.join(d, "chains.json")
        if not os.path.isdir(d) or not os.path.exists(cjp):
            continue
        if pdb not in SITE:
            print(f"■ {pdb:14} 건너뜀 — 붙는자리 분류 없음(2026-07-27 제외 결정)")
            continue
        cj = json.load(open(cjp))
        seqs = chain_seqs(cj)
        ag, H, L, scored = pick_roles(cj, seqs)
        out.append(dict(pdb=pdb, dir=d, name=f"{pdb}_{H}{L}", ag=ag, H=H, L=L,
                        seqs=seqs, scored=scored))
    return out


def find_a3m(d, pdb, want_seq):
    """질의서열이 want_seq와 같은 a3m 중 서열 수가 가장 많은 것."""
    cands = sorted(set(glob.glob(os.path.join(d, f"msa_{pdb}", "**", "*.a3m"), recursive=True)
                       + glob.glob(os.path.join(d, "**", "*.a3m"), recursive=True)))
    best = None
    for p in cands:
        try:
            hs, ss = read_a3m(p)
        except Exception:
            continue
        if ss and norm(ss[0]) == want_seq and (best is None or len(ss) > best[1]):
            best = (p, len(ss))
    return best


# ── 단계 ① prep용 CSV + structures ──────────────────────────────────────
def stage_csv(a, tg):
    need = ["pdb", "Hchain", "Lchain", "antigen_chain", "antigen", "AB", "label"]
    print(f"\n[① prep용 CSV] → {a.csv_out} · 구조 → {a.struct}/\n")
    print(f"{'pdb':6}{'타깃':12}{'항원':8}{'중쇄':6}{'경쇄':6}{'표식':5}  중/경쇄 판정근거")
    print("-" * 92)
    rows = []
    for t in tg:
        site, ab = SITE[t["pdb"]]
        why = " ".join(f"{c}({'W' if h else '-'}{'F' if l else '-'}{n})" for c, h, l, n in t["scored"])
        print(f"{t['pdb']:6}{t['name']:12}{'|'.join(t['ag']):8}{t['H']:6}{t['L']:6}{ab:5}  {why}")
        rows.append({"pdb": t["pdb"], "Hchain": t["H"], "Lchain": t["L"],
                     "antigen_chain": "|".join(t["ag"]), "antigen": GROUP,
                     "AB": ab, "label": site})
    print("\n판정근거 표기: 사슬(W=중쇄모티프 WGxG · F=경쇄모티프 FGxG · 길이)")

    if not a.apply:
        print(f"\n→ 실제로 쓰려면 --apply. 그다음:\n"
              f"   python prep_targets.py --csv {a.csv_out} --struct {a.struct} --outdir {a.targets_dir}")
        return
    with open(a.csv_out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=need); w.writeheader(); w.writerows(rows)
    print(f"\n→ {a.csv_out} ({len(rows)}행)")

    os.makedirs(a.struct, exist_ok=True)
    for t in tg:
        dst = os.path.join(a.struct, f"{t['pdb']}.cif")
        if os.path.exists(dst):
            print(f"   {t['pdb']}.cif 이미 있음"); continue
        # 앞선 실행이 targets/에 받아둔 원본이 있으면 재사용(중복 내려받기 방지)
        prev = os.path.join(a.targets_dir, t["name"], "native.cif")
        if os.path.isfile(prev) and not os.path.islink(prev):
            shutil.copy2(prev, dst)
            print(f"   {t['pdb']}.cif ← 앞선 실행분 재사용 ({os.path.getsize(dst)//1024} KB)")
            continue
        try:
            urllib.request.urlretrieve(f"https://files.rcsb.org/download/{t['pdb'].upper()}.cif", dst)
            print(f"   {t['pdb']}.cif 받음 ({os.path.getsize(dst)//1024} KB)")
        except Exception as e:
            print(f"   !! {t['pdb']}.cif 실패: {e}")
    print(f"\n다음: python prep_targets.py --csv {a.csv_out} --struct {a.struct} --outdir {a.targets_dir}")


# ── 단계 ② a3m을 개명된 사슬에 붙이고 목록에 등록 ────────────────────────
def stage_msa(a, tg):
    print(f"\n[② a3m 연결 + 목록 등록]  사다리 → {a.data}/ladders/<타깃>/<개명ID>/rung0.a3m\n")
    rows = list(csv.DictReader(open(a.list))) if os.path.exists(a.list) else []
    fields = list(rows[0].keys()) if rows else ["target", "pdb", "group", "ab",
                                                "dirtype", "ag_chains", "label"]
    have = {r["target"] for r in rows}
    dts = sorted({r.get("dirtype", "") for r in rows if r.get("dirtype")})
    dtype = dts[0] if len(dts) == 1 else "targets"
    new, plans = [], []

    for t in tg:
        name = t["name"]
        cjp = os.path.join(a.targets_dir, name, "chains.json")
        if not os.path.exists(cjp):
            print(f"■ {name:12} ❌ targets/{name}/chains.json 없음 — prep_targets.py 를 먼저 돌릴 것")
            continue
        cj = json.load(open(cjp))
        if "chains" not in cj:
            print(f"■ {name:12} ❌ chains.json이 prep 형식이 아님(키: {sorted(cj)[:6]}) — prep 재실행 필요")
            continue
        agc = [c for c in cj["chains"] if c.get("role") == "antigen"]
        site, ab = SITE[t["pdb"]]
        print(f"■ {name:12} 표식={ab}  항원 사슬 {[(c['id'], c.get('src'), len(c['seq']), c.get('crop')) for c in agc]}")
        okall = True
        for c in agc:
            want = norm(c["seq"])
            best = find_a3m(t["dir"], t["pdb"], want)
            if not best:
                print(f"    {c['id']}(원본 {c.get('src')}): ❌ 질의서열이 일치하는 a3m 없음"
                      f" — 크롭됐거나 서열이 달라졌을 수 있음(길이 {len(want)})")
                okall = False
                continue
            dst = os.path.join(a.data, "ladders", name, c["id"], "rung0.a3m")
            print(f"    {c['id']}(원본 {c.get('src')}): {best[1]:>5}개 서열 ← {os.path.relpath(best[0], t['dir'])}"
                  + ("  [이미 있음 — 건드리지 않음]" if os.path.exists(dst) else ""))
            plans.append((dst, best[0], want))
        if not okall:
            continue
        if name in have:
            print(f"    목록에 이미 등록됨")
            continue
        row = {k: "" for k in fields}
        for k, v in (("target", name), ("pdb", t["pdb"]), ("group", GROUP), ("ab", ab),
                     ("dirtype", dtype), ("ag_chains", "|".join(c["id"] for c in agc)),
                     ("label", site)):
            if k in row:
                row[k] = v
        miss = [k for k in fields if not row[k]]
        if miss:
            print(f"    ⚠️ 값을 못 채운 열: {miss}")
        new.append(row)

    print("\n" + "=" * 78)
    print(f"a3m 연결 {len(plans)}건 · 목록 추가 {len(new)}건")
    if not a.apply:
        print("\n→ 실제로 쓰려면 --apply.")
        return
    for dst, srcp, want in plans:
        if os.path.exists(dst):
            continue
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        hs, ss = read_a3m(srcp)              # 읽는 순간 머리말 오염 제거
        assert norm(ss[0]) == want, "질의서열 불일치 — 기록 중단"
        with open(dst, "w") as fh:
            for h, s in zip(hs, ss):
                fh.write(f">{h}\n{s}\n")
        print(f"   {dst} ({len(ss)}개 서열)")
    if new:
        with open(a.list, "a", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=fields)
            for r in new:
                w.writerow(r)
        print(f"   {a.list} 에 {len(new)}행 추가")
    print("\n다음: python prep_a3m_check_match.py   (반드시 게이트로 확인)  →  python prep_ladder_neff.py")


# ── 단계 ③ 잘못된 형식으로 만들어졌던 산출물 정리 ───────────────────────
def stage_undo(a, tg):
    """2026-07-27 첫 시도가 남긴 것만 지운다. 우리가 만든 이름 외에는 절대 손대지 않는다."""
    names = {t["name"] for t in tg}
    print(f"\n[③ 정리]  대상 타깃 {len(names)}개: {', '.join(sorted(names))}\n")
    hits = []
    for n in sorted(names):
        td = os.path.join(a.targets_dir, n)
        cjp = os.path.join(td, "chains.json")
        if os.path.isdir(td):
            bad = ""
            if os.path.exists(cjp):
                try:
                    bad = "" if "chains" in json.load(open(cjp)) else "  ← 형식 잘못됨"
                except Exception:
                    bad = "  ← 읽기 실패"
            hits.append(("폴더", td, bad))
        ld = os.path.join(a.data, "ladders", n)
        if os.path.isdir(ld):
            subs = sorted(os.listdir(ld))
            hits.append(("사다리", ld, f"  하위 {subs}"))
    rows = list(csv.DictReader(open(a.list))) if os.path.exists(a.list) else []
    fields = list(rows[0].keys()) if rows else []
    drop = [r for r in rows if r["target"] in names]
    for kind, p, note in hits:
        print(f"  [{kind}] {p}{note}")
    print(f"  [목록] {a.list} 에서 {len(drop)}행 제거: {[r['target'] for r in drop]}")
    if not hits and not drop:
        print("  정리할 것 없음.")
        return
    if not a.apply:
        print("\n→ 실제로 지우려면 --apply. (여기 나열된 경로 외에는 건드리지 않음)")
        return
    for kind, p, _ in hits:
        shutil.rmtree(p)
        print(f"  삭제 {p}")
    if drop:
        keep = [r for r in rows if r["target"] not in names]
        with open(a.list, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=fields); w.writeheader(); w.writerows(keep)
        print(f"  {a.list} → {len(keep)}행 (10행 제거)")
    print("\n정리 완료. 다음: python prep_import_rbd.py --stage csv --apply")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--src", default=os.path.expanduser(
        "~/projects/bk21-antibody-ml/pipeline/runs_rbd"))
    ap.add_argument("--stage", choices=["plan", "csv", "msa", "undo"], default="plan")
    ap.add_argument("--targets-dir", default="targets")
    ap.add_argument("--struct", default="structures")
    ap.add_argument("--csv-out", default="rbd_offhot.csv")
    ap.add_argument("--data", default=os.environ.get("DATA", "/mnt/data/admuser/msadepth"))
    ap.add_argument("--list", default="sweep_targets.csv")
    ap.add_argument("--apply", action="store_true")
    a = ap.parse_args()

    if not os.path.isdir(a.src):
        sys.exit(f"!! 원본 폴더 없음: {a.src}")
    print(("[실제 기록]" if a.apply else "[dry-run — 아무것도 쓰지 않음]")
          + f"  단계={a.stage}  원본={a.src}\n")
    tg = src_targets(a.src)
    print(f"\n대상 {len(tg)}개 (B군 {sum(1 for t in tg if SITE[t['pdb']][1]=='B')}"
          f" · A군 {sum(1 for t in tg if SITE[t['pdb']][1]=='A')})")

    if a.stage == "plan":
        print("\n단계를 골라 실행:  --stage undo → --stage csv → (prep_targets.py) → --stage msa")
    elif a.stage == "csv":
        stage_csv(a, tg)
    elif a.stage == "msa":
        stage_msa(a, tg)
    elif a.stage == "undo":
        stage_undo(a, tg)


if __name__ == "__main__":
    main()
