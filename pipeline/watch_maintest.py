#!/usr/bin/env python3
"""[진행 확인] 본 검정이 잘 돌고 있는지 본다. 읽기만 하고 아무것도 바꾸지 않는다.

GPU를 쓰지 않으며, 돌고 있는 작업의 파일을 건드리지 않는다(열어 읽기만 한다).
따라서 본 검정과 동시에 아무 때나 돌려도 안전하다.

무엇을 보여주나
  ① 타깃별 진행     끝난 실행 / 목표 실행, 지금 몇 번째가 도는 중인지
  ② 전체 진행과 예상 남은 시간   최근 완료 간격의 중앙값으로 추정(초반엔 부정확)
  ③ 조용한 실패 감시  run.log의 MSA 경고 — 모델이 MSA를 버리고 단일서열로 도는 경우.
                      이건 끝까지 가도 오류가 안 나므로 도중에 봐야 한다.
  ④ 자세 수 이상     한 실행에 5개가 정상. 0개로 끝난 실행은 실패다.

사용 (DockQ가 있는 환경이면 --score 도 가능)
  python watch_maintest.py                 # 한 번 출력
  python watch_maintest.py --loop          # 60초마다 갱신 (Ctrl+C 종료)
  python watch_maintest.py --loop --every 300
  python watch_maintest.py --score 8t4d_OQ # 지금까지 나온 것만 채점(CPU, 몇 분)

--score 는 dump_seedrep_full.py 를 그대로 부른다. 아직 안 끝난 타깃도 그 시점까지의
자세로 채점하므로 "지금 어떤 값이 나오고 있나"를 미리 볼 수 있다. 다만 이건 중간 확인일
뿐이고, 판정은 전부 끝난 뒤 6.4절 기준으로만 한다.
"""
import argparse, csv, glob, os, re, statistics, subprocess, sys, time

WARN_PAT = re.compile(
    r"does not match input sequence|creating dummy|msa .*mismatch|query.*size mismatch",
    re.I)


def rows(csv_path):
    with open(csv_path) as fh:
        return [r for r in csv.DictReader(fh) if r.get("status") == "run"]


def scan(base):
    """실행 폴더별 (자세 수, 마지막 수정 시각)."""
    out = []
    for d in sorted(glob.glob(os.path.join(base, "d*", "seed*_r*"))):
        cifs = glob.glob(os.path.join(d, "results", "**", "*sample*.cif"), recursive=True)
        mt = max((os.path.getmtime(c) for c in cifs), default=os.path.getmtime(d))
        out.append((d, len(cifs), mt))
    return out


def hms(sec):
    if sec is None or sec != sec or sec < 0:
        return "?"
    sec = int(sec)
    return f"{sec//3600}시간 {sec%3600//60}분" if sec >= 3600 else f"{sec//60}분"


def bar(done, want, w=22):
    n = 0 if want <= 0 else min(w, int(w * done / want))
    return "█" * n + "·" * (w - n)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", default="maintest.csv")
    ap.add_argument("--data", default=os.environ.get("DATA", "/mnt/data/admuser/msadepth"))
    ap.add_argument("--loop", action="store_true")
    ap.add_argument("--every", type=int, default=60)
    ap.add_argument("--score", default="", help="이 타깃을 지금까지 나온 자세로 채점")
    a = ap.parse_args()

    if a.score:
        print(f"[채점] {a.score} — 지금까지 나온 자세만. 아직 진행 중이면 값이 더 늘어난다.\n")
        cmd = ["python", "dump_seedrep_full.py", "--data", os.path.join(a.data, "compreps"),
               "--only", a.score, "--csv-out", f"results/compreps_{a.score}.csv"]
        raise SystemExit(subprocess.call(cmd))

    if not os.path.exists(a.csv):
        raise SystemExit(f"!! {a.csv} 없음. pipeline 폴더에서 실행할 것.")

    while True:
        rs = rows(a.csv)
        print("\033[2J\033[H" if a.loop else "", end="")
        print(f"[{time.strftime('%m-%d %H:%M:%S')}] 본 검정 진행 상황  ({a.csv})\n")
        print(f"{'타깃':13}{'군':5}{'진행':24}{'실행':>10}  {'자세':>6}  상태")
        print("-" * 78)

        tot_done = tot_want = 0
        finish_times, warns, empties, running = [], [], [], []
        for r in rs:
            t = r["target"]; model = r.get("model") or "protenix"
            want = int(r.get("n_comp") or 6) * int(r.get("n_reps") or 4) + int(r.get("n_full") or 8)
            base = os.path.join(a.data, "compreps", "seedrep_cand", model, t)
            runs = scan(base) if os.path.isdir(base) else []
            full = [x for x in runs if x[1] >= 5]
            part = [x for x in runs if 1 <= x[1] < 5]
            zero = [x for x in runs if x[1] == 0]
            npose = sum(x[1] for x in runs)
            tot_done += len(full); tot_want += want
            finish_times += [x[2] for x in full]

            last = max((x[2] for x in runs), default=0)
            fresh = last and (time.time() - last) < 1800   # 30분 안에 뭔가 쓰였나
            if len(full) >= want:
                st = "완료"
            elif not runs:
                st = "시작 전"
            elif part or fresh:
                st = "진행 중"; running.append(t)
            else:
                st = f"멈춤/대기 ({hms(time.time()-last)} 전)"
            if zero:
                empties.append((t, len(zero)))
                st += f" · 빈 실행 {len(zero)}개"

            for d, n, _ in runs:
                lg = os.path.join(d, "run.log")
                if os.path.exists(lg):
                    try:
                        with open(lg, errors="replace") as fh:
                            if WARN_PAT.search(fh.read()):
                                warns.append(d)
                    except OSError:
                        pass

            print(f"{t:13}{r.get('group',''):5}{bar(len(full), want):24}"
                  f"{len(full):>4}/{want:<5}{npose:>6}  {st}")

        pct = 100 * tot_done / tot_want if tot_want else 0
        print("-" * 78)
        print(f"{'전체':13}{'':5}{bar(tot_done, tot_want):24}{tot_done:>4}/{tot_want:<5}"
              f"{'':6}  {pct:.1f}%")

        # 최근 완료 간격으로 남은 시간 추정
        finish_times.sort()
        recent = finish_times[-15:]
        gaps = [b - a_ for a_, b in zip(recent, recent[1:]) if 0 < b - a_ < 7200]
        if len(gaps) >= 3:
            med = statistics.median(gaps)
            left = tot_want - tot_done
            print(f"      최근 실행 1회당 중앙값 {med/60:.1f}분 · 남은 {left}회 → 약 {hms(med*left)}")
            since = time.time() - finish_times[-1]
            if since > 3 * med and since > 1800:
                print(f"      ⚠️ 마지막 완료가 {hms(since)} 전이다 — 멈췄는지 확인할 것"
                      f" (tmux attach -t maintest · nvidia-smi)")
        else:
            print("      (완료가 아직 적어 남은 시간 추정 불가)")

        if warns:
            print(f"\n⚠️⚠️ MSA 경고가 있는 실행 {len(warns)}개 — 모델이 MSA를 버리고 단일서열로 "
                  f"돌았을 수 있다. 그러면 그 실행은 무효다:")
            for d in warns[:5]:
                print("   " + d)
            if len(warns) > 5:
                print(f"   … 외 {len(warns)-5}개")
            print("   → 즉시 확인: grep -i 'does not match\\|dummy' <위 폴더>/run.log")
        if empties:
            print(f"\n⚠️ 자세가 0개인 실행이 있는 타깃: "
                  + ", ".join(f"{t}({n})" for t, n in empties[:8]))
            print("   → 해당 폴더의 run.log 확인. 일시적이면 다시 돌리면 채워진다.")
        if not warns and not empties:
            print("\n✅ MSA 경고 없음 · 빈 실행 없음")
        if running:
            print(f"진행 중인 타깃: {', '.join(running)}")

        print("\n채점을 미리 보려면:  python watch_maintest.py --score <타깃>   (DockQ 환경)")
        if not a.loop:
            return
        time.sleep(a.every)


if __name__ == "__main__":
    main()
