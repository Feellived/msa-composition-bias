#!/usr/bin/env bash
# 원천 데이터 다운로드 — SAbDab2 전체 summary (항체-항원 구조 메타데이터 = 매니페스트 구축의 출발점).
# ⚠️ SAbDab은 live DB — 재다운로드하면 세트가 달라질 수 있다. 확정 세트(pilot_lean_full.csv)는
#    committed manifest를 그대로 쓰고, 이 스크립트는 "어떻게 가져왔나"(provenance) 재현·갱신용.
# 포맷: SAbDab2 = 콤마 CSV(~44컬럼), GET만 허용. 구조 날짜는 'date' 컬럼(SABDABdepo_date=DB재구축일이라 부적합).
# 사용: bash fetch_sabdab.sh [출력경로]
set -uo pipefail
OUT="${1:-sabdab2_summary.csv}"
URL="https://sabdab.opig.stats.ox.ac.uk/api/download/all-summary"
echo "[sabdab] $URL"
curl -sL "$URL" -o "$OUT" || { echo "!! 다운로드 실패"; exit 1; }
n=$(wc -l < "$OUT")
echo "[done] $OUT — $n 행"
echo "  헤더: $(head -1 "$OUT" | cut -c1-100)..."
[ "$n" -lt 1000 ] && echo "  ⚠️ 행 수가 적음 — SPA 셸(HTML)을 받았을 수 있음. URL/포맷 확인."
echo "  다음: build_manifest.py(세트 재구축·live라 세트 변동) 또는 committed manifest 그대로 사용."
