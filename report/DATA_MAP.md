# 원자료 지도 — MSA 조성 실험 (2026-07-30 확정)

발표·보고서·그림 작업에서 **어떤 주장이 어느 파일에서 나오는지**를 한곳에 모은 문서.
`발표흐름_MSA조성_20260731.md`의 절 번호를 기준으로 한다.

---

## 0. 저장 위치 세 곳

| 위치 | 내용 | 접근 |
|---|---|---|
| `~/projects/msa-composition-bias/pipeline/results/` | 본 검정·사다리 채점 결과 (**CSV 196 · json 30 · txt 59 · `analysis/` · `.bak` 2**) | 대부분 레포 커밋됨 |
| `~/projects/epitope-guided-docking/pipeline/` | 유도 재도킹(데모) · 물리 도킹 · 모델 비교 | `results/*.csv` · `dataset/*.csv` 커밋됨 |
| `/mnt/data/msadepth/` | 예측 원본(구조·MSA·pose) · iDist | **레포 밖.** 용량 큼 |

⚠️ `_stale_0727_1505` · `_stale_0727_1516` 계열 폴더는 **a3m 사고로 밀어둔 무효 자료**다. 채점기가
긁지 않는지 확인할 것.

---

## 1. 절별 원자료 지도

| 절 | 주장 | 파일 | 열·규모 |
|---|---|---|---|
| 1.2 | 타깃 내 ipTM–DockQ 상관 중앙값 +0.04 | `results/pose_features.csv` | 5,761행 · `iptm,ptm,plddt` |
| 2.1 | 5모델 × 5복합체 DockQ | `pipeline/results/dockq_cofolder.csv` · `summary.csv` | — |
| 2.1 | 항원별 공개 구조 수 (9K6J 2,212 / 776) | RCSB 조회 (문서 표에 기록) | 저장 파일 없음 |
| **2.2** | A/B/C recall 0.72 / 0.42 / 0.30 | `results/epitope_recall.csv` · `report/epitope_recall.csv` | 505행 · 47복합체 × 칸 × `best/mean/min_recall·prec·f1·mcc·auprc` |
| **2.2** | A/B/C DockQ 성공 9/20 · 3/18 · 1/9 | `results/analysis/success_rates.csv` | 9행 · 모델 × scope × `full_23/49/80` · `orc_23/49/80` |
| **5.1** | 깊이 이득 · 53종 | `results/dockq_sweep.csv` | 1,150행 · 59타깃 × 모델 × rung × `neff80,best_dockq` |
| **5.1** | 평탄 31/47 · 37/44 | `results/analysis/shape_counts.csv` | 10행 |
| 5.1 | medium 도달 2종(8ulr·9azr) | `results/analysis/rescue_candidates.csv` | 10행 · `shape,kind,best_rung,best_neff` |
| 5.1 | 타깃별 이득 요약 | `results/screen_candidates.csv` | 45행 · `tier,full,peak,gain,peak_rung,peak_neff` |
| **5.2.1** | 8ulr_HL 조성별 성공 · 흔들림 | `results/compreps_8ulr_HL.csv` | **311행 = 62실행 × 5자세** |
| **5.2.2** | 8k3k_D 18/24 대 0/20 | `results/compreps_8k3k_D.csv` | **221행 = 44실행 × 5자세** |
| 5.2 | 성공/실패 접촉면·precision | `results/site_repro_<T>.csv` | 타깃별 · `within,between,ratio,perm_p,n_res,true_covered,precision` |
| **5.3 · 부록 A** | 세 기준 29종 | `results/criteria.csv` | 30행 · `het_p,perm_p,n_cand,best_cover,full_cover,dq_succ` |
| **5.4 · 부록 B** | 방향 5 대 5 | `results/anchor_tests_recall.csv` · `_dockq.csv` | 4행 · `comp_rate,full_rate,direction,per_group` |
| **5.5** | 후보별 잔기·수렴 조성 | `results/sites_<T>.json` (30개) | `candidates[{cand,n_comp,comps,residues}]` · `perm_p,within,between` |
| 5.5 | 실행별 자카드 | `results/epitope_cluster_<T>.csv` (29개) | `jac_to_true,jac_to_popular,jac_to_succ/fail_consensus` |
| **5.7 · 6.2** | 데모 네 조건 최고·중앙 | `pipeline/results/demo_dockq.csv` | `arm,n_pose,dockq_max,dockq_med,recall_max,recall_med` |
| **6.1** | 928회 중 DockQ ≥ 0.49는 4회 | `results/maintest_poses.csv` | 349KB · 자세 단위 전수 |
| **6.2** | ipTM 선택 10/47 대 천장 23/47 | `results/pose_features.csv` + `analysis/success_rates.csv` | 두 파일이 서로 검증 |
| 6.4 | Protenix full acceptable 1/44 | `analysis/success_rates.csv` | `protenix,ALL,44,1,0,0,6,2,0` |
| 4.1 | 명단 59 · 판독 54 · run 30 / no_response 24 | `pipeline/maintest.csv`(55행) · `sweep_targets.csv`(60행) · `set3_targets.csv` | `status,neff_pick,stratum,recall_max,dockq_max` |
| 4.1 | 명단의 원본 956 후보 | `pipeline/dataset/manifest_labeled.csv` | 956행 · `antigen,epitope_class,n_epi,frac,epitope_residues` |
| 4.4 | 칸별 Neff80 실측 | `/mnt/data/.../ladders/<T>/<ch>/neff.tsv` | 레포 밖 |

---

## 2. CSV 계열 전수 (results/, 헤더 기준)

| 개수 | 행 합 | 헤더 | 무엇인가 |
|---|---|---|---|
| 83 | 166 | `target,n_full,n_red,med_full,med_red,succ49_*,succ23_*,p_fisher*,p_ranktest,p_heterogeneity,verdict` | `summary_<T>_{dockq,recall,overrep}.csv` — 타깃×지표 1행 |
| 31 | 9,891 | `target,model,depth,seed,pose,dockq,recall,overrep,n_contact` | **`compreps_<T>.csv` — 자세 단위.** `seed` = `seed<조성>_r<반복>`, `seedfull` = 원래 MSA |
| 30 | 104 | `target,model,depth,within,between,ratio,perm_p,n_cand,cand,n_comp,n_res,true_covered,precision,comps` | `site_repro_<T>.csv` — 자리(후보) 단위 |
| 29 | 957 | `run,dockq,recall,n_pred,jac_to_true,jac_to_fail_consensus,jac_to_succ_consensus,jac_to_popular,success` | `epitope_cluster_<T>.csv` — 실행 단위 + 자카드 |
| 3 | 1,690 | `target,group,ab,model,rung,neff80,best_dockq,n_pose` | `dockq_sweep{,_boltz,_protenix}.csv` |
| 3 | 178 | `target,model,depth,seed,pose,dockq,recall` | `seedrep_poses.csv` 등 **구 자료(boltz, a3m 사고분)** |
| 2 | 4 | `...,comp_rate,full_rate,direction,per_group` | `anchor_tests_{dockq,recall}.csv` |
| 1 | 5,761 | `target,group,ab,label,model,rung,neff80,pose,dockq,recall,n_contact,overrep,true_rank,pop_rank,dcc_true,dcc_pop,iptm,ptm,plddt` | **`pose_features.csv` — 유일하게 ipTM 포함** |
| 1 | 529 | `...,native_overrep,mean_overrep,oracle_overrep,mean_recall,oracle_recall,mean_true_rank,...` | `epitope_shift_protenix.csv` |
| 1 | 505 | `...,n_true,n_ag,n_pose,best_recall,mean_recall,min_recall,best_prec,best_f1,best_mcc,best_auprc` | `epitope_recall.csv` — **2.2의 원본** |
| 1 | 45 | `target,tier,full,peak,gain,peak_rung,peak_neff,n_rungs_better,n_rungs_succ` | `screen_candidates.csv` |
| 1 | 41 | `target,family,ab,n_ref,min_dist,nearest,mean_knn,frac_ndup,n_0.02~0.06` | `overrep_idist.csv` — iDist 중복도 |
| 1 | 36 | `target,grp,site,n_rung,full_recall,red_recall,d_recall,full_overrep,red_overrep,d_overrep` | `shift_direction.csv` — 자리 단위 이동 |
| 1 | 30 | `target,group,stratum,het_p,perm_p,n_cand,best_cover,full_cover,dq_succ,dq_het,scored` | **`criteria.csv` — 부록 A** |
| 1 | 17 | `...,dq_full,dq_peak,dq_gain,rec_*,over_*,coincide,responsive` | `crosscheck_depth.csv` |
| 1 | 17 | `target,grp,site,antigen,full_recall,full_overrep,full_dockq,best_rung,...,shift,n_rung_ok` | 자리 단위 이동(다른 판) |
| 1 | 83 | 위 83개와 같은 헤더 | `maintest_summary.csv` — 병합판 |
| 1 | 5 / 4 / 2 | `dq_mean,dq_sd…` / `ctl_*,seedrep_*,p_one_sided` / `len_a3m_raw,len_a3m_clean,verdict` | `seedrep_cand_scored` · `fullmsa_ctl_scored` · a3m 점검 |

**`.bak` 2개는 구버전**(`pose_features` 5,161 < 5,761 · `dockq_sweep` 1,033 < 1,150). 현재 파일이 상위집합.

---

## 3. `results/analysis/` — 2026-07-24 에피토프 이동 분석

| 파일 | 내용 |
|---|---|
| **`success_rates.csv`** | 모델 × {ALL,A,B,C} × `full_23/49/80` · `orc_23/49/80` — **2.2·6.2·6.4의 근거** |
| **`shape_counts.csv`** | 깊이 곡선 모양 분포 (flat / spiky / mid-peak / deep-better / shallow-better) |
| `rescue_candidates.csv` | 10행 · `model_rescue` 6 + `strong_rescue` 4 · `shape,best_rung,best_neff` |
| `family_summary.csv` | 계열 × ab × `mean_full,mean_best,delta_lowfull,mean_spearman,spearman_pos` |
| `epitope_cross_model.csv` | `biased_boltz,biased_protenix,both_biased` — 편향의 모델 특이성 |
| `epitope_family_summary.csv` · `epitope_biased_filter.csv` · `epitope_biased_trajectory.csv` · `per_target_summary.csv` | 보조 집계 |
| `epitope_summary.txt` · `summary.txt` | 사람이 읽는 요약 |
| **`fig_depth_sensitivity.png` · `fig_rescue_candidates.png`** | **이미 그려 둔 그림 2장** |

---

## 4. epitope-guided-docking/pipeline

| 파일 | 내용 |
|---|---|
| **`results/demo_dockq.csv`** | 데모 8종 × 다섯 조건 · `dockq_max/med · recall_max/med` (blind는 protenix 40자세) |
| `results/dockq_cofolder.csv` · `summary.csv` | 모델 비교(2.1) |
| `results/zdock_dockq.csv` · `zdock_poses_<pdb>.csv` ×10 | 물리 도킹(강체) |
| `results/haddock_dockq_patch.csv` · `haddock_score_rank.csv`(+cache) | 물리 도킹(정제) |
| `results/guided_*_{zdock,haddock}*.csv` ×8 · `dockq_{boltz,protenix,tfold}_guided.csv` | 물리 힌트 유도 |
| **`dataset/manifest_labeled.csv`** | **956행 — 명단의 원본.** `epitope_class,n_epi,frac,epitope_residues` |
| `figures/` | `heatmap_dockq_fab_fv.png` · `boxplot_dockq_fab_vs_fv.png` · `cases_5model/` · `cases_diverse/` |
| `results/diagnose_*.txt` ×10 · `where_docked_*.txt` ×11 | **미추적.** Phase 0 진단 로그 |
| `runs_*/` 폴더 13개 · `haddock*/` · `zdock/` | **미추적.** 구조 원본, 용량 큼 |

---

## 5. 그림 후보와 원자료

| 그림 | 원자료 | 무엇이 보이나 |
|---|---|---|
| **H1** 타깃 29 × 조성 성공률 히트맵 | `compreps_<T>.csv` 30개 | 조성이 자리를 정한다 — 한 장 |
| **H2** 타깃 47 × 깊이칸 recall 히트맵 | `epitope_recall.csv` | 5.1의 평탄함을 격자로 |
| H3 데모 8종 × 조건 히트맵 | `demo_dockq.csv` | 1종만 통과 |
| **S1** 928 자세 DockQ 대 recall | `maintest_poses.csv` | 위치 축과 자세 축의 분리 |
| S2 조성 대 원래 MSA 성공률 | `anchor_tests_recall.csv` | 대각선 위/아래 5 대 5 |
| S3 원래 MSA 덮음 대 최고 후보 덮음 | `criteria.csv` | 후보 생성 4종 |
| **S4** ipTM 대 DockQ (타깃별 회귀선 54개) | `pose_features.csv` | **부호가 절반씩 갈림 = 선택 불가** |
| R1 Neff80 대 best_dockq 회귀 | `dockq_sweep.csv` | 문헌 주장을 우리 데이터로 검정 |
| R2 후보 수 대 최고 후보 덮음 | `criteria.csv` | 5.5 각주의 정량화 |

---

## 6. 취급 주의

- **`eval_dump_seedrep.py`는 `--data $DATA/compreps` 필수.** 빼면 옛 `seedrep_cand`(무효 boltz 3건)를
  읽고 **CSV를 덮어쓴다**(2026-07-30 실제 사고). 예측 CIF는 무사하므로 재채점으로 복구된다.
- **DockQ가 없는 env에서 돌리면** 파일을 쓰지 않고 멈춘다(정상 동작). `conda activate boltz`.
- **Boltz 사다리 자료는 a3m 사고분**이다. 깊이 축 해석에 쓸 수 없다. 단 한 실행 안의 pose–confidence
  관계나 rung0 성적은 유효하다.
- `criteria.csv`의 `stratum`(층) 열은 **주장에 쓰지 않기로 합의**되었다.
- 데모는 `MODELS=boltz`로 제한한다. Protenix는 제약이 반영되지 않는다.
- 원격 MSA 서버(ColabFold MMseqs2)는 큐 대기가 길다. `PENDING` 반복은 정상.
