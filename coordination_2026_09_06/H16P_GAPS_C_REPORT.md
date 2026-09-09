# H16P — completeness pass over the 248 branch-only files

Prepared 2026-09-05. Scope: every path in
`/Users/scottg/Claude_all/H16P_GAPS_branch_paths.txt` (248 files, all present,
all opened), read against `H16P_SUMMARY.md` §0, §2(a) Engines 4–5, §2(c),
§4 R12–R14 and §7. Worktrees are read-only; nothing was modified.

**Global negative result, established mechanically over all 248 files.** A
recursive scan of every `total`, `count`, `new_total`, `max_total`,
`origin_roots`, `roots`, `nests` and `stab` field in every `.json`/`.jsonl`
here returns **exactly zero records with a value ≥ 5** (the only ≥5 hit is
`q4/sixth/boundary_diagnostic.json:"count": 40`, a row count). A `grep` for
`CANDIDATE|TRIGGER` over all 248 files matches **four files only**, all of
them D1 trigger *machinery* or the trigger's own rejection — never a data
record. Maximum cycle count anywhere in this set is **4**, always as 3 origin
+ 1 remote.

---

## 1. COVERAGE MANIFEST

One line per file. No file is UNREAD; five `.npy` arrays are marked
"UNREAD numerically" because this Mac has no `numpy` — their shapes and
contents are nevertheless pinned by the generating script and the paired log,
both of which were read in full.

- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/.gitignore` — gitignore (see NEW FACT G1: it ignores data/queue.log and data/QUEUE_DONE)
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-d1-2026-09-05/FASTRA_D1_COUNTER_DISCREPANCY_2026_09_05.md` — THE COUNTER-DISCREPANCY REPORT. 5 fields / 10 binary128-confirmed roots the production counter missed; exact rational vectors; counter settings table; per-field forensics.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/astra_afternoon_2026_09_05/.gitignore` — gitignore (see NEW FACT G1: it ignores data/queue.log and data/QUEUE_DONE)
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/astra_afternoon_2026_09_05/audit_proposition_a.py` — Independent sympy elimination + reflection certificates for D2 Proposition A; builds proposition_a_exact.json.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/astra_afternoon_2026_09_05/center_sign_map.csv` — 24-row table: K, exact c, inner/outer root radii, outside D, edge D, signs, comparison, full rational vector.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/astra_afternoon_2026_09_05/center_sign_map.jsonl` — 24 records, one per K on the centerward positive sheet. Each: binary128 fold, 21-point local profile, EXACTLY 2 roots, outside/edge returns, 12-point tail to log-r 40. All 816 return points OK_NUMERICAL, all comparisons 'agree'.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/astra_afternoon_2026_09_05/center_sign_map.log` — stdout of run_center_sign_map.py: 24 lines, roots ~5.5304 and ~8.2379, both signs +1, edge_kind configured_grid_end.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/astra_afternoon_2026_09_05/domain_edge_map.csv` — 24 rows: last success / first failure log-radius (112.5 to 629.3), failure status ANGULAR_CHART_UNRESOLVED, edge log-displacement ~94-106, all agree.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/astra_afternoon_2026_09_05/domain_edge_map.jsonl` — 24 records with 338 binary128 probes (202 OK, 136 ANGULAR_CHART_UNRESOLVED) bisecting the return-domain edge.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/astra_afternoon_2026_09_05/domain_edge_map.log` — stdout of extend_domain_edge.py; edge brackets per K.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/astra_afternoon_2026_09_05/domain_endpoint_checks.jsonl` — 2 records (K=1e-10 and K=0.001953125) probing log-radius 60/100/200/500/1000.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/astra_afternoon_2026_09_05/example_request.json` — Sample full_return128 request (KKL c=969/1000, m=37).
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/astra_afternoon_2026_09_05/example_result.json` — Its result: 5 OK_NUMERICAL points at log-r -20,1.9,2.3,37,44.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/astra_afternoon_2026_09_05/extend_domain_edge.py` — Driver that bisects the binary128 return-domain edge for the same 24 fields.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/astra_afternoon_2026_09_05/full_return128.py` — Python wrapper: rational/string-only one-call binary128 full-turn return; rejects float input; auto y_scale; compiles+caches the C++ binary.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/astra_afternoon_2026_09_05/full_return128_validation.json` — 6 analytic controls (see NEW FACT A4) with absolute errors 0 to 4.6e-31.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/astra_afternoon_2026_09_05/inherited_center_events.jsonl` — 10 ACCEPTED_NUMERICAL_FOLD events (positive_center_000..009), K 1.953125e-3 down to 1e-10, fold r~6.7579, c~0.96862; seeds for the sign map.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/astra_afternoon_2026_09_05/inherited_half_quad.cpp` — binary128 two-half-passage matcher (beta-free variant of half_beta_quad) used to re-correct folds.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/astra_afternoon_2026_09_05/proposition_a_exact.json` — Exact verdict: literal Proposition A GAP (false converse, missing H=0); corrected implication VERIFIED incl. b=0 and l=-1; counterexample vector and its three equilibrium traces -1,0,+1.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/astra_afternoon_2026_09_05/run_center_sign_map.py` — The centerward sign-map driver (24 K values, Newton fold correction, 21-point local profile, 27-step bisection per root, tail to e^40, edge bisection).
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/astra_afternoon_2026_09_05/summarize_and_check.py` — Re-validates the 24 records and writes summary.json; asserts 24 K, 48 root brackets, 24+24 sign agreements, 0 disagreements.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/astra_afternoon_2026_09_05/validate_full_return128.py` — The 6 analytic controls plus the float-input rejection test.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/audit/fable_d2_C_signmap.py` — D2 statement (C): scans 161x141 (a,l) on the order-three Shi stratum, sign(sigma*eta3) over all saddles.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/audit/fable_d2_loop_div.py` — D2 step 1: integrates the saddle separatrix loop, records divergence structure and tests 5 Dulac candidates on the interior.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/audit/fable_d2_theoremA.py` — Exact sympy proof of the one-way Proposition A implication (div=0 off-origin iff a(b+2l)=0 or C3=0).
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-d1-2026-09-05/audit/fable_engine/.gitignore` — gitignore (see NEW FACT G1: it ignores data/queue.log and data/QUEUE_DONE)
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/audit/fable_engine/data/D2_C_signmap.npy` — 22539x7 float array (a,l,x,y,sigma,eta3,sign) of saddles on the order-three stratum. Binary; not decoded here (numpy unavailable on this Mac) - UNREAD numerically, but generated by fable_d2_C_signmap.py above.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/audit/fable_engine/data/D2_cross_a-1.5.log` — D2 neutral-crossing refinement at a=-1.5: 4 bisections then 'lost branch at b=1.265625' and a TypeError crash (NoneType not subscriptable) in fable_d2_crossing_refine.py line 20.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/audit/fable_engine/data/D2_cross_a-2.log` — Same crash at a=-2, also 'lost branch at b=1.265625'.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/audit/fable_engine/data/D2_loop_div.json` — 22 focus-type homoclinic loops: sigma, eta2, integral of div, div sign changes (always 2), winding 1.000, Dulac positive-fraction for k in {m/a, -(2l+b)/b, 1, -1, 2}, and the range of 1+ax+by on the interior.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/audit/fable_engine/data/F11.log` — F11 loop scan: 160 lines, 'DONE loops found: 159'; each line a=..,b=..,l=..,saddle,trace,eta2.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/audit/fable_engine/data/F11_branch_a-2.jsonl` — 31 branch records at a=-2 (a,b,l,m,which,saddle,trace,eta2,sig_lo,sig_hi).
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/audit/fable_engine/data/F11_branch_a-2.log` — Its stdout, tagged CENTER/FOCUS per row.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/audit/fable_engine/data/F11_branch_a-3.jsonl` — 27 branch records at a=-3.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/audit/fable_engine/data/F11_branch_a-3.log` — Its stdout.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/audit/fable_engine/data/F11_branch_a2.jsonl` — 2 records at a=+2 (mirror of a=-2 with all signs flipped).
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/audit/fable_engine/data/F11_branch_a2.log` — Its 2-line stdout.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/audit/fable_engine/data/F11_loops.jsonl` — 159 loop records - the master F11 ledger consumed by fable_d2_loop_div.py.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/audit/fable_engine/data/F11_neutral_a-3.jsonl` — EMPTY (0 bytes) - no neutral point was written.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/audit/fable_engine/data/F11_neutral_a-3.log` — 3 coarse rows + DONE; trace crosses zero between b=1.2 and b=1.4 at a=-3 but eta2 crosses too.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/audit/fable_engine/data/F11_refine_a-3.log` — 20 rows refining b in [1.22,1.36]; wherever trace hits 0 exactly, eta2 collapses to ~1e-8 (a center), and 'no loop found' at b=1.27,1.28.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/audit/fable_engine/data/F12.log` — KKL double centre c*=0.968620633553494, alpha*=-37.136414809497211; 25-radius displacement profile |D|<=4.4e-9; 45 ovals; singular values show span dim 3; zero histogram over 300k random directions {0:207401,1:92245,2:342,3:12}; max 3.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/audit/fable_engine/data/F12b.log` — Same centre out to x0=3000: dim 3, histogram {0:273151,1:124439,2:2405,3:5}; five 3-zero directions with zero locations and min|M|/max|M| 5e-10..1e-8.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/audit/fable_engine/data/F13_origin.jsonl` — 19 Yu-Han origin-annulus records (65 ovals each): dim 3 always, max_zeros 2 (one record 1), phi3_zeros {0:13,1:4,2:2}, phi2_zeros 0 in 17 of 19.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/audit/fable_engine/data/F13_origin.log` — Its stdout - the same dicts including full phi3 profiles (65 points each).
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/audit/fable_engine/data/F13_second.jsonl` — 19 records for the second centre: dim 3, max_zeros 2 (14) or 1 (5), phi3_zeros 1 in 18 of 19, phi2_zeros 0 in all 19.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/audit/fable_engine/data/F13_second.log` — Its stdout.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/audit/fable_engine/data/F14_kklx_seeds_evolve.jsonl` — 7 elite records, total=3 throughout, score 3.0 -> 3.935.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/audit/fable_engine/data/F14_kklx_seeds_evolve.log` — 5 generations, best score 3.9991, per-generation totals max 3.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/audit/fable_engine/data/F14_kklx_seeds_evolve2.jsonl` — EMPTY (0 bytes).
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/audit/fable_engine/data/F14_kklx_seeds_evolve2.log` — 4 generations, best score 4.8042 (score, not count), totals max 3.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/audit/fable_engine/data/F14_kklx_seeds_evolve3.jsonl` — 77 elite records, EVERY one total=4 score=4.0 - 60 generations with no improvement past four.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/audit/fable_engine/data/F14_kklx_seeds_evolve3.log` — 53 generations; population totals reach 4 in 1-6 members per generation, never 5.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/audit/fable_engine/data/F15_kklstar_sweep.jsonl` — THE 24 FOUR-CYCLE KKL* FIELDS. c in [0.9491,0.9883], K in [-0.0574,0.0410], al in [-38.60,-35.77], |beta|<=1.9e-4, p~-10, q~2.2. Every one is 3 origin (SUS/USU) + 1 remote at radius 1e6-1e9.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/audit/fable_engine/data/F15_mvneutral_sweep.jsonl` — EMPTY (0 bytes) - the --store=3 filter never fired.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/audit/fable_engine/data/F15_mvneutral_sweep.log` — 313 progress lines; final '20032 sets 5906s hist {0: 11574, 1: 6453, 2: 1988, 3: 17}'. Max 3.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/audit/fable_engine/data/F15_targeted.jsonl` — 270 records; total hist {0:213,1:52,2:4,3:1}.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/audit/fable_engine/data/F15_targeted.log` — 271 lines ending 'DONE {0: 213, 1: 52, 2: 4, 3: 1}'.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/audit/fable_engine/data/F16_b0.5.log` — Dulac-coefficient rank at b=0.5: base centre |D|<=2.44e-14; da and db give pure noise (~1e-9), e0 gives O(1), e1/e2 give O(y).
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/audit/fable_engine/data/F16_b1.00.log` — Same at b=1.0 plus the full 3-window least-squares fits and leading-coefficient ranks.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/audit/fable_engine/data/F16_b1.5.log` — Same at b=1.5.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/audit/fable_engine/data/F17.log` — Second-order Melnikov at (a,b)=(-1,1) along the first-order-null direction: upper annulus max|D/eps| ~1e-11 (null confirmed) and M2 has ZERO zeros at every eps; lower annulus max|D/eps|=1.078 (not null) and M2 is strictly negative.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/audit/fable_engine/data/F18_third_cycle.jsonl` — EMPTY (0 bytes) - the HIT filter (max nest >=3 or total >=4) never fired in 3600 fields.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/audit/fable_engine/data/F18_third_cycle.log` — 145 lines, 'DONE {0: 1206, 1: 1506, 2: 642, 3: 246}' over 3600 fields. No HIT line anywhere.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/audit/fable_engine/data/F18b_alien_points.jsonl` — EMPTY (0 bytes).
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/audit/fable_engine/data/F18b_alien_points.log` — 'DONE {0: 750, 1: 1122, 2: 474, 3: 54}' over 2400 fields at a0=-0.5. No HIT.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/audit/fable_engine/data/F18c_scaled.jsonl` — EMPTY (0 bytes).
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/audit/fable_engine/data/F18c_scaled.log` — 'DONE {0: 630, 1: 2016, 2: 1638, 3: 252}' over 4536 fields. No HIT.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/audit/fable_engine/data/F19_q4_alien.jsonl` — 126 records. first_order_zeros {3:90, 2:36}; ACTUAL origin cycles {3:54, 2:42, 1:28, 0:2}; stabilities SUS/US/SU/S; 27 distinct directions x 3 eps. Total across all nests never exceeds 3.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/audit/fable_engine/data/F19_q4_alien.log` — The FIRST (pre-fix) F19 run: 70 ovals, span dim 4, random-direction histogram {0:294422,1:5531,2:47} - only 2-zero directions found, so only 2-zero picks were shot.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/audit/fable_engine/data/F19_q4_alien_fixed.log` — The fixed run with the targeted-null-vector construction: 'targeted: 3-zero directions 15  4-zero directions 0'; 45 shots at 3-zero directions; no 'FOUR ZEROS' line, no 'ALIEN?' flag.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/audit/fable_engine/data/F19b_q4_alien.log` — A re-run of the fixed lane, same 15 directions, differing only in the double-precision origin counts of a few borderline rows (edge 3.68e15 vs 2.01e17).
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/audit/fable_engine/data/F20_kklstar_evolve.jsonl` — 1 record: total 4, the elite.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/audit/fable_engine/data/F20_kklstar_evolve.log` — 42 lines; 24 seeds, best SCORE 5.50 held for all 40 generations, 'DONE best_total 4'. The 5.50 is a score, not a cycle count.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/audit/fable_engine/data/F21_shi_compact.jsonl` — EMPTY (0 bytes) - the store filter never fired in 30016 fields.
- `/Users/scottg/Claude_all/H16P_branches/fable_compute-f3-lam0/audit/fable_engine/data/F3_lam0_L2_worker.jsonl` — 2873 cloud-worker records, count hist {1:2854, 2:19}.
- `/Users/scottg/Claude_all/H16P_branches/fable_compute-f3-lam0/audit/fable_engine/data/F3_lam0_L2_worker.log` — 46 lines, last '94208 sets, 769s, hist {0: 153852, 1: 2854, 2: 19}' - never reached DONE.
- `/Users/scottg/Claude_all/H16P_branches/fable_compute-f3-lam0/audit/fable_engine/data/F3_lam0_L8_worker.log` — 391 lines, 'DONE 800000 6390.46213722229 {0: 1307502, 1: 23009, 2: 133}'. The largest single sweep in the repository.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/audit/fable_engine/data/F4_kkl_evolve.jsonl` — 351 elite records, every one total=4; score max 4.998215940519204.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/audit/fable_engine/data/F4_kkl_evolve.log` — Seed score 4.952 with nests ([6050.4],U) + ([0.938,3.021,23.393],SUS); 30 generations; per-generation total-4 counts 3-17 of 96.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/audit/fable_engine/data/F4_kkl_evolve2.log` — 60 generations, 'DONE best_total 4 elite score 4.998177894069219'; up to 26 of 96 population members at total 4.
- `/Users/scottg/Claude_all/H16P_branches/fable_compute-f5-shi/audit/fable_engine/data/F5_shi_full_L4_worker.jsonl` — 18 worker records, all count=2.
- `/Users/scottg/Claude_all/H16P_branches/fable_compute-f5-shi/audit/fable_engine/data/F5_shi_full_L4_worker.log` — 129 lines, last '264192 sets, 2357s, hist {0: 408772, 1: 30835, 2: 18}' - never reached DONE.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/audit/fable_engine/data/RECOUNT_fixed_counter.jsonl` — 24 records, source data/F15_kklstar_sweep.jsonl, old_total->new_total is (4,4) for all 24. Confirms the noise-floor fix changed no verdict.
- `/Users/scottg/Claude_all/H16P_branches/fable_compute-evolve/audit/fable_engine/data/W_evolve_kkl_all.jsonl` — 34 elites, all total=3.
- `/Users/scottg/Claude_all/H16P_branches/fable_compute-evolve/audit/fable_engine/data/W_evolve_kkl_all.log` — seed 3 (SUS origin only, no remote), 'DONE best_total 3 elite score 3.981348964145437'.
- `/Users/scottg/Claude_all/H16P_branches/fable_compute-evolve/audit/fable_engine/data/W_evolve_kkl_wide.jsonl` — 2 elites, total=3.
- `/Users/scottg/Claude_all/H16P_branches/fable_compute-evolve/audit/fable_engine/data/W_evolve_kkl_wide.log` — sigma 0.15, 'DONE best_total 3 elite score 3.790945811233903'.
- `/Users/scottg/Claude_all/H16P_branches/fable_compute-evolve/audit/fable_engine/data/W_evolve_shi_all.jsonl` — 1 elite, total=4: remote ([7.0968],S) + origin ([0.00414,0.00613,0.06434],USU).
- `/Users/scottg/Claude_all/H16P_branches/fable_compute-evolve/audit/fable_engine/data/W_evolve_shi_all.log` — 'DONE best_total 4 elite score 4.0' - 60 generations, no improvement.
- `/Users/scottg/Claude_all/H16P_branches/fable_compute-evolve/audit/fable_engine/data/W_evolve_shi_wide.jsonl` — EMPTY (0 bytes).
- `/Users/scottg/Claude_all/H16P_branches/fable_compute-evolve/audit/fable_engine/data/W_evolve_shi_wide.log` — 59 generations at sigma 0.15, best score 4.0000, never finished.
- `/Users/scottg/Claude_all/H16P_branches/fable_compute-evolve/audit/fable_engine/data/W_evolve_yz_all.jsonl` — 1 elite, total=4: remote ([241.534],S) + origin ([0.0294,0.0643,0.0952],USU).
- `/Users/scottg/Claude_all/H16P_branches/fable_compute-evolve/audit/fable_engine/data/W_evolve_yz_all.log` — 'DONE best_total 4 elite score 4.846413315779456'.
- `/Users/scottg/Claude_all/H16P_branches/fable_compute-evolve/audit/fable_engine/data/W_evolve_yz_wide.jsonl` — 1 elite, identical to yz_all.
- `/Users/scottg/Claude_all/H16P_branches/fable_compute-evolve/audit/fable_engine/data/W_evolve_yz_wide.log` — 'DONE best_total 4 elite score 4.846413315779456' - the wide search found nothing the narrow one did not.
- `/Users/scottg/Claude_all/H16P_branches/fable_compute-pert/audit/fable_engine/data/W_kklx_L3.jsonl` — 1141 worker records; count hist {2:1115, 3:26}. Root pairs span radii 0.013 to 494.
- `/Users/scottg/Claude_all/H16P_branches/fable_compute-pert/audit/fable_engine/data/W_kklx_L3.log` — 28 lines, last '57344 sets, 4015s, hist {0: 63135, 1: 27938, 2: 1085, 3: 24}' - partial, no DONE.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/audit/fable_engine/data/W_kklx_L3_hits_reeval.jsonl` — 11 fields re-evaluated with the fixed counter: total {3:7, 4:4}. All four 4s are 3-origin (S,U,S) + 1 remote U at radius 1.3e3 to 2.9e4.
- `/Users/scottg/Claude_all/H16P_branches/fable_compute-pert/audit/fable_engine/data/W_mv_L3.jsonl` — 1 record, count 2.
- `/Users/scottg/Claude_all/H16P_branches/fable_compute-pert/audit/fable_engine/data/W_mv_L3.log` — 'DONE 100096 610.85 {0: 201865, 1: 602, 2: 1}'. Max 2 in 100k sets.
- `/Users/scottg/Claude_all/H16P_branches/fable_compute-pert/audit/fable_engine/data/W_mvpert_e2.jsonl` — 113 records, all count=2.
- `/Users/scottg/Claude_all/H16P_branches/fable_compute-pert/audit/fable_engine/data/W_mvpert_e2.log` — 'DONE 60160 2055.21 {0: 100660, 1: 19544, 2: 113}'. Max 2.
- `/Users/scottg/Claude_all/H16P_branches/fable_compute-pert/audit/fable_engine/data/W_q3rpert_e2.jsonl` — EMPTY (0 bytes).
- `/Users/scottg/Claude_all/H16P_branches/fable_compute-pert/audit/fable_engine/data/W_q3rpert_e2.log` — 'DONE 60160 117.85 {0: 58678, 1: 1482}'. MAX 1.
- `/Users/scottg/Claude_all/H16P_branches/fable_compute-pert/audit/fable_engine/data/W_q3rpert_e4.jsonl` — EMPTY (0 bytes).
- `/Users/scottg/Claude_all/H16P_branches/fable_compute-pert/audit/fable_engine/data/W_q3rpert_e4.log` — 'DONE 60160 116.98 {0: 58694, 1: 1466}'. MAX 1.
- `/Users/scottg/Claude_all/H16P_branches/fable_compute-pert/audit/fable_engine/data/W_q4pert_e2.jsonl` — EMPTY (0 bytes).
- `/Users/scottg/Claude_all/H16P_branches/fable_compute-pert/audit/fable_engine/data/W_q4pert_e2.log` — 'DONE 60160 145.84 {0: 119550, 1: 770}'. MAX 1.
- `/Users/scottg/Claude_all/H16P_branches/fable_compute-pert/audit/fable_engine/data/W_q4pert_e4.jsonl` — EMPTY (0 bytes).
- `/Users/scottg/Claude_all/H16P_branches/fable_compute-pert/audit/fable_engine/data/W_q4pert_e4.log` — 'DONE 60160 138.24 {0: 119517, 1: 803}'. MAX 1.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/audit/fable_engine/data/astra_missed_roots.jsonl` — Fable's copy of Astra's counter_check/missed_roots.jsonl - byte-identical (86409 bytes each). 5 fields, 10 roots.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/audit/fable_engine/data/f16_G_b0.50.npy` — 6x97 array (Y grid + 5 direction derivatives) for b=0.5. Binary; UNREAD numerically (no numpy on this host); content described by F16_b0.5.log.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/audit/fable_engine/data/f16_G_b1.00.npy` — Same for b=1.00. UNREAD numerically.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/audit/fable_engine/data/f16_G_b1.50.npy` — Same for b=1.50. UNREAD numerically.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/audit/fable_engine/data/f17_M2_lower.npy` — 2x120 array (Y, D/eps^2) lower annulus. UNREAD numerically; summarised by F17.log (strictly negative, no zeros).
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/audit/fable_engine/data/f17_M2_upper.npy` — 2x120 array upper annulus. UNREAD numerically; F17.log records zero sign changes.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/audit/fable_engine/data/f19_directions.json` — The picked directions: hist {0:19630, 1:369, 2:1}; picks['2'] holds the single 2-zero direction found by RANDOM sampling; picks['3'] and picks['4'] are EMPTY LISTS.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/audit/fable_engine/data/f19_q4_basis.npy` — 70x13 Melnikov basis (x0 + 12 monomial integrals) at the Q4 resonant graphic. UNREAD numerically.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/audit/fable_engine/data/kklx_seeds.jsonl` — 11 seeds for F14: count {3:7, 4:4}, all in the extended KKL chart with free x^2/xy coefficients.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/audit/fable_engine/data/queue2.log` — 0 bytes by construction (stdout redirected per job).
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/audit/fable_engine/data/queue3.log` — 0 bytes.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/audit/fable_engine/data/queue4.log` — 0 bytes.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/audit/fable_engine/data/queue5.log` — 0 bytes.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/audit/fable_engine/data/queue6.log` — 0 bytes.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/audit/fable_engine/data/queue7.log` — 0 bytes.
- `/Users/scottg/Claude_all/H16P_branches/fable_compute-f3-lam0/audit/fable_engine/data/worker_driver.log` — 1 line: 'STAGE1_DONE rc=0'.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/audit/fable_engine/f15_targeted.py` — F15 targeted lane: 3 b-values x 5 da x 3 magnitudes x 6 random directions on the MV neutral two-centre family. Prints 'FIVE OR MORE' if total>=5 - never printed.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/audit/fable_engine/f16_dulac_rank.py` — F16: numerical Dulac-coefficient rank at the neutral hemicycle; fits G_i(u) ~ c0 + c1 u w + c2 w + c3 u w^2 + ... in w=e^-u and reports the rank of the leading-coefficient matrix over 5 unfolding directions.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/audit/fable_engine/f17_second_order.py` — F17: second-order Melnikov along the first-order-null direction e1=1, e2=-2 at the holomorphic point.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/audit/fable_engine/f18_third_cycle.py` — F18: codim-2 hemicycle unfolding, 3 b x 6 e1 x 8 da x 5 residual x 5 e0 = 3600 fields. HIT filter mx>=3 or tot>=4.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/audit/fable_engine/f18b_alien_points.py` — F18b: identical grid with a0=-0.5 instead of -1 (off the neutral line), 2400 fields.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/audit/fable_engine/f18c_scaled.py` — F18c: the SCALED unfolding e0 ~ delta^2, da ~ delta implied by D ~ pi e0 - pi delta y + c2 y^2 + da delta y log y; 4536 fields.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/audit/fable_engine/f19_q4_alien.py` — F19: the Q4 alien-cycle lane. Q4 centre X'=Y+6X^2+4XY-2Y^2, Y'=-X+2X^2+8XY-2Y^2, annulus boundary x*=0.2272111321. Includes the TARGETED 3-zero null-vector construction and a 'FOUR ZEROS at first order' print that never fired.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/audit/fable_engine/push_fold.py` — THE FOLD-PUSH LANE (not described anywhere in the summary). Finite-difference gradient descent on 7 free coefficients driving the interior minimum of D/r between the 2nd and 3rd origin roots toward zero, to birth a 4th and 5th origin cycle. Prints 'CROSSED ZERO' on success. No data/push_fold_result.json exists on the branch - the lane produced no committed output.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/audit/fable_engine/queue2.sh` — Night-watch queue 2: F4 kkl evolve, F8 yz/shi evolve, F6 q3rpert, F9 mvpert, F7 q4pert, F3b shi lam0, F3c evolve, F9 mv, F5 kklx (L=1.5).
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/audit/fable_engine/queue3.sh` — Queue 3: adds F14 kklx-seeds evolve2; F5 kklx at L=3.0.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/audit/fable_engine/queue4.sh` — Queue 4: starts with F14 kklx-seeds evolve3; otherwise queue 3 minus the first job.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/audit/fable_engine/queue5.sh` — Queue 5: sweep_log mvneutral 20000 sets seed 41, then kklstar 10000 sets seed 43, both --store=3; writes data/F15_QUEUE_DONE.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/audit/fable_f11_branch.py` — F11 branch tracker at fixed a over b, 36 l-values per b, 3-way Pool.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/audit/fable_f11_neutral.py` — F11 neutral-point finder: brentq on saddle trace along the focus-type loop branch.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/audit/fable_f12_basis.npy` — 45x13 Melnikov basis at the KKL double centre. UNREAD numerically; F12.log gives the singular values and histogram.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/audit/fable_f12_double_center.py` — F12: KKL double centre (beta=0, K=0, J(c)=305+634c-11c^2-1000c^3=0); verifies the origin is a centre and computes the 12-direction Melnikov basis.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/audit/fable_f12b_basis.npy` — 90x13 basis out to x0=3000. UNREAD numerically.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/audit/fable_f12b_double_center_far.py` — F12b: the same basis pushed to large amplitude toward the neutral two-saddle infinity graphic.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-afternoon-2026-09-05/audit/fable_f13_yuhan_melnikov.py` — F13: Yu-Han reversible two-centre family x'=y(1+a1 x), y'=-x+x^2+a4 y^2; Melnikov span, phi2/phi3 constrained elements, Taylor-rank test on the Yu-Han curve a4=(a1-5)/3.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-d1-2026-09-05/fastra_d1_2026_09_05/.gitignore` — gitignore (see NEW FACT G1: it ignores data/queue.log and data/QUEUE_DONE)
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-d1-2026-09-05/fastra_d1_2026_09_05/accepted.jsonl` — 44 accepted D1 folds with full history. Labels listed in NEW FACT D6.
- `/Users/scottg/Claude_all/H16P_branches/codex_fastra-d1-fold-counts-2026-09-05/fastra_d1_2026_09_05/artifact_audit.json` — codex-branch record audit: PASS, 303 exact vectors validated, fields.jsonl has 266 records (vs 176 on the astra branch), checks include 'raw trigger rejection' and 'credential scan'; interval_certification false.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-d1-2026-09-05/fastra_d1_2026_09_05/audit.json` — astra-branch D1 audit: PASS, accepted 44, field_records 176, dense_profiles 88, 132 independent two-half sign brackets, 4 cross-engine brackets, fable_sources_unchanged true, exhaustive_root_isolation FALSE.
- `/Users/scottg/Claude_all/H16P_branches/codex_fastra-d1-fold-counts-2026-09-05/fastra_d1_2026_09_05/audit_artifacts.py` — The codex audit script; asserts verified_precursor has exactly 4 roots with origin USU and remote S, each with a sign bracket, and fold residuals |F|<=1e-19, |G|<=1e-15.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-d1-2026-09-05/fastra_d1_2026_09_05/audit_records.py` — The astra audit script; asserts 44 accepted, 176 fields, 88 dense, the exact vector shape [0,0,1,1,1,0,0,-m,beta,-10,11/5,c], and that four Fable sources are byte-identical to git 'fable/current'.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-d1-2026-09-05/fastra_d1_2026_09_05/build_engines.sh` — Builds libretmap.so, libretmap_log.so (gcc -O3 -fopenmp) and /tmp/d1_half, /tmp/d1_full (g++ -O2 -lquadmath).
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-d1-2026-09-05/fastra_d1_2026_09_05/build_report.py` — Builds the scope-honest D1 ledger (sign_map.csv, root_ledger.csv, summary) from saved records.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-d1-2026-09-05/fastra_d1_2026_09_05/calibration.json` — 5 calibration fields with their counter output. See NEW FACT D7.
- `/Users/scottg/Claude_all/H16P_branches/codex_fastra-d1-fold-counts-2026-09-05/fastra_d1_2026_09_05/center_precision.log` — 15 lines: positive_center_005..009 pair SU, hopf_1 USU, hopf_0.1 USU after binary128 repair.
- `/Users/scottg/Claude_all/H16P_branches/codex_fastra-d1-fold-counts-2026-09-05/fastra_d1_2026_09_05/center_precision_final.log` — 3 lines: positive_center_004 pair SU, hopf USU.
- `/Users/scottg/Claude_all/H16P_branches/codex_fastra-d1-fold-counts-2026-09-05/fastra_d1_2026_09_05/check_sign_disagreements.py` — Re-runs matching_profile on the labels summary.json lists as sign_disagreements.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-d1-2026-09-05/fastra_d1_2026_09_05/continuation.jsonl` — 44 rows, completed=True for all 44.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-d1-2026-09-05/fastra_d1_2026_09_05/counter_check/audit.json` — {status PASS, fields 5, full_return_pair_roots 10, same_ray_grid_endpoints 20, fable_production_sources_unchanged true, no_new_sheets_or_sweeps true, D1 OPEN}.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-d1-2026-09-05/fastra_d1_2026_09_05/counter_check/build_report.py` — Builds the discrepancy handoff markdown + discrepancy_ledger.json from calls/missed_roots/refined_coordinates.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-d1-2026-09-05/fastra_d1_2026_09_05/counter_check/calls.jsonl` — 178 binary128 full-return calls, ALL NUMERICAL_ONLY, 388.5 s total. 88 negative_extension_1, 70 center_0.0000000001, 20 positive_extension_1.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-d1-2026-09-05/fastra_d1_2026_09_05/counter_check/center_tolerance.json` — The 4-level double-precision tolerance ladder at the C0 grid points. See NEW FACT D2.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-d1-2026-09-05/fastra_d1_2026_09_05/counter_check/check_misses.py` — The audit driver: transfers the section at 3 log-offsets, brackets, secant-refines in binary128, then replays the two adjacent DEFAULT grid points with the unmodified Fable counter.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-d1-2026-09-05/fastra_d1_2026_09_05/counter_check/coefficient_rounding.jsonl` — 6 records comparing binary128 on the exact rational field vs on the double-rounded field. See NEW FACT D4.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-d1-2026-09-05/fastra_d1_2026_09_05/counter_check/first_step.c` — A 12-line instrument that #includes the unmodified retmap_log.c and prints the FIRST Dormand-Prince step's error and the original rejection test.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-d1-2026-09-05/fastra_d1_2026_09_05/counter_check/first_step.json` — Its output for two fields. See NEW FACT D3 - the NaN-acceptance smoking gun.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-d1-2026-09-05/fastra_d1_2026_09_05/counter_check/full_ray_float_coeff_quad.cpp` — binary128 full-turn shooter with q_xy passed in (so the double-rounded 11/5 can be used).
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-d1-2026-09-05/fastra_d1_2026_09_05/counter_check/full_ray_quad.cpp` — binary128 full-turn angular shooter, modified-midpoint + Richardson, 9 variational states, tol-scaled step control, emits L, L_z, L_zz, L_c, L_m, L_zc, L_zm, multiplier, period. Also the clearest statement of the D1 chart: P = y + x^2 + xy, Q = -m x + beta y - 10 x^2 + (11/5) xy + c y^2.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-d1-2026-09-05/fastra_d1_2026_09_05/counter_check/grid_return_values.csv` — 20 rows - THE CORE TABLE. For each of the 10 missed roots, the two adjacent default grid points with binary128 vs native double log-return, displacement, status and whether the original run visited the point.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-d1-2026-09-05/fastra_d1_2026_09_05/counter_check/half_float_coeff_quad.cpp` — binary128 two-half matcher on the double-rounded coefficients.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-d1-2026-09-05/fastra_d1_2026_09_05/counter_check/missed_roots.jsonl` — 5 records, the full forensic payload. See NEW FACT D1.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-d1-2026-09-05/fastra_d1_2026_09_05/counter_check/precision.log` — 6 lines: the binary128 status/displacement at the six probed grid points (2 UNRESOLVED, both positive_extension_1).
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-d1-2026-09-05/fastra_d1_2026_09_05/counter_check/precision_check.py` — Separates double COEFFICIENT rounding from double TRAJECTORY error.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-d1-2026-09-05/fastra_d1_2026_09_05/counter_check/refine.log` — 8 lines: the safeguarded-Newton refined roots and residuals down to 5.6e-31.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-d1-2026-09-05/fastra_d1_2026_09_05/counter_check/refine_coordinates.py` — Safeguarded full-return Newton refinement inside recorded grid brackets; skips positive_extension_1 deliberately.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-d1-2026-09-05/fastra_d1_2026_09_05/counter_check/refined_coordinates.jsonl` — 8 refined roots with residuals 8.2e-30 to 3.9e-20 and 4-9 Newton iterations each.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-d1-2026-09-05/fastra_d1_2026_09_05/counter_check/rounded_half.jsonl` — 3 half-map results on the double-rounded near-fold field at log-offsets -0.8, 0, +0.8: F = +8.44e-18 at all three. See NEW FACT D5.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-d1-2026-09-05/fastra_d1_2026_09_05/counter_check/rounding_half.log` — Its 3-line stdout.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-d1-2026-09-05/fastra_d1_2026_09_05/counter_check/run.log` — 10 lines - the printed root/grid summary for all 10 missed roots.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-d1-2026-09-05/fastra_d1_2026_09_05/dense.jsonl` — 88 dense double-precision origin profiles on log-r [-25,46] at spacing 0.125 and two tolerances. Stabilities: SU 20, US 16, USU 14, '' 14, SUS 14, S 10. MAX 3 brackets. Zero sign_mismatch rows.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-d1-2026-09-05/fastra_d1_2026_09_05/dense.log` — 88 lines label/kind/stability/edge_kind (54 scan_cap, 34 integration_failure).
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-d1-2026-09-05/fastra_d1_2026_09_05/dense_profiles.py` — The dense-profile driver; uses an empirical two-tolerance resolution check e=max(5e-12, 10|d0-d1|) and marks unresolved signs 0.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-d1-2026-09-05/fastra_d1_2026_09_05/endpoint.log` — 3 lines: center_0.0000000001 pair full returns at offsets -0.8, 0, +0.8 - a genuine sign change (+9.9e-14, -9.2e-14, +5.5e-13).
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-d1-2026-09-05/fastra_d1_2026_09_05/endpoint_finish.log` — 11 lines: the hopf twin's identical bracket plus a tolerance control at 2e-26 vs 2e-27 on the positive endpoint.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-d1-2026-09-05/fastra_d1_2026_09_05/endpoint_full.jsonl` — 24 binary128 full-return replays (16 NUMERICAL_ONLY, 8 UNRESOLVED - all 8 are 25-second timeouts on positive_extension_1).
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-d1-2026-09-05/fastra_d1_2026_09_05/endpoint_replay.py` — Runs those replays at offsets -0.8, 0, 0.8, 2 with tol 2e-29.
- `/Users/scottg/Claude_all/H16P_branches/codex_fastra-d1-fold-counts-2026-09-05/fastra_d1_2026_09_05/events_negative.jsonl` — 30 events, 29 ACCEPTED_NUMERICAL_FOLD + 1 CORRECTOR_UNRESOLVED. K from -1e-4 to -2.598e12 (chart switches K->m at index 15), fold r from 6.7579 to 4.05e10.
- `/Users/scottg/Claude_all/H16P_branches/codex_fastra-d1-fold-counts-2026-09-05/fastra_d1_2026_09_05/events_positive_center.jsonl` — 10 accepted, K 1.953125e-3 down to 1e-10, fold r 6.75940 -> 6.75794.
- `/Users/scottg/Claude_all/H16P_branches/codex_fastra-d1-fold-counts-2026-09-05/fastra_d1_2026_09_05/events_positive_infinity.jsonl` — 26 events, 25 accepted + 1 CORRECTOR_UNRESOLVED at target log_r 42. K 0.00390625 -> 7.0694, fold r 6.7609 -> 1.739e18, c 0.9692 -> 1.5937.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-d1-2026-09-05/fastra_d1_2026_09_05/example.log` — 1 JSON line: the branch-B exact rational 3+1 field replayed in log-polar AND Cartesian, all four cycles bracketed both ways.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-d1-2026-09-05/fastra_d1_2026_09_05/fields.jsonl` — 176 field records (44 labels x 4 kinds fold/pair/hopf_small/hopf). Origin-nest counts {0:60, 3:46, 2:40, 1:30}; TOTAL roots {1:52, 3:40, 0:36, 4:32, 2:16} - 32 of 176 double evaluations reported four. Nest edge kinds: 200 integration_failure_bracket, 152 scan_cap.
- `/Users/scottg/Claude_all/H16P_branches/codex_fastra-d1-fold-counts-2026-09-05/fastra_d1_2026_09_05/finish_positive.log` — 6 lines ending 'ACCEPTED FINAL positive_infinity_024_polished c=1.5933653820754588 K=7.06316120411702 r=2.3539e17'.
- `/Users/scottg/Claude_all/H16P_branches/codex_fastra-d1-fold-counts-2026-09-05/fastra_d1_2026_09_05/finish_positive.py` — Far-positive endpoint finisher using the repaired controller (.matching_quad_v3.so); raises 'STOP FOUR ORIGIN' if a profile ever has 4 roots. Never raised.
- `/Users/scottg/Claude_all/H16P_branches/codex_fastra-d1-fold-counts-2026-09-05/fastra_d1_2026_09_05/finish_positive_extension.log` — 6 lines ending 'ACCEPTED FINAL positive_infinity_025_polished c=1.5937029026273091 K=7.06938707146707 r=1.7393e18'.
- `/Users/scottg/Claude_all/H16P_branches/codex_fastra-d1-fold-counts-2026-09-05/fastra_d1_2026_09_05/finish_positive_extension.py` — Same for index 025, using finite-difference curvature instead of the archived derivatives.
- `/Users/scottg/Claude_all/H16P_branches/codex_fastra-d1-fold-counts-2026-09-05/fastra_d1_2026_09_05/float64_field_collapses.json` — 5 labels (positive_infinity_021,022,023,024_polished,025_polished) where the exact FOLD vector and the exact PAIR vector round to the IDENTICAL float64 vector. See NEW FACT D8.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-d1-2026-09-05/fastra_d1_2026_09_05/full_beta_quad.cpp` — binary128 full-turn shooter with beta, used as /tmp/d1_full.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-d1-2026-09-05/fastra_d1_2026_09_05/half_beta_quad.cpp` — binary128 two-half matcher with beta, used as /tmp/d1_half.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-d1-2026-09-05/fastra_d1_2026_09_05/hp_calls.jsonl` — 639 binary128 calls: 608 NUMERICAL_TWO_HALF_PASSAGES, 23 NUMERICAL_ONLY, 8 UNRESOLVED (all 25-second /tmp/d1_full timeouts). Purposes: 132 pair, 132 hopf_small, 132 hopf pair brackets, 116 fold corrector, 44 pair-side centre tests, 44 beta calibrations, 20 endpoint replays, 8 fixed-radius, 6 tolerance controls, 5 missing-record replays.
- `/Users/scottg/Claude_all/H16P_branches/codex_fastra-d1-fold-counts-2026-09-05/fastra_d1_2026_09_05/make_manifest.py` — Writes MANIFEST.json; records fable_input_commit 4cece20 and astra_fold_input_commit 7db8597cb7d9bb34e119e85bec3f229270eaf1aa.
- `/Users/scottg/Claude_all/H16P_branches/codex_fastra-d1-fold-counts-2026-09-05/fastra_d1_2026_09_05/matching.py` — ctypes wrapper for .matching_quad.so plus matching_profile(): u-grid = arange(umin,umax,du) union {-8..1} union a dense +-0.65 window at spacing 0.065 around the fold radius; roots with min|F| below `noise` go to uncertain_sign_changes instead of roots.
- `/Users/scottg/Claude_all/H16P_branches/codex_fastra-d1-fold-counts-2026-09-05/fastra_d1_2026_09_05/matching_quad.cpp` — binary128 two-sided angular matcher with y=sqrt(m)Y conditioning AND the sensitivity-tolerance floor max(tol,1e-24) (see line 38) - the fix that unblocked the far-positive endpoint.
- `/Users/scottg/Claude_all/H16P_branches/codex_fastra-d1-fold-counts-2026-09-05/fastra_d1_2026_09_05/matching_quad_initial.cpp` — The pre-fix version: single tol for both components, no x0/y0/th0/sense fields. This is the one that exhausted the step guard.
- `/Users/scottg/Claude_all/H16P_branches/codex_fastra-d1-fold-counts-2026-09-05/fastra_d1_2026_09_05/matching_quad_remote.cpp` — The remote-equilibrium variant (arbitrary x0,y0,th0,sense) used by verify_precursor.py for the remote cycle; still uses the single tol.
- `/Users/scottg/Claude_all/H16P_branches/codex_fastra-d1-fold-counts-2026-09-05/fastra_d1_2026_09_05/negative.log` — 28 lines; each negative_* field is 1 remote U + origin US (pair) / SUS (hopf) = total 4. ACCEPTED rows give exact r, c, K.
- `/Users/scottg/Claude_all/H16P_branches/codex_fastra-d1-fold-counts-2026-09-05/fastra_d1_2026_09_05/negative_resume.log` — 1 line: 'FAILED negative -20'.
- `/Users/scottg/Claude_all/H16P_branches/codex_fastra-d1-fold-counts-2026-09-05/fastra_d1_2026_09_05/negative_resume2.log` — 85 lines; ends in a Traceback.
- `/Users/scottg/Claude_all/H16P_branches/codex_fastra-d1-fold-counts-2026-09-05/fastra_d1_2026_09_05/negative_resume3.log` — 63 lines; the successful continuation into the m chart.
- `/Users/scottg/Claude_all/H16P_branches/codex_fastra-d1-fold-counts-2026-09-05/fastra_d1_2026_09_05/positive_center.log` — 50 lines. Indices 000-005 give remote 1S + origin SU(pair)/USU(hopf) = total 4; from index 006 on the ORIGIN COUNT COLLAPSES TO 0 - these are exactly the fields the counter_check later showed have 2 real roots.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-d1-2026-09-05/fastra_d1_2026_09_05/positive_full_control.jsonl` — 6 binary128 full-return controls on positive_extension_1 at tolerances 2e-26 and 2e-27, three offsets each; displacements agree to 3 digits and keep the +,-,+ bracket.
- `/Users/scottg/Claude_all/H16P_branches/codex_fastra-d1-fold-counts-2026-09-05/fastra_d1_2026_09_05/positive_infinity.log` — 104 lines. Indices 001-018 give remote 1S + origin SU/USU = total 4. Then the FOUR_ORIGIN_TRIGGER fires: RuntimeError('STOP: four origin roots require hostile reproduction') from engine.py line 116.
- `/Users/scottg/Claude_all/H16P_branches/codex_fastra-d1-fold-counts-2026-09-05/fastra_d1_2026_09_05/positive_infinity_resume.log` — 52 lines. Every field re-done with 'binary128 full-grid recheck'. Indices 019-024: the REMOTE nest is gone (0 roots) and the origin gives SU/US - so total drops to 2, not 4.
- `/Users/scottg/Claude_all/H16P_branches/codex_fastra-d1-fold-counts-2026-09-05/fastra_d1_2026_09_05/precision_controls.jsonl` — 3 rows for positive_infinity_024_pair at tol 2e-27/2e-29/2e-30: F values agree to 4 digits (2.75e-25, -1.50e-25, 1.14e-25), all status 0, engine matching_quad_v3.
- `/Users/scottg/Claude_all/H16P_branches/codex_fastra-d1-fold-counts-2026-09-05/fastra_d1_2026_09_05/precision_repairs.jsonl` — 28 binary128 re-profiles. Stabilities USU 17, SU 9, SUS 2; root counts {3:19, 2:9}. NEVER 4.
- `/Users/scottg/Claude_all/H16P_branches/codex_fastra-d1-fold-counts-2026-09-05/fastra_d1_2026_09_05/preserved_hopf.jsonl` — 6 records driving beta down to keep the large fold pair while adding the innermost Hopf cycle.
- `/Users/scottg/Claude_all/H16P_branches/codex_fastra-d1-fold-counts-2026-09-05/fastra_d1_2026_09_05/preserved_hopf.log` — 4 lines, beta -3.56e-13, -2.23e-15, -1.48e-17, -1.02e-19, all USU.
- `/Users/scottg/Claude_all/H16P_branches/codex_fastra-d1-fold-counts-2026-09-05/fastra_d1_2026_09_05/preserved_hopf_final.log` — 2 lines, beta -7.25e-22 and -5.24e-24, still USU. See NEW FACT D9.
- `/Users/scottg/Claude_all/H16P_branches/codex_fastra-d1-fold-counts-2026-09-05/fastra_d1_2026_09_05/repair_counts.py` — The targeted precision-repair driver. Contains a SECOND four-origin trigger: writes FOUR_ORIGIN_TRIGGER_VALIDATED.json and raises 'STOP FOUR' if a repaired profile has >=4 roots. Never fired.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-d1-2026-09-05/fastra_d1_2026_09_05/replay_example.py` — Reproduces the branch-B rounded rational 3+1 field with BOTH log-polar and Cartesian returns; writes rational_3_plus_1.json with status NUMERICAL_3_PLUS_1_NOT_A_COUNTEREXAMPLE.
- `/Users/scottg/Claude_all/H16P_branches/codex_fastra-d1-fold-counts-2026-09-05/fastra_d1_2026_09_05/resume.py` — Resumes run_sheets.py by exec'ing its definition prologue, preserving the original driver as the experiment plan.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-d1-2026-09-05/fastra_d1_2026_09_05/run.log` — 132 lines of the initial D1 continuation: per field c, m, beta, origin stability list and edge. Shows pair SU -> hopf USU throughout the positive_K family, and origin [] for both negative_extension fields.
- `/Users/scottg/Claude_all/H16P_branches/astra_fastra-d1-2026-09-05/fastra_d1_2026_09_05/run_d1_initial.py` — The original D1 driver. Contains the trigger: appends to triggers.jsonl and raises RuntimeError('FOUR_ORIGIN_TRIGGER_REQUIRES_REVIEW') when the origin nest has >=4 roots. Uses beta trial -sign(K)*min(1e-6, |K|*1e-5).
- `/Users/scottg/Claude_all/H16P_branches/codex_fastra-d1-fold-counts-2026-09-05/fastra_d1_2026_09_05/run_sheets.py` — The codex-branch three-sheet continuation driver. Gives the exact K ladders per sheet and the beta choice b = -sign(K)*min(margin*sqrt(m)/(20 pi), 5e-4 sqrt(m)).
- `/Users/scottg/Claude_all/H16P_branches/codex_fastra-d1-fold-counts-2026-09-05/fastra_d1_2026_09_05/shooting.jsonl` — 274 half-map shots: 268 NUMERICAL_TWO_HALF_PASSAGES, 6 UNRESOLVED (3+1 12-second timeouts, 2 focus/section gate). 223 s wall total.
- `/Users/scottg/Claude_all/H16P_branches/codex_fastra-d1-fold-counts-2026-09-05/fastra_d1_2026_09_05/sign_disagreement_checks.log` — 5 lines - the 5 disagreeing labels re-profiled in binary128; all resolve consistently (3 positive_center USU, 2 negative SUS).
- `/Users/scottg/Claude_all/H16P_branches/codex_fastra-d1-fold-counts-2026-09-05/fastra_d1_2026_09_05/start_fold.json` — The seed fold for all three codex sheets: r=6.75939506, c=0.96888848, K=0.001953125, alpha=-37.1188, F=1.26e-22, multiplier 1.0000000000000000000018.
- `/Users/scottg/Claude_all/H16P_branches/codex_fastra-d1-fold-counts-2026-09-05/fastra_d1_2026_09_05/summarize.py` — Builds sign_map.csv, root_ledger.csv, summary.json AND float64_field_collapses.json; reports max_selected_origin_count, max_beta_zero_pair_count, max_hopf_count, raw_double_profile_max, raw_baseline_max.
- `/Users/scottg/Claude_all/H16P_branches/codex_fastra-d1-fold-counts-2026-09-05/fastra_d1_2026_09_05/trigger_resolution.json` — The trigger verdict: REJECTED_NUMERICAL_FALSE_POSITIVE, label positive_infinity_019_fold, raw_apparent_count 6, suspect_brackets 6, all rechecked matching endpoints positive.
- `/Users/scottg/Claude_all/H16P_branches/codex_fastra-d1-fold-counts-2026-09-05/fastra_d1_2026_09_05/verify_precursor.log` — The 4 exact-rational 3+1 roots with brackets and multipliers, then a Traceback (numpy ufunc exp on an mpf) that aborted the field_record supplement.
- `/Users/scottg/Claude_all/H16P_branches/codex_fastra-d1-fold-counts-2026-09-05/fastra_d1_2026_09_05/verify_precursor.py` — Builds the exact rational precursor c=9688912553490597/10^16, m via K=1/512, beta=-1/10^7, and brentq-brackets its four cycles in binary128 using .matching_quad_v2.so.
- `/Users/scottg/Claude_all/H16P_branches/codex_fastra-d1-fold-counts-2026-09-05/fastra_d1_2026_09_05/write_report.py` — Renders FASTRA_D1_REPORT from summary.json and verified_precursor.json.
- `/Users/scottg/Claude_all/H16P_branches/astra_q4-determinant-2026-09-05/q4/seventh/a1_determinants.json` — The a=1 exact face determinants. Confluent face: determinant_M numerator Q(M,t) degree 4 in M; endpoint face: 27 F^3 t^2 (M-1)(t-1)(M^2 t + 4Mt - 3M - 2)/1600. Positivity polynomial in (v,z) has 29 terms, ALL COEFFICIENTS POSITIVE.
- `/Users/scottg/Claude_all/H16P_branches/astra_q4-determinant-2026-09-05/q4/seventh/a1_determinants.py` — The sympy script proving it, with the substitution t=z/(1+z), M=(1/6 + v(1+z)/(6+z))/(1+v) and the assertion that all 29 coefficients are positive.
- `/Users/scottg/Claude_all/H16P_branches/astra_q4-determinant-2026-09-05/q4/seventh/explore_upper.py` — Bounded first-Duhamel upper-bound diagnostics (r in .825/.99/.99999, both faces, a in .875/.99/1). NUMERICAL_ONLY; its output json is NOT on the branch.
- `/Users/scottg/Claude_all/H16P_branches/astra_q4-determinant-2026-09-05/q4/seventh/finite_basis.json` — The 3-dimensional closed finite-a subspace: functions F, Ft, Gt(1-t) with exact coefficient rows, and basis determinant 2263040(486a^3-441a^2-486a+236)/81.
- `/Users/scottg/Claude_all/H16P_branches/astra_q4-determinant-2026-09-05/q4/seventh/finite_basis.py` — Exact sympy verification of the forcing identities, six centre data, and that determinant; asserts P(1) = -205.
- `/Users/scottg/Claude_all/H16P_branches/astra_q4-determinant-2026-09-05/q4/seventh/parameter_connection.json` — EXACT_IDENTITY_ONLY. The connection U, the scalar c = (1458a^4-2871a^3+1809a^2-222a+236)/(2a(a-1)(486a^3-441a^2-486a+236)) and h = (t-1)/(a(1-a)).
- `/Users/scottg/Claude_all/H16P_branches/astra_q4-determinant-2026-09-05/q4/sixth/boundary_diagnostic.json` — 40 frozen boundary determinants: nonnegative_D_count = 0. Independent 45-digit quadrature control at r=.95, a=1, endpoint face agrees with the ODE to 7.1e-18.
- `/Users/scottg/Claude_all/H16P_branches/astra_q4-determinant-2026-09-05/q4/sixth/boundary_diagnostic.py` — The diagnostic script; integrates the determinant in log time rather than reconstructing it by subtracting large products.
- `/Users/scottg/Claude_all/H16P_branches/astra_q4-determinant-2026-09-05/q4/sixth/determinant_exploration.json` — 60 two-anchor determinant samples over r in {.25,.6,.9,.99,.9999} x gap {.01,.5,.99} x a in {.5,.875,.95,.999}. first_gate PASS COUNT = 0. K ranges -1.017e-5 to -6.58e-6 - NEVER POSITIVE. P_B>0 in 42 of 60.
- `/Users/scottg/Claude_all/H16P_branches/astra_q4-determinant-2026-09-05/q4/sixth/environment.json` — Python 3.12.13, Linux 6.18.35, mpmath 1.3.0, numpy 2.3.5, scipy 1.17.0, sympy 1.14.0, base_commit 3b94f34, numeric_status 'diagnostic only'.
- `/Users/scottg/Claude_all/H16P_branches/astra_q4-determinant-2026-09-05/q4/sixth/exact_checks.json` — q4/sixth exact identities + rational interval certificates. See NEW FACT Q1.
- `/Users/scottg/Claude_all/H16P_branches/astra_q4-determinant-2026-09-05/q4/sixth/explore_determinant.py` — The floating exploration script; prints 'PASS FIRST GATE' if P_B>0 and det>0 - never printed.
- `/Users/scottg/Claude_all/H16P_branches/astra_resonant-joint-2026-09-05/resonant/check_resonant.py` — The resonant exact-identity script: verifies the base Hamiltonian H=(x^2+y^2+1/4)/y, the weighted divergence, the resonant jet 2(b-1), F3 = -4 pi (b-1) e1/sqrt(2-b-e1^2), the endpoint division D_1, the logarithmic generator J(h)=2pi[(h-1)-2log((h+1)/2)], its s->0 limits, and a 70-digit quadrature of three moments.
- `/Users/scottg/Claude_all/H16P_branches/astra_resonant-joint-2026-09-05/resonant/data/exact_checks.json` — Its output. Quadrature errors 0 to 4.7e-64. The (p,q,k)=(1,3,20) two-root Melnikov control: u roots 0.20848 and 2.00752, h roots 1.41697 and 5.01505, F signs at u=0,1,10 = -1, +1.137, -25.796. Limitations: no interval cycle certificate, no global bound for H(2), pure endpoint cyclicity NOT computed.
- `/Users/scottg/Claude_all/H16P_branches/astra_resonant-joint-2026-09-05/resonant/data/shoot_control.json` — 18 ORIGINAL-FIELD shooting rows. See NEW FACT R1.
- `/Users/scottg/Claude_all/H16P_branches/astra_resonant-joint-2026-09-05/resonant/shoot_control.py` — The original-field control: family a=-1, b=1-20t, e0=t^2, e1=-t, e2=2t+6t^2; four compact h-brackets 1.2/1.8/4/6 and two lower-endpoint starts s/t=0.5,1.5; refuses to proceed if the sign changes between rtol 2e-11 and 2e-13.

---

## 2. NEW FACTS — not in H16P_SUMMARY.md

### A. The centerward positive-sheet sign map (`astra_afternoon_2026_09_05/`)

The summary's §2(a) Engine 5 describes `full_return128.cpp` as "the best-engineered
numerical artifact in the repository" but never says **what it was used for**. It
was used for a completed 24-field experiment that closes the centerward sheet
numerically. This is the campaign's most systematic exclusion of a third origin cycle.

**A1 — 24 exact rational KKL fields, each with EXACTLY two origin cycles, and
nothing else out to the chart boundary.** `run_center_sign_map.py` walks
K = 1e-10, 2e-10, 5e-10, … 0.0012, 0.001953125 (24 values), Newton-corrects the
binary128 fold at each, offsets `c` by `-0.02·exp(backward_log_sensitivity)·G_z/F_c`
to open the pair, then profiles 21 points at spacing 0.05 in log-radius around the
fold and bisects each bracket **27 times**. Result in `center_sign_map.jsonl`:

- `nroots` is **2 for all 24 fields** (`Counter({2: 24})`), never 1 and never 3;
- all **816** return evaluations across local profiles and tails are
  `OK_NUMERICAL` — zero failures;
- inner radius runs 5.5304134744 → 5.5316222081 and outer 8.2378533486 →
  8.2396594286 as K rises over five orders of magnitude — a remarkably rigid pair;
- outside the outer root the displacement is **positive at every one of the 24
  fields** and stays positive through the whole tail grid `[3,4,6,8,12,16,20,24,28,32,36,40]`.

`summarize_and_check.py` re-asserts this independently: 24 K values, 48 tight root
brackets, 24 outside/grid-edge agreements, 24 outside/resolved-edge agreements,
`"sign_disagreements": 0`, `"unresolved_sign_comparisons_after_budget_rechecks": 0`.

**A2 — The return domain reaches log-radius 112 to 629, not 40.** The summary
records the engine's `umax=45` / `e^36` domain caps as a hard limit. They are a
limit of the *double* engine only. `domain_edge_map.csv` bisects the true
binary128 edge for the same 24 fields:

```
0,0.0000000001,1,1,627.7109375,628.7109375,629.296875,ANGULAR_CHART_UNRESOLVED,99.441763230575766230318523701746767,agree
23,0.001953125,1,1,111.5390625,112.5390625,112.65625,ANGULAR_CHART_UNRESOLVED,94.723244936362276287635491914210307,agree
```

The last-success log-radius decreases monotonically from **628.71** at K=1e-10 to
**112.54** at K=0.001953125. The displacement at the edge is ~+94 to +106 in log
units — enormous and unambiguously positive. The first failure is
`ANGULAR_CHART_UNRESOLVED` in **all 24 cases** (never an evaluation limit), and
`extend_domain_edge.py` deliberately backs off one log unit before reading the
sign, with the comment *"numerically resolved returning endpoint; failure beyond
it is not a proof of nonreturn"*. Across the 24 edge searches, 202 probes returned
`OK_NUMERICAL` and 136 `ANGULAR_CHART_UNRESOLVED`.

**Consequence for a fifth-cycle hunt: on this sheet there is no third origin
cycle anywhere below log-radius 112, and no "the engine just could not see out
there" excuse remains.** `domain_endpoint_checks.jsonl` pushes two of them
further still: K=1e-10 returns cleanly at log-radius 60, 100, 200 and **500**
(log-displacement +0.4768, radial displacement 8.6e216) and only hits
`EVALUATION_LIMIT` at 1000.

**A3 — The exact rational vectors are published.** `center_sign_map.csv` carries
the full 12-vector for each K. Example (index 23, K = 1/512):

```
["0","0","1","1","1","0","0","-2100097656250000000000000000000000000000000000/56578038088396572996137790707919064106286049","0","-10","11/5","9688912553490597545103435518901733100571459/10000000000000000000000000000000000000000000"]
```

Note this is the **same c** as the branch-A exact 3+1 field
(`c = 9688912553490597/10^16`) carried to 43 digits: the D1 3+1 seed and the
sign-map's K=1/512 row are the same point of the fold surface, and the sign map
proves that at that point there are exactly two β=0 origin cycles out to log-radius 112.

**A4 — The six analytic controls, with numbers.** `full_return128_validation.json`:

| control | target | max absolute error |
|---|---|---|
| `near_integrable_1e_minus17` (ε=1e-17, log-r −20/0/3) | 6.283185307179586e-17 | 4.47e-32 |
| `beyond_exp36` (log-r 37/44/60, ray 1.2345678901234568) | 6.283185307179586e-3 | 4.57e-31 |
| `large_coefficient_1e14` (M=1e14, scale 1e7) | 6.283185307179586e-10 | 3.01e-30 |
| `exact_center_scaled_nonzero_ray` (M=100, ray 1.7) | 0 | **exactly 0** |
| `nonlinear_exact_reversible_center` | 0 | 3.06e-28 |
| `chart_failure_has_null_displacement` | — | returns `ANGULAR_CHART_UNRESOLVED`, displacement `null` |

The reversible-centre control is the **Proposition-A counterexample field itself**
(`["0","0","-1","0","1","1","0","1","0","1","1","0"]`), i.e. the D2 audit's
counterexample doubles as the engine's nonlinear-centre control.

**A5 — The exact Proposition A output.** `proposition_a_exact.json` gives the
pieces the summary's R12 paraphrases but does not quote:

- candidate equilibrium `x = -a/(a²-b(l+1))`, `y = (l+1)/(a²-b(l+1))`;
- `P` at that point `= (-a²(b+2l+1) + (b+1)(l+1)²)/(a²-b(l+1))²`;
- the **uncancelled** numerator `-b²(l+1)(a²(b+2l+1)-(b+1)(l+1)²)` — kept
  separately so the cancellation is not mistaken for the invariant;
- the `l = -1` degenerate polynomial `N = b³ - 3b² - b m² + m² + 4`, with its own
  reflection certificate;
- the counterexample's three equilibria `(-1,0), (0,0), (0,1)` with traces
  `-1, 0, +1` — i.e. the neutral equilibrium is the *origin itself*, which is why
  the "existential" reading of the converse fails.

### B. Lanes and parameter regions the summary does not name

**B1 — `push_fold.py` is an entire undocumented lane.** It is the only piece of
code in the repository that attacks the fifth cycle *directly* rather than by
sampling. It loads the best `total ≥ 4` record from an evolve ledger, then runs
40 steps of finite-difference gradient descent on the seven free coefficients
`DIMS=[3,4,7,8,9,10,11]`, minimising the **interior minimum of D/r strictly
between the 2nd and 3rd origin roots**, with a line search over step
`0.02 → 0.02/16`, aiming to push that minimum through zero and birth a fourth and
fifth origin cycle:

```python
    if m<0:
        print("CROSSED ZERO: interior minimum negative -> new pair"); break
```

**There is no `data/push_fold_result.json` on the branch and no log.** The lane
either never ran to completion or its output was never committed. For anyone
hunting a fifth cycle this is the single most directly relevant unfinished piece
of code in the repository.

**B2 — Cloud-worker sweep sizes exceed everything the summary reports.** §7.2
describes the four `fable/compute-*` branches only as "raw sweep `.jsonl` output".
Their terminal histograms:

| lane | file | sets | histogram | max |
|---|---|---|---|---|
| F3 λ=0 Shi, L=8 | `F3_lam0_L8_worker.log` | **800000** (DONE, 6390 s) | `{0: 1307502, 1: 23009, 2: 133}` | 2 |
| F3 λ=0 Shi, L=2 | `F3_lam0_L2_worker.log` | 94208 (partial) | `{0: 153852, 1: 2854, 2: 19}` | 2 |
| F5 Shi full, L=4 | `F5_shi_full_L4_worker.log` | 264192 (partial) | `{0: 408772, 1: 30835, 2: 18}` | 2 |
| kklx L=3 | `W_kklx_L3.log` | 57344 (partial) | `{0: 63135, 1: 27938, 2: 1085, 3: 24}` | **3** |
| mv L=1 | `W_mv_L3.log` | **100096** (DONE) | `{0: 201865, 1: 602, 2: 1}` | 2 |
| mv-pert ε=1e-3 | `W_mvpert_e2.log` | **60160** (DONE) | `{0: 100660, 1: 19544, 2: 113}` | 2 |
| q3r-pert ε=1e-3 | `W_q3rpert_e2.log` | **60160** (DONE) | `{0: 58678, 1: 1482}` | **1** |
| q3r-pert ε=1e-3 (seed 2) | `W_q3rpert_e4.log` | 60160 (DONE) | `{0: 58694, 1: 1466}` | 1 |
| q4-pert ε=1e-3 | `W_q4pert_e2.log` | 60160 (DONE) | `{0: 119550, 1: 770}` | 1 |
| q4-pert ε=1e-3 (seed 2) | `W_q4pert_e4.log` | 60160 (DONE) | `{0: 119517, 1: 803}` | 1 |

The F3 λ=0 L=8 sweep at 800000 sets is **2.7× larger than the `F3_lam0_L4`
run of 300032 sets** the summary quotes, and it still tops out at 2. The four
`pert` lanes (q3r and q4, two seeds each, 240640 sets total) top out at **one**
cycle — the perturbative attack on the Q3-reversible and Q4 centres produces
essentially nothing, a fact absent from the summary.

**B3 — The evolve lanes and their ceilings, per family.** `DONE` lines:

- `W_evolve_kkl_all` / `_wide`: **best_total 3**, elite scores 3.981 / 3.791.
  The kkl seed here is a *pure origin* `[0.626, 2.029, 14.399] SUS` with **no
  remote cycle**, and 60 generations at both σ=0.03 and σ=0.15 never added one.
- `W_evolve_shi_all`: best_total 4, elite score exactly 4.0. Seed
  `([7.0968] S) + ([0.00414, 0.00613, 0.06434] USU)`.
- `W_evolve_yz_all` / `_wide`: best_total 4, elite score 4.846413315779456,
  identical in both. Seed `([241.534] S) + ([0.0294, 0.0643, 0.0952] USU)`.
- `F4_kkl_evolve` (30 gens) / `F4_kkl_evolve2` (60 gens): both stall at
  `best_total 4`, scores 4.998215940519204 and 4.998177894069219. The elite
  ledger `F4_kkl_evolve.jsonl` has **351 records, every single one total=4** —
  the descent reached four immediately and then spent 30 generations moving
  within the four-cycle stratum.
- `F14_kklx_seeds_evolve3.jsonl`: **77 records, every one `total=4, score=4.0`**
  across 53 generations. The score never moved at all.
- `F20_kklstar_evolve`: 24 seeds, best **score** 5.50 held constant for all 40
  generations, `DONE best_total 4`. **The 5.50 is a score, not a cycle count** —
  a reader skimming the log will misread it.

**B4 — F15 targeted, the fine structure.** `F15_targeted.jsonl` is 270 records
over `b ∈ {0.5,1.0,1.5} × da ∈ {−0.02,−0.005,0,0.005,0.02} × mag ∈ {1e-4,1e-3,1e-2}`
× 6 random unit directions in `(e₀,e₁,e₂)`, on the MV neutral two-centre family
`[(b−2)/4, e₁, 1−b, −1+da, e₂, b, e₀, 0,0,0, −2, 0]`. Totals
`{0: 213, 1: 52, 2: 4, 3: 1}`. **79% of this targeted region has no cycle at all.**

**B5 — F18/F18b/F18c: the "no single nest ever held three" claim is provable from
the file sizes.** All three lanes write a record only when
`max nest ≥ 3 or total ≥ 4`. `F18_third_cycle.jsonl`, `F18b_alien_points.jsonl`
and `F18c_scaled.jsonl` are all **0 bytes**, and there is not one `HIT` line in
any of the three logs, over 3600 + 2400 + 4536 = **10536 fields**. The
`total: 3` entries in the histograms are therefore all 2+1 or 1+1+1 splits across
two nests. F18c is the important one: it scans the *scaling* the theory demands
(`e₀ ~ δ²`, `da ~ δ`, from `D ≈ πe₀ − πδy + c₂y² + da·δ·y·log y`) at
`a₀ ∈ {−1, −1/2}`, `b ∈ {0.5,1,1.5}`, `e₁ = ±0.01`, `δ/e₁ ∈ {±0.1,±0.01,±0.001}`,
`k₀ ∈ {0,±0.3,±1,±3,±10}`, `k_d ∈ {0,±0.3,±1,±3}` — and still never gets three
in one nest.

**B6 — F19: two runs, and the targeted construction the summary conflates.**
There are two distinct F19 runs with different search strategies.

- `F19_q4_alien.log` (first run): 70 ovals, span **dim 4**, singular values
  `[20.276069, 6.585868, 1.508192, 0.087302, 0, 0]`. Random sampling of 300000
  directions gives `{0: 294422, 1: 5531, 2: 47}` — it **never found a 3-zero
  direction at all**, so only 2-zero directions were shot. Every shot returned
  exactly 2 actual origin cycles.
- `F19_q4_alien_fixed.log` / `F19b_q4_alien.log` (second run): adds a *targeted*
  construction — pick 3 oval indices, take the null vector of the 3×4 basis
  submatrix — and reports
  `targeted: 3-zero directions 15  4-zero directions 0` over 20000 triples.
  The script prints `FOUR ZEROS at first order:` whenever `z ≥ 4`. **That line
  appears nowhere.**
- The 45 shots at the 15 three-zero directions give actual origin counts
  `{3: 54, 2: 42, 1: 28, 0: 2}` across 126 records — i.e. the actual count is
  **at most, and often below, the first-order count**. The `ALIEN?` flag
  (printed when actual > first-order) never appears. `f19_directions.json`
  makes it explicit: `"picks": {"2": [one direction], "3": [], "4": []}`.
- The Q4 field and its annulus are pinned in `f19_q4_alien.py`:
  `X' = Y + 6X² + 4XY − 2Y²`, `Y' = −X + 2X² + 8XY − 2Y²`, boundary
  `x* = 0.2272111321`, ovals sampled to `x₀ = 0.22718841098679`.

**B7 — F12/F12b: the KKL double-centre Melnikov ceiling is 3, and it is rare.**
`c* = 0.968620633553494`, `α* = −37.136414809497211` (root of
`J(c) = 305 + 634c − 11c² − 1000c³`). Origin displacement is ≤ 4.45e-9 in
magnitude out to r = 200 (a numerically verified centre). Span dimension **3**.
Zero histograms: near field (45 ovals to x₀=60) `{0: 207401, 1: 92245, 2: 342, 3: 12}`
— **3 zeros in 12 of 300000 directions**; far field (90 ovals to x₀=3000)
`{0: 273151, 1: 124439, 2: 2405, 3: 5}` — **5 of 400000**. The five far 3-zero
directions and their zero locations are listed, e.g.
`[-0.2185,-0.765,0.9037] → x₀ = [0.8586, 1.8026, 13.0283]`, with
`min|M|/max|M|` between 5e-10 and 1e-8 — i.e. these are near-degenerate,
not robust, 3-zero directions.

**B8 — F13: the Yu-Han family, both annuli.** 19 `(a₁,a₄)` points on and near the
Yu-Han curve `a₄ = (a₁−5)/3`, with `a₁ ∈ {−1.1, −1.25, −1.5, −2, −2.5, −3, −30/7,
−5, −6, −8, −10, −15, −25}` plus ±0.5 offsets at `a₁ ∈ {−30/7, −2, −8}`, 65 ovals
each. Origin annulus: `dim = 3` everywhere, `max_zeros` 2 (18 of 19) or 1;
`phi3_zeros` `{0: 13, 1: 4, 2: 2}`; `phi2_zeros = 0` in 17 of 19. Second annulus:
`dim = 3`, `max_zeros` 2 (14) or 1 (5), `phi3_zeros = 1` in 18 of 19,
`phi2_zeros = 0` in **all 19**. Taylor rank ratios reach 2.0e-16, i.e. the 3×3
Taylor matrix really is singular on the curve — and the resulting constrained
element still has at most one interior zero.

**B9 — F16 and F17 numbers.** F16 (Dulac rank at the neutral hemicycle, three
`b` values): the base centre's residual is 1.8e-14 to 2.4e-14, and the two
*chart* directions `da`, `db` produce **pure noise** (|G| ≤ 1.1e-9 at every
sampled y) while `e₀` is O(1) (2.57 at b=0.5, 4.40 at b=1.5) and `e₁`, `e₂` are
O(y). So only three of the five unfolding directions do anything at first order —
which is exactly why the hemicycle emits two cycles and not three. F17: on the
upper annulus `max|D/ε| ≈ 1e-11` (the direction really is first-order null) and
`M₂` has **zero sign changes at all three ε** — its values are uniformly positive
(+6.5e-9 to +3.3e-7). On the *lower* annulus `max|D/ε| = 1.078` at every ε — the
direction is **not** null there — and `M₂` is uniformly negative. The
center-preserving verdict is therefore an upper-annulus statement only.

**B10 — D2 loop divergence, quantified.** `D2_loop_div.json` holds 22 focus-type
homoclinic loops with `a ∈ {−3,−2,−1.5,−1,−0.5}`, `b ∈ {−4,−3,−2,−1,0,4}`. In
every one of the 22: winding number is exactly 1.000, the divergence has exactly
**2 sign changes** along the loop, and `1 + ax + by` on the loop interior spans a
range that always **straddles values far from 1** (as low as 0.001 and as high as
4.665) — so the natural Dulac factor is nowhere sign-definite. The five tested
Dulac exponents give positive fractions ranging 0.17 to 0.74; the closest to
definite is `k = −(2l+b)/b ≈ −2.15` at `a=−1.5, b=−4` with 0.919 positive — still
8% negative. `sign(σ·η₂) < 0` holds in all 22 (σ and η₂ always have opposite
signs). The two refinement runs `D2_cross_a-1.5.log` and `D2_cross_a-2.log` both
**crashed** with `TypeError: 'NoneType' object is not subscriptable` after
"lost branch at b= 1.265625" — the neutral crossing was never located, at either `a`.

**B11 — F11 refinement shows why: trace-zero and η₂-zero coincide.**
`F11_refine_a-3.log` scans `b ∈ [1.22, 1.36]` at `a = −3`. Every row where
`trace = +0.0000000` has `eta2` of order **1e-8 or smaller** (`−2.77e-8`,
`−3.04e-8`, `−5.02e-8`, `+1.03e-8`, `+4.07e-9`, `+1.96e-9`, `+9.11e-10`,
`−2.66e-10`, `−6.25e-10`) — i.e. wherever the saddle becomes neutral, the focus
becomes a **centre**. At `b = 1.27` and `b = 1.28` the loop is simply lost
("no loop found"). This is the numerical shadow of the rigidity R12 is chasing.

### C. Facts that correct or sharpen the summary

**C1 — `data/QUEUE_DONE` is `.gitignore`d, so its absence proves nothing.**
§7.2 says: *"`data/QUEUE_DONE` does not exist — queues 1–4 never finished."*
`H16P_branches/astra_fastra-d1-2026-09-05/audit/fable_engine/.gitignore` reads:

```
audit/fable_engine/libretmap.so
audit/fable_engine/data/queue.log
audit/fable_engine/data/QUEUE_DONE
```

Both the sentinel and `queue.log` are deliberately untracked. The queues may well
have finished; the repository simply cannot say. (`data/queue2.log` … `queue7.log`
*are* tracked and *are* 0 bytes, as the summary says, because each job redirects
its own stdout.) The same file also ignores `libretmap_log.so` — so the engine
binaries are never committed and every replay must rebuild via `build_engines.sh`.

**C2 — There is a fourth "genuine bug", and it is not in `REVIEW_engine.md`'s
list.** §2(a) lists A1, A2, A3, B1, B2, B3, D1, E4, and separately notes the NaN
step-control acceptance. `counter_check/first_step.c` **instruments it directly**,
by `#include`-ing the unmodified `retmap_log.c` and printing the *first* step's
error and rejection decision. `first_step.json` for `negative_extension_1`:

```json
"first_step": {"initial_du": 27755027833338.312, "initial_dtheta": -27755027833333.39,
 "h": 0.001, "error_finite": false, "next_u_finite": false, "next_theta_finite": false,
 "original_reject_condition": false}
```

At the very first sample radius ρ = 0.001 the rates are ~2.8e13, the DP error
estimate is **non-finite**, the proposed next state is **non-finite**, and because
`NaN > 1` is `false` the original rejection test **accepts the invalid step**.
The run then grinds to the time cap and reports status 2. The control field
`center_0.0000000001` at the same radius has `error_finite: true` — so this is a
large-`m` conditioning failure, not a universal one. **`status 2` in the D1
ledgers is therefore not evidence of a dynamical return-domain edge.**

### D. The D1 counter-discrepancy ledger — contents

This is the material the task specifically asks for. It lives on
`astra/fastra-d1-2026-09-05` under `fastra_d1_2026_09_05/counter_check/` plus the
report `FASTRA_D1_COUNTER_DISCREPANCY_2026_09_05.md`. Fable's copy of the payload
(`audit/fable_engine/data/astra_missed_roots.jsonl`) is **byte-identical** to
`counter_check/missed_roots.jsonl` (86409 bytes each).

**D1 — Five fields, ten cycles the production counter missed.** All confirmed by
archived binary128 full returns and same-ray full-return brackets.

| ID | field label / kind | θ | missed ρ and stability | counter's own output |
|---|---|---:|---|---|
| C0 | `center_0.0000000001` / pair | −0.9605113314536193 | **7.9251590983 S**, **17.961187925 U** | `roots: [], stab: [], k:188, redge: 2.01e17, edge_kind: scan_cap` — **zero roots** |
| Cβ | `center_0.0000000001` / hopf | −0.9605113314536193 | 7.8749259617 S, 18.030201261 U | zero roots |
| P0 | `positive_extension_1` / pair | −1.0660159994339242 | [3.145e18, 4.038e18] S, [6.657e18, 8.548e18] U | `roots: [4454785377.297]`, one root only |
| N0 | `negative_extension_1` / pair | −0.785398163397555 | 1.9201578544e11 U, 4.2884355619e11 S | `roots: [], k: 0, edge_kind: integration_failure_bracket` |
| Nβ | `negative_extension_1` / hopf | −0.785398163397555 | 1.9201578524e11 U, 4.2884355641e11 S | zero roots |

The refined roots (`refined_coordinates.jsonl`, safeguarded full-return Newton,
4–9 iterations) with residuals:

```
center_0.0000000001 pair 1 S r= 7.92515909833610608667505136940302379  D=  8.221409746600232e-30
center_0.0000000001 pair 2 U r= 17.9611879246041517157314437686938594  D=  1.511434537318522e-26
center_0.0000000001 hopf 1 S r= 7.87492596169087785558732325911677685  D= -5.627567297499503e-31
center_0.0000000001 hopf 2 U r= 18.0302012605018682847737971510579245  D=  6.809081021992351e-27
negative_extension_1 pair 1 U r= 192015785438.105991370760811697672029 D= -1.036096357168234e-21
negative_extension_1 pair 2 S r= 428843556189.896092709279665453903628 D= -1.189642639319015e-20
negative_extension_1 hopf 1 U r= 192015785244.521468801129113012312139 D=  1.440472842199731e-20
negative_extension_1 hopf 2 S r= 428843556414.283897953356069459834715 D= -3.941720850015825e-20
```

**The consequence is stated in the report and it is the single most important
methodological fact on any branch:** *"Treat sweep negatives as lower bounds."*
Every histogram in §R14 is a lower bound with a demonstrated failure mode.

**D2 — Why C0 was missed: the noise filter is innocent; double precision is not.**
The two adjacent default grid points bracketing the S root at ρ≈7.925 were both
**actually visited** by the original run (`visited_in_original_run = True`), and
they do bracket the sign change:

| u | ρ | binary128 D | double D | status |
|---:|---:|---:|---:|---:|
| 1.8422447210178632 | 6.310688108089 | **+6.250426e-14** | −1.846456e-11 | 0 |
| 2.092244721017863 | 8.103083927575 | **−6.355235e-15** | −2.267742e-11 | 0 |

The true displacement is ~1e-14, but the double integration's own error is ~2e-11,
**three orders larger**, and it reports the *wrong sign* at the left endpoint.
The report's conclusion is exact: *"Lowering the noise filter alone cannot fix
this."* `center_tolerance.json` proves it with a four-rung ladder at the same four
grid points: at `rtol=1e-12` D is `[−1.85e-11, −2.27e-11, −3.83e-11, −4.44e-11]`
(all negative, no bracket); at 1e-13 `[−1.88e-12, −2.28e-12, −3.81e-12, −4.30e-12]`;
at 1e-14 `[−7.42e-14, −2.76e-13, −6.43e-13, −2.85e-13]`; at `rtol=2e-15`
`[+7.55e-14, +5.15e-14, −4.63e-13, +1.54e-13]` — **the signs are still not
consistent with the binary128 answer even at the tightest double tolerance.**

**D3 — Why N0 was missed: the NaN first step (see C2).** The counter's first
chunk had `status 2` at all eight radii and stopped there, so it never reached
ρ≈1.9e11. The two bracketing grid points are *inside* the default grid
(`inside_default_grid = True`) but `visited_in_original_run = False`. Direct
replay at those points gives binary128 D of `−3.186e-4` and `+2.291e-4` against
double `−1.439e-2` and `−2.318e-2` — **the root signal exceeds the 1e-10 noise
filter by six orders of magnitude. This is not a noise-floor rejection.**

**D4 — Coefficient rounding is NOT the explanation (for C0 and N0).**
`coefficient_rounding.jsonl` re-runs binary128 on the *double-rounded* coefficient
vector:

| field | u | binary128 on rounded field | binary128 on exact field |
|---|---|---|---|
| center_0.0000000001 | 1.8422447210178632 | +6.256704007953282e-14 | +6.250426306771985e-14 |
| center_0.0000000001 | 2.092244721017863 | −6.272441549360039e-15 | −6.355235294190589e-15 |
| negative_extension_1 | 25.842244721017863 | −3.186336431472478e-4 | −3.186336431472028e-4 |
| negative_extension_1 | 26.092244721017863 | +2.291210510379819e-4 | +2.291210510380339e-4 |
| positive_extension_1 | 42.59224472101786 | **UNRESOLVED** | +7.001135533782749e-5 |
| positive_extension_1 | 42.84224472101786 | **UNRESOLVED** | −3.559553123465363e-4 |

The sign brackets survive rounding for C0 and N0 to ~5e-17 in log displacement.

**D5 — P0 is different, and the report says so.** The exact field's transverse pair
offset in `c` is ≈**6.1e-25**, while conversion to double perturbs `c` at the
**1e-17** scale. *"the double vector does not faithfully specify this near-fold
rational field."* `rounded_half.jsonl` confirms it: on the rounded field the
half-map residual is **positive at all three original bracket points**:
`F = 8.44003705954692e-18`, `8.44003676144864e-18`, `8.44003684548544e-18`
(versus the exact field's splitting near 1e-26). The saved sign brackets are
destroyed by rounding. Both P0 roots also lie beyond the `u < 40` scan cap
(at u ≈ 42.6–43.6), so no genuine default grid point brackets them, and both
binary128 checks on the rounded field hit the step-resolution guard — recorded as
**unresolved, explicitly not a no-return proof**.

**D6 — The counter settings actually used** (from the report's own table, all
defaults, unchanged from Fable commit `afbcdd419309e30222e494e075c45b3049350020`):
initial grid radius `1e-3·scale` with scale capped at 1; u-grid
`np.arange(log(0.001), 40, 0.25)`; last actual u `39.84224472101786`
(ρ ≈ 2.010328497117133e17); `rtol 1e-12`, C `atol 1e-15`; sign threshold
`min(|D_left|,|D_right|) > 1e-10`; `umax=45`; `Smax=2000`; `maxsteps=300000`;
initial `h=1e-3`; 8 radii per chunk, **stop at the first nonzero status**;
8 bisections of edge refinement, only when a previous valid point exists; ray
chosen away from the nearest other equilibrium. The older supplemental D1
profiles used a *different* configuration — spacing 0.125 on log-radius
`[−25, 46]`, horizontal ray, rtol 1e-12 **and** 1e-13, `umax=60`, `Smax=10000`,
`maxsteps=1000000` — and 32 further evaluations used rational coordinate scaling.

**D7 — Reproduction is a two-liner** (verbatim from the report):

```bash
g++ -O2 -std=c++17 -fext-numeric-literals fastra_d1_2026_09_05/counter_check/full_ray_quad.cpp -o /tmp/d1_full_ray -lquadmath
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 python fastra_d1_2026_09_05/counter_check/check_misses.py
```

The audit cost 178 binary128 full-return calls, **all `NUMERICAL_ONLY`**, 388.5
CPU-seconds total (`calls.jsonl`).

**D8 — Distinct exact fields collapse to the same float64 vector.**
`float64_field_collapses.json` lists **five** labels — `positive_infinity_021`,
`022`, `023`, `024_polished`, `025_polished` — where the exact *fold* vector and
the exact *pair* vector, which differ by the deliberate `c` offset that opens the
cycle pair, round to the **identical** binary64 vector. Anyone re-running the
double engine on the "pair" field at these points is silently integrating the
**fold** field instead. This is generated automatically at the end of
`summarize.py` and has no analogue anywhere in the summary.

**D9 — The β ladder: the innermost Hopf cycle survives to β = −5.24e-24.**
`preserved_hopf.log` + `preserved_hopf_final.log`, on the far-positive sheet:

```
positive_infinity_018_hopf_preserved -3.563520345169391e-13 USU
positive_infinity_019_hopf_preserved -2.2323678901476827e-15 USU
positive_infinity_020_hopf_preserved -1.4832349336726276e-17 USU
positive_infinity_021_hopf_preserved -1.0235510170637541e-19 USU
positive_infinity_022_hopf_preserved -7.252680515205544e-22 USU
positive_infinity_023_hopf_preserved -5.240197753710807e-24 USU
```

β is chosen as `−0.1·|F₀/F_β|` from a measured β-derivative of the *bounded
matching residual*, which is why it can be pushed 11 orders of magnitude while the
stability word stays `USU`. But `preserved_hopf.jsonl` also records each row's
double-precision baseline: for 018 it is `[(0,''), (2,'SU')]` and for 019–023 it
is `[(0,''), (1,'S')]` — **the remote cycle is gone at these K**. So the far-positive
sheet buys a third origin cycle at the price of the remote one: still 3, never 5.

**D10 — Three more four-origin triggers exist, and none fired.** §7.5 names four
live triggers. There are at least three additional ones in this file set, all
independent:

- `run_d1_initial.py`: appends to `triggers.jsonl` and raises
  `RuntimeError('FOUR_ORIGIN_TRIGGER_REQUIRES_REVIEW')`;
- `repair_counts.py`: writes `FOUR_ORIGIN_TRIGGER_VALIDATED.json` and raises
  `RuntimeError('STOP FOUR')` if a *repaired binary128* profile has ≥4 roots;
- `finish_positive.py` and `finish_positive_extension.py`:
  `if len(p['roots'])>=4: raise RuntimeError('STOP FOUR ORIGIN')`.

The only one that ever fired is the `engine.py:116` trigger visible as a
traceback at the end of `positive_infinity.log`, and `trigger_resolution.json`
records the verdict:
`{"status": "REJECTED_NUMERICAL_FALSE_POSITIVE", "label": "positive_infinity_019_fold",
"raw_apparent_count": 6, "suspect_brackets": 6, "all_rechecked_matching_endpoints_positive": true}`.

**D11 — The three continuation sheets, with their exact endpoints.**

| sheet | accepted folds | K range | fold radius range | endpoint c |
|---|---|---|---|---|
| `positive_center` | 10 | 1.953125e-3 → **1e-10** | 6.75939506 → 6.75794352 | 0.968620633567205247249 |
| `positive_infinity` | 25 (+1 `CORRECTOR_UNRESOLVED` at log-r 42) | 0.00390625 → **7.069387071467069** | 6.76085 → **1.7393e18** | 1.593702902627309087 |
| `negative` | 29 (+1 `CORRECTOR_UNRESOLVED`) | −1e-4 → **−2.5975e12** (chart switches K→m at index 15) | 6.75787 → 4.0489e10 | — |

The positive-infinity ladder is worth quoting in full because it is the campaign's
longest continuation and the summary does not give it: K = 0.00390625, 0.015625,
0.0625, 0.125, 0.2242157839673 (**where c crosses exactly 1**), 0.5, 1, 1.5, 2, 3,
4, 5, 5.5, 6, then by log-radius z = 8, 10, 12, 16, 20, 24, 28, 32, 36, 40, 42.
Fold radius grows 6.76 → 1.74e18; c grows 0.9692 → 1.5937 and appears to be
converging to a limit near 1.5937–1.5938 as K → ~7.07.

**D12 — Double vs binary128 disagreement on the sheets, in aggregate.** From the
astra branch `fields.jsonl` (176 double-precision counter evaluations, 44 labels ×
4 kinds): origin-nest counts `{0: 60, 3: 46, 2: 40, 1: 30}`; **totals across all
nests `{1: 52, 3: 40, 0: 36, 4: 32, 2: 16}` — 32 of 176 evaluations reported
four**. Nest edge kinds: **200 `integration_failure_bracket` vs 152 `scan_cap`** —
i.e. more than half of all nests ended in an integration failure, not a scan cap.
The 88 dense re-profiles (`dense.jsonl`, spacing 0.125 on log-r `[−25,46]`, two
tolerances) top out at **3 brackets**, with stabilities
`SU 20, US 16, USU 14, '' 14, SUS 14, S 10`, edge kinds 54 `scan_cap` / 34
`integration_failure`, and **zero `sign_mismatch` rows**. On the codex branch,
`precision_repairs.jsonl` (28 binary128 re-profiles) gives root counts
`{3: 19, 2: 9}` and stabilities `USU 17, SU 9, SUS 2` — **never four**.

**D13 — Cross-engine agreement on the branch-B 3+1 field.** `example.log` (one
JSON line) replays every one of the four cycles of the branch-B exact rational
field in **both** log-polar and Cartesian coordinates, at radii bracketing each
root by ±0.001 in log:

```
remote  r = 19913343.5355  S   log_D ±1.0994e-5 / ±1.0996e-5, cartesian_D_over_r ±1.098e-5 / ±1.101e-5
origin  r = 0.0695518418   U   log_D ∓4.021e-11 / ±4.025e-11, cartesian agrees to 1%
origin  r = 7.8802462384   S   log_D ±1.1233e-8, cartesian agrees to 4 digits
origin  r = 18.0424042756  U   log_D ∓2.630e-8 / ±2.641e-8, cartesian agrees to 4 digits
```

with `cartesian_status = 0` at all eight points. `audit_records.py` asserts
`log_D * cartesian_D_over_r > 0` at every point.

**D14 — Calibration fields.** `calibration.json` holds five reference points the
D1 run used to check the counter, including the **incumbent** four-cycle field
(exact vector with `m = 363889/5000`, `β = 3/2000`, `c = 7/10`) giving
`remote [6107.33] U` + `origin [0.9678, 3.3779, 24.9593] SUS` = 4; and — importantly —
the `kkl-table` field (`m = 8403125/226052`, `c = 9683/10000`, β = 0) giving only
`[1.211e10] U` + `[12929.2] S` = **2**, i.e. the published KKL table point does
not itself reproduce four cycles under this counter. The `events_logm` calibration
returns `k: 0, redge: null` for both nests — a total counting failure at
`m ≈ 2.766e12`.

**D15 — A third sandbox path.** §7.1 lists Astra's sandbox as
`/workspace/scratch/97757d9f13f6/H16P`. The codex branch tracebacks in
`positive_infinity.log`, `verify_precursor.log` and the `shooting.jsonl` timeout
errors all name a **different** one: `/workspace/scratch/16d9227a7ce0/H16P`.
`make_manifest.py` pins the provenance: `fable_input_commit '4cece20'`,
`astra_fold_input_commit '7db8597cb7d9bb34e119e85bec3f229270eaf1aa'`.

**D16 — The codex branch has 266 field records, not 176.** `artifact_audit.json`:
`{"fields.jsonl": 266, "shooting.jsonl": 274, "events_negative.jsonl": 30,
"events_positive_infinity.jsonl": 26, "precision_repairs.jsonl": 28,
"events_positive_center.jsonl": 10, "preserved_hopf.jsonl": 6,
"precision_controls.jsonl": 3}`, `"exact_vectors_validated": 303`. The astra
branch's `audit.json` reports 176 field records and 88 dense profiles. **The two
D1 branches are not the same experiment at two sizes; they are two different runs
with different drivers, and the codex one is 50% larger.**

**D17 — The controller fix, located precisely.** §2(a) mentions the `1e-24`
sensitivity floor. It is on **line 38 of `matching_quad.cpp`**:

```cpp
err=std::max(err,fabsq(tab[0][v]-tab[1][v])/((v==0?tol:std::max(tol,1e-24Q))*(1+fabsq(tab[0][v]))));
```

`matching_quad_initial.cpp` (the pre-fix version) has the same line **without**
the `std::max(tol,1e-24Q)` — a single shared `tol` for both position and
sensitivity. `matching_quad_remote.cpp`, used by `verify_precursor.py` for the
remote cycle of the exact rational 3+1 field, **also still lacks the fix**; it
compensates with the extra `x0,y0,th0,sense` fields for an arbitrary section.
`precision_controls.jsonl` is the acceptance test: `positive_infinity_024_pair`
at `tol` 2e-27 / 2e-29 / 2e-30 gives F values agreeing to four digits
(2.7203e-25, 2.7524e-25, 2.7524e-25 at the first grid point), all status 0.

### E. q4/sixth and q4/seventh — the determinant reduction

**Q1 — `q4/sixth/exact_checks.json` contains the only rational *interval*
certificates in this entire file set.** Status
`EXACT_IDENTITIES_AND_RATIONAL_INTERVAL_CERTIFICATES`. The exact identities:
corner `(A, B, η) = (11843/9623, −833/9623, 13320/9623)`, corner
`Y₀ = −81/19246`, **global infimum `q = 185/108`**, first bootstrap
`r = 988/1331`, final lower bound `r ≥ 1 − (7/22)^{3/2}`, lift polynomial
(ascending) `[−135, 440, 2584, 3936, 2816, 1024]`, root enclosures

```
u  ∈ [368493931857/2500000000000, 1473975727429/10000000000000]
κ  ∈ [2.899241080973277432530648, 2.899241080974989225893177]
```

and two confluent anchor certificates at `r = 7/10` and `r = 4/5`, each summing
128 terms, giving `q = η/(−192 Y₀)` enclosed in
`[1.858596725263160627569580, 1.858596725263160627696723]` (strictly exceeding
167/90, margin ≥ 0.00304) and
`[1.901881476015943662190798, 1.901881476023078681745372]` (strictly exceeding
19/10, margin ≥ 0.00188). The stated scope is explicit and limiting:
*"The global sign implications use the analytic notes; this script does not
certify the remaining determinant faces."*

**Q2 — The two-anchor determinant never has the sign the Q4 argument needs.**
`determinant_exploration.json`: 60 samples over `r ∈ {.25,.6,.9,.99,.9999}`,
gap `∈ {.01,.5,.99}`, `a ∈ {.5,.875,.95,.999}`. The gate is `P_B > 0 AND det > 0`.
`P_B > 0` in **42 of 60**, but **`det = K` is negative in all 60**, ranging
`−1.0168e-5` to `−6.579e-6`. `first_gate` PASS COUNT = **0**. `explore_determinant.py`
prints `PASS FIRST GATE` on success — that line appears nowhere.

**Q3 — The frozen boundary determinant is also uniformly negative, with an
independent 45-digit control.** `boundary_diagnostic.json`: 40 samples over
`r ∈ {.825,.95,.999,.99999}` × 2 faces × `a ∈ {.66,.875,.99,.99999,1.}`;
`"nonnegative_D_count": 0`. The determinant is integrated **in logarithmic time**
rather than reconstructed by subtracting large products, and the difference
between the two routes is exactly 0.0 on the sampled rows. The independent check
is a direct 45-digit `mpmath.quad` of P and Φ at `r=.95, a=1`, endpoint face:
`quadrature_D = −0.000011886801474563956523658974025596068` vs
`ODE_D = −1.1886801474556853e-05`, **absolute difference 7.10e-18**. Its own
warning is retained: *"Finite samples do not determine the sign between samples,
at other lifts, or at unsampled joint limits."*

**Q4 — `q4/seventh` proves an exact positivity statement at `a = 1` only.**
`a1_determinants.py` reconstructs the four solutions `Y_j` exactly and computes
two 4×4 determinants:

- **endpoint face**:
  `det = 27 F³ t² (M−1)(t−1)(M² t + 4M t − 3M − 2)/1600`, with the identity
  `Q_e = (M+2)(M−1) − (1−t)(M²+4M)` verified;
- **confluent face**: `det = −81 F⁴ t³ · Q(M,t)/2620618000` with `Q` an explicit
  quartic in `M`.

The payoff is a positivity certificate: under `t = z/(1+z)`,
`M = (1/6 + v(1+z)/(6+z))/(1+v)` with `z, v > 0`, the numerator becomes a
polynomial in `(v,z)` with **29 terms and every coefficient positive** (the script
asserts both `all(co>0)` and `len(P.terms())==29`). So the confluent face has a
fixed sign on the whole Stieltjes region — **at `a = 1`**. The JSON's own `scope`
field is blunt: *"a=1 only; no conclusion for all finite a follows by continuity."*

**Q5 — The finite-`a` obstruction is a cubic.** `finite_basis.json`: the closed
3-dimensional subspace is spanned by `F`, `tF`, `t(1−t)G`, and its basis
determinant is

```
2263040·(486a³ − 441a² − 486a + 236)/81
```

with `P(1) = −205`. `parameter_connection.json` shows the same cubic in **every
denominator**: `c = (1458a⁴ − 2871a³ + 1809a² − 222a + 236)/(2a(a−1)(486a³ − 441a² − 486a + 236))`,
and the connection `U` divided by `61600·a(a−1)(486a³ − 441a² − 486a + 236)`.
**The whole finite-`a` reduction degenerates at `a = 0`, `a = 1`, and at the three
real roots of `486a³ − 441a² − 486a + 236`** — and `a = 1` is precisely the value
at which the only exact positivity result was obtained. `finite_basis.json`'s
scope: *"Three-dimensional subspace only; not the complete finite-a family."*

**Q6 — `explore_upper.py` ran and its output is not on the branch.** It writes
`explore_upper.json` (first-Duhamel upper-bound diagnostics at
`r ∈ {.825,.99,.99999}` × 2 faces × `a ∈ {.875,.99,1.}`, DOP853 rtol 3e-12). That
file is **not in the branch tree**. The lane is unfinished, and it is the one that
would have turned Q4's negative determinant samples into a bound.

**Q7 — The q4 environment.** `environment.json`: Python 3.12.13 / Clang 22.1.3,
Linux 6.18.35, mpmath 1.3.0, numpy 2.3.5, scipy 1.17.0, sympy 1.14.0,
`base_commit 3b94f34`, `arithmetic_proof: "Fraction interval arithmetic and
symbolic identities"`, `numeric_status: "diagnostic only"`.

### F. The resonant obstruction (`resonant/`)

**R1 — There is a working original-field NUM control for 2 compact + 1 endpoint
cycle, and its numbers are stable across three ε.** `shoot_control.py` uses the
family `a = −1`, `b = 1 − 20t`, `ε₀ = t²`, `ε₁ = −t`, `ε₂ = 2t + 6t²` and shoots
half-returns in the **original quadratic field** (not the Melnikov model), at
`t ∈ {1e-3, 3e-4, 1e-4}` and two tolerances (2e-11 and 2e-13), refusing to
proceed if the sign changes between them. `shoot_control.json`, normalised by `t²`:

| start | t=1e-3 | t=3e-4 | t=1e-4 |
|---|---:|---:|---:|
| compact h=1.2 (y₀=0.931662) | **+0.2313** | +0.2338 | +0.2345 |
| compact h=1.8 (y₀=1.64833) | **−0.3062** | −0.2974 | −0.2949 |
| compact h=4.0 (y₀=3.93649) | **−0.2923** | −0.2571 | −0.2470 |
| compact h=6.0 (y₀=5.95804) | **+0.1618** | +0.2079 | +0.2211 |
| lower s/t=0.5 | **−1.4734** | −1.5419 | −1.5612 |
| lower s/t=1.5 | **+1.7979** | +1.6374 | +1.5929 |

Precision-change magnitudes are 2.9e-19 to 5.6e-16, i.e. the signs are far above
the numerical floor. The sign word on the compact side is `+ − − +` → **two
compact sign changes**; on the lower endpoint side `− +` → **one**. That is
exactly the "two compact plus one endpoint" configuration, realised in the
original field at three values of t, with `D/t²` converging as t → 0.

**R2 — The exact obstruction, in the file's own symbols.** `data/exact_checks.json`:

- persistent resonance `F₃ = −4π e₁ (b−1)/√(2 − b − e₁²)` — **it vanishes
  identically at `b = 1`**, which is the holomorphic point;
- first normal endpoint division
  `D = −π(s−2)(e₀s² + 4e₀s + 4e₀ − 2e₁s − e₂s)/(4(s+2))`;
- mixed generator
  `D_mixed = π s (s(s − 8·log((s(s+4)+4)/(8s)) − 4) + 4)/(2(s² − 4))`, with
  `lim_{s→0⁺} D_mixed/s = −π/2` and `lim_{s→0⁺} D_mixed = 0` (both asserted);
- the third compact generator `J(h) = 2π[(h−1) − 2 log((h+1)/2)]`;
- 70-digit quadrature of three moments at `h ∈ {1.02, 1.5, 3, 10, 50}` with
  errors from **exactly 0.0** to 4.75e-64, and the mixed-derivative identity
  matching to ≤ 4.34e-66.

**R3 — The positive control's exact roots.** `(p,q,k) = (1,3,20)` with
`F(z) = 3 − 4(1+z) + 20(1 − log(1+z)/z)`:
`u ∈ {0.2084839277105037018877434993403395755138909505749561593706871433598084,
2.007524695383949113253606375735464955775506168659705564428855142266509}`,
`h = 1 + 2u ∈ {1.416967855421007…, 5.015049390767898…}`, with
`F(0) = −1`, `F(1) = +1.137056…`, `F(10) = −25.795790…`.

**R4 — The stated limitations, verbatim.** `exact_checks.json` ends with
`"limitations": ["No interval cycle certificate.", "No global bound for H(2).",
"Pure endpoint cyclicity is not computed by these checks."]` and
`shoot_control.json` is labelled `"NUM only; no interval certificate"`.
**Pure endpoint cyclicity — the case of five cycles all tending to infinity that
§R13 explicitly leaves open — is not merely unresolved; it is not computed by any
file on this branch.**

### G. Practical notes for a fifth-cycle hunt

1. **The single most reusable asset here is the D1 chart, written out explicitly**
   in `counter_check/full_ray_quad.cpp`:
   `P = y + x² + xy`, `Q = −m·x + β·y − 10x² + (11/5)xy + c·y²`, with
   `α = −m` and `K = m(11c−5)/5 − 42`. Every exact rational vector in the D1
   ledgers is `[0,0,1,1,1,0,0,−m,β,−10,11/5,c]`; `audit_records.py` asserts this
   shape for all 176 records.
2. **Never trust a sweep negative from the double engine.** D1–D5 above give a
   demonstrated, reproduced, five-field failure mode with three distinct causes
   (double trajectory error swamping a 1e-14 signal; a NaN-accepted first step
   producing a fake status-2 "edge"; and a near-fold field whose splitting is
   1e-25 and therefore invisible to a float64 vector).
3. **Never trust a "pair" field's double vector on the far-positive sheet** —
   see D8, five labels where the pair and fold vectors are the same float64.
4. **The two most-searched regions with a hard ceiling of four** are (i) the KKL*
   fold surface near `c ≈ 0.95–0.99`, `K ≈ ±0.05`, `α ≈ −38…−36`, `|β| ≤ 2e-4`
   — 10048 fields, 24 fours, every one 3+1, all recount to 4 after the noise fix
   (`RECOUNT_fixed_counter.jsonl`, `old→new = (4,4)` × 24); and (ii) descent from
   four-cycle seeds, where 351 + 77 + 34 elite records never once moved past four.
5. **The regions that never even reached two** are the q3r-pert and q4-pert lanes
   (240640 sets, max 1) and the mv lane (100096 sets, max 2). If there is a fifth
   cycle it is not a small perturbation of the Q3-reversible or Q4 centres.
6. **Three unfinished threads** that a fifth-cycle hunt would want first:
   `push_fold.py` (no committed output), `q4/seventh/explore_upper.py` (output
   json missing from the branch), and `fable_d2_crossing_refine.py` (crashes at
   `b = 1.265625` for both `a = −1.5` and `a = −2`).
