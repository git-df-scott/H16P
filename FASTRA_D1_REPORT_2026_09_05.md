# FASTRA D1 — Joint fold work, 2026-09-05

**D1 STATUS: OPEN. No four-origin field and no five-cycle candidate.**
The deterministic continuation and beta tests are completed for the grid below.
The requested global sign-map dichotomy is not justified by these computations.

This continues Astra `7db8597` and uses the unchanged Fable sources from
`afbcdd419309e30222e494e075c45b3049350020` on
`claude/conjecture-progress-report-ixsmgv`. No sweep, descent, literature audit,
or re-derivation of the section theorem was performed. No Fable-hostile-review
claim is made; the four-origin trigger never occurred.

## Requested report

FOLD SHEETS CONTINUED: positive K toward center / positive K toward infinity /
center-selected negative K. Numerical checkpoints, all on the positive horizontal
ray (x=r,y=0):

| Direction | Endpoint c | Endpoint K | Origin fold radius |
|---|---:|---:|---:|
| positive, centerward | 0.9686206335672053 | 1e-10 | 6.757943516002 |
| positive, outward | 1.593739327603797 | 7.070058783171657 | 2.1862836088e+18 |
| negative, outward | 0.3362310484231757 | -14463255200474.76 | 207137206639.2 |

These are continuation checkpoints, not proven mathematical endpoints. Every
endpoint's exact rational coefficient vector and exact K derived from that
vector appear in `summary.json`. Fold positions are numerical; rounding a fold
field to a rational vector does not make its double root exact.

MAX ORIGIN COUNT ON ANY SHEET (compactified): **3**, reproduced at the rational
field below, including its remote cycle. This is the maximum observed count,
not an upper bound throughout the sheets.

FOUR-IN-ORIGIN-NEST FIELD: **NONE found**.

SIGN-MAP VERDICT: **UNRESOLVED globally.** No extra-pair signature was found in
the resolved sampled profiles. Neither an everywhere-negative conclusion nor
a newly located second-pair region is established.

D1 STATUS: **OPEN**.

## Exact rational 3+1 field and independent replay

Coefficient order is `(P_1,P_x,P_y,P_xx,P_xy,P_yy,
Q_1,Q_x,Q_y,Q_xx,Q_xy,Q_yy)`, with ordinary monomial coefficients, not derivatives.
The reproducible field is

```
["0", "0", "1", "1", "1", "0", "0", "-37101199745401/1000000000000", "-1/25600000", "-10", "11/5", "242288563571/250000000000"]
```

Equivalently,

\[
\dot x=y+x^2+xy,\qquad
\dot y=-\frac{37101199745401}{10^{12}}x-\frac{1}{25600000}y
-10x^2+\frac{11}{5}xy+\frac{242288563571}{250000000000}y^2.
\]

Its exact parameters are
`c = 242288563571/250000000000`,
`K = 4660642062301237256681/1250000000000000000000000`,
`beta = -1/25600000`.

| Nest | Ray angle | Radius | Stability |
|---|---:|---:|---|
| remote | 2.180914664056 | 19913343.53554 | S |
| origin | -0.9606779895342 | 0.06955184177825 | U |
| origin | -0.9606779895342 | 7.880246238408 | S |
| origin | -0.9606779895342 | 18.04240427559 | U |

Both Fable engines return with matching displacement signs on both sides of
each of these four roots. `rational_3_plus_1.json` records the coefficients,
focus coordinates, rays, root estimates, endpoint signs, return statuses, and
Cartesian periods. This is numerical cross-verification, not an interval
certificate. It is a 3+1 field, not the requested 4+1 field.

## Work completed and count limitations

- 44 accepted fold checkpoints; 176 exact-rational field records.
- At every accepted checkpoint, the unchanged `sweep_log.evaluate(coef)` was called at the fold, the pair-present side, and two beta amplitudes: 176 evaluations.
- 32 additional evaluations used the exact rational conjugacy `x=X, y=sY, tau=s*t` to balance large m. These coefficient vectors are also saved.
- 88 supplemental full origin profiles used the horizontal ray, log-radius spacing 0.125 on [-25,46], two tolerances, and explicit integration-failure records. Maximum sampled count: 3.
- 30 binary128 full-return endpoint checks; 8 unresolved. All successes and failures are retained.

At a fold, an even-multiplicity zero has no sign change. Fable's unchanged
counter has spacing 0.25 in log radius, starts at radius 0.001 times its scale,
and rejects sign changes whose endpoints have displacement at most 1e-10.
Its reported roots are geometric bracket midpoints, not refined roots.
Consequently, its count at an exact fold is not a count including the double
cycle, and it can miss an entire narrow pair or a small Hopf cycle. Near the
center endpoint, binary128 full returns recover the pair that the unchanged
counter reports as absent.

At the large positive endpoint, rational coefficient rounding to double and
integration error can swamp the splitting of the two cycles. The inherited
fold radius already exceeds exp(40); the final extension exceeds it further.
The observed single stable root near 1e9 at those fields is not evidence that
the inherited pair disappeared. On the large-m negative branch the unchanged
counter can fail before its first return. Rational balancing improves some
profiles but does not repair all large-m cases. Zero returned roots after an
integration failure are **unresolved**, not evidence of no cycles.

All 132 pair-side and beta fields retain independent two-half sign brackets
for the finite pair. At the final positive endpoint, strict-tolerance full
returns hit the 25-second wall fuse. Follow-up full returns at tolerances
2e-26 and 2e-27 both recover signs +,-,+ at fold log-radius offsets
-0.8,0,+0.8; their displacement differences are at most 6.2e-6, compared with
minimum absolute displacement 4.65e-4. The final negative pair and beta field
both reproduce full-return signs -,+,-. `endpoint_full.jsonl` and
`positive_full_control.jsonl` preserve these results and the earlier failures.
This repairs the local pair verification, not global root coverage.

`edge_kind=scan_cap` means the last sampled radius, not a return-domain
boundary. `integration_failure_bracket` locates a numerical loss of return
under the stated caps, not a certified dynamical edge. No actual mathematical
return-domain edge is certified here.

## Sign map and the stopping rule

The full per-field map is `field_ledger.csv`: exact rational vector and
parameters, origin roots, stability string, section angle, edge radius and
kind, edge displacement, displacement just outside the last sampled root,
and remote roots. Raw profiles and failure codes remain in `fields.jsonl`;
tolerance comparisons are in `dense.jsonl`.

There are 0 opposite-sign comparisons between the outside-last-root sample and the recorded edge in the unchanged counter's usable records.


This does not prove that a second pair cannot be born. A saddle-node pair can
appear inside an interval while the displacement signs at both interval
endpoints stay unchanged. Thus endpoint-sign agreement cannot exclude an
unresolved even pair. It also cannot establish coverage of all sheets from a
finite parameter grid. A closure needs zero isolation/enclosures over the
relevant intervals and parameter ranges, or an additional theorem controlling
the number of extrema/zeros. Neither is supplied by the present engines.

## Continuation grid

The following table uses horizontal fold radii. Each row's exact rational
field vector is in `accepted.jsonl`, with all corrector iterations; every
pair/beta field has its own vector in `field_ledger.csv` and `fields.jsonl`.

| Checkpoint | c | K | Fold r |
|---|---:|---:|---:|
| positive_K_0.00390625 | 0.9691564318214 | 0.00390625 | 6.760850755245 |
| positive_K_0.015625 | 0.9707663765577 | 0.015625 | 6.769672344523 |
| positive_K_0.0625 | 0.9772434072956 | 0.0625 | 6.806479770357 |
| positive_K_0.15 | 0.9894828683482 | 0.15 | 6.881931021479 |
| positive_K_0.3 | 1.010849341764 | 0.3 | 7.033103093933 |
| positive_K_0.65 | 1.061831565959 | 0.65 | 7.507605837277 |
| positive_K_1.2 | 1.141141504509 | 1.2 | 8.676022186059 |
| positive_K_2.0 | 1.244351416199 | 2 | 11.65248081 |
| positive_K_3.0 | 1.346236366682 | 3 | 19.12457732358 |
| positive_K_4.0 | 1.423468746521 | 4 | 37.25591027918 |
| positive_K_6.0 | 1.537759687273 | 6 | 784.7258650092 |
| center_0.001953125 | 0.9688884792671 | 0.001953125 | 6.759395061611 |
| center_0.0001 | 0.968634344655 | 0.0001 | 6.758017733793 |
| center_0.00001 | 0.968622004651 | 1e-05 | 6.757950936625 |
| center_0.000001 | 0.9686207706631 | 1e-06 | 6.757944257392 |
| center_0.00000001 | 0.9686206349246 | 1e-08 | 6.757943522681 |
| center_0.0000000001 | 0.9686206335672 | 1e-10 | 6.757943516002 |
| events_half_0 | 1.576461970531 | 6.747476031852 | 556117.4511128 |
| events_half_3 | 1.581698645612 | 6.846096786303 | 11169917.59801 |
| events_half_7 | 1.585881135147 | 6.924332020286 | 609856836.8899 |
| events_quad_0 | 1.587328151261 | 6.95128974397 | 4506266379.957 |
| events_quad_3 | 1.590307012951 | 7.006609697245 | 1817957608823 |
| events_quad_6 | 1.592151344142 | 7.04074234868 | 7.334164447483e+14 |
| events_quad_9 | 1.593405805278 | 7.063907004368 | 2.958813114325e+17 |
| events_negative_0 | 0.9686069227328 | -0.0001000000000007 | 6.757869307583 |
| events_negative_3 | 0.9645201293489 | -0.03 | 6.736164517535 |
| events_negative_6 | 0.9418020694474 | -0.2 | 6.630177401517 |
| events_negative_9 | 0.8491037262464 | -1 | 6.437207763353 |
| events_negative_12 | 0.709609736937 | -2.999999999972 | 6.927364243725 |
| events_negative_15 | 0.5828776925127 | -7.999999999752 | 8.958660131026 |
| events_negative_17 | 0.5219418685599 | -14.99999999995 | 11.56010881593 |
| events_negative_19 | 0.4728769348207 | -29.999999994 | 16.27191083631 |
| events_m_0 | 0.46 | -37.76781889053 | 18.44210222157 |
| events_m_2 | 0.44 | -57.522580478 | 23.48466887586 |
| events_m_4 | 0.4 | -201.4531461075 | 52.44163698332 |
| events_m_5 | 0.35 | -43472.13498924 | 3523.033376474 |
| events_logm_0 | 0.3461511059136 | -201848.1452104 | 12987.04961277 |
| events_logm_3 | 0.3405161512919 | -19110424.78475 | 706258.3315338 |
| events_logm_6 | 0.3381881113001 | -1755383011.567 | 43260344.81057 |
| events_logm_10 | 0.3366860086554 | -717314104495.9 | 11970315900.32 |
| positive_extension_0 | 1.593576894534 | 7.067063091065 | 8.042887922477e+17 |
| positive_extension_1 | 1.593739327604 | 7.070058783172 | 2.1862836088e+18 |
| negative_extension_0 | 0.3364423725987 | -3221424291117 | 49665051896.66 |
| negative_extension_1 | 0.3362310484232 | -1.446325520047e+13 | 207137206639.2 |

## Replay

Run from the repository root:

```bash
bash fastra_d1_2026_09_05/build_engines.sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 python fastra_d1_2026_09_05/replay_example.py
python fastra_d1_2026_09_05/audit_records.py
```

`run_d1_initial.py` is the driver used for the archived continuation;
`run_d1.py` tightens the positive-extension corrector's residual threshold
relative to fold curvature for future continuations. No historical budget
ledger was edited. This authorized D1 run is recorded separately.
`half_beta_quad.cpp` and `full_beta_quad.cpp` extend the archived binary128
engines by the beta terms in radial/angular rates and divergence. The Fable
sources are unchanged. `source_manifest.json` supplies source hashes.

Unfinished: an exhaustive origin-root count in every accepted field; true
return-domain edge enclosures; coverage of the full fold sheets; any second
pair-birth region; a four-origin candidate; rigorous exclusion/certification.
No work is left running after this report is finalized.
