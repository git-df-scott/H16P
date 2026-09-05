"""Build a scope-honest D1 ledger from saved numerical evidence."""
import json,csv,hashlib,collections,math
from pathlib import Path
from fractions import Fraction as Q
HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
def load(name):
    p=HERE/name
    return [json.loads(s) for s in p.read_text().splitlines()] if p.exists() else []
def origin(pack):return next((n for n in pack['nests'] if abs(n['pt'][0])+abs(n['pt'][1])<1e-7),{})
fields=load('fields.jsonl');accepted=load('accepted.jsonl');dense=load('dense.jsonl');full=load('endpoint_full.jsonl')+load('positive_full_control.jsonl')
example=json.loads((HERE/'rational_3_plus_1.json').read_text())
rows=[];mismatch=[];missed=[]
for row in fields:
    n=origin(row['fable']);v=[Q(x) for x in row['coefficient_vector']]
    outside=n.get('outside_outer');edge=n.get('edge_D')
    disagree=bool(outside and edge is not None and outside['D']*edge<0)
    if disagree:mismatch.append((row['label'],row['kind']))
    checks=row.get('half_pair_checks',[])
    signs=[float(x['result']['F']) for x in checks if 'F' in x['result']]
    pair_supported=len(signs)==3 and signs[0]*signs[1]<0 and signs[1]*signs[2]<0
    if pair_supported and len(n.get('roots',[]))<2:missed.append((row['label'],row['kind']))
    rows.append(dict(label=row['label'],kind=row['kind'],c_exact=str(v[11]),m_exact=str(-v[7]),K_exact=str(-v[7]*(11*v[11]-5)/5-42),beta_exact=str(v[8]),coefficient_vector=json.dumps(row['coefficient_vector']),origin_count=len(n.get('roots',[])),origin_stability=''.join(n.get('stab',[])),origin_roots=json.dumps(n.get('roots',[])),section_angle=n.get('theta'),return_edge=n.get('redge'),edge_kind=n.get('edge_kind'),edge_D=edge,outside_outer_D=outside['D'] if outside else None,sign_mismatch=disagree,known_pair_supported=pair_supported,remote_roots=json.dumps([x['roots'] for x in row['fable']['nests'] if abs(x['pt'][0])+abs(x['pt'][1])>1e-7])))
with (HERE/'field_ledger.csv').open('w') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator="\n");w.writeheader();w.writerows(rows)
endpoints={label:next(r['fold'] for r in accepted if r['label']==label) for label in ['center_0.0000000001','positive_extension_1','negative_extension_1']}
for label,a in endpoints.items():
    a['K_exact_from_vector']=str(-Q(a['coefficient_vector'][7])*(11*Q(a['coefficient_vector'][11])-5)/5-42)
summary=dict(accepted_fold_checkpoints=len(accepted),field_records=len(fields),unchanged_fable_evaluations=len(fields),rationally_scaled_fable_evaluations=sum('scaled_fable' in r for r in fields),dense_profiles=len(dense),max_observed_origin_count=max(r['origin_count'] for r in rows),max_dense_sampled_count=max(len(r['brackets']) for r in dense),four_origin_field=None,endpoint_records=endpoints,raw_sign_mismatches=mismatch,known_pair_missed_by_unmodified_counter=missed,binary128_full_checks=len(full),binary128_full_unresolved=sum(r['result']['status']!='NUMERICAL_ONLY' for r in full),status='OPEN',sign_map_verdict='No additional-pair signature detected in resolved sampled profiles; global existence unresolved.',exhaustive=False)
(HERE/'summary.json').write_text(json.dumps(summary,indent=2))

text='''# FASTRA D1 — Joint fold work, 2026-09-05

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
'''
for name,label in [('positive, centerward','center_0.0000000001'),('positive, outward','positive_extension_1'),('negative, outward','negative_extension_1')]:
    a=endpoints[label];text+=f"| {name} | {float(a['c']):.16g} | {float(a['K']):.16g} | {float(a['r']):.13g} |\n"
text+='''
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
'''+json.dumps(example['coefficient_vector'])+'''
```

Equivalently,

\\[
\\dot x=y+x^2+xy,\\qquad
\\dot y=-\\frac{37101199745401}{10^{12}}x-\\frac{1}{25600000}y
-10x^2+\\frac{11}{5}xy+\\frac{242288563571}{250000000000}y^2.
\\]

Its exact parameters are
`c = 242288563571/250000000000`,
`K = 4660642062301237256681/1250000000000000000000000`,
`beta = -1/25600000`.

| Nest | Ray angle | Radius | Stability |
|---|---:|---:|---|
'''
for r in example['replays']:
    nest='origin' if abs(r['focus'][0])<1e-7 else 'remote'
    text+=f"| {nest} | {r['theta']:.13g} | {r['r']:.13g} | {r['stability']} |\n"
text+='''
Both Fable engines return with matching displacement signs on both sides of
each of these four roots. `rational_3_plus_1.json` records the coefficients,
focus coordinates, rays, root estimates, endpoint signs, return statuses, and
Cartesian periods. This is numerical cross-verification, not an interval
certificate. It is a 3+1 field, not the requested 4+1 field.

## Work completed and count limitations

'''
text+=f"- {len(accepted)} accepted fold checkpoints; {len(fields)} exact-rational field records.\n"
text+=f"- At every accepted checkpoint, the unchanged `sweep_log.evaluate(coef)` was called at the fold, the pair-present side, and two beta amplitudes: {len(fields)} evaluations.\n"
text+=f"- {summary['rationally_scaled_fable_evaluations']} additional evaluations used the exact rational conjugacy `x=X, y=sY, tau=s*t` to balance large m. These coefficient vectors are also saved.\n"
text+=f"- {len(dense)} supplemental full origin profiles used the horizontal ray, log-radius spacing 0.125 on [-25,46], two tolerances, and explicit integration-failure records. Maximum sampled count: {summary['max_dense_sampled_count']}.\n"
text+=f"- {len(full)} binary128 full-return endpoint checks; {summary['binary128_full_unresolved']} unresolved. All successes and failures are retained.\n"
text+='''
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

'''
text+=f"There are {len(mismatch)} opposite-sign comparisons between the outside-last-root sample and the recorded edge in the unchanged counter's usable records.\n\n"
text+='''
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
'''
for r in accepted:
    a=r['fold'];text+=f"| {r['label']} | {float(a['c']):.13g} | {float(a['K']):.13g} | {float(a['r']):.13g} |\n"
text+='''
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
'''
(ROOT/'FASTRA_D1_REPORT_2026_09_05.md').write_text(text)
