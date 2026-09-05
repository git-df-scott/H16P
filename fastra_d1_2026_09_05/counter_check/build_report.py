"""Build Fable's same-ray full-return discrepancy handoff."""
import json,csv,hashlib
from pathlib import Path
import mpmath as mp
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent.parent
mp.mp.dps=65
def load(f):
    p=HERE/f;return [json.loads(x) for x in p.read_text().splitlines()] if p.exists() else []
misses=load('missed_roots.jsonl');refs={(r['label'],r['kind'],r['root_index']):r for r in load('refined_coordinates.jsonl')}
ids={('center_0.0000000001','pair'):'C0',('center_0.0000000001','hopf'):'Cbeta',('positive_extension_1','pair'):'P0',('negative_extension_1','pair'):'N0',('negative_extension_1','hopf'):'Nbeta'}
records=[];gridrows=[]
for r in misses:
    fid=ids[(r['label'],r['kind'])]
    rec=dict(id=fid,label=r['label'],kind=r['kind'],coefficient_vector=r['coefficient_vector'],parameters=r['parameters'],original_counter=r['original_fable'],roots=[])
    for x in r['roots']:
        refined=refs.get((r['label'],r['kind'],x['index']))
        root=dict(index=x['index'],stability=x['stability'],section_theta=x['section_theta'],section='(x,y)=rho*(cos(theta),sin(theta)), origin focus',approx_r=refined['r'] if refined else None,root_residual=refined['D'] if refined else None,numerical_sign_bracket=[g['native']['r'] for g in x['grid']],grid=x['grid'])
        root['classification']='finite-radius pair at center parameter limit; NOT small-radius roots' if fid.startswith('C') else 'large-radius root'
        for side,g in enumerate(x['grid']):
            assert g['quad']['status']=='NUMERICAL_ONLY'
            u=mp.mpf(g['u_decimal']);L=mp.mpf(g['quad']['L'])
            gridrows.append(dict(field=fid,root=x['index'],side=side,theta=x['section_theta'],u=g['u_decimal'],rho=mp.nstr(mp.exp(u),36),quad_log_return=mp.nstr(u+L,36),quad_radial_return=g['quad']['return_coordinate'],quad_log_displacement=g['quad']['L'],native_log_return=g['native']['u_return'],native_log_displacement=g['native']['D'],native_status=g['native']['status'],inside_default_grid=g['inside_default_grid'],visited_in_original_run=g['visited_in_original_run']))
        assert mp.mpf(x['grid'][0]['quad']['L'])*mp.mpf(x['grid'][1]['quad']['L'])<0
        rec['roots'].append(root)
    records.append(rec)
(HERE/'discrepancy_ledger.json').write_text(json.dumps(records,indent=2))
with (HERE/'grid_return_values.csv').open('w') as f:
    w=csv.DictWriter(f,fieldnames=list(gridrows[0]),lineterminator='\n');w.writeheader();w.writerows(gridrows)
txt='''# FASTRA D1 — Counter discrepancy audit, 2026-09-05

**D1 stays OPEN. No sheets were extended and no sweeps were repeated.**
The original D1 snapshot was published on `astra/fastra-d1-2026-09-05`
as `6f3ad357930623d664843ff9c254d2865c6a2bb2`, with an identical tree to
local commit `6f744d5997249047e1bfeaa16b2a513100f48686`.

## Scope correction

There are **five fields, ten missed pair roots** supported by successful
archived binary128 **full returns**, and now by same-ray full-return brackets.
The other D1 half-return checks must not be represented as full-return
cross-verification. In particular, `positive_extension_1/hopf` had full-return
timeouts and is **not** in this list.

My phrase "near the center" meant the **center stratum in parameter space**,
not the small-radius end of the return domain. The demonstrated C0/Cbeta misses
are at finite radii about 8 and 18 on Fable's ray. This record does **not**
establish a binary128-verified missing small-radius Hopf cycle.

Sweep negatives remain lower bounds. These examples do not demonstrate that
an overnight sampled field has five cycles, or identify which overnight runs
need repetition. Check each run's parameter range, failure records and precision.

## 1. Every full-return-supported missed field

Section: `(x,y)=rho*(cos(theta),sin(theta))` about the origin, using the ray
selected by `sweep_log.evaluate`. Decimal theta strings preserve its recorded
angle; binary128 uses that decimal angle, Fable its double representation.
Coordinates and sign brackets are numerical, not interval certificates.

| ID | `field_ledger.csv` label / kind | theta | Missed rho and stability | Location |
|---|---|---:|---|---|
'''
for r in records:
    coords=[]
    for x in r['roots']:
        pos=f"{float(x['approx_r']):.11g}" if x['approx_r'] else '['+', '.join(f'{v:.9g}' for v in x['numerical_sign_bracket'])+']'
        coords.append(pos+' '+x['stability'])
    loc='finite radius; center parameter limit' if r['id'].startswith('C') else 'large radius'
    txt+=f"| {r['id']} | `{r['label']}` / `{r['kind']}` | {r['roots'][0]['section_theta']} | {'; '.join(coords)} | {loc} |\n"
txt+='''
P0 is reported by adjacent extended-grid brackets. Its last secant iterate was
not promoted to a precise root estimate. The other coordinates were refined
by safeguarded full-return Newton iteration; residuals remain in the JSON.

Exact rational vectors follow. Order:
`(P_constant,P_x,P_y,P_x2,P_xy,P_y2,Q_constant,Q_x,Q_y,Q_x2,Q_xy,Q_y2)`.

'''
for r in records:txt+=f"**{r['id']}**\n\n```json\n{json.dumps(r['coefficient_vector'])}\n```\n\n"
txt+='''The standalone `counter_check/discrepancy_ledger.json` includes every exact
vector, original counter output, section, stability, coordinate and bracket.

## 2. Counter settings: unchanged defaults

The original 176 `sweep_log.evaluate(coef)` calls used unchanged Fable sources
from `afbcdd419309e30222e494e075c45b3049350020`.

| Setting | Value actually used |
|---|---|
| Initial grid radius | `1e-3*scale`; scale is capped at 1 and equals 1 for these origin nests |
| u-grid | `np.arange(log(0.001),40,0.25)` |
| Last actual u / radius | `39.84224472101786` / approximately `2.010328497117133e17` |
| Grid step | `0.25` |
| rtol / C atol | `1e-12` / `1e-15` |
| Sign threshold | `min(abs(D_left),abs(D_right)) > 1e-10`, D=`u_return-u` |
| Integrator escape cap | `umax=45`, because the counter passes grid cap + 5 |
| Rescaled-time cap | `Smax=2000` |
| Step cap | `300000` |
| Initial integrator step | `h=1e-3` |
| Chunk / failure handling | 8 radii per chunk; stop at first nonzero status |
| Edge refinement | 8 bisections, only when a previous valid point exists |
| Ray | Away from nearest other equilibrium |

None was changed in the original counter evaluations. The older supplemental
D1 profiles were separate: spacing 0.125 on log-radius [-25,46], horizontal
ray, rtol 1e-12 and 1e-13, umax 60, Smax 10000, maxsteps 1000000. Another 32
evaluations used rational coordinate scaling. These do not replace the original
counter outputs.

For this check, direct Fable calls use the exact default settings above.
Unvisited and extended-grid samples are labeled. Binary128 uses independent
modified-midpoint extrapolation, tolerances 2e-25 for C/N and 2e-27 for P,
recorded per call. The separate center tolerance ladder uses rtol 1e-13,
1e-14 and 2e-15; it is not a new counter or sweep result.

## 3. Full returns at adjacent grid points

`counter_check/grid_return_values.csv` supplies all 20 endpoints: radial input,
radial full return, log full return, displacement, status, and whether the
original counter visited the point.

### C0: stable root near rho=7.9251591

Two actual, visited, adjacent default grid points:

| u | rho | Binary128 full log return | Binary128 D | Double D | Status |
|---:|---:|---:|---:|---:|---:|
'''
def add_table(fid):
    return ''.join(f"| {g['u']} | {float(g['rho']):.14g} | {g['quad_log_return']} | {float(g['quad_log_displacement']):.12e} | {g['native_log_displacement']:.12e} | {g['native_status']} |\n" for g in gridrows if g['field']==fid and g['root']==1)
txt+=add_table('C0')
txt+='''
The grid brackets the sign change; edge logic is not involved. True D is much
smaller than 1e-10. Default double integration additionally has error about
2e-11 and reports the wrong sign at the left endpoint. Lowering the noise
filter alone cannot fix this. Even tighter double tolerances do not consistently
resolve every sign in the saved tolerance ladder. Binary128 on the exact
double-rounded coefficient vector retains this sign change, ruling out
coefficient rounding as its cause.

### N0: unstable root near rho=1.920157854e11

Adjacent points inside the default grid, but not visited by `evaluate`, which
already stopped at its first sample. Double columns are direct replay calls:

| u | rho | Binary128 full log return | Binary128 D | Double D | Direct status |
|---:|---:|---:|---:|---:|---:|
'''
txt+=add_table('N0')
txt+='''
The original first chunk has status 2 at all eight radii. Instrumentation of
the unchanged first step at rho=0.001 produces non-finite error, new u and
new theta. Since `NaN > 1` is false, the rejection test accepts the invalid
step; execution eventually reports the time cap. This status 2 is not evidence
of a dynamical return edge. `first_step.c` and `first_step.json` reproduce the
failure without changing production code. Direct large-radius calls also show
inaccurate displacements despite status 0. The root signal exceeds the 1e-10
filter by orders of magnitude: this is not noise-floor rejection. Binary128
full returns on the exact double-rounded coefficients preserve the same
negative/positive sign bracket to about 5e-17 in log displacement; coefficient
rounding is not the explanation for these N0 grid discrepancies.

### P0: scan cap and coefficient representation

Both roots lie outside u<40. There are no actual default grid points bracketing
them. The CSV supplies adjacent points on a hypothetical extension of the
same grid, labeled `inside_default_grid=false` and `visited_in_original_run=false`.

The exact field's transverse pair offset is approximately 6.1e-25 in c, while
conversion to double changes c at about 1e-17 scale and changes m too. Thus
the double vector does not faithfully specify this near-fold rational field.
Direct double discrepancies beyond the cap do not establish missed roots of
the double-rounded field itself. Binary128 full returns on the rounded field
reach the step-resolution guard; these are unresolved, not no-return proofs.
The exact and rounded vectors and results are in `coefficient_rounding.jsonl`.
The rounded-field half-map residual is positive at all three original pair
bracket points (about 8.4400e-18, versus the exact-field splitting near 1e-26).
This destroys those saved sign brackets; it is not a global zero-exclusion
proof for the rounded field. This half-map test is separate from full-return evidence.

## Consequences and reproduction

Treat sweep negatives as lower bounds. Reject or flag non-finite stages;
distinguish early failures from a zero count. Flat profiles need a precision
check, not simply a smaller noise filter. Cap changes must be explicit, and
near-fold input coefficients must preserve the splitting being tested.

No counter patch, overnight rerun, sheet extension, or new D1 closure is
claimed. The centerward positive-sheet sign map remains next after this
handoff and Fable's reproduction.

From the repository root, build the original engines with
`bash fastra_d1_2026_09_05/build_engines.sh`, then:

```bash
g++ -O2 -std=c++17 -fext-numeric-literals fastra_d1_2026_09_05/counter_check/full_ray_quad.cpp -o /tmp/d1_full_ray -lquadmath
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 python fastra_d1_2026_09_05/counter_check/check_misses.py
```

The driver skips completed field records. `calls.jsonl` distinguishes full
returns from partial passages used solely to transfer sections.
`refine_coordinates.py` refines C/N roots inside full-return grid brackets;
P0 is reported by brackets. Every result is numerical, not interval-certified.
'''
(ROOT/'FASTRA_D1_COUNTER_DISCREPANCY_2026_09_05.md').write_text(txt)
manifest={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in HERE.iterdir() if p.suffix in ['.py','.cpp','.c']}
(HERE/'source_manifest.json').write_text(json.dumps(manifest,indent=2))
print(json.dumps(dict(fields=len(records),missed_pair_roots=sum(len(r['roots']) for r in records),grid_endpoints=len(gridrows),status='OPEN')))
