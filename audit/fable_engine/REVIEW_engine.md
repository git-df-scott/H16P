# Hostile correctness review — Fable return-map engine

Files reviewed: `retmap.c`, `retmap.py`, `sweep_shi.py`, `sweep_family.py`, `evolve.py`
(read in full; nothing in them was modified). All experiments were run against a
private copy of `retmap.c` compiled in a scratch directory
(`gcc -O2 -fopenmp -shared -fPIC`), `OMP_NUM_THREADS=1` unless noted, each test < 60 s.

**Bottom line.** The core algorithm (winding-gated ray return, DP45, sign changes of
`D = R - r`) is sound and reproduces KKL exactly (`0.683210, 2.183700, 15.962784`).
But the *production settings* used by the sweeps are far too loose to resolve the
regime the sweeps are searching, and the reported multi-cycle hits are dominated by
integration noise. Two headline results:

* **The Yu-Zeng four-cycle field, run through `sweep_shi.py`'s exact settings
  (`NR=40`, `RMIN=1e-3·scale`, `rtol=1e-9`, `Tmax=2e3`, `maxsteps=3e5`), reports 2
  cycles, not 4.** With `NR=200, rtol=1e-12, Tmax=1e5` the same code reports 4.
* **3 of the 4 `count=2` records I re-checked in `data/F3_lam0_L4.jsonl` are actually
  `count=1`**: the inner sign change is at `|D| ~ 1e-14`, at/below the integrator's
  own `atol`, and disappears when `rtol` is tightened from `1e-9` to `1e-12`.

Severity key: **C** critical (wrong published counts), **H** high, **M** medium, **L** low.

---

## A. False cycles

### A1 [C] `atol = 1e-14*(1+r)` is an absolute floor: every sign change with `|D| < ~1e-13` is noise
`retmap.c:61` sets `atol = 1e-14*(1.0+r)`, and `dp_step` scales by
`sx = atol + rtol*max(|x|,|xn|)`. For `r ≲ atol/rtol` the step controller is
absolute-tolerance-limited, so `D` cannot be resolved below ~`1e-13` no matter how
small `r` is. `sweep_shi.py` uses `rtol=1e-9` and `RMIN=1e-3·scale`, landing squarely
in that regime, and `--lam0` mode forces `tr = 0` exactly (a linear centre), where
`D(r) ~ V₃r³` is *by construction* tiny.

Failing example (the engine's own output, `data/F3_lam0_L4.jsonl`, Shi chart `lam=0`,
`par = (l,m,n,a,b)`; re-run on the record's own stored radius grid):

| record | reported | `min|D|` at inner bracket | `rtol=1e-9` | `rtol=1e-12` |
|---|---|---|---|---|
| 0 (`par=[0.2037,-3.3932,0.9532,2.7724]`, l=-1.9608) | count=2 | 6.4e-15 | 2 sign changes | **1** |
| 5 | count=2 | 6.7e-15 | 2 | **1** |
| 11 | count=2 | 2.1e-14 | 2 | **1** |
| 23 | count=2 | 5.0e-03 | 2 | 2 (genuine) |

In record 11 the innermost `D` is `-2.06e-14` at `rtol=1e-9` and `+1.27e-13` at
`rtol=1e-12` — a pure sign flip of noise. Every one of the 25 `count=2` records in
that file except one has inner-bracket `min|D|` between `6e-15` and `2e-13`.

**Minimal fix:** `atol = 1e-16*r` (scale the absolute tolerance to the orbit, not to 1),
and raise the sweep tolerance to `rtol ≤ 1e-12`. Then apply the noise rule in §F.

### A2 [C] Tolerance is scaled by `|x|`, i.e. by the *focus offset*, not by the orbit radius
`sx = atol + rtol*max(|x|,|xn|)` uses the absolute coordinate. When the antisaddle sits
far from the origin, the per-step tolerance is `rtol·|focus|` while the orbit radius is
`r ≪ |focus|`, so the effective relative accuracy in the local coordinate degrades by
`|focus|/r`. Below `r ~ |focus|·10⁻¹¹` the *roundoff* in `x` itself swamps the orbit.

Failing example (translate KKL so the focus sits at `(X0,0)` — identical dynamics;
`be = 1e-13` makes the origin a numerical centre, so true `D/r` is smooth and monotone;
grid = 20 radii on `[1e-5,1e-2]`, `rtol` as shown):

| `rtol` | `X0=0` | `X0=1e2` | `X0=1e4` | `X0=1e6` |
|---|---|---|---|---|
| 1e-8  | 0 sign changes | 0 | **12** | 3 (+11 failures) |
| 1e-10 | 0 | **3** | **5** | 20 failures |
| 1e-12 | 0 | **5** | **10** | 20 failures |

Tightening `rtol` does **not** help (5 → 10 spurious cycles at `X0=1e4`): the error is
roundoff in the global coordinate. At `X0=1e3` with `be=1e-13`, a 40-point grid on
`[1e-4,1e-1]` produces **10 spurious limit cycles** where the untranslated system
produces 0.

**Minimal fix (verified):** re-expand the quadratic field about the focus once, at the
top of `full_return`, and integrate in local coordinates. Six lines:

```c
double cl[12];
for (int k=0;k<12;k+=6){
  double a0=c[k],a1=c[k+1],a2=c[k+2],a3=c[k+3],a4=c[k+4],a5=c[k+5];
  cl[k+0]=a0+a1*fx0+a2*fy0+a3*fx0*fx0+a4*fx0*fy0+a5*fy0*fy0;
  cl[k+1]=a1+2*a3*fx0+a4*fy0;
  cl[k+2]=a2+a4*fx0+2*a5*fy0;
  cl[k+3]=a3; cl[k+4]=a4; cl[k+5]=a5;
}
c = cl; fx0 = fy0 = 0.0;
double x = r*dx, y = r*dy, t = 0;
```
With this plus `atol = 1e-16*r`, the table above becomes **0 spurious sign changes at
every `X0` up to `1e8`, and `max|D/r|` is identical (3.390e-07) at all offsets**. KKL still
returns `0.683210, 2.183700, 15.962784`. This is the single highest-value fix in the file.

### A3 [M] The ray-crossing direction filter is dead code
`retmap.c:80` — the third disjunct `(steps>0 && gprev*gnew<=0)` logically *subsumes*
the two sense-checked clauses before it, so a crossing of the ray in the **wrong**
transverse direction is accepted. A/B test (512 Sobol Shi sets, 856 nests, 44 240 radii,
production settings): removing the third clause leaves all statuses identical and the
total sign-change count unchanged (64 vs 64), but **22 of the 10 699 successful returns
report a different `R`**. In all 22 the as-is code returns `R < 0` (e.g. `R = -9.75e-11`
at `r = 1.94e-3`), i.e. a crossing on the *opposite* ray — the `dot > 0` test is applied
to the pre-bisection endpoint only and is not re-checked at the bisected root.

**Minimal fix:** delete `|| (steps>0 && gprev*gnew<=0)`, and after the bisection reject
the return with a new status if `*Rout <= 0`.

### A4 [L] Near-focus angle wrapping is possible but was not observed
`dang = atan2(...)` is confined to `(-π,π]`, so a step turning more than `π` about the
focus silently loses `2π`. There is no angular step limiter — `h` is controlled by the
local error in `(x,y)` only. Instrumented runs on KKL (origin and second antisaddle,
`r` up to 1e5) recorded `max|dang| = 0.31`, so this never fired here, but it is one
unlucky field away.
**Minimal fix:** after computing `dang`, if `fabs(dang) > 0.5` reject the step and set
`h *= 0.5/fabs(dang)` (cheap; only bites near the focus). Also guard `relx=rely=0`.

### A5 [L] `if (fabs(sp) < 1e-300) return 5` is an absolute test
The rotation sense is taken from the initial transverse velocity with an absolute
threshold. It should be relative to the speed (`fabs(sp) < 1e-12*speed`), otherwise a
ray that is nearly tangent to the flow at the start point yields a sense determined by
roundoff. Note also that the ray direction is hard-wired to `(1,0)` everywhere in the
drivers and its transversality (`J21(focus) ≠ 0`) is never checked.

---

## B. Missed cycles

### B1 [C] The sweep radius grid is far too coarse — 1.5× between consecutive points
`NR=40` over `[1e-3·scale, 3e3·max(scale,1)]` is 6.5 decades, i.e. a ratio of
`(3e6)^(1/39) ≈ 1.50` per interval. Two limit cycles closer than 1.5× cancel and are
counted as **zero**. Yu-Zeng's three inner cycles are at `r = 0.0228, 0.0589, 0.1013`
— ratios 2.59 and 1.72, i.e. 1.3–2.3 grid intervals apart. Combined with A1 this loses
two of them outright:

```
Yu-Zeng, sweep_shi.py production settings   -> total 2  (nest at origin: 1 of 3)
Yu-Zeng, NR=200, RMIN=1e-5, rtol=1e-12      -> total 4  (0.0232, 0.0561, 0.1010 + 234.2)
```
An even number of cycles inside one interval is invisible, and the `near`/`gaps`
"near-miss" heuristic cannot see them either — it only inspects `D` at grid points.

**Minimal fix:** `NR ≥ 160` for the same span (ratio 1.09), or two-pass: coarse grid,
then locally refine every interval where `D/r` has a local extremum or where
`|D/r|` drops below 10× its neighbours. Also add a semi-stable detector: a local
extremum of `D` with `min|D|` below the noise floor of §F is an unresolved
double root, not a "near miss".

### B2 [H] Any single failed radius truncates the whole nest
`count_nest` (`retmap.py:110`) and both sweeps do `k = argmin(ok)` and discard
*everything beyond the first non-zero status*. Statuses 2 (`Tmax`) and 3 (`maxsteps`)
mean "I gave up", not "nest boundary". In the 512-set batch only **10 699 of 44 240
radii (24 %) returned successfully**: 17 508 escaped (st 1), 4 657 hit `Tmax` (st 2),
1 356 stalled (st 4). Yu-Zeng's origin nest truncates at grid index 15 of 40 under
production settings (`Tmax=2e3`); KKL's second antisaddle truncates at 17 of 40.

**Minimal fix:** keep the full grid, mark failed radii as `NaN`, and only bracket sign
changes between *adjacent successful* radii; treat status 1 and 4 as a genuine boundary
(escape / graphic) and statuses 2 and 3 as "unknown" — retry those with a larger budget
before declaring a boundary.

### B3 [H] `Tmax = 2e3` and `maxsteps = 3e5` in the sweeps are below what real nests need
Return times diverge logarithmically as the outer boundary of a nest approaches a
graphic, so the outermost (and often the decisive) cycle is exactly where `Tmax` bites.
`count_nest`'s own default is `Tmax=1e5`; the sweeps use `2e3`, 50× smaller.
Also `RMAX = 3e3` in both sweeps is **smaller than the KKL remote cycle at r ≈ 3711**
that the engine is claimed to validate against — anything at that radius is off the grid
unless `scale > 1.24`.

### B4 [M] I could not reproduce the "remote cycle near 3711"
With `dir=(1,0)` from the origin, `D(r)` is smooth and strictly negative for all
`r ∈ [16, 1e5]` (60 points, `rtol=1e-10`, `Rmax=1e8`, `Tmax=1e6`, all status 0). I
verified the engine's returns independently with `scipy.solve_ivp` (`rtol=1e-12`,
`max_step=1e-4`): at `r=3711` scipy gives `T=0.52175, R=277.65`, matching the engine to
5 digits — so the engine is right and there is no root there *on this ray*. With
`dir=(-1,0)` every radius from 100 to 1e5 returns status 3. Around the second antisaddle
every `r ≥ 5.96` returns status 3 (orbit leaves that antisaddle's winding region after
`t ≈ 370` with `cum ≈ 0.6`, having burnt the whole step budget).
**This validation claim should be re-derived or retired**; as stated it does not hold
for the default ray.

---

## C. Stability classification

### C1 [M] `stab = 'S' if Dlo > 0 else 'U'` is correct only for a transversal simple root
`retmap.py:127`. The convention is right (`D>0` inside, `D<0` outside ⇒ stable). Two
failure modes:
* If the refine loop aborts (`if sm[0,0] != 0: break`), `Dlo` is whatever survived from
  a partial bisection and `root = 0.5*(lo+hi)` is an unrefined bracket midpoint — the
  label is kept with no flag that refinement failed.
* Semi-stable cycles and even-multiplicity roots are never labelled at all (B1).
Also `root = 0.5*(lo+hi)` is an **arithmetic** midpoint of a *geometric* bracket; with
`refine=False` (and in `sweep_shi.py`, which reports `rad[i,j]`, the left grid point)
the reported radius can be 13–20 % low — KKL under sweep settings reports
`0.65315, 2.0571, 13.9205` for true `0.68321, 2.18370, 15.96278`.
**Minimal fix:** return a `refined: bool` flag; use `sqrt(lo*hi)`; require the
multiplier `dD/dr` at the root to exceed the noise floor before assigning S/U.

### C2 [L] `antisaddles` counts nodes as nests
`det J > 0` admits nodes; the `focus` flag is computed but never used by any driver.
Node nests can only burn budget (they never wind 2π and always end in status 2/3).
Skipping `focus == False` would cut sweep cost measurably.

---

## D. Equilibria

### D1 [H] `c5 == 0 and c11 == 0` ⇒ *all* equilibria silently lost
`equilibria` builds the resultant of two quadratics in `y`. With `A1 = A2 = 0` the
expression `(A1C2-A2C1)² - (A1B2-A2B1)(B1C2-B2C1)` is identically zero, the
`np.all(|co| < 1e-14)` guard fires, and the function returns an empty array.

```
c = [0,0,1, 1,1,0,  0,-1,0, 1,1,0]      # P = y+x²+xy, Q = -x+x²+xy
rm.equilibria(c) -> array([], shape=(0,2))     # the origin is an equilibrium
rm.shi_coef(0.1,1.0,2.0,0.5,0.3, n=0.0) -> []  # Shi chart with n = 0
```
The Shi chart has `c11 ≡ 0` *always*, so this is a codimension-1 wall at `n = 0` inside
the very family `sweep_shi.py` samples. It also creates a blind spot for `evolve.py`
`--dims=...,11,...` (KKL seed has `c5 = 0`, and `c11` starts at 0.7 with step
`0.05·0.7`): as `c11 → 0` a spurious-looking equilibrium runs off to `y ≈ 2.2/c11`
(`(-1, 2.198e13)` at `c11 = 1e-13`) and at `c11 = 0` the count drops to **zero**, so any
candidate that lands on that wall scores `total = 0` and is discarded.

**Minimal fix:** if `|c5| + |c11| < tol·scale`, both `P` and `Q` are linear in `y`;
solve `B1 y + C1 = 0`, `B2 y + C2 = 0` by eliminating `y` (a cubic in `x`:
`B2 C1 - B1 C2 = 0`) and handle `B1 ≡ B2 ≡ 0` separately. Similarly, if only one of
`c5, c11` vanishes, drop the degree explicitly rather than relying on `np.poly1d`
stripping near-zero leading coefficients.

### D2 [M] `scale = min(others)` can exceed `rmax`, silently producing a descending grid
`cycles_all_nests` passes `rmin = rmin_frac*scale` with a fixed `rmax` (default 1e6).
When a runaway equilibrium sets `scale > rmax/rmin_frac`, `np.geomspace(rmin, rmax, nr)`
returns a **decreasing** grid; the sign-change scan then walks inward, inverting every
S/U label and reporting roots in reverse. Nothing guards `rmin < rmax`. `count_nest`
can also return a zero-length `rad`, which `cycles_all_nests` callers index without
checking.
**Minimal fix:** `rmin = min(rmin, rmax/1e3)`; assert `nr >= 2` after truncation.

### D3 [M] The equilibrium acceptance gate is quadratically loose and Newton is not re-validated
`abs(P)+abs(Q) < 1e-6*(1+|x|+|y|)**2` is `1e-2` at `|x| ~ 100` and `1e12` at
`|x| ~ 1e9`. The subsequent 5 Newton steps have no convergence test and the point is
appended **unconditionally**, `linalg.solve` failures included (`break`, then append).
I did not manage to make it emit a false equilibrium — a 4096-set Shi audit gave
13 648 equilibria, all with residual `< 1e-10` relative, and a 256-set comparison
against `fsolve` from a 17×17 grid found **zero** missed equilibria — so this is a
latent risk rather than an observed one. Still: re-test the residual *after* Newton
and drop the point if `|P|+|Q| > 1e-12*(1+|x|+|y|)`.

### D4 [L] Double roots survive, triple roots would not
At a saddle-node the quartic's double root comes out with `Im = 0` exactly (verified by
bisecting to `m* = 1.7526110771` in the Shi chart: at `d = 0` the roots are
`-4.10602, -4.10602, 0, 0`, all real; at `d = -1e-8` they become a conjugate pair with
`|Im| = 4.1e-4`, correctly rejected). The `|Im| > 1e-7*(1+|Re|)` filter is therefore
adequate for double roots but would reject a triple root, whose numerical imaginary
part scales as `ε^{1/3} ≈ 6e-6`. Low priority; widen to `1e-4·scale` and let the
post-Newton residual test do the filtering.

---

## E. Bisection and step-size controller

### E1 [M] The `g = 0` bisection has no bracket validity check and mishandles `gnew == 0`
`retmap.c:82-91`. The invariant assumed is `sign(g(0)) = sign(gprev) ≠ sign(gnew)`. Two
holes: (i) `if ((gm>0) == (gnew>0)) hi = mid;` — when `gnew == 0` exactly, the test is
`(gm>0) == false`, so `hi` collapses toward 0 and the routine returns essentially the
*starting* point of the step as the crossing; (ii) when the crossing was accepted via
the wrong-direction clause (A3), `gprev` and `gnew` may not bracket the root that the
winding gate intended. Also `*Tout = t + 0.5*(lo+hi)` uses the final bracket midpoint
while `(xm,ym)` was evaluated at the *previous* midpoint (`O(1e-15·h)`, harmless).
**Minimal fix:** bail with a failure status if `gprev*gnew > 0`; use
`if (gm == 0) { lo = hi = mid; break; }` and compare `gm` against `gprev`'s sign.

### E2 [M] No lower bound on `h`, no cap on the step relative to the crossing
On a rejected step `h *= fmax(0.2, 0.9*pow(err,-0.2))` with no floor: a stiff patch
drives `h → 0` and the run silently burns the whole `maxsteps` budget and returns
status 3 rather than status 4. This is the dominant cost of B2. The accept branch
`h *= fmin(5.0, 0.9*pow(fmax(err,1e-10),-0.2))` is fine (the `1e-10` clamp only
matters when `err == 0`).
**Minimal fix:** `if (h < 1e-13*(1.0+t)) return 4;` and cap growth once
`sense*cum > 2π - 0.2` so the crossing step is not overshot.

### E3 [L] `dir` is never normalised
`R = (xm-fx0)*dx + (ym-fy0)*dy` is a projection, so a non-unit `dir` silently rescales
`R` while `r` is not rescaled. All current callers pass `(1,0)`; add
`double nn = hypot(dx,dy); dx/=nn; dy/=nn;` to make the API safe.

### E4 [L] `evolve.py` near-miss gate is a broken expression
`evolve.py:72`: `rad[i,j] > 3e-4*rad[i,-1]/1e8*1e8` reduces to `rad[i,j] > 3e-4*rad[i,-1]`
and, since `rad = geomspace(1e-4·scale, 1e4·scale)`, to `rad > 3·scale`. The near-miss
term of the fitness score therefore **ignores every radius inside 3× the distance to
the nearest equilibrium** — exactly where small-amplitude nested cycles live. The
`/1e8*1e8` looks like leftover algebra. This mis-steers the whole evolutionary search.

---

## F. Proposed noise-floor rule for rejecting sign changes

Two error sources, with different scalings — a sign change must clear both.

**1. Roundoff floor (irreducible, does not shrink with `rtol`).**

```
sigma_rd(r) = 10 * u * (|focus| + r) * sqrt(N)        u = 2.22e-16, N = accepted steps
```
Calibration: at `|focus| = 1e4`, `r = 1e-5`, `N ≈ 1e3` this predicts `7e-10`; the
measured `D` noise is `5e-10` (§A2 table), and at `|focus| = 0` it predicts `7e-19`,
consistent with the observed absence of any spurious crossing. `N` is not currently
returned — expose the loop counter alongside `T` (one extra `long*` argument), or use
`N ≈ 1e3` as a stand-in.

**2. Truncation floor (differential, the reliable test).** Recompute `D` at the two
bracket endpoints with `rtol' = rtol/100` and `atol' = atol/100`:

```
sigma_tol(r) = |D_rtol(r) - D_rtol'(r)|
noise(r)     = max(sigma_rd(r), sigma_tol(r))
```

**Acceptance rule.** For a bracket `[r_i, r_{i+1}]` with `D_i·D_{i+1} < 0`, classify:

```
m = min(|D_i|, |D_{i+1}|) ;  z = max(noise(r_i), noise(r_{i+1}))
m > 8 z   -> ACCEPT   (real cycle)
m < 2 z   -> REJECT   (noise)
otherwise -> UNRESOLVED: tighten rtol by 1e-2 and re-test; give up after 3 rounds
             and report as unresolved rather than as a cycle.
```
Additionally require the sign change to **survive** the refinement (not merely to stay
above threshold): if the sign change vanishes at the tighter tolerance, reject
unconditionally — this is what kills the F3 records in §A1.

**Validation of the rule on six known points** (all with `focus` at the origin unless
noted, so `sigma_rd` is negligible and `sigma_tol` binds):

| case | `min|D|` | `noise` | rule | truth |
|---|---|---|---|---|
| KKL r≈0.575 | 5.3e-05 | 3.0e-10 | ACCEPT | real |
| KKL r≈2.13  | 1.2e-04 | 1.1e-09 | ACCEPT | real |
| KKL r≈13.4  | 6.3e-02 | 7.0e-09 | ACCEPT | real |
| Yu-Zeng r≈0.0228 @rtol 1e-12 | 7.9e-13 | 1.5e-14 | ACCEPT | real |
| Yu-Zeng r≈0.0228 @rtol 1e-10 | 6.1e-13 | 1.9e-13 | UNRESOLVED → refine → ACCEPT | real |
| F3 record 0 inner root | 6.4e-15 | 2.2e-13 | REJECT | spurious |
| shifted KKL, `X0=1e4`, r≈1e-5 | 5e-10 | 7e-10 | REJECT | spurious |

**Cheap static prefilter** for the small-`r` end, applied before any integration: for a
focus with trace `τ` and imaginary part `ω`, `D(r) ≈ (e^{πτ/ω} - 1) r + O(r²)`. Any sign
change at `r` below

```
r_lin = 100 * noise(r) / |e^{π*tr/ω} - 1|
```
is spurious, and if `|e^{πτ/ω} - 1| · r` is itself below `noise(r)` for the whole inner
part of the grid (the `--lam0` case, where `τ = 0` exactly), then **no** sign change in
that range can be trusted at the current tolerance — refine or drop the range.

**Orthogonal cross-check (cheap, catches section artefacts rather than noise):** recompute
the candidate root with the ray rotated by 90° (`dir = (0,1)`). A genuine limit cycle is a
property of the orbit and must appear on both sections at the corresponding radii; a
crossing artefact (§A3) will not.

---

## G. Ranked fix list

1. **A2 + A1** — integrate in focus-centred coordinates, `atol = 1e-16*r`. Verified: removes
   all spurious cycles up to `|focus| = 1e8`, preserves KKL and Yu-Zeng exactly. ~10 lines.
2. **§F noise rule** in `count_nest` / both sweeps; re-audit every stored `count ≥ 2` record.
3. **B1** — `NR ≥ 160`, or coarse+refine. **B3** — `rtol ≤ 1e-12`, `Tmax ≥ 1e5`,
   `maxsteps ≥ 3e6`, `RMAX ≥ 1e5` in the sweeps.
4. **B2** — stop truncating a nest on statuses 2/3; retry with a larger budget.
5. **D1** — handle `c5 = c11 = 0`; **D2** — guard `rmin < rmax`.
6. **A3/E1** — delete the wrong-direction crossing clause, reject `R ≤ 0`, harden the bisection.
7. **E4** — fix `evolve.py`'s near-miss radius gate; **C1** — geometric root midpoint and a
   `refined` flag; **A4/A5/E2/E3** — angle limiter, relative sense test, `h` floor, unit `dir`.
8. **B4** — re-derive or retire the "remote cycle near 3711" validation claim.
