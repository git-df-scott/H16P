# KKL origin-cycle exclusion: the weighted-orbit comparison, carried out

Groundwork continuation of `45f4ea9`. Zero ODE calls; the shared 4096-call
ledger is untouched. Replay: [`uniqueness_2026_09_07/`](uniqueness_2026_09_07/).

**Outcome.** The comparison that
[`theory_obstructions.md`](fold_closure_2026_09_05/theory_obstructions.md) §3
listed as the remaining proof obligation is now carried out, and it yields the
first cycle-count theorem for this family: for every `1<=c<8/5` with
`Theta>0` on `x>0` — an explicitly decidable condition that holds exactly on
`K>=K*(c)` — the beta-zero KKL field has **no periodic orbit at all in the
half-plane `x>-1`**. The residual search window for the K1 target
(three origin cycles) is therefore the bounded set `0<K<K*(c)`. In the
residual window the same certificate gives an exact lower bound on the
amplitude of any stable, semistable or multiplier-one cycle, quantifying the
unquantified `delta` of [`kkl/notes_lienard.md`](kkl/notes_lienard.md) §8.

**This is not a result about H(2).** It concerns one two-parameter slice of
one family. It supplies no five-cycle candidate, no exclusion of any other
family, and it does not contradict or exclude any archived candidate: every
archived point in the strip has `K<=7.0643`, and all 124 distinct archived
points lie strictly inside the residual window.

## 1. Setting

The field, its Liénard reduction and the polynomials `W`, `N` are inherited
verbatim from [`kkl/notes_lienard.md`](kkl/notes_lienard.md):

    xdot = (1+x)y+x^2,  ydot = -m x-10x^2+(11/5)xy+c y^2,
    u = 1+x,  e = 11c-5,  d = 16-10c,  m = 5(K+42)/e,
    W = m+(2m+10)x+(m+111/5)x^2+(61/5-c)x^3,
    N = {d u+(c+1)(21+d x)}W - u(21+d x)W'.

On `x>-1` the field is `zdot=w-F(z)`, `wdot=-g(z)` with `f=F'`, and

    f(phi(x)) = -x(21+dx)/(5u),   g(phi(x)) = x W(x) u^{-c-2}.

Throughout `1<=c<=8/5, K>0` the inherited facts are: `m>50/3`, `W>0` on
`u>0` (so `g(z)z>0` for `z!=0` and the origin is the only equilibrium in
`x>-1`), `21+dx>=5+10c>0` on `x>-1`, and — the fold-closure theorem — `N>0`
on `-1<=x<=0`.

Define

    psi(x) = -f/g = (21+d x) u^{c+1} / (5 W),   lambda0 = psi(0) = 21/(5m),
    Theta(x) = 5 m W (psi(x) - lambda0) = m(21+d x) u^{c+1} - 21 W(x).

`Theta(0)=0`, and `Theta` has the sign of `psi-lambda0` because `m,W>0`.

## 2. The comparison

Two exact moments are available on any closed orbit `Gamma` in `x>-1`:

    (M0)  oint g dt = 0            (g = -wdot, w periodic),
    (M1)  log M = -oint f dt       (divergence of the Liénard field is -f).

Hence for **every** constant `lambda`,

    log M = oint [-f - lambda g] dt = oint g (psi - lambda) dt.

`g` has the sign of `x`. So the integrand is positive off `x=0` exactly when
`psi-lambda` has the sign of `x` there. Because `psi` is continuous at `0`
with value `lambda0`, and a closed orbit's `x`-range is an interval
containing `0` in its interior, **`lambda` is pinned to `lambda0`**: no other
constant can work. This is why the earlier constant-combination search in
`notes_lienard.md` §5 found no global certificate — it was searching a
one-parameter family whose only admissible member is `lambda0`.

**Lemma 1 (left side is free).** `d psi/dx = u^c N/(5W^2)`, so on
`1<=c<=8/5, K>0` the inherited `N>0` on `[-1,0]` makes `psi` strictly
increasing there; hence `psi<lambda0`, i.e. `Theta<0`, on all of `(-1,0)`.
The left half of the certificate therefore needs no new hypothesis.

**Lemma 2 (sharpness).** For any `C^1` function `P` and any orbit,
`oint P'(z) zdot^2 dt = oint P g dt` (integrate `d/dt[P(z)zdot]` and use
`oint P f dz = 0`). So

    log M = oint [-f - P g] dt + oint P' zdot^2 dt

for every `P`. If a certificate is sought by making both integrands
pointwise nonnegative — the Cherkas/Zhang form — then continuity at `0`
forces `P(0)=lambda0`, and `P'>=0` forces `psi>=P>=lambda0` on the right and
`psi<=P<=lambda0` on the left. **The whole nondecreasing-multiplier family
reduces to the sign condition on `Theta`.** Nothing is lost by using the
constant `lambda0`.

**Theorem 1 (nonexistence).** Fix `1<=c<8/5` and `K>0` and suppose

    Theta(x) > 0   for every x > 0.                                    (T)

Then the field has **no periodic orbit whatsoever in the half-plane `x>-1`**.

*Proof.* Let `Gamma` be a periodic orbit in `x>-1`. It cannot cross `x=-1`,
where `xdot=1>0` forces one-way crossing, and it must enclose an
equilibrium; `W>0` leaves only the origin. Its `x`-range is
`[p,q]` with `-1<p<0<q`. By Lemma 1 and (T), `g(psi-lambda0)>0` off `x=0`
on that range. The orbit meets `z=0` transversally (there `zdot=w`, and
`w=0` only at the equilibrium), so at finitely many instants; the integrand
is therefore positive off a finite set and
`log M = oint g(psi-lambda0) dt > 0`. Every periodic orbit in `x>-1` is thus
hyperbolic and repelling.

The origin is a repelling weak focus: in the harmonic coordinate
`s=sgn(z)sqrt(2V)` of `notes_lienard.md` §7 the system is
`s'=w-F(z(s)), w'=-s`, whose damping `h(s)=dF/ds=s R(z(s))` with `R=f/g`
has `s^2` coefficient `h_2 = R_x(0)/sqrt(m) = -K/m^{5/2} < 0`; averaging
gives `drho/dtau = -h_2 rho^3/8 > 0`.

Now let `D` be the interior of `Gamma` minus a small disc around the origin
crossed strictly outward by the flow. `D` is compact, forward invariant
(orbits leave the disc into `D`, and leave the repelling `Gamma` into `D`),
and contains no equilibrium. Every periodic orbit in `D` is hyperbolic hence
isolated, and an accumulation of them would be a non-isolated periodic orbit
or would contain an equilibrium; so `D` carries finitely many. Pick `p` in
`D` on none of them. By Poincaré–Bendixson `omega(p)` is a periodic orbit
`Gamma'` in `D` with `p` not on `Gamma'` — impossible, since `Gamma'`
repels. []

**Theorem 2 (amplitude gate).** Fix `1<=c<=8/5`, `K>0`, `X>0`, and suppose
`Theta>0` on `(0,X]`. Then every periodic orbit with `x`-range inside
`(-1,X]` is hyperbolic repelling. Equivalently, **every stable, semistable
or multiplier-one origin cycle has `x_max > X`.**

Theorem 2 is the quantified form of the local statement of
`notes_lienard.md` §8, with the unspecified `delta` replaced by an
explicit `X` and the neighbourhood of the origin replaced by the whole
range `(-1,X]`. It also strictly contains the inherited "must reach the
negative band" restriction: `X` always exceeds the first positive root
of `N`.

## 3. Where (T) holds

`Theta` is nondecreasing in `m` at fixed `c` on `x>=0`, since
`dTheta/dm = u^2[(21+dx)u^{c-1}-21] >= 0` for `c>=1, d>=0`, strictly for
`x>0`; and `m` is increasing in `K`. So (T) holds on a half-line
`K>=K*(c)`, and `K*(c)` is well defined.

**Exact endpoint.** At `c=1`, `u^{c+1}=u^2` makes `Theta` a cubic:

    Theta = (6m-210)x + (12m-2331/5)x^2 + (6m-1176/5)x^3.

All three coefficients are nonnegative exactly when `m>=196/5`, i.e.
`K>=126/25`, and then `Theta>0` for `x>0`. Below that the cubic coefficient
is negative and `Theta<0` for large `x`. Hence

    K*(1) = 126/25 = 5.04   exactly.

**Exact decision at rational `c`.** With `u=t^q`, `q` the denominator of
`c+1`, `Theta` becomes a polynomial in `t`, and `x>0` is `t>1`; a Sturm
count over `Q` decides (T) exactly. Bisection on rational `K` gives

| `c` | `K*(c)` is in |
|---|---|
| `1` | `= 126/25` (exact) |
| `11/10` | `(26855/4096, 53715/8192]` = `(6.55640, 6.55701]` |
| `6/5` | `(137505/16384, 8595/1024]` = `(8.39264, 8.39355]` |
| `13/10` | `(43355/4096, 1355/128]` = `(10.58472, 10.58594]` |
| `7/5` | `(57075/4096, 7135/512]` = `(13.93433, 13.93555]` |
| `3/2` | `(88525/4096, 177075/8192]` = `(21.61255, 21.61560]` |

At `c=8/5` the coefficient `d` vanishes and `Theta ~ m x^{c+1} - 21(61/5-c)x^3`
is negative for large `x` at every `K`: `K*(8/5)=infinity`. `K*` increases
with `c` on the sampled grid and diverges as `c -> 8/5`; monotonicity in `c`
is observed, not proved.

**A closed-form sufficient condition, all real `c` in the strip.** For
`1<=c<=2` and `x>=0`, Bernoulli with exponent `2-c` gives
`(1+x)^{c-1} >= (1+x)/(1+(2-c)x)`, so `Theta >= Theta_3/(1+(2-c)x)` with
`Theta_3 = m(21+dx)(1+x)^3 - 21W(1+(2-c)x)`. Its coefficients are

    x:    m e - 210 = 5K,
    x^2:  3(20cm+10m+350c-1477)/5,
    x^3:  3(45m-15cm+812c-1981)/5,
    x^4:  (80m-50cm+1491c-105c^2-2562)/5,

so **`K>0` together with**

    m >= max{ 7(211-50c)/(10(2c+1)),  7(116c-283)/(15(c-3)),
              21(-5c^2+71c-122)/(10(5c-8)) }

**implies (T)**, hence Theorem 1. This bound is exactly sharp at `c=1`
(it returns `m>=196/5`, `K>=126/25`) and degrades away from it: it gives
`K>=17.58` at `c=11/10` against the true `6.5564`, `K>=216.4` at `c=3/2`
against the true `21.613`, and nothing usable as `c -> 8/5`. Use the Sturm
decision when a sharp threshold matters.

## 4. Residual window and exact amplitude gates

The residual window is `{1<=c<8/5, 0<K<K*(c)}` — bounded in `K` for every
`c` bounded away from `8/5`. Exact gates there (any stable, semistable or
multiplier-one origin cycle must have `x_max` strictly greater than the
listed dyadic-derived value; all decided by exact Sturm counts):

| `c` | `K` | `x_max >` |
|---|---|---|
| `1` | `1` | `0.12866` |
| `1` | `3` | `0.65527` |
| `1` | `5` | `24.171` |
| `11/10` | `1` | `0.094093` |
| `11/10` | `6` | `1.7864` |
| `6/5` | `4` | `0.44674` |
| `6/5` | `8` | `2.2850` |
| `13/10` | `10` | `2.5481` |
| `7/5` | `13` | `3.2825` |
| `3/2` | `7` | `0.64661` |
| `3/2` | `21` | `7.8512` |
| `8/5` | `7` | `0.57754` |

The gate diverges as `K` approaches `K*(c)` from below, joining Theorem 1
continuously. It is **weak deep inside the window**: at the largest saved
pair field, `c=1.59340580527813710990835865677884849`,
`K=7.06390700436779910773804298664181037`, the gate is only
`x_max > 0.5889` (NUM-class evaluation, since that `c` has denominator
`10^35`), whereas the archived returns for that field are at section radii
near `3e17`. The gate therefore does **not** explain the observed amplitudes
and must not be presented as doing so.

## 5. Falsification check against the archive

[`archived_point_check.py`](uniqueness_2026_09_07/archived_point_check.py)
re-reads every saved `(c,K)` record in `fold_surface_2026_09_05/`,
`fold_closure_2026_09_05/`, `staged_2026_09_05/` and `reversible_reseed/`.
52 files are scanned. The largest archived `K` in the strip is `7.0643`;
of the 124 distinct archived points with `1<=c<8/5`, **none** lies in the
certified region.
No archived cycle, fold or two-cycle pair is contradicted, and none is newly
excluded: the entire 3297-call fold continuation ran inside the residual
window. This is a consistency test that the theorem passes, not evidence
for it.

## 6. What this closes and what it does not

Closed, within `1<=c<8/5`:

- The `K>=K*(c)` part of the strip contains no origin-surrounding limit
  cycle, hence no `3+1` configuration with three origin cycles, hence no K1
  candidate. This is the family's first cycle-count theorem; the previous
  status line was "No at-most-two origin-cycle theorem has been obtained".
- The nondecreasing-multiplier (Cherkas/Zhang) certificate family is
  exhausted by the `Theta` sign condition (Lemma 2): improving on this
  requires a genuinely different comparison, not a better multiplier.

Not closed:

- The residual window `0<K<K*(c)`, which is where every archived candidate
  and the whole numerical fold component lives. Theorem 1 does not touch it.
- Any cycle in `x<-1`, including the remote focus and its cycles.
- Monotonicity of `K*` in `c`, and `K*` for irrational or high-denominator
  `c` (only a `c`-uniform sufficient condition is proved there).
- Everything outside this family: no statement about `H(2)`, no five-cycle
  candidate, no other stratum.

**Next task, well posed.** Theorem 2 gives a lower bound `X(c,K)` on the
amplitude of any stable or fold cycle. An *upper* amplitude bound
`X_max(c,K)` for origin cycles in the residual window — from the energy
function `E(r)` of `notes_lienard.md` §3 and the behaviour at the two
boundaries — would close every parameter point with `X(c,K) >= X_max(c,K)`.
That is a bounded, zero-ODE analytic task, and it is the natural way to
attack the residual window without spending returns.

## 7. Replay

    python3 uniqueness_2026_09_07/uniqueness_exact.py      # identities, thresholds, Sturm
    python3 uniqueness_2026_09_07/amplitude_gate.py        # exact amplitude gates
    python3 uniqueness_2026_09_07/archived_point_check.py  # archive consistency

Requires SymPy (exact rational arithmetic throughout) and mpmath for the
NUM-class screening in the archive check only. No ODE integration, no
parameter sweep, no charged returns: the shared ledger stays at 4096/4096.
