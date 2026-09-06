# `Q_3^R` rank repair: the claimed eight-dimensional blind spot does not exist

**Date 2026-09-06. Branch `opus/rank-repair-cusp-compatibility-2026-09-06`.
Source commits: `main` 45f4ea9, `opus/degeneracy-collision` ca4fe90,
`fable/lane2-cusp` e3a6174, `astra/fastra-afternoon` 55537ae.**

This corrects [Q3R_FIRST_ORDER.md](Q3R_FIRST_ORDER.md), which is **wrong** in
its central structural claim. The error was found by Astra.

## 0. The correction, and the attribution

`Q3R_FIRST_ORDER.md` asserted that `T_{a-2}, T_{a-1}, T_a, U` span a
*four*-dimensional first-order space, that the first-order kernel is therefore
eight-dimensional, and that the tiny 4x4 determinants it measured at 70 digits
were an ill-conditioning problem. All three statements are false. The four
functions satisfy the exact relation

\[
 (a+2)\,U \;+\; 3\Bigl[\tfrac{b-2}{4}T_{a-2} + (1-b)\,T_{a-1} + b\,T_a\Bigr] \;=\; 0 .
\]

**This is not new.** It is equation (11) of
[REVERSIBLE_RESEED_2026_09_05.md](REVERSIBLE_RESEED_2026_09_05.md), written
there in area moments,

\[
 (a+2)J + bI_a + (1-b)I_{a-1} + \tfrac{b-2}{4}I_{a-2}=0,
\]

with `I_j = \iint_{D_h} y^j dx dy` and `J = \iint_{D_h} x^2y^{a-2}dxdy`.
Carrying out the inner integrals gives `I_j = 2T_j` and `J = (2/3)U`, and
substituting turns (11) into the boxed identity exactly. That report also
states in terms the earlier document's own words that the repaired direction
does not supply "a fourth independent first-order function". The new work
regressed behind a result the repository already held, and this note is a
reconciliation, not a discovery.

### Proof of the identity

The constants are exactly what makes the following work:

\[
 yR'(y) + aR(y) = -\Bigl[\tfrac{b-2}{4} + (1-b)y + by^2\Bigr]
\]

(verified symbolically: the difference is identically `0`). Hence

\[
 \frac{d}{dy}\bigl(y^{a-1}R^{3/2}\bigr)
 = -\tfrac12\,y^{a-2}\sqrt R\;\Bigl[(a+2)R + 3\bigl(\tfrac{b-2}{4}+(1-b)y+by^2\bigr)\Bigr],
\]

using `2(a-1)R + 3yR' = -(a+2)R - 3\Phi`. Integrating between the turning
points, where `R` vanishes, kills the primitive and gives the identity.

### Verified against the numbers this repository actually produced

Evaluating the relation on the generators computed by `melnikov.basis` — the
same code whose determinants were misread — gives residuals at relative
`8e-18` to `7e-16` across ten `(a,b,h)` points. **The 4x4 determinants are
identically zero.** No amount of precision could have rescued that test; the
70-digit run in `Q3R_FIRST_ORDER.md` was measuring nothing.

## 1. Corrected rank

With the relation eliminating `U`, the first-order image is spanned by
`T_{a-2}, T_{a-1}, T_a`, so its dimension is **at most three**.

Numerically it is exactly three, and better than that: over 12 samples across
the annulus at 60 digits, all 220 of the `3x3` determinants have the **same
sign** at `(a,b) = (-1/2,1), (-1/4,9/10), (-5/2,6/5)`, with magnitudes up to
`1e-4` — an extended Chebyshev system on the sampled range, hence **at most two
zeros** of `M_1` there.

This *reproduces* Fable's F6 result ("span three, at most two zeros") by an
independent route. It does not prove a two-zero bound off the sample: three
generators alone would not imply it, and the Chebyshev property is numerical
here.

### Exceptional loci

`A, B, C` require `a \notin \{0,-1,-2\}`; the elimination of `U` requires
`a \ne -2`. At `a = 1` the coefficients of both `T_{a-2}` and `U` carry the
factor `(a-1)` and vanish, so the image drops to `span\{T_{a-1}, T_a\}` and the
rank is 2. The two-centre region of interest has `a<0`, so `a=1` is outside it,
but the drop is recorded rather than extrapolated through.

## 2. Corrected kernel

Eliminating `U`, `M_1 \equiv 0` iff

\[
 q_{00}=\frac{(b-2)\,q_{20}}{4(a+2)},\qquad
 p_{10}=-a\,q_{01}-\frac{(a-1)(b-1)}{a+2}q_{20},\qquad
 p_{11}=-(a+1)q_{02}+\frac{(a-1)b}{a+2}q_{20}.
\]

Three conditions on the six visible coefficients, plus the six that never enter
`M_1`: the kernel has dimension **nine**, not eight. The old conditions
`q_{00}=q_{20}=0` were an artefact of treating `U` as independent.

## 3. The kernel is exactly the geometrically trivial span

Nine directions are trivial for geometric reasons: six infinitesimal affine
coordinate changes (`\delta X = -[V,X_0]`, `V` affine), time rescaling
(`\delta X = X_0`), and the two directions along the reversible centre family
(`\partial_a X_0`, `\partial_b X_0`).

Computed symbolically, that `9 x 12` matrix has **rank exactly 9**, and **all
nine directions satisfy the three corrected kernel conditions identically**.

> **The corrected first-order kernel is precisely the span of the trivial
> directions.** There is no blind spot of first-order-invisible mechanisms.
> Everything invisible to first order is, to first order, a coordinate change,
> a time rescaling, or motion along the centre family.

This is the exact opposite of what `Q3R_FIRST_ORDER.md` claimed.

### What this does *not* establish

A perturbation `X_0 + \varepsilon\delta` with `\delta` in the kernel is
conjugate to a member of the centre family **only to first order in
`\varepsilon`**. It may leave the centre variety at second order and produce
cycles. This note therefore does **not** close higher-order reversible
perturbations; it removes a spurious reason for expecting a large reservoir of
them. Any claimed higher-order opportunity must first be separated from
coordinate changes and motion within the integrable family — which is exactly
the quotient computed above.

## 4. Mandatory regression test

The affine shear `y \mapsto y + \varepsilon x` induces

\[
 \delta P = (b-1)x - 2bxy,\qquad
 \delta Q = (a+2)x^2 + \tfrac{b-2}{4} + (1-b)y + by^2 .
\]

Derived independently here: `\delta P = -x\,\partial_yP_0`,
`\delta Q = -x\,\partial_yQ_0 + P_0`. Its `M_1` must vanish. Under the old
claimed conditions it could not, since `q_{20}=a+2 \ne 0` and
`q_{00}=(b-2)/4 \ne 0`. Under the corrected conditions all three residuals are
identically zero, and the coefficient of every generator vanishes.

The mechanism is worth stating: substituting the shear into the *uncorrected*
four-term expression gives `\tfrac{2}{3}(a-1)` times exactly the left side of
the moment identity. **The shear's `M_1` vanishes precisely because the
relation holds** — the coordinate-change direction is the geometric reason for
the dependency.

`quotient.py` and `corrected_rank.py` run this test symbolically.

## 5. Replay

```bash
python3 q3r_2026_09_06/rank_correction.py   # the exact identity y R' + a R = -Phi
python3 q3r_2026_09_06/rank_numeric.py      # relation vs this repo's own generators
python3 q3r_2026_09_06/corrected_rank.py    # rank 3, corrected kernel, shear test
python3 q3r_2026_09_06/quotient.py          # trivial span = kernel, rank 9
```

## 6. Status of the reversible route

Unchanged and open. The correction removes a false reason to expect a large
first-order-invisible reservoir; it does not close the route, and it supplies
no five-cycle candidate. `cheb_hp.py` and `zeros.py` in `q3r_2026_09_06/`
remain in the tree as the record of the failed four-function test, now
labelled.
