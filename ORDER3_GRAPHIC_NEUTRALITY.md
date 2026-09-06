# The order-3 stratum's boundary graphic: neutrality is exactly `eta_3 = 0`

**Date: 2026-09-06.  Evidence class: exact algebra (no ODE returns charged),
with a 50-digit numerical confirmation and one bounded NUM scan.**

This answers questions 2, 3 and 4 of section 4b of
[CLAUDE_THOUGHT_SESSION.md](CLAUDE_THOUGHT_SESSION.md), the item recorded there
as "route 1 in its only surviving form".  The answer closes that route on this
stratum.  **No five-cycle candidate was produced and none is claimed.**

## Setting

Shi chart, third-order weak-focus stratum (`lambda=0`, `m=5a`, `b=3l+5`):

\[
\dot x=-y+\ell x^2+5axy+y^2,\qquad \dot y=x+ax^2+(3\ell+5)xy .
\]

Poincare compactification.  With `s=(s1,s2,s3)`, `x=s1/s3`, `y=s2/s3`, and
`Pbar,Qbar` the degree-2 homogenisations, the sphere field rescaled by `s3` is

\[
 s'=\bigl(\bar P-s_1W,\;\bar Q-s_2W,\;-s_3W\bigr),\qquad W=s_1\bar P+s_2\bar Q .
\]

Infinite singularities lie at the real roots of

\[
 G(u)=a+(2\ell+5)u-5au^2-u^3 ,\qquad u=y/x,
\]

with eigenvalues along and transverse to the equator

\[
 \lambda_{\rm eq}(u)=-3u^2-10au+2\ell+5,\qquad
 \lambda_{\rm tr}(u)=-(u^2+5au+\ell).
\]

The coefficient of `y^3` in `xQ_2-yP_2` is `-1`, so the direction `x=0` is never
singular and every infinite singularity is captured by `G`.

## Result 1 (general, exact): antipodal infinite saddles are reciprocal

For **every** quadratic field the sphere field above satisfies `F(-s)=F(s)`,
and the antipodal map has `d(iota)=-I`, hence

\[
 DF(-s)=-DF(s).
\]

Both identities are verified symbolically in twelve free coefficients by
`neutrality_2026_09_06/antipodal_general.py`.  So the linearisations at an
antipodal pair of infinite singularities are exact negatives, and a graphic

\[
 A\;\xrightarrow{\ \text{plane}\ }\;B\;\xrightarrow{\ \text{equator arc}\ }\;A
\]

through an antipodal pair `A,B` enters `A` along the equator and leaves it
transversally, and enters `B` transversally and leaves it along the equator.
Therefore

\[
 r_A=\Bigl|\tfrac{\lambda_{\rm eq}}{\lambda_{\rm tr}}\Bigr|(u^*),\qquad
 r_B=1/r_A,\qquad r(\Gamma)=r_Ar_B\equiv 1 .
\]

**The first stability coefficient of such a graphic is identically neutral.**
It is not a curve in parameter space; it is everything.  The premise of
question 3 of section 4b — "is the neutrality curve nonempty" — is therefore
the wrong question whenever the graphic's two saddles are antipodes.

## Result 2 (NUM): in that case the graphic does not exist on the stratum

The Attack-2 box `l in [-12,-8]`, `a in [4/5,6/5]` has exactly **one** real
infinite direction, so a graphic there is forced to use an antipodal pair.
Its splitting function — the separation, on the section `{x=0}`, between the
transverse separatrix leaving `A` and the one entering `B` — was evaluated at
268 stratum points covering `l in [-30,6]`, `a in [0.2,4]` restricted to the
single-direction region:

| | |
|---|---|
| splitting range | `[-9.8578, -0.7063]` |
| sign changes | none |
| closest approach to zero | `-7.06e-1` at `l=-30, a=0.2` |

The splitting is bounded away from zero throughout.  **The connection never
occurs there**, so the identically neutral graphic is not realised, and the
boundary of the origin nest in that region is not this graphic.  This is
numerical evidence over a scanned window, not a theorem of absence.

## Result 3 (exact): where the graphic has two non-antipodal saddles,
## neutrality is exactly `eta_3 = 0`

Where `G` has three real roots, two of them are saddles that are *not*
antipodes, and the graphic

\[
 \text{equator arc}\to S_1\to\text{plane orbit}\to S_2\to\text{equator arc}
\]

has first stability coefficient
`r=|\lambda_{\rm eq}/\lambda_{\rm tr}|(u_1)\cdot|\lambda_{\rm tr}/\lambda_{\rm eq}|(u_2)`.
Writing `lambda_eq = Au^2+Bu+C`, `lambda_tr = Du^2+Eu+F`, the neutrality
condition factors as

\[
 \lambda_{\rm eq}(u_1)\lambda_{\rm tr}(u_2)-\lambda_{\rm tr}(u_1)\lambda_{\rm eq}(u_2)
 =(u_1-u_2)\bigl[(AE-BD)e_2+(AF-CD)e_1+(BF-CE)\bigr]
\]

with `e_1=u_1+u_2`, `e_2=u_1u_2`.  Eliminating the third root through the Vieta
relations of `u^3+5au^2-(2\ell+5)u-a` gives the resultant

\[
 -250\,a\,(2a^2+\ell+2)\,(5a^2\ell+6a^2-3\ell^3-12\ell^2-15\ell-6)\;=\;640\,\eta_3 .
\]

**The neutrality set is exactly the zero set of `eta_3`.**  Both factors of
`eta_3` appear, with no extra component.  Confirmed independently at 50 digits
by `check_curve.py`: on `l=-2-2a^2` the computed `r-1` is `O(10^{-51})` at
`a=3,3.5,4,4.5,5,7,10`, and it crosses zero transversally in `l` off the curve
(`r-1 = -0.00330` at `l=-34.05`, `+0.00331` at `l=-33.95`, for `a=4`).

## Consequence: route 1 of section 4b is closed on this stratum

Section 4b asked whether the neutrality curve is nonempty, and (question 4)
whether the two graphic unfolding directions are independent of the three focus
directions `lambda, delta, epsilon`.  The identity answers both:

- At every point of the stratum with a genuine third-order weak focus
  (`eta_3 != 0`) the boundary graphic is **hyperbolic**, `r != 1`, with a
  definite sign.  It cannot emit a pair of cycles, so the fourth and fifth
  cycles of the origin nest cannot come from it.  This is exactly the sign
  agreement Li--Cherkas forces, now with the mechanism identified.
- The graphic-stability direction is **not** independent of the focus
  directions: it is the `eta_3` direction.  The hoped-for codimension-five
  point, where three focus directions and two graphic directions unfold
  independently to give `3+2=5`, does not exist by this mechanism — the two
  conditions are one equation.

So on the third-order weak-focus stratum the origin nest can only ever give
three plus one, which is the alternative section 4b itself named as the
route's death.

## Scope and what is not claimed

- Results 1 and 3 are exact algebra.  Result 2 is NUM over a scanned window.
- Result 3 identifies the neutrality set of the two-root condition.  The
  elimination does not by itself certify which pair of roots is the saddle
  pair; that identification is numerical (`neutrality_set.py`), and the
  50-digit check is on the saddle pair.
- Whether the two-saddle graphic is actually *realised* on `eta_3=0` (the
  required connections exist) was not computed.  It does not affect the
  conclusion: the route needs neutrality **and** `eta_3 != 0`, and those are
  incompatible.
- Cyclicity statements use the standard hyperbolic-polycycle theory for a
  graphic with regular transitions; the transitions were not re-verified here.
- Nothing here bears on the other quadratic center classes, on higher-order
  Melnikov functions, or on any other attack in `ATTACK_MATRIX.md`.

## Replay

```bash
python3 neutrality_2026_09_06/inf_singularities.py     # exact chart quantities
python3 neutrality_2026_09_06/antipodal_general.py     # Result 1, 12 free coefficients
python3 neutrality_2026_09_06/classify_infinity.py     # infinite singularity types
python3 neutrality_2026_09_06/splitting.py             # Result 2, Attack-2 box
python3 neutrality_2026_09_06/splitting_wide.py        # Result 2, wide scan
python3 neutrality_2026_09_06/neutrality_set.py        # locates r=1 numerically
python3 neutrality_2026_09_06/neutrality_proof.py      # Result 3, the resultant
python3 neutrality_2026_09_06/check_curve.py           # 50-digit confirmation
```

Dependencies: NumPy, SciPy, SymPy, mpmath.  The ODE work uses the
drift-stabilised sphere field in `sphere3.py`; `|s|-1` stays below `4e-12`
over the reported integrations.  No interval certificate is claimed.
