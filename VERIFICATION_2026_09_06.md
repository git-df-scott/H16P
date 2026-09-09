# Verification audit of the 2026-09-06 structural results

Every load-bearing claim in [DEGENERACY_COLLISION.md](DEGENERACY_COLLISION.md)
and [ORDER3_GRAPHIC_NEUTRALITY.md](ORDER3_GRAPHIC_NEUTRALITY.md), re-derived by
hand or re-checked by an independent computation. Two errors were found; both
were in checking scripts, and one omission was found in the write-up. They are
recorded below with their corrections.

## A. Hand-derived, independent of any script

**A1. The Poincare-sphere field.** With `s=(x,y,1)/rho`, `rho=sqrt(1+x^2+y^2)`,
`s3=1/rho`, differentiating gives `s1' = s3[xdot - s1(s1 xdot + s2 ydot)]` and
`s3' = -s3^2(s1 xdot + s2 ydot)`. Substituting `xdot = Pbar/s3^2` with
`Pbar = s3^2 P(s1/s3, s2/s3)` and rescaling time by `s3` (positive on `s3>0`,
so orientation is preserved) yields exactly

```
s1' = Pbar - s1 W,   s2' = Qbar - s2 W,   s3' = -s3 W,   W = s1 Pbar + s2 Qbar.
```

**A2. Eigenvalues at an infinite singularity.** In chart U1 (`x=1/z`, `y=u/z`):
`udot = z(Q-uP)` and `zdot = -z^2 P`. Expanding by homogeneous parts and
rescaling by `z` gives `udot = G(1,u) + O(z)`, `zdot = -z P_2(1,u) + O(z^2)`.
The Jacobian at `(u*,0)` is upper triangular, so the eigenvalues are

```
lam_eq = G'(u*)  (along the equator),   lam_tr = -P_2(1,u*)  (transverse).
```

For the Shi chart `P_2(1,u)=l+mu+u^2`, `Q_2(1,u)=a+bu`, so
`G(1,u)=a+(b-l)u-mu^2-u^3`, `lam_eq=(b-l)-2mu-3u^2`, `lam_tr=-(l+mu+u^2)`.

**A3. The antipodal lemma.** `Pbar` and `Qbar` are homogeneous of degree 2, so
even; `W` is homogeneous of degree 3, so odd. Hence `F(-s)=F(s)`. With
`g(e)=F(-s+ev)=F(s-ev)` one gets `DF(-s)v = -DF(s)v`, so `DF(-s) = -DF(s)`:
antipodal infinite singularities carry negated linearisations.

**A4. The neutrality factorisation.** With `lam_eq=Au^2+Bu+C`,
`lam_tr=Du^2+Eu+F`, expanding `lam_eq(u1)lam_tr(u2)-lam_tr(u1)lam_eq(u2)`
term by term leaves exactly

```
(u1-u2) * [ (AE-BD) e2 + (AF-CD) e1 + (BF-CE) ],   e1=u1+u2, e2=u1u2,
```

and for the Shi chart `AE-BD = m`, `AF-CD = 2l+b`, `BF-CE = m(l+b)`.
These match the symbolic output.

**A5. OMISSION IN THE ORIGINAL WRITE-UP (now corrected).** `r=1` means
`lam_eq(u1)lam_tr(u2) = +/- lam_tr(u1)lam_eq(u2)`, and only the `+` branch was
used. This is legitimate but the reason was not stated: at a saddle `lam_eq`
and `lam_tr` have opposite signs, so the *signed* ratio is negative at both
saddles, and the `-` branch would require them to have opposite signs. It is
unreachable for two saddles. The derivation is complete; the argument is now
recorded.

**A6. `N = 640 eta_3` checked by hand.** At Shi's seed `l=-10, a=1` (so `m=5`,
`b=-25`), summing the fourteen terms of `N`:

```
125 + 15625 + 37500 + 30000 + 18750 + 15625 + 8000 + 15000 + 6250
+ 1953125 + 2343750 - 781250 - 500000 - 312500  =  2,850,000
```

and `eta_3 = 285000/64 = 4453.125`, so `640*eta_3 = 2,850,000`. Exact match.

**A7. The collision point is exactly a centre.** At `l=-10, a=2, m=10, b=-25`:
`eta_1 = ab+2al-lm-m = -50-40+100-10 = 0`; the `eta_2` bracket sums to
`-600-480-1500-20-15625+6250-2500+15000-2750+2500+125-600+100+100 = 0`; and
`2a^2+l+2 = 8-10+2 = 0`, so `eta_3 = 0`. Three vanishing Lyapunov quantities,
so a centre by Bautin — no numerics required.

**A8. `(0,1)` is a node.** `J = [[m,1],[1+b,0]] = [[10,1],[-24,0]]`: `det=24`,
`tr=10`, `tr^2-4det=4>0`. Matches the symbolic `det=-3l-6`,
`tr=(5/2)sqrt(-2l-4)`, `tr^2-4det = -(l+2)/2 = a^2 > 0` along the whole curve.

## B. Independent numerical checks

**B1. The Lyapunov chain, against the actual return map.** Measuring
`d(x) ~ C x^(2k+1)` by direct integration:

| case | predicted exponent | measured |
|---|---|---|
| generic, `eta_1 != 0` | 3 | 3.0104, 3.0058 |
| `eta_1 = 0`, `eta_2 != 0` | 5 | 5.0315, 5.0108, 5.0212, 5.0130 |
| order-3 stratum | 7 | 7.3817, 7.4689, 7.2996 |

(The order-7 fits are the least clean, as expected from higher-order
contamination at the sample radii used.)

**B1a. ERROR FOUND, in the checking script.** The first run reported sign
mismatches. Cause: one test point had `V_4 = 0` exactly, so it was not the
generic case it was meant to test; and my formal-integral convention is the
negative of the repository's inherited `eta_3`. Resolved exactly:
`V_8 = -eta_3` on the order-3 stratum, verified symbolically with ratio `-1`.
With that single global sign, **all** valid cases agree: `sign(d) = -sign(V_k)`
throughout. Since the resultant identity is a statement about zero sets, it is
unaffected. As a by-product `V_4` and `V_6` vanish identically on `m=5a`,
`b=3l+5`, independently confirming that stratum.

**B2. Neutrality from the sphere Jacobian**, independent of the chart formulas.
At `l=-10, a=2`: `|lam_eq/lam_tr| = 2.0000000002` and `2.0000000001` at the two
saddle directions, so `r = 1.000000000020`.

**B2a. ERROR FOUND, in the checking script.** The first version paired a saddle
with the *antipode* of the other and reported `r = 4`. The graphic runs
`equator -> S1 -> plane -> S2 -> equator` and so uses the two saddles on the
same side. Corrected in `verify_geometry.py`; the result is `r = 1` as the
chart formula predicts.

**B3. The splitting zero is genuine, not a discontinuity.** This was the
weakest point, since the phase portrait changes qualitatively at `a_deg`. A
21-point scan across `a_deg = 2` at `l = -10` shows the splitting varying
smoothly and monotonically from `-7.57e-4` to `+6.09e-4`, crossing zero at
`a = 2.000000` with value `2.08e-16`. The saddle directions move smoothly, both
hyperbolicity ratios stay equal, and no branch selection flips. It is a
transversal crossing.

**B4. The graphic really does bound the nest.** At the collision point the
origin is a centre, so it has a genuine period annulus; its outer boundary
crosses `y=0` at `x* = 0.06392`, and the boundary orbit approaches **both**
infinite saddles (`2.08e-3` and `8.51e-3`) while staying far from the node
(`3.37e-1`). The annulus boundary is the two-saddle graphic, as assumed.

**B5. `a* = a_deg` to 14 digits** at `l = -20, -14, -10, -7, -5, -3`
(differences `6.6e-14` down to `1.4e-14`).

## C. What remains unproved

1. **`a* = a_deg` is numerical, not symbolic.** This is the one load-bearing
   claim without an exact proof. It cannot be settled by resultants: the
   splitting is a transcendental function of the flow, not an algebraic one.
   The direction "centre implies the graphic exists" is easy — a centre has a
   period annulus and its boundary is the graphic. The converse, "the
   connection forces a centre", is the surprising half and is currently
   supported only by the table in B5 and the transversality in B3.
2. **"Cyclicity adds" is a bookkeeping heuristic**, not a theorem about
   simultaneous unfolding.
3. **Scope.** One chart (Shi), one graphic type (elementary, through two
   non-antipodal infinite saddles). Graphics through finite saddles, nilpotent
   points and saddle-nodes — the hard end of the DRR list — are untouched.
4. **The centre class of the collision curve is not identified.** It is not
   Hamiltonian, has no invariant straight line, is codimension four, and its
   second finite singularity is a node — all consistent with `Q_4`, which the
   repository's Theorem N already caps at four zeros. Not confirmed.

## D. Replay

```bash
python3 structure_2026_09_06/verify_lyapunov.py   # B1: return-map exponents and signs
python3 structure_2026_09_06/verify_chain.py      # B1a: V_8 = -eta_3 on the stratum
python3 structure_2026_09_06/verify_geometry.py   # B2, B4: sphere-Jacobian r, annulus boundary
python3 structure_2026_09_06/verify_splitting.py  # B3: transversality across a_deg
python3 structure_2026_09_06/collision.py         # B5: a* vs a_deg
```
