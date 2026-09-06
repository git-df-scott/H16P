# The degeneracy-collision obstruction: why the fifth cycle keeps not being there

**Date: 2026-09-06.  Exact algebra plus high-precision NUM.  This is not a
proof of `H(2)=4` and produces no counterexample.  It is a structural
explanation of why every route this campaign has tried has failed in the same
way, and a statement of exactly which routes it does not cover.**

## 0. The reframe

`H(2) < infinity` is equivalent to every limit periodic set of the compactified
quadratic family having finite cyclicity — this is the Dumortier–Roussarie–Rousseau
program, which reduced it to 121 graphics in 1994 and is still open.  The same
catalogue governs `H(2)=4`, because five cycles have to be *born* somewhere and
Bautin caps a focus or centre at three.  So a counterexample needs

\[
 \underbrace{\text{cyclicity of the focus}}_{\le 3\ \text{(Bautin)}}
 \;+\;\underbrace{\text{cyclicity of the nest's boundary graphic}}_{\ge 2}
 \;\ge\; 4\ \text{in one nest},\qquad +\,1\ \text{remote}.
\]

The campaign has been searching parameters.  The question is really a finite
enumeration over limit periodic sets, and a **codimension budget**: each cycle
costs one independent degeneracy condition, and the chart has five parameters.

The known four-cycle systems spend exactly four: `lambda=0`, `eta_1=0`,
`eta_2=0` (three cycles from an order-3 weak focus) plus one Hopf condition at
the remote focus.  A fifth cycle needs a **fifth independent condition**.  This
note shows that for the boundary graphics reachable in the Shi chart, that
fifth condition is not independent — it collides with the ones already spent.

## 1. Endpoint theorems (literature, verified)

| | |
|---|---|
| Li (1986), Cherkas (1986) | **no** limit cycle surrounds a third-order weak focus in any quadratic system |
| same line of work | **at most one** limit cycle surrounds a second-order weak focus, hyperbolic if it exists |
| Zhang / Zegeling | with two nests, one of them holds at most one cycle |
| Dumortier–Roussarie–Rousseau 1994 | finiteness reduces to 121 graphics; still open |
| Rousseau et al. | elementary graphics surrounding a focus or centre have finite cyclicity (value not pinned) |

**Parity consequence.**  On a transversal from the focus to a hyperbolic
boundary graphic `Gamma` with ratio `r != 1`, the number of cycles has parity
fixed by whether the focus's stability and `Gamma`'s stability agree.  Li–Cherkas
forces agreement at an order-3 weak focus, so the count there is **even** beyond
the three Bautin cycles: a lone fourth cycle in that nest is parity-forbidden,
and a *pair* requires `r = 1`.  That is why neutrality is the whole question.

## 2. The neutrality condition, exactly

Shi chart `xdot = lam x - y + l x^2 + m x y + y^2`, `ydot = x + a x^2 + b x y`.
Infinite directions are the roots of `G(u) = a + (b-l)u - m u^2 - u^3`, with

```
lam_eq(u) = -3u^2 - 2m u + (b-l)        (along the equator)
lam_tr(u) = -(u^2 + m u + l)            (transverse, into the plane)
```

Neither involves `lam`: neutrality lives on `(l,m,a,b)` alone.

**Exact lemma (every quadratic field).**  The Poincaré-sphere field is even,
`F(-s) = F(s)`, and `d(iota) = -I`, so `DF(-s) = -DF(s)`.  Antipodal infinite
singularities carry negated linearisations, hence reciprocal ratios, so a
graphic through an *antipodal* pair is neutral identically — no information.
The informative case is two **non-antipodal** infinite saddles, where

\[
 r(\Gamma)=\Bigl|\tfrac{\lambda_{eq}}{\lambda_{tr}}\Bigr|(u_1)\cdot
           \Bigl|\tfrac{\lambda_{tr}}{\lambda_{eq}}\Bigr|(u_2),
\]

i.e. neutrality says the two saddles have **equal hyperbolicity ratios**.
(`r=1` allows `lam_eq(u1)lam_tr(u2) = +/- lam_tr(u1)lam_eq(u2)`; only `+` is
reachable, because at a saddle `lam_eq` and `lam_tr` have opposite signs, so the
signed ratio is negative at both and the `-` branch would need them to differ.)
Eliminating the roots gives the exact neutrality polynomial

```
N(l,m,a,b) = a^2 m^3 - a b^3 - 6 a b^2 l - 12 a b l^2 + 3 a b l m^2 - a b m^4
           - 8 a l^3 + 6 a l^2 m^2 - a l m^4 + b^4 m + 3 b^3 l m + b^2 l m^3
           - 4 b l^3 m + b l^2 m^3
```

which does not factor.

## 3. Collision at order 3 (exact)

On the order-3 stratum `m=5a`, `b=3l+5`:

\[
 N\big|_{\text{stratum}} \;=\; 640\,\eta_3 .
\]

**Neutrality is exactly `eta_3 = 0`.**  So a cyclicity-2 boundary graphic and a
genuine third-order weak focus cannot coexist: the fifth condition is the third
one over again.  (Full derivation: `ORDER3_GRAPHIC_NEUTRALITY.md`.)

## 4. Collision at order 2 (exact family, numerical connection)

Drop to the cheaper target: an **order-2** weak focus (2 Bautin cycles) plus a
neutral graphic (2 more) — four in one nest at codimension 4 rather than 5.
On `{eta_1 = 0}` the neutrality polynomial splits,

```
N   |_{eta_1=0} = m^3 * S * W ,      S = a^3 - a^2 m + 3 a l + a - l m - m
eta_2|_{eta_1=0} = m * (5a - m) * W
```

and the branch `S = 0` gives neutrality **without** `eta_2 = 0`.  Solving it:

\[
 m=\frac{a\,(a^{2}+3\ell+1)}{a^{2}+\ell+1},\qquad
 b=-2\ell+\frac{m(\ell+1)}{a},
\]

an explicit two-parameter family on which `eta_1 = 0`, `r(Gamma) = 1` to 12
digits, `eta_2 != 0`, the origin is a weak focus of order exactly two, and
`(0,1)` is a genuine second focus.  13 of 40 grid points pass every
non-degeneracy check.  **The order-2 route survives every algebraic test.**

It dies on the last one.  The graphic also has to *exist* — its two transverse
separatrices must connect.  Bisecting that splitting along the family:

| `l` | `a_deg = sqrt(-(l+2)/2)` (where `eta_2=0`) | `a*` (splitting zero) | difference |
|---|---|---|---|
| −20 | 3.000000000000 | 3.000000000000 | `+6.6e-14` |
| −14 | 2.449489742783 | 2.449489742783 | `−1.4e-14` |
| −10 | 2.000000000000 | 2.000000000000 | `−2.6e-14` |
| −7 | 1.581138830084 | 1.581138830084 | `−3.1e-14` |
| −5 | 1.224744871392 | 1.224744871392 | `−3.4e-14` |
| −3 | 0.707106781187 | 0.707106781187 | `−3.2e-14` |

**The connection is realised only where `eta_2 = 0`.**  Same collision, one
level down.

## 5. What sits at the collision

`lam=0`, `a^2 = -(l+2)/2`, `m = 5a`, `b = 3l+5` — a curve on which
`eta_1 = eta_2 = eta_3 = 0`, so by Bautin the origin is a **centre**.  Confirmed:
the return map is the identity to `1.4e-17` at `l=-3, a=1/sqrt(2), m=5/sqrt(2),
b=-4`.  At `l=-10` it is the rational field

\[
 \dot x=-y-10x^{2}+10xy+y^{2},\qquad \dot y=x+2x^{2}-25xy .
\]

The second finite singularity `(0,1)` has, symbolically along the whole curve,

\[
 \det = -3\ell-6,\qquad \operatorname{tr}=\tfrac{5}{2}\sqrt{-2\ell-4},\qquad
 \operatorname{tr}^2-4\det = -\tfrac{\ell+2}{2}=a^{2}>0,
\]

so it is **always a node** — no Hopf, so the remote `+1` is gone too.  The
repository's own audit records a node as the second finite singularity of every
`Q_4` centre (`audit/claude_center_identify.py`), and this curve is a
codimension-four centre stratum, not Hamiltonian and with no invariant straight
line.  If it is `Q_4`, Theorem N already caps it at four zeros and the picture
closes on itself.

## 6. The obstruction, stated

> Within the Shi chart, for a nest bounded by an elementary graphic through two
> non-antipodal infinite saddles, the degeneracy conditions of the focus and of
> the graphic are **not independent**.  Imposing graphic neutrality on the
> order-3 stratum forces `eta_3 = 0`; imposing neutrality plus the connection on
> the order-2 stratum forces `eta_2 = 0`; the common solution is a centre whose
> second singularity is a node.  The five independent conditions a fifth cycle
> needs collapse to four — which is what the known four-cycle systems already
> spend.

This is the first mechanism this campaign has produced that *explains* the
repeated near-misses instead of recording another one.

## 7. What it does not cover — the live surface

The obstruction is proved for one chart and one graphic type.  Everything below
is untouched and is where a counterexample would still have to live:

1. **Non-elementary graphics.**  Through nilpotent points, saddle-nodes, and
   degenerate graphics — precisely the hard end of the DRR 121 list.  Nothing
   here applies to them.
2. **Graphics through finite saddles** (homoclinic loops, heteroclinic
   polycycles).  The collision was tested only at infinity.
3. **Charts outside Shi.**  `SEARCH_SPACE.md` already records that no single
   five-parameter chart covers every affine class.
4. **The `a* = a_deg` dependency is numerical** — 14 digits at six values of
   `l`, not eliminated symbolically.  A resultant computation should settle it
   and would upgrade section 4 to a theorem.
5. **"Cyclicity adds" is a heuristic**, an upper-bound bookkeeping device, not
   a theorem about simultaneous unfolding.
6. **Higher-order Melnikov / the Bautin ideal** — the eight-dimensional
   first-order blind spot of `Q3R_FIRST_ORDER.md` is a different hole and this
   note does not touch it.

## 8. Verification

Every claim above was re-derived by hand or re-checked independently; the audit,
including two errors found in checking scripts and one omission in this
write-up, is [VERIFICATION_2026_09_06.md](VERIFICATION_2026_09_06.md). The one
load-bearing claim still without an exact proof is `a* = a_deg` (section 4): it
is transcendental, not algebraic, so resultants cannot settle it.

## 9. Replay

```bash
python3 structure_2026_09_06/general_neutrality.py   # N(l,m,a,b), exact
python3 structure_2026_09_06/lyapunov.py             # eta_1, eta_2, eta_3 in the Shi chart
python3 structure_2026_09_06/compare.py              # the splitting of N on {eta_1=0}
python3 structure_2026_09_06/order2_neutral.py       # the explicit order-2 + neutral family
python3 structure_2026_09_06/splitting_gen.py        # graphic splitting
python3 structure_2026_09_06/collision.py            # a* vs a_deg, and the centre test
python3 structure_2026_09_06/second_sing.py          # (0,1) is a node along the curve
```

Dependencies: NumPy, SciPy, SymPy.  No interval certificate is claimed and no
ODE budget was charged to the KKL/Shi ledger.
