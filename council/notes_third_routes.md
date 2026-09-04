# Post-Q4 council: third-route eligibility and endpoint corrections

2026-09-04. Bounded source and algebra audit only; no return-map search,
continuation, new cycle certificate, or full construction attack was run.
The recovered frontier packet is the audit-only publication at `de39ea7`.
This note also checks `FASTRA_COUNCIL_2026_09_04.md` as published at
`2136896`. It does not edit the Claude/Fable records.

## 1. Binding state and corrections to the new council

**PROVED IN THE LOCAL CANON:** [Theorem N](../Q4_THEOREM_N.md) excludes
five distinct interior zeros of the fixed generic Q4 generating integral
for finite `kappa>1`. The bound is three on the strict primitive-lobe
region and four globally. A four-zero original integral is not certified.

**OPEN:** cyclicity of the original Q4 endpoint and its simultaneous
compatibility with interior cycles. The old one-saddle endpoint closure
is withdrawn. The original boundary has two distinct saddles at infinity;
the elliptic covering is singular there. Gavrilov–Iliev's Q4 paper uses a
double ramified covering and states its result for the open annulus.
The essential-perturbation/Nash-space paper explicitly excludes the
closed-annulus/polycycle displacement from its scope. Those results do
not supply a uniform original endpoint theorem.
[GI 2009, section 2](https://arxiv.org/pdf/0811.4602);
[Françoise–Gavrilov–Xiao, section 1.1](https://arxiv.org/html/1610.07582).

The new council correctly reopens this question but needs these repairs:

1. **Linear dependence is neither a kill test nor an alien-cycle test.**
   At a fixed generic center, four independent first-order interior
   coefficient functionals form coordinates on the four-dimensional
   normal-control space. Every scalar first-order endpoint functional on
   that space is therefore their linear combination. This algebraic fact
   supplies no zero count for the full displacement. A uniform bound
   needs a proved displacement normal form, controlled remainders and a
   real-zero theorem in all relevant parameter scalings. Conversely an
   extra independent scalar coefficient would not by itself prove two
   endpoint cycles or compatibility with three interior cycles.
2. **Alien occurrence is an existence result with hypotheses.** The
   Caubergh–Dumortier–Roussarie construction treats particular Hamiltonian
   two-saddle unfoldings, including a persistent connection; it does not
   assert an alien cycle at every two-saddle graphic. Their one-saddle
   exclusion concerns the Hamiltonian setting of that comparison with
   Abelian integrals. It is not a license to identify arbitrary original
   Q4 graphics through a singular covering.
   [CDR 2005](https://comptes-rendus.academie-sciences.fr/mathematique/item/10.1016/j.crma.2005.03.009.pdf).
   In the genuinely quadratic Hamiltonian example, two actual cycles can
   coexist with only one nearby Melnikov zero; the paper separately proves
   a closed-annulus bound of three. The discrepancy concerns the full
   return map, not a count of independent parameter functionals.
   [Gavrilov–Iliev 2015, Theorems 1, 3 and Appendix A](https://www.math.univ-toulouse.fr/~gavrilov/publications/50.pdf).
3. **The stated base is generic, not symmetric.** At the council's
   `lambda=(0,1,3,-10,0,1)`, the generic Darboux nonvanishing condition is
   satisfied. The symmetric Darboux stratum requires `lambda2=0`.
   In the family below it is the degenerate limit `r=0`, outside the
   proposed `r in [1/2,2]` box. Essential order two at that intersection
   does not mean five essential Melnikov coefficients: the cited theorem
   gives four. A symmetric-point attack would require its own exact base,
   covering, graphic and essential unfolding.
   [Buică–Giné–Grau, Lemma 5(iv),(vii), Theorem 6(iv),(vii)](https://arxiv.org/html/1406.7612).
4. **Primitive roots are not original Melnikov roots.** The certified
   three-root boxes in [Q4_LOBE_REGION.md](../Q4_LOBE_REGION.md) concern
   the auxiliary primitive `H`. They do not establish the proposed three
   original interior cycles. Those original integral roots must be
   independently produced and certified in the same parameter direction.
5. **Parameter counting is not a cyclicity theorem.** Even the analytic
   one-parameter displacement `D(s,e)=e product_j(s-s_j)` can have any
   prescribed finite number of simple zeros. This is not a quadratic
   counterexample; it exposes the missing logical implication. Center
   codimension, essential order and dimension of a function space cannot
   by themselves rule out five cycles. In particular neither “three
   transverse directions forbids five” nor “only the symmetric Q4 point
   has enough parameters” is an established exclusion.
6. **Retain numerical-only status for the old infinity samples.** Positive
   splitting at sampled Shi-stratum parameters is not a uniform proof
   that the splitting never vanishes at finite parameters. The frozen
   route is the previously tested/stated construction, not an unproved
   global theorem.

The council's proposed dependent/independent-coefficient cross-review
trigger should be replaced by review of an actual original displacement
normal form or a proved obstruction to the required simultaneous roots.

## 2. Strongest third-route eligibility test: original generic Q4 endpoint

This is a concrete **OPEN compatibility problem**, suitable for a bounded
analytic preflight. It is not a five-cycle mechanism already established.

Use the actual quadratic family

```text
X' = tau X - Y - (2+r^2) X^2 + [2(r+u)+w] XY + Y^2,
Y' = X + tau Y + (r+u) X^2 + [-1-3r^2+v] XY - (r+u) Y^2.
```

The first exact base is `r=1`, `tau=u=v=w=0`, namely

```text
X' = -Y - 3X^2 + 2XY + Y^2,
Y' = X + X^2 - 4XY - Y^2.
```

If a later bounded modulus variation is justified, retain
`r in [1/2,2]`, `max(|tau|,|u|,|v|,|w|) <= 2^-8`.
These are experiment limits, not a global exclusion region.

**EXACT SYMBOLIC.** Put `d=1+r^2`. Its Bautin parameters are
`(tau,r+u,2+r^2,-5d+v,w,1)`. In the normal-control order
`(tau,w,v,u)`, the derivative of the four focus generators at the center is

```text
diag(1,d,-5 r d^2,10 r^2 d^3),
det = -50 r^3 d^6 != 0  for r>0.
```

At `r=1` this is `diag(1,2,-20,80)`. This proves normal-control rank;
it does not prove rank or versality of the original endpoint unfolding.
At `r=0` this same derivative has rank two, illustrating why the generic
and symmetric targets cannot silently be exchanged.

For the base field, in the slope chart `z=Y/X`,

```text
Q2(1,z)-z P2(1,z) = -(z+r)(z^2+2rz-1).
```

The two saddle directions are `z=-r +/- sqrt(d)`, with desingularized
radial/angular eigenvalues `(d,-2d)`. The remaining direction `z=-r`
is a node with `(2d,d)`. Antipodal time orientation changes do not alter
these types. In GI's original coordinates the linear map is
`(x_G,y_G)=(-d X/4,d Y/4)`, so the saddle slopes become
`y_G/x_G=r +/- sqrt(d)`, as in the recovered audit. At `r=1`, use the
algebraic directions `-1 +/- sqrt(2)` directly; no approximate saddle
location is needed.

The compactification and normal-rank identities were independently
checked by tiny symbolic calculations; the root council's
[exact checker](check_council_algebra.py) reproduces them. This is not
verification of the connecting itinerary or a Dulac expansion.

**Conditional five-cycle mechanism.** One admissible small perturbation
arc would need three simple original Melnikov zeros on a compact
subannulus, together with two additional simple return-map roots tending
to the original infinity graphic. The latter two may involve an alien
cycle only if the actual two-saddle displacement proves it. The five must
coexist along the same arc. The alternative four-interior-plus-one-endpoint
route additionally lacks its four-zero original integral premise.

**Bounded first experiment, in dependency order:**

1. At exact `r=1`, identify the original annulus boundary's saddle
   branches, finite connecting orbit and infinity connection. Fix oriented
   transverse sections and write the full return as the appropriate
   original Dulac and regular transitions. Check which connection persists
   under the permitted quadratic perturbations.
2. Derive the genuine separation and saddle-ratio functions in
   `(tau,u,v,w)`, and their transport to the four generating coefficients.
   Record which exact ranks are needed for the proposed unfolding. The
   focus-rank calculation above cannot substitute for this calculation.
   Separation-derivative methods are available, but their hypotheses must
   be verified on this graphic.
   [Marín–Villadelprat 2025](https://link.springer.com/article/10.1007/s12346-025-01379-8).
3. Attempt original Melnikov anchors `t=1/4,1/2,3/4` at that same base:
   verify rank three of the original `3 x 4` evaluation matrix, a nonzero
   null direction, and three simple zeros. No anchored-`H` theorem is being
   applied here. Failure is a failure of this seed, not all Q4.
4. Only if these gates pass, allow the stated modulus box and original
   roots in `[1/8,3/8]`, `[3/8,5/8]`, `[5/8,7/8]`, with gaps at least
   `1/16`. Test a common connection-degeneracy condition while retaining
   all three original roots. Then derive a parameter-dependent endpoint
   return expansion with a remainder strong enough to prove or obstruct
   two further roots. Stop at the finite box if unsuccessful.

A rigorous negative first result could be an obstruction to coexistence
of the interior-root condition and the endpoint two-root sector, or a
uniform closed-annulus zero theorem. Merely finding dependent leading
coefficients, or seeing no extra root in samples, cannot give that result.

**Certificate if a candidate is found.** Freeze one exact rational
quadratic field. Provide five disjoint transverse intervals, validated
complete first returns with the stated compactified itinerary, opposite
displacement signs at interval endpoints, and derivative separation from
one or interval Newton for each root. Prove distinctness and nesting and
independently replay all five at the same parameter vector. An analytic
existence argument must likewise prove simultaneous roots and a parameter
range before it can replace a finite-field certificate.

## 3. Resonant infinity backup and rejected alternatives

The exact resonant-family starting point is

```text
x' = (b-2)/4 + epsilon1 x + (1-b)y + a x^2 + epsilon2 xy + b y^2,
y' = epsilon0 - 2xy,
```

at `(a,b)=(-1,1)`, `epsilon=0`. The 2025 hemicycle theorems give
individual cyclicity two off `a=-1`; their simultaneous count is two or
three in stated nonresonant regions. On `a=-1` they give lower bounds,
not the corresponding upper bounds. The coincidence of the secondary
exponents is a precise missing estimate, not evidence of five cycles.
[Marín–Villadelprat, Theorems B–D and Lemma 3.1](https://arxiv.org/html/2501.16924v1).

At `(-1,1)` the compact double-center cyclicity is two; the theorem
explicitly leaves infinity cycles outside that count.
[Françoise–Gavrilov 2022, Theorem 11 and conclusion](https://arxiv.org/html/2011.08316).
Thus a putative mechanism is two compact cycles plus three endpoint
cycles in the same unfolding, with an allowed distribution such as
`(4,1)`. The bounded first task would be a joint resonant displacement
expansion at this exact base and a comparison of its degeneracy conditions
with the compact double-Bautin conditions. Neither the required three
endpoint cycles nor compatibility with two compact cycles is supplied
here. This remains a backup research question, below the original-Q4
preflight in present readiness.

Do not substitute the naive three-small-plus-two-endpoint proposal at
`a=-1,b=1/5`: the packet's exact first-normal calculation forces
`epsilon0=0`, `epsilon2=-(2/3)epsilon1` for the leading endpoint degeneracy,
whereas the center trace then equals `(2/3)epsilon1`. It obstructs that
common nonzero first direction; it is not an exclusion of every higher
arc or blow-up.

Other suggestions fail a required gate:

- The withdrawn reversible finite-loop-plus-external-focus proposal lacks
  the claimed external focus: its off-axis Jacobians are triangular and
  have real eigenvalues. For the printed samples `k>0,p<0,q>0`, their
  determinant `-2 k q Y_*^2` is negative, contradicting the description
  of those off-axis points as antisaddles. The separate four-real-equilibria
  configuration theorem also obstructs its proposed `(4,1)` neighborhood.
  [The reversible council audit](notes_reversible.md) supplies the detailed
  calculation. A new
  reversible `(5,0)` proposal would need a new common five-cycle mechanism;
  parameter counting alone is not its exclusion.
- The finite Hamiltonian two-saddle alien example is a useful validation
  control, but its applicable closed-annulus theorem prevents five in that
  local construction. Local maxima cannot be added.
- The Swirszcz reversible infinity bicycle is a different graphic. The
  retrieved primary abstract reports cyclicity two or three, but this
  audit did not recover the original curve and complete hypotheses.
  There is no explicit common five-cycle unfolding to execute.
  [Swirszcz 1999, primary record](https://www.sciencedirect.com/science/article/pii/S0022039698935380).
- Continuing the certified Galias–Tucker four-cycle field supplies a
  rigorous incumbent, not a fifth-cycle mechanism or a certified open
  parameter box. Its extreme separated scales increase continuation and
  certification cost relative to the finite-distance KKL proposal.
  [Galias–Tucker 2022](https://www.sciencedirect.com/science/article/pii/S009630032100775X).

## 4. Relative certificate burden

These are ordinal engineering judgments, **not probabilities or existence
evidence**. Here the cost scale is `0 = cheapest`, `5 = hardest`; it is the
reverse of the council's favorable-score convention.

| Route, conditional on a numerical candidate | Cost | Reason |
|---|---:|---|
| KKL Hopf completion | 3/5 | Finite-distance returns and an explicit Hopf completion; large remote return and weak multipliers still need validation. |
| Original generic Q4, three interior plus two endpoint | 4/5 | Original compactification, two saddle passages, small endpoint scales and common-parameter compatibility must all be certified. |
| Resonant infinity, compact two plus endpoint three | 5/5 | Joint resonance remainder and compatibility are missing before a candidate can even reach certification. |
| Reversible loop plus external focus as proposed | Ineligible | The stated equilibrium and configuration mechanism fails before certification. |
| A hypothetical repaired reversible `(5,0)` mechanism | Unrated | No explicit eligible five-cycle target has been supplied. |

The original-coordinate Q4 preflight is preferable to the withdrawn
reversible proposal because it names an actual unclosed graphic,
admissible normal controls and a falsifiable common-unfolding question.
It is not yet preferable to KKL on certificate readiness. No route audited
here supplies an accepted fifth cycle, and no probability of success is
assigned.
