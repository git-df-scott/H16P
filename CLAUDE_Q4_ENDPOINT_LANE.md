# Fable lane: Q4 closed annulus at the two-saddle infinity graphic

2026-09-04, after the FASTRA council. Scripts: `audit/claude_q4_original_infinity.py`,
`claude_q4_graphic_itinerary.py`, `claude_q4_boundary_level.py`,
`claude_q4_endpoint_c0.py`, `claude_q4_crossing_count.py`,
`claude_q4_tworoot_scan.py`, `claude_q4_triple_center_face.py`,
`claude_q4_triple_center_face2.py`, `claude_q4_triple_center_point.py`.
Everything below is HIGH-PRECISION NUMERICAL unless marked EXACT.

## 1. The boundary graphic in original coordinates (rho = 1, kappa = 2)

- EXACT: the rational first integral `Hcal=[8y(1+Yv)-\frac23(1+\kappa Yv^3)]^2/[1-8y+\kappa Yv^2]^3`,
  `Yv=cx-(2+b)y`, satisfies `dHcal/dt\equiv0`; `Hcal(0,0)=4/9`; the loop
  value is `4/(9\kappa)=h_{loop}^2`. The level set `Hcal=4/(9\kappa)` is a
  quartic whose top-degree part factors as `768(v-1)^2(v^2/2-v-1/2)`: real
  asymptotic slopes `1` (double, the infinity node) and `1\pm\sqrt2` (the two
  infinity saddles, eigenvalues `(8,-4)`, hyperbolicity ratios `2` and `1/2`
  in the traversal sense, product one, each `1:2` resonant).
- The annulus boundary on the positive `x`-axis is `x^*=0.2272111321`
  (bisection, agreeing with the quartic's root `0.22721128717`). Tracing the
  level-curve component through it: one end goes to infinity along slope
  `-0.414`, the other along slope `2.414`; it crosses the negative axis at
  `-0.1586` and comes within `0.078` of the origin. Itinerary confirmed:
  **finite heteroclinic orbit between the two infinity saddles around the
  center, closed by the equator arc between them** (the arc not containing
  the node). Direct separatrix integration is numerically unreliable here
  (the connection is codimension-one and integrability-forced); the level
  curve is the right tool.

## 2. First-order compatibility face `c_0=0`

`I(\text{loop})=-\frac{aC}{2}\sqrt{1-a}\,X(1)`, `X(1)=\int_0^1Y/(1-au)^{3/2}`,
so the first-order Melnikov value at the graphic vanishes exactly on the
face `X(1)=0`. Since `Y_0<0` on the lobe region, `I>0` near the center, hence
`X(1)<0` gives an even interior count and `X(1)>0` an odd one.

- On the lobe region `X(1)` **does** change sign as a function of `a`: for
  nine of ten anchor triples the face is crossed at `a^*\in(0.87,0.99)`
  (e.g. `(0.46,0.6,0.8)`: `a^*=0.907`; `(0.5,0.75,0.9)`: `0.928`;
  `(0.9,0.95,0.99)`: `0.985`). The late-root triples `(0.99,0.999,0.9999)`
  and the corner point keep `X(1)<0` for all `a\le0.99`.
- Across the face the interior count goes **0 to 1**: the zero enters from
  the graphic (at `t\approx0.996` just past `a^*`). Consistent with the
  Dulac form `I\approx c_0+c_1w\log w`, `c_1\propto H(1)<0`: a near-boundary
  zero exists exactly on the side `c_0c_1>0`.
- Parity consequence: on the lobe region, `X(1)<0` forces at most two
  interior zeros, and every three-zero configuration would need `a>a^*`.

## 3. Closed-interval count never reaches four (first order)

- Lobe region: sampled counts are 0 (before the face) and 1 (after).
- Two-root region (Astra's fibre), 546 samples over `r,s,\eta,a`: interior
  counts 0 (315) and 1 (231); never 2, 3 or 4.
- Near-center hierarchy on the face (`Y_0<0` tiny, `\eta<0` small,
  `X(1)=0`): exactly **two** small interior zeros plus the boundary zero,
  closed count three, at every `a` and every hierarchy tested; never three
  interior.
- The unique triple-center point `Y_0=Y_1=Y_2=0` is
  `(A,B,\eta)=(1,-17/12,0)` (EXACT, for all `a`). There `X(1;a)>0` for all
  `a\in[0.05,0.995]` (from `1.4\cdot10^{-4}` to `0.97`), so **no lift makes
  the triple center zero coexist with a boundary zero**. The
  codimension-four "three small plus one at the graphic" configuration does
  not exist in Q4.

## 4. Consequence for the (5,0) endpoint target

Five cycles from the closed annulus need first-order closed count `n_1`
plus alien cycles `n_a` at the two-saddle graphic with `n_1+n_a=5`. Theorem
N gives `n_1\le4` (four interior, or three interior plus a boundary zero).
The computations above find `n_1\le3` in every configuration reachable
(generic sampling gives `\le1`; the best hierarchical construction gives
`2+1`). Hence the target needs **at least two alien cycles** at the graphic,
one more than the only known alien mechanism produces (Dumortier–Roussarie:
cyclicity two where the Melnikov count is one, and Gavrilov–Iliev 2015 cap
the Hamiltonian two-saddle closed annulus at three, equal to the parameter
count). No theorem forbids two aliens here; nothing suggests them either.

Lane status: **first-order compatibility gate failed at every tested
configuration; the route survives only as a two-alien hypothesis.**
This is numerical evidence over explicit samples, not a theorem.

## 5. What would revive the lane

1. A four-zero Q4 integral (Astra's two-root determinant criterion): then
   `n_1=4` and a single alien would suffice. The scans here found no
   two-root point with more than one interior zero.
2. A proof or disproof of two-alien cyclicity at a `1:2`-resonant
   two-saddle infinity graphic under a four-parameter unfolding. This is
   the analytic problem the council named; it is now the only content of
   the lane.
