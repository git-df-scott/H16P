# D2: the order-two loop law (toward a theorem)

2026-09-05. Status: CONJECTURE with numerical evidence; proof program below. Nothing here is proved.

## 1. Statement

Work in the Shi chart at zero trace,

    x' = -y + l x^2 + m x y + y^2,     y' = x + a x^2 + b x y,

with the first focal value zero, i.e. m = a(b+2l)/(l+1). The second focal value is (exact)

    eta_2 = a (b+2l) (b-3l-5) [ a^2 (b+2l+1) - (b+1)(l+1)^2 ] / (48 (l+1)^2).

Its three factors are the three center strata through the order-two locus: a(b+2l) = 0 (reversible),
b = 3l+5 (the Shi line, order three or center), and the third factor (the reversible center family C3 met
in F11).

**Conjecture D2 (order-two loop law).** Suppose the origin is a weak focus of order two (eta_2 != 0) and
there is a homoclinic loop Gamma through a hyperbolic saddle S with the origin the only singularity inside
Gamma. Let sigma = div X(S) be the saddle quantity. Then

    sigma * eta_2 < 0.

Consequences. (i) Gamma has the opposite stability to the focus, so the open region between them is
filled with spirals and contains no limit cycle: the codimension-three configuration "order-two focus plus
loop" is cycle-free. (ii) A neutral loop (sigma = 0) around an order-two focus is impossible unless eta_2 = 0,
i.e. the origin is a center; this is exactly what F11 observed. (iii) Together with Zhang (1999, at most one
hyperbolic cycle around an order-two focus) and Li-Cherkas (order three: no cycle), it is the missing
middle rung of the "3 - k" ladder and closes the neighbour of the K1 question.

## 2. Evidence (NUMERICAL)

46 focus-type loops located by separatrix splitting in F11 (`audit/fable_f11_neutral_loop.py`), a in
[-3, 3], b in [-4, 4]: sigma * eta_2 < 0 in every case. The ratio sigma/eta_2 ranges from -0.027 (a = -3)
to -33 (a = -0.5) and is not constant, so the law is a sign law, not a proportionality. It is exactly
antisymmetric under a -> -a, as the symmetry (x, t) -> (-x, -t) requires. Along the a = -3 branch the ratio
was constant to four digits and the branch crossed the C3 center stratum at b = 1.275 with sigma and eta_2
vanishing together.

## 3. Exact ingredients (SYMBOLIC, checked with sympy in this session)

The saddle S lies on the line 1 + a x + b y = 0 (the y' = 0 factor). Its coordinates are algebraic in
(a, b, l) with the discriminant

    Delta = a^2 (b + 3l + 1)^2 - 4 l (l+1)^2 (b+1)      (after expansion of the sympy output),

and its trace factors as

    sigma = -(b+2l) * B_pm(a, b, l, sqrt(Delta)) / ( 2 (l+1) * [a^2 (b+2l+1)(b-1) ... ] )

with the common factor (b+2l) in front. Since eta_2 also carries a(b+2l), the sign law is equivalent to a
sign statement about a * (b-3l-5) * [third factor] * B_pm / [denominator] on the loop locus. The loop locus
itself is transcendental, so the law cannot be a pure algebraic identity; it must use the loop.

## 4. Proof program

1. **Divergence along the loop.** Compute div X along Gamma for the 46 loops (`audit/fable_d2_loop_div.py`).
   If div changes sign exactly once on Gamma and the loop integral of div has the sign of sigma, the loop's
   stability is governed by sigma alone and the problem reduces to a statement about where the loop can be.
2. **Dulac route (Ye's method).** Look for B = (1 + a x + b y)^k or a polynomial multiple with div(BX)
   sign-definite on the interior of Gamma minus a neighbourhood of the origin. Sign-definiteness proves the
   region is cycle-free (Dulac) and, evaluated near S, fixes the sign of sigma relative to the focus. The
   classical exponent that kills the linear terms of div(BX) is incompatible with the order-two condition
   except on the reversible stratum, so a two-parameter family B = (1+ax+by)^k (1 + c x)^j is the first
   candidate; test numerically on the 46 loops before attempting algebra.
3. **Liénard route (Zhang's method).** Reduce the order-two system to a generalized Liénard form and express
   sigma and eta_2 through the Liénard functions; Zhang's uniqueness proof for the order-two focus already
   controls the monotone quantities that decide loop stability.
4. **Melnikov route.** At sigma = 0 the loop's stability is the integral of div along Gamma; show this
   integral vanishes only on the center strata (F11 says it does), which is the statement that the
   integrable configurations are the only neutral ones.

## 5. What a proof would buy

Not a counterexample. It would explain, rigorously, why every attempt to stack a degenerate boundary on a
degenerate focus collapsed onto a center, and it would reduce the five-cycle question to the closed annuli
of the codimension-three centers, which is where the campaign's remaining open questions already sit.

## 6. Step-1 results (NUMERICAL, `audit/fable_d2_loop_div.py`, 23 loops with a < 0; the a > 0 loops are mirror images)

- The interior of every loop lies in {1 + a x + b y > 0}, touching the line only at the saddle (minimum of the
  line factor on the interior 0.001 to 0.019). The rescaling by 1/(1+ax+by) is therefore admissible inside.
- The divergence changes sign exactly twice on every loop (its zero set is a line through the focus); the loop
  integral of the divergence has the sign of sigma in all cases.
- Dulac candidates (1+ax+by)^k for k in {m/a, -(2l+b)/b, 1, -1, 2} are not sign-definite on any interior
  (positive fraction 0.17 to 0.92). An order-two focus forces div(BX) to vanish to fourth order at the origin,
  so no polynomial-times-power B of degree below four can work; Ye's route needs a quartic factor or a
  different method.

## 7. Continuity architecture (the route now being pursued)

On the loop locus L (a two-dimensional transcendental surface in (a, b, l)), sigma * eta_2 can change sign
only through sigma = 0 or eta_2 = 0. Zeros of eta_2 on L are center strata (the loop persists in the center
families). If zeros of sigma on L are also centers, i.e. "a homoclinic loop through a weak saddle around an
order-two weak focus forces a center", then sign(sigma eta_2) is constant on each component of L minus the
center strata, and D2 reduces to one evaluation per component plus a description of the components.
F11 observed exactly this at a = -3: the neutral loop is the crossing with the C3 center stratum.

## 8. Crossing search results (NUMERICAL, `audit/fable_d2_center_crossings.py`)

a = -2 and a = -2.5, 58 focus-type loops on b in [-4, 4]: sigma * eta_2 < 0 at every point (lemma evidence
now 104 loops in total). At a = -2.5 the focus branch has a neutral loop at b = 1.28125, l = -1.13719 with
sigma = -1.3e-9 and eta_2 = +2.0e-8; the three center factors there are a(b+2l) = 2.48, (b-3l-5) = -0.31 and
C3 = -2.3e-8. The neutral loop lies on the C3 center stratum to eight digits: lemma (A) holds at a = -2.5 as
it did at a = -3. At a = -2 sigma and eta_2 both change sign between b = 1.0 and 1.5 (bisection running).
