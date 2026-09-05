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

## 9. Proposition A (EXACT, verified by polynomial algebra in `audit/fable_d2_theoremA.py`)

**Proposition A.** In the Shi chart at zero trace with eta_1 = 0, let E be any equilibrium other than the
origin. Then div X(E) = 0 if and only if the origin is a center; precisely, iff a(b+2l) = 0 (E = (0,1), where
div = m) or C3 = 0 (E on the line 1+ax+by = 0), the chart degeneracies b = 0 and l = -1 aside.

Proof. The non-origin equilibria are (0,1) and the points of the line 1+ax+by = 0 where P = 0. At (0,1)
the divergence equals m = a(b+2l)/(l+1). On the line, the divergence is the affine function
-(b+2l)(a^2 x + a - b(l+1) x)/(b(l+1)) of x, which vanishes at the single point x0 = -a/(a^2 - b(l+1)) when
b+2l != 0 (and identically when b+2l = 0). An equilibrium sits at x0 iff P(x0, -(1+a x0)/b) = 0, and the
numerator of that expression factors exactly as -b^2 (l+1) C3 (sympy: remainder modulo C3 is zero). Both
a(b+2l) and C3 are factors of eta_2, so each condition places the origin on a center stratum. QED.

**Corollary A1 (no neutral loops).** On the order-two stratum, a homoclinic loop around the origin through a
saddle of zero divergence exists only when the origin is a center. This is the exact form of what F11 saw.

**Corollary A2 (structure of the sign law).** Along the focus sheet L of the loop locus, sigma vanishes only
where C3 = 0 or b+2l = 0, and there eta_2 vanishes too, each to first order generically (sigma is
proportional to x_S - x0, which is a simple zero of C3 when the saddle is a simple root of P on the line).
Hence sigma * eta_2 is proportional to C3^2 or (b+2l)^2 near those crossings and does not change sign. The
only other zero of eta_2 is the Shi line b = 3l+5, where sigma has no factor; so Conjecture D2 reduces to:

  (C) the focus sheet L meets the Shi line b = 3l+5 only at points where (b+2l) also vanishes (the chart
      point (b,l) = (2,-1)), i.e. there is no homoclinic loop around an order-three weak focus;
  (D) one sign evaluation per connected component of L minus these curves.

(C) is exactly Lane C's numerical finding on the order-three stratum (loops only at centers, where the
center subset of the Shi line is C3 intersected with the Shi line plus 2a^2 + l + 2 = 0 and a = 0; the
quintic factor of eta_3 on the Shi line is C3 restricted to it, checked by hand: a^2(5l+6) = 3(l+1)^2(l+2)).
(D) has been evaluated at 104 loops on the components sampled: negative in every case.

Status: Proposition A is proved. Conjecture D2 is reduced to (C), an order-three loop statement one rung
above Li-Cherkas, plus a finite component check (D). No counterexample content; this is the rigidity result.

## 10. Lemma (B) at the crossings (NUMERICAL)

Bisection for the neutral loop loses the focus branch at b = 1.2656 for both a = -2 and a = -1.5: at the C3
crossing the focus and center sheets of the loop locus meet, the splitting function has a double zero in l,
and sign-based detection fails there. Interpolating the brackets instead: at a = -2, sigma vanishes at
b = 1.276 and eta_2 at b = 1.273 (slopes +0.58 and -4.3 per unit b); at a = -1.5, at b = 1.281 and 1.277
(slopes +0.75 and -1.9). Both vanish at the same point with opposite-signed slopes, as Proposition A
predicts (both are simple zeros of C3 along the sheet). Combined with a = -2.5 (neutral loop on C3 to eight
digits) and a = -3 (F11), lemma (B) is confirmed at every crossing found.

## 11. Attempt at (C): loops around an order-three weak focus (NUMERICAL + one exact identity)

Parity plus Li-Cherkas: a homoclinic loop around a genuine order-three focus encloses no limit cycle, so the
loop's stability must oppose the focus's, sigma * eta_3 < 0. Sign map over 48,254 saddles on the Shi line
(b = 3l+5, m = 5a, a in [-4,4], l in [-4,3]): 22 percent are loop-compatible (7,544 at the saddle (0,1),
3,170 on the line 1+ax+by = 0). So parity alone does not exclude loops. Exact identity at (0,1):
sigma = 5a and eta_3 = -25 a (2a^2+l+2) q / 64 with q the quintic factor (q = C3 restricted to the Shi line),
so sigma * eta_3 = -125 a^2 (2a^2+l+2) q / 64 and (0,1) is loop-compatible exactly where q < 0 (the factor
2a^2+l+2 is positive on the whole sampled region). For saddles on the line no such factor law holds (43
percent agreement with the q-sign rule): the sign depends on which root of P on the line is the saddle.
Conclusion: (C) does not follow from sign algebra. It stands as Lane C's numerical statement (loops on the
order-three stratum occur only on the center curve q = 0) unless the Llibre-Schlomiuk classification of
third-order weak-focus portraits already contains it (being checked).

## 12. Statement (C) is a theorem in the literature

Llibre and Schlomiuk, "The geometry of quadratic differential systems with a weak focus of third order",
Canad. J. Math. 56 (2004) 310-343, Theorem 16: among the eighteen phase portraits of the class QW3, those
with a graphic are exactly W13, W15 and W18, and in each "the graphic is unique, surrounding a strong focus";
all graphics have their singular points at infinity, and the bifurcation set (Theorem 12) contains no finite
saddle-loop bifurcation. Limit cycles around the third-order weak focus are excluded by citing Li (1986).
Hence no homoclinic loop surrounds a third-order weak focus in any quadratic system: the focus sheet L of the
order-two loop locus never meets the Shi line b = 3l+5 at a non-center point.

## 13. Status of Conjecture D2

Proved: Proposition A (section 9, exact) and statement (C) (Llibre-Schlomiuk 2004). Consequence: on every
connected component of L minus the center strata, sigma and eta_2 are continuous and nonzero, so
sign(sigma * eta_2) is constant per component. Remaining: (D), a description of the components of L minus the
center strata and one sign evaluation on each. Evaluated: 104 loops over a in [-3, 3], b in [-4, 4], every
component met, all negative. D2 is therefore a theorem modulo the finiteness and enumeration of those
components, which is a computable but not yet completed step. The Zhang uniqueness theorem (1999) adds that
a hypothetical component with positive sign would carry exactly one hyperbolic cycle inside the loop.
