# `Q_3^R`: the first-order generating space, and where it goes blind

> **RETRACTED IN PART, 2026-09-06.** The central claim below — a
> four-dimensional generating space, an eight-dimensional kernel, and tiny
> determinants explained as ill-conditioning — is **false**. The four
> generators satisfy an exact linear relation, already present as equation (11)
> of `REVERSIBLE_RESEED_2026_09_05.md`; the rank is three, the kernel is nine
> and equals the span of the geometrically trivial directions, and the
> determinants are identically zero. Error found by Astra. See
> [Q3R_RANK_CORRECTION.md](Q3R_RANK_CORRECTION.md), which supersedes
> Findings 1 and 2 below. The exact first integral, integrating factor and
> flow-verified reduction in section "Exact structure" remain correct.

**Date: 2026-09-06.  Evidence class: exact algebra plus NUM.
No five-cycle candidate was found and none is claimed.**

This addresses the gap named in section 4a of
[CLAUDE_THOUGHT_SESSION.md](CLAUDE_THOUGHT_SESSION.md) — the reversible class
as "the largest unexplored integrable component" — and the observation there
that on reversible perturbations the first-order Melnikov function may vanish
identically, so higher order is where the cycles live.

## Setting

The two-centre reversible family of
[REVERSIBLE_RESEED_2026_09_05.md](REVERSIBLE_RESEED_2026_09_05.md), equation (1):

\[
 \dot x=\tfrac{b-2}{4}+(1-b)y+ax^2+by^2,\qquad \dot y=-2xy,
\]

with two centres exactly when `0<b<2`, `a<=0` (and `a` outside `{0,-1,-2}`).

## Exact structure

`mu = y^{a-1}` is an integrating factor — the unique one of the form `y^s`,
since `div(mu F) = -2x y^s (s+1-a)` — and

\[
 H(x,y)=y^{a}\Bigl(x^{2}+\frac{b-2}{4a}+\frac{1-b}{a+1}\,y+\frac{b}{a+2}\,y^{2}\Bigr),
\]

verified symbolically (`H_x = -mu Q_0`, `H_y = mu P_0`).  The ovals around the
upper centre `(0,1/2)` are `x^2 = R(y) = h y^{-a} - A - By - Cy^2`.

Under `x -> -x` the ovals are symmetric, so of the twelve quadratic
perturbation coefficients only `q00,q01,q02,q20` (through `dx`) and `p10,p11`
(through `dy`) reach `M_1`.  Integrating by parts,

\[
 M_1\in\operatorname{span}\{\,T_{a-2},\;T_{a-1},\;T_{a},\;U\,\},\qquad
 T_s=\!\int\!\sqrt{R}\,y^{s}dy,\quad U=\!\int\! R^{3/2}y^{a-2}dy,
\]

with the explicit coefficients

| generator | coefficient |
|---|---|
| `T_{a-2}` | `2(a-1) q00` |
| `T_{a-1}` | `2(a q01 + p10)` |
| `T_a` | `2((a+1) q02 + p11)` |
| `U` | `(2/3)(a-1) q20` |

**The generating space is four dimensional** — the same dimension as `Q_4`.

The reduction is not assumed: `verify_flow.py` integrates the unperturbed field
around an oval for one full period and accumulates the six loop moments
directly, agreeing with the reduced expressions to **1.3e-12 relative** at four
`(a,b,h)` points.

## Finding 1 (exact): an eight-dimensional blind spot

`M_1 == 0` exactly when

\[
 q_{00}=q_{20}=0,\qquad p_{10}=-a\,q_{01},\qquad p_{11}=-(a+1)\,q_{02}.
\]

Six coefficients never enter `M_1` at all, and four conditions cut the
remaining six down by four.  So **an eight-dimensional space of quadratic
perturbations, out of twelve, has `M_1` identically zero.**  Every limit cycle
born there is invisible to first order.

This is the precise size of the hole.  The repository's reversible work — the
64-sample moment search of the re-seed report — is a first-order search, so it
cannot see any of that eight-dimensional set, and the repository contains no
second-order Melnikov or Bautin-ideal machinery to search it with.

## Finding 2 (undecided): the Chebyshev property is not reachable by quadrature

Whether `M_1` can have four zeros on the annulus decides the route: three zeros
give `3+1 = 4` cycles, four would give a `4+1` five-cycle lead.  For a
four-dimensional space this is exactly the extended-Chebyshev question, tested
by the sign of `det[f_i(h_j)]` over increasing four-point subsets.

That test does not resolve here.  With rows and columns positively normalised
(which cannot change any sign), at 22 samples across the annulus:

| `(a,b)` | decided subsets | undecided | `|det|` range |
|---|---:|---:|---|
| `(-1/2, 1)` | 4354 | 2961 | `1.0e-55 .. 8.1e-42` |
| `(-1/4, 9/10)` | 2874 | 4441 | `1.0e-55 .. 4.3e-40` |
| `(-5/2, 6/5)` | 0 | 7315 | all below tolerance |

At **70 decimal digits** the determinants sit at or below the working
precision: the four generators are so close to linearly dependent on the
annulus that direct quadrature cannot separate them.  Double precision is far
worse — a first attempt there returned "up to 12 sign changes", which is pure
quadrature noise and is recorded here only as the negative result it is.

So this is not a bound of three and not a lead of four; it is a statement that
the question needs the Picard--Fuchs / Wronskian machinery the repository
already built for `Q_4`, applied to these four generators, rather than more
compute.

## What this does and does not settle

- The generating space, the first integral, the integrating factor and the
  coefficient table are exact and independently checked.
- The eight-dimensional `M_1 == 0` set is exact.
- The maximum number of zeros of `M_1` on the annulus is **open**.
- Nothing here computes a second-order Melnikov function; Finding 1 says why
  one is needed, not what it is.
- No perturbation with four annulus zeros was found, and no five-cycle field
  was produced.

## Replay

```bash
python3 q3r_2026_09_06/derive.py           # integrating factor and first integral
python3 q3r_2026_09_06/verify_flow.py      # reduction checked against the flow
python3 q3r_2026_09_06/kernel.py           # Finding 1, exact
python3 q3r_2026_09_06/cheb_hp.py          # Finding 2, 70-digit Chebyshev attempt
python3 q3r_2026_09_06/zeros.py            # the double-precision attempt (noise)
```

Dependencies: NumPy, SciPy, SymPy, mpmath.  No interval certificate is claimed.
