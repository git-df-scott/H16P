# KKL fold closure follow-up

Continues `7db8597cb7d9bb34e119e85bec3f229270eaf1aa`. The user resumed the
previously stopped strike and asked for a candidate or obstruction.

**Outcome: a proved obstruction to the proposed scalar quartic Dulac strategy,
and a new exact multiplier-band restriction. No three-origin-cycle candidate
or component-wide exclusion has been obtained.** These conclusions must not
be reported as a K1 kill or a solution of H(2).

## New exact results

For the original beta-zero KKL family, throughout
`1 <= c <= 8/5, K > 0`, the multiplier polynomial N is strictly positive on
`-1 <= x <= 0` and has at most two positive roots, counting multiplicities.
Its negative set on the right is therefore at most one interval. The restoring
polynomial W is positive everywhere on the origin side `x>-1`.

Consequently every attracting or multiplier-one origin cycle in this strip
must visit the same possible right-side negative band. The entire left side
contributes positively to the multiplier identity. The proof covers all
positive K in the strip; it is not an inference from the saved parameter grid.
Exact Bernstein coefficients and Descartes sign exclusions are reproduced by
[`theory_exact.py`](fold_closure_2026_09_05/theory_exact.py).

The proposed quartic certificate in the Lienard coordinates is

    Psi(z,T) = T^4 + C3 T^3 + A2(z) T^2 + A1(z) T + A0(z),
    X(Psi) + 4 kappa f Psi = Phi(z),   kappa=c/(2c+1), div X=-f.

At a multiplier-one periodic orbit, an integrating-factor identity forces
any one-sign continuous Phi to vanish throughout the orbit's z-projection.
If Phi is analytic, it must vanish identically on the connected interval.
Allowing isolated extra zeros cannot rescue the certificate.

There is a further obstruction to this degenerate escape. At a first-order
weak focus, a nonzero analytic solution of
`X(Psi)-4 kappa (div X)Psi=0` must have vanishing order `n=16 kappa`.
A monic quartic in T has vanishing order at most four. For c>1/2,
`16 kappa>4`, a contradiction. The statement uses analytic coefficients near
the origin and a connected domain containing both the focus and the cycle.

Thus **no certificate of this specified analytic monic-quartic scalar-residual
form can have one sign across a true fold in the stated regime**, even with
Phi identically zero. This is a complete obstruction to that strategy under
its hypotheses. It does not exclude a two-variable residual, off-fold
certificates on domains omitting the focus, nonanalytic/singular coefficients,
a different degree or multiplier exponent, or another orbit-comparison proof.

The complete proofs and hypotheses are in
[`theory_obstructions.md`](fold_closure_2026_09_05/theory_obstructions.md).
The number of limit cycles is still not bounded by the number of sign bands:
different cycles supply different positive weights in the multiplier integral.
The missing comparison between those weights remains a real proof obligation.

## Specific construction test: a possible triple-cycle degeneration

A regular finite fold creates at most two nearby return roots. A triple root
would supply a different local route to three origin cycles. To test that
possibility without repeating the parameter scans, one additional coefficient
was released:

    xdot = (1+x)y+x^2,
    ydot = -m x-10x^2+Bxy+c y^2.

The original family is B=11/5. The exact first-focus parameter is

    K=m(Bc-1)-10(B+2),
    l1=K/(8 m^(3/2))

in the normalized convention recorded in `generalized_exact.py`. K=1/512 was
held fixed, so the origin was not deliberately moved onto the center stratum.
For the reviewed two-sided matching map the augmented equations were
`F=0, G=0, G_z=0` in `(log r,c,B)`. A successful solution would still require
cubic nondegeneracy, unfolding rank, topology and complete-return checks.

Four damped Newton steps, using sixteen charged half-passage evaluations,
**did not converge**. The normalized residual norm increased from approximately
0.01146 to approximately 0.05922. The iterates are not accepted fold points,
not a new continued branch, and not triple-cycle candidates. Their motion
toward larger radii and smaller B does not prove an endpoint or an exclusion.
Every request, derivative estimate and failed correction is archived in
`fold_closure_2026_09_05/`.

An exact reversible-center parity polynomial for this enlarged family is

    H(B,c)=B^2 c^2+B^2 c-2Bc^2-Bc-B+40c^3-28c-10.

At c=1 it is `2(B-1)^2`. However at B=c=1 every finite m gives K=-30,
so this projected point cannot be accepted as a finite positive-K weak-focus
candidate. This is a coordinate/algebraic restriction, not a proof that the
failed Newton sequence approaches that point.

## Largest saved pair: previously unfinished replay completed

The last three available calls were complete clockwise angular returns at
one common field, with

    c=1.59340580527813710990835865677884849
    K=7.06390700436779910773804298664181037
    m=5(K+42)/(11c-5).

The following logarithmic displacement signs were obtained:

| Positive horizontal section radius | log(R/r) |
|---:|---:|
| 2.19194266668433e17 | +2.13829603794768e-4 |
| 2.95881311432548e17 | -5.12233717944016e-5 |
| 3.99397994234363e17 | +2.43937901467020e-4 |

All three complete returns succeeded at the identical saved coefficient
strings. Their angular variational/flux discrepancies are below 5e-23.
The signs reproduce two numerical root brackets at the largest saved pair
field. They do not show a third root or prove interval existence/uniqueness.
Full decimal data are in `outermost_check.json`. The earlier report's
unexecuted-outermost-replay limitation is superseded by this follow-up.

## Accounting and unresolved target

| Recorded work | Charged evaluations |
|---|---:|
| Historical KKL/Shi | 756 |
| Previous finite-fold strike | 3297 |
| Parallel merged reversible re-seed | 24 |
| This follow-up | 19 |
| **Known shared total** | **4096** |

The shared budget is exhausted. Exact algebra and saved-data checks above
add zero ODE evaluations. The cusp and outermost-return calls are separately
identified in the append-only follow-up ledger.

**Still unfinished:** the complete connected fold component, exhaustive
origin-root coverage, the weighted-orbit comparison, endpoint certificates,
a three-origin-cycle field, 3+1 coexistence, five-cycle completion, hostile
reproduction of a >=5-cycle candidate, and interval certification. No
mathematical stopping condition A-D from the construction request has been
established. The exact obstruction proved here closes the specified
certificate strategy; it does not close the counterexample route.
