# Q4 candidate search

## Design

The screen in q4/q4_search.py never integrates a perturbed quadratic field. It
evaluates the exact four-function target (Q4-I) numerically and ranks coefficient
directions. Its stages are:

1. choose \(1<\kappa<85/23\) and build a sine-clustered \(s\)-grid;
2. precompute the four basis functions by fixed Gauss quadrature;
3. choose three separated grid locations and obtain the projective
   \(\mu\)-direction in the nullspace of their \(3\times4\) basis matrix;
4. apply Zhao's exact \(\beta_0\)-strip and \(P_2(\beta_0)\) filters;
5. count robust sign changes, refine brackets at higher quadrature order, and
   record root slopes and conditioning.

Arbitrary precision is available through q4_integrals.basis_mp and is the
promotion path for any leader. The code canonicalizes projective directions,
so \(\mu\) and \(-\mu\) are not duplicated.

## Adaptive promotion pipeline

A retained leader must pass, in order:

- resample every bracket at doubled Gauss order;
- evaluate the same points at 80, 160, and 320 decimal digits;
- recursively subdivide wherever an interval contains a sign change or where
  \(|I|/\|\mu\cdot B\|\) is below the conditioning threshold;
- isolate roots with Brent only for ranking, then interval Newton for proof;
- estimate \(I'\), reject a slope that trends to zero under precision
  doubling, and report distance from both endpoints;
- perturb each of \(\kappa\) and three local projective coordinates, tracking
  roots by predictor/corrector continuation;
- hash the canonical tuple to deduplicate continuation branches.

Surrogates may interpolate the four basis functions and rank boxes. They may
not reject a box unless an outward error enclosure makes the rejection
rigorous.

## Tiny strike performed

The only exploratory run in this audit was:

    python q4/q4_search.py --mode smoke --cpu-hours 0.02 \
      --candidate-mode triple --kappa-count 9 --samples-per-kappa 256 \
      --grid-points 121 --quad-order 64 --output q4/data/smoke.json

It consumed 0.50 CPU seconds (0.40 wall seconds), evaluated 2,060 projective
directions, and retained three after the necessary \(\beta\)-strip. Each
retained direction showed three sampled crossings; none showed four or five.
Two leaders refined to three roots, with two roots close to the homoclinic
endpoint. The third lost its brackets when quadrature order changed, exposing
the expected ill-conditioning failure mode.

This sample is far too small and too deliberately structured to estimate a
frequency of five-zero functions. It validates throughput, pruning,
deduplication, and promotion failure handling. It is a lead generator, not an
exclusion and not a search result about the conjecture.

## Continuation strategy

The first real campaign should continue the manifold of triple-zero
directions in coordinates \((\kappa,r_1,r_2,r_3)\), watching for either:

- two additional sign-change pairs;
- a new double root \(I=I_s=0\), across which root count may change by two;
- a root entering from a controlled endpoint expansion.

This is mathematically sharper than random sampling in \(\mathbb{RP}^3\).
Near rank-deficient basis matrices, use multiprecision SVD and do not infer a
root count from binary64 values.

## Stop rules

Stop a tranche on its CPU fuse, any five-root lead, or proof that the retained
component has at most three roots. Freeze all leaders and failure labels.
Never roll unused time into \(\kappa\ge85/23\), generic quadratic coefficients,
or an endpoint with no asymptotic chart.
