# Astra handoff: Q4

## Decision

# YELLOW — Q4 is valid but remains computationally broad

Q4 can logically yield five quadratic limit cycles: five distinct simple zeros
of a realizable first nonzero Q4 generating function persist as five cycles for
all sufficiently small nonzero perturbations. No theorem bounds the relevant
integral below five. However, no accepted four- or five-zero Q4 example exists,
Żołądek's conjectural bound is three, and an exhaustive endpoint-aware cover is
still a continuous four-dimensional problem.

## Exact target

Search \(1<\kappa<85/23\) and \([\mu]\in\mathbb{RP}^3\), subject to Zhao's
beta-zero strip and \(P_2(\beta_0)\) inequality, for five distinct simple
zeros of

\[
I_\mu=\mu_1hI_{00}+\mu_2I_{10}+\mu_3I_{01}
+\mu_4(2I_{-1,0}+3\kappa hI_{-1,1})
\]

on \(h\in(-2/3,-2/(3\sqrt\kappa))\). Use the normalized root coordinate
\(r=(s-1)/(\kappa-1)\), with \(h=-2\sqrt{s/\kappa}/3\).

## First Astra command

Run exactly this first 24-CPU-hour-fused tranche:

    OPENBLAS_NUM_THREADS=1 python q4/q4_search.py \
      --mode astra --cpu-hours 24 --max-cpu-hours 24 \
      --candidate-mode triple --seed 160926 \
      --kappa-min 1.01 --kappa-max 3.69 --kappa-count 193 \
      --samples-per-kappa 50000 --grid-points 257 --quad-order 96 \
      --output q4/data/astra_tranche_001.json

The process stops earlier if it refines five sign-change roots. It must not be
expanded beyond this command without reviewing the retained branches.

## Promotion gates

1. NUM-LEAD: five roots persist under grid/order/precision doubling and endpoint
   changes.
2. Q4-ZERO: q4/zero_certificate.schema.json replays five interval-Newton
   inclusions with \(0\notin I'(S_i)\).
3. Q4-REALIZED: exact original-coordinate quadratic coefficients, all lower
   Melnikov functions zero, target \(M_k\) verified.
4. CAP-CANDIDATE: five validated return-map fixed points and disjoint flow
   tubes for one rational nonzero \(\varepsilon\).
5. COUNTEREXAMPLE: independent replay and human-readable proof.

Do not call a five-zero floating-point hit a counterexample.

## Stop conditions

Stop on the command's CPU fuse, a five-root lead, or a certified at-most-three
bound on the retained component. Freeze data and label every branch as
promoted, rejected, or unresolved. Do not search \(\kappa\ge85/23\), arbitrary
quadratic coefficients, or an endpoint without its analytic expansion.

## Expected cost

The old broad estimate was about 500 CPU-hours. Precomputed basis tables,
Zhao's exact pruning, and triple-zero coordinates reduce a serious numerical
campaign to approximately 20--60 CPU-hours. A well-conditioned Abelian-zero
certificate is estimated at another 2--20 CPU-hours. Original-field
realization and CAPD validation are deferred until a five-zero lead and may
cost 20--200 CPU-hours.

## Current evidence

- published Q4 lower construction: three zeros/cycles, asymptotic coefficients;
- finite numerical control reproduced here: three zeros at \(\kappa=4\);
- tiny structured strike: 2,060 directions, three analytic-filter survivors,
  maximum three sampled crossings, no four/five lead;
- audit verdict: YELLOW, unchanged.

The tiny strike is only an implementation validation and has no negative
mathematical force.
