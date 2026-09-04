# Q4 cost model

## Why the old estimate was about 500 CPU-hours

The feasibility audit budgeted roughly 500 CPU-hours for repeated
128--512-bit elliptic-integral evaluation over a generic four-dimensional
\((\kappa,[\mu])\) cover, followed by adaptive root isolation. That estimate
also included eight broad \(\kappa\)-charts, most of which are now proved
irrelevant to a five-zero result.

## Measured kernel

On the audit host, nine tables of 121 nodes plus 2,060 structured projective
directions took 0.50 CPU seconds. The expensive part is the basis table; dot
products for many \(\mu\) values at fixed \(\kappa\) are cheap. A 60-digit
single-point adaptive quadrature took about 0.33 seconds; the binary64 Gauss
version took about 0.006 seconds and an independent orbit evaluation about
0.024 seconds.

These are engineering measurements, not portable benchmarks.

## Pruning factors

| Reduction | Effect |
|---|---|
| Zhao strip | removes all \(\kappa\ge85/23\) and all \(\beta_1=0\) cases |
| \(P_2(\beta_0)\) bound | rejects another analytic region |
| basis precomputation | changes repeated integral calls into matrix products |
| triple-zero coordinates | samples only functions already capable of rich structure |
| continuation | reuses roots and periods between nearby parameters |
| interval boxes | rejects sign/variation boxes in parallel |
| endpoint charts | avoids arbitrary global precision escalation |

In the random-direction smoke control, 64 of 2,304 directions (2.8%) survived
both Zhao tests. In the triple-zero smoke test, 3 of 2,060 survived the strip.
These percentages are sampling diagnostics, not rigorous volume estimates.

## Revised estimates

| Campaign | CPU estimate | Meaning |
|---|---:|---|
| old direct grid | about 500 h | broad high-precision four-dimensional scan |
| pruned numerical discovery | 20--60 h | tables, triple-zero continuation, adaptive promotion |
| first Astra tranche | at most 24 h | bounded candidate mapping and near-miss census |
| interval proof per good zero lead | 2--20 h | five root boxes, depending on conditioning |
| explicit-field/CAPD backend | 20--200 h | only after realization; saddle proximity dominates |

Surrogate ranking and parallel matrix products could lower wall time, but do
not change logical coverage. A credible first attack is now below 500
CPU-hours; an exhaustive exclusion of all \(\mathcal D_\delta\), including
endpoint limits, is not. The continuous four-dimensional manifold, thin
multiple-root loci, conjectural bound three, and missing coefficient-realization
implementation keep the verdict YELLOW.

No more than 0.0002 CPU-hours was spent on the exploratory strike in this
audit, far below the five-hour ceiling.
