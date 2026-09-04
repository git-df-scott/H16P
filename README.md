# H16 post-Q4 frontier audit — September 4, 2026

This branch contains only the work produced by this frontier audit: six reports, two current-state updates, and the supporting source analysis, calculations, data and replay instructions. Its history starts with this standalone audit package.

Start with the [full frontier report](H16_POST_Q4_FRONTIER_2026_09_04.md) and [Astra Strike #5 handoff](ASTRA_FIFTH_STRIKE_HANDOFF.md).

| Report | Contents |
|---|---|
| [Global frontier audit](H16_POST_Q4_FRONTIER_2026_09_04.md) | Accepted status, findings, route ranking and coverage limits |
| [Five-cycle mechanisms](FIVE_CYCLE_MECHANISMS.md) | Configuration restrictions, scores and remaining routes |
| [Four-cycle seed ledger](FOUR_CYCLE_SEED_LEDGER.md) | Exact families, cycle reproductions and continuation probes |
| [Historical claims](HISTORICAL_FIVE_CYCLE_CLAIMS.md) | Five-plus and upper-bound claims, corrections and source gaps |
| [Certification plan](H16_CERTIFICATION_PLAN.md) | Required validated return-map checks for five cycles |
| [Strike #5 handoff](ASTRA_FIFTH_STRIKE_HANDOFF.md) | Exact KKL family, bounded experiment and certification gate |

Current updates: [attack matrix](ATTACK_MATRIX.md) and [canonical state](CANONICAL_STATE.md).

Supporting material: [source audits, independent checks, data and replay](frontier_2026_09_04/README.md).

The audit found no accepted five-cycle quadratic example. The accepted lower bound remains four, and uniform finiteness remains unknown. New computations here are numerical evidence, not interval proofs. The recommended KKL precursor remains unproved.

Earlier campaigns, Q4 code and Theorem N remain in the [source repository history at 5bcfe11](https://github.com/git-df-scott/H16P/tree/5bcfe1172c124cb9a162a2e88c26aa5be26b40c1). References to that work are links; those files are not included in this branch.
