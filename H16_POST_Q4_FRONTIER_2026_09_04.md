# H16 post-Q4 global frontier audit — 2026-09-04

Source baseline: [`5bcfe11`](https://github.com/git-df-scott/H16P/tree/5bcfe1172c124cb9a162a2e88c26aa5be26b40c1), verified against remote main during the audit. This independent branch, `audit/post-q4-frontier-2026-09-04`, contains only this audit's reports, updates and evidence; earlier campaign files are linked from the source history.

## Decision

**No accepted five-cycle smooth real planar quadratic example was found. The accepted lower bound remains H(2) >= 4; uniform finiteness of H(2) remains unknown.** Each fixed quadratic vector field has finitely many limit cycles, by [Bamón (1986)](https://www.numdam.org/item/PMIHES_1986__64__111_0.pdf); that is a different assertion from a uniform degree-two bound.

This status is supported independently of this repository by Gasull–Santana, *Proc. AMS* 153(2) (2025), 669–677, [DOI 10.1090/proc/17116](https://doi.org/10.1090/proc/17116), and especially Artés–Cairó–Llibre's paper published **2026-08-07**, whose introduction still states the quadratic maximum is unknown and four is the established example count: [QTDS 25:145](https://doi.org/10.1007/s12346-026-01563-4). The audit found no accepted later change by September 4.

**Recommended Strike #5: finite-distance continuation to a first-order Hopf precursor in the KKL family.** Find three nonzero hyperbolic cycles around an order-one weak focus, together with the remote singleton; then one ordinary Hopf cycle yields (4,1). The exact two-parameter problem and five-return certification gate are in [ASTRA_FIFTH_STRIKE_HANDOFF.md](ASTRA_FIFTH_STRIKE_HANDOFF.md). No precursor or fifth cycle has been found. This is a selected research problem, not a promised construction.

## What the audit changed

1. **Only (5,0) and (4,1) survive as five-cycle distributions.** Zegeling's 2024 Theorem 1.2 supplies the repaired (n,0)/(n,1) theorem; Theorem 5.4 further restricts four-real-equilibrium systems to (n,0) or (1,1). [Primary paper](https://doi.org/10.1515/anona-2024-0012). A finite saddle cannot simply be added to a robust (3,1) seed while preserving its four cycles.
2. **The infinitesimal Shi finite-saddle proposal is closed.** Exact focus and equilibrium algebra shows that its second-focus condition forces the remaining finite equilibria to be nonreal. Li's order-three theorem and Llibre–Schlomiuk's QW3 classification also exclude a preexisting cycle or graphic around the maximal weak focus. Some published neighborhood assertions in the latter paper are explicitly numerical; this audit does not promote them to theorems. [Classification](https://doi.org/10.4153/CJM-2004-015-2).
3. **Q4's five-interior-zero route stays CLOSED.** Theorem N gives <=3 distinct interior zeros in the strict lobe and <=4 globally for finite kappa>1. The [proof](https://github.com/git-df-scott/H16P/blob/5bcfe1172c124cb9a162a2e88c26aa5be26b40c1/Q4_THEOREM_N.md) and earlier Q4 files remain intact in the source history; the audit-only branch does not copy them. Its status is a repository mathematical result; this audit does not represent it as a peer-reviewed publication or rerun an independent full proof audit.
4. **The older Q4 endpoint closure was unjustified.** The elliptic-chart saddle's two branches map to two different saddles at infinity in the original quadratic field. The covering is singular; the one-saddle alien exclusion does not transfer. The old Lane B note also uses the wrong global-three bound and infers too much from sampled coefficient ranks. Endpoint status is **UNKNOWN**. [Original Q4 reduction](https://arxiv.org/abs/0811.4602), [exact audit](frontier_2026_09_04/Q4_GRAPHICS_AUDIT.md).
5. **M1 identically zero is not a new generic-Q4 basis.** The four-function generating space already represents the first nonzero Melnikov function. Higher-order multiplicity splitting, exceptional center intersections and endpoint limits require separate work. [Buica–Gine–Grau, Theorems 4–6](https://arxiv.org/abs/1406.7612).
6. **The original five-cycle claim is recovered.** Shi's December 1978 attempt combined three Bautin cycles with two Sommerfeld cycles. His 1990 retrospective supplies the family and inequalities, and explains why the corrected V3 sign leaves four. The source says V3, correcting the [earlier README's fifth-focus-quantity description](https://github.com/git-df-scott/H16P/blob/5bcfe1172c124cb9a162a2e88c26aa5be26b40c1/README.md#L80). [Primary retrospective, pp.202–203](https://www.labmath.uqam.ca/~annales/volumes/14-2/PDF/193-206.pdf).
7. **Recent six/seven claims must be logged, not accepted silently.** The 2026 Montisnigri six-cycle proof adds local counts at different values of lambda; its simultaneous construction fails an explicit Jacobian check. The August 2026 Hernandez Rosales seven-cycle preprint has no accessible manuscript or coefficients in the record located. It remains unverified, not refuted here. [Claims ledger](HISTORICAL_FIVE_CYCLE_CLAIMS.md).

## Four-cycle mechanisms and actual reproduction

The distinct construction classes located are Shi's order-three weak focus plus remote trapping cycle; Chen–Wang's order-two focus plus a preexisting intermediate and remote cycle; the Leonov/Kuznetsov small/large-cycle constructions; and reversible near-integrable two-center perturbations of Yu–Han. Galias–Tucker rigorously certify a fixed Shi/Songling instance; they do not supply an independent fifth mechanism. “Songling” and “Shi” name the same example. “Two large plus two small” describes scales, not a (2,2) configuration. No accepted four-around-one seed was located.

This audit independently reproduced all four returns of the published visual Shi, Chen–Wang, KKL and Yu–Zeng rational fields. For exact Galias–Tucker coefficients it reproduced the remote cycle and opposite return signs bracketing each of the three tiny cycles using 900-bit MPFR Taylor integration at orders 112 and 128. All six tiny-cycle signs agree to 35 printed significant digits across orders. These new calculations are **NUMERICAL ONLY**, with no validated remainder enclosure. The rigorous four-cycle theorem remains the published Galias–Tucker result. [Seed ledger](FOUR_CYCLE_SEED_LEDGER.md).

The exact Galias–Tucker point has no published surrounding parameter box in the inspected paper: its coefficient data are a singleton of exact rationals, with interval enclosures for cycle positions. We additionally tested the exact continuation path delta=-s, epsilon=-s^4, lambda=-10^8 s^16 from s=10^-13 toward 10^-2. Three sign-bracketed local cycles persist at the sampled intermediate values, and the remote return persists at the relaxed endpoint. This conditions a benchmark; it neither proves persistence over the whole path nor produces a fifth cycle.

## Modern developments that matter

- Yeung's criticism of an assertion in Ilyashenko's 1991 finiteness proof was published in January 2025. This proof discussion does not supply a quadratic counterexample or invalidate Bamón's earlier quadratic fixed-field result. The claimed gap, publication status and scope are recorded in the historical ledger. [Dulac's Theorem Revisited](https://doi.org/10.1007/s12346-025-01220-2).

- Marín–Villadelprat, *JDE* 433 (2025), 113281, obtain individual hemicycle cyclicity two and **simultaneous** cyclicity two or three off a resonance. Their a=-1 line has an omitted upper bound and alien examples. Independent compatibility algebra blocks the naive maximal-center-plus-degenerate-endpoint addition in a natural slice. [Primary manuscript](https://arxiv.org/abs/2501.16924).
- Their 2025 separation-function derivative work provides practical tools for actual infinity connection calculations. [QTDS 24:227](https://doi.org/10.1007/s12346-025-01379-8).
- Haibo Lu's active August 2026 H^3_14 manuscript claims existential local finite cyclicity of one semihyperbolic graphic. It is an unverified preprint and does not establish H(2) finite. [Version record](https://arxiv.org/abs/2607.13785).
- New quadratic phase-portrait classifications in 2025–26 organize restricted families; classification modulo cycles is not a global cycle-count theorem. Higher-degree replication, three-dimensional examples and piecewise quadratic examples do not improve H(2).

## Route ranking

Scores in [FIVE_CYCLE_MECHANISMS.md](FIVE_CYCLE_MECHANISMS.md) are judgments on research value, not probabilities or proof. The top three are:

1. **Finite-distance first-order Hopf completion**, using KKL as a tractable seed and the full return map as the target.
2. **Original-coordinate Q4 endpoint compatibility**, using its true two-saddle infinity boundary, with simultaneous interior and endpoint constraints.
3. **Resonant infinity hemicycles**, treating the 2025 a=-1 gap and center compatibility jointly.

Raw Galias–Tucker continuation ranks below the better-conditioned KKL representative for discovery, but above it as a published certificate benchmark. The old Shi finite-saddle route and the generic-Q4 larger-M2-space story are killed in their stated forms.

## Coverage and remaining obligations

| Requested audit area | Result / limitation |
|---|---|
| Current global status, 2025–26 changes, upper claims | Primary version/publication checks; accepted status unchanged; no universal proof inferred from a claimed theorem |
| Four-cycle mechanisms and coefficients | Four substantive construction classes, exact fields and reproducible evidence; no enumeration of every affine-equivalent field or every point in a parameter continuum |
| Historical five-plus claims | All located claims recorded, including surviving coefficients; unavailable Qin counterexample and Hernandez Rosales text explicitly retained as retrieval gaps |
| Shi/focus/codimension/topology | Exact focus quantities and equilibrium/infinity gates; same-stratum finite-saddle architecture excluded |
| Graphics, Q4 endpoint, higher order | Hypothesis-specific bounds and explicit unknowns; no invented maximum for every saddle-node graphic |
| Numerical seeds and continuation | Four visual fields, exact GT brackets, finite/infinite singularities, short probes, one GT homotopy; no complete global bifurcation atlas or nearest-surface minimization |
| Certification | Exact five-return gate, MPFR/interval precision plan, distinctness and itinerary rules; no new interval certificate claimed |
| Best strike / hostile check | Complete bounded specification, necessary inequalities, kill gates and explicit failure interpretation |

This is a substantial primary-source frontier audit, **not a guarantee that no historical paper was missed**. It did not complete an exhaustive MathSciNet/zbMATH review, every Chinese/Russian journal archive, or every inaccessible full text. Publication alone is not expert acceptance; lack of a retrieved rebuttal is not a proof. Exact unresolved retrieval and mathematical tasks are listed in the companion ledgers. No authors were contacted and no external messages were sent.

## Final status vocabulary

“NO” below means *not known / no accepted example located as of the cutoff*, not a theorem that five cycles are impossible or that H(2)=infinity. “LIVE” means a well-defined route not excluded by the audited results; it does not mean a candidate exists.

```text
ACCEPTED FIVE-CYCLE QUADRATIC EXAMPLE EXISTS: NO
H(2) KNOWN FINITE: NO
BEST ACCEPTED LOWER BOUND: 4
Q4 FIVE-INTERIOR-ZERO ROUTE: CLOSED
Q4 ENDPOINT ROUTE: UNKNOWN
SHI 4+1 ROUTE: NEEDS REFORMULATION
GALIAS-TUCKER CONTINUATION: LIVE
HIGHER-ORDER MELNIKOV ROUTE: UNKNOWN
BEST NEXT MECHANISM: Finite-distance first-order Hopf completion in the KKL family
BEST NEXT TARGET: Find three nonzero hyperbolic cycles around the beta=0, K>0 origin and one remote hyperbolic cycle in the exact KKL slice, then certify a fifth after beta becomes negative.
```
