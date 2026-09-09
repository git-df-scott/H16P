# Higher-order focus candidates and coexistence test

**Search-direction correction:** Published focus nonexistence/uniqueness results obstruct the proposed exact-focus precursor. The numerical record below is retained, but its suggested continuation is superseded by [the correction](H16P-focus-route-correction.md).

No counterexample found. Two numerical points on the cubic-degeneracy surface also have nearly vanishing fifth-order return coefficient and clearly nonzero seventh-order coefficient. They are plausible third-order weak-focus candidates, but neither has the additional sampled outer cycle needed for the proposed five-cycle construction.

The calculation uses the same exact-trace-zero quadratic family and rational parameterization as the preceding Bautin precursor search. In normalized local coordinates the angular radial equation is
\[
\frac{dr}{d\theta}=\frac{r^2R(\theta)}{1+rT(\theta)}.
\]
I integrated its formal power series through order seven. The cubic output agrees with the previously independently checked analytic cubic Hopf formula to 3.03e−15 on three nondegenerate controls.

Fifth-order coefficients were first evaluated on all 27 admitted near-cubic-degenerate fields. Five sign-changing parameter brackets were then solved numerically while remaining on a rational approximation to the cubic-degeneracy branch.

| Candidate | a | b | k | Seventh return coefficient | Sampled upper outer crossings | Sampled lower crossings |
|---|---:|---:|---:|---:|---:|---:|
| 1 | −1.1843004252363332 | 1/3 | 1/10 | +0.03207823 | 0 | 1 |
| 2 | −1.980999263255256 | 1 | 2/5 | −0.09690122 | 0 | 0 |

Three other numerical quintic zeros also had seventh coefficients near zero. Those may lie near centers and are not classified as third-order weak foci here. Vanishing several numerical coefficients is not a center proof.

All five points received complete two-sided passage profiles on 49 upper and 33 lower section radii: **820 half-passages**, all resolved. Candidate 1's lower crossing lies between radii 177.82794 and 316.22777. No upper outer crossing appeared from radius 0.03 to 10^7. This sampled absence is not a global exclusion.

For the two nonzero-seventh candidates, eight additional returns were integrated directly in the original Cartesian field. As the local radius decreased from 0.12 to 0.04, displacement/radius^7 moved toward the formal coefficients:

- Candidate 1: 0.03473059 → 0.03294471, predicted 0.03207823.
- Candidate 2: −0.11711926 → −0.10261800, predicted −0.09690122.

This independent finite-radius check supports the signs and magnitudes. It does not rigorously locate simultaneous coefficient zeros. The fields use rational approximations to numerical roots; interval enclosures, exact degeneracy, and generic unfolding conditions remain unverified. Three small cycles have not been constructed.

The coexistence gap is concrete: candidate 1 has evidence for the remote cycle, but lacks the additional outer cycle around the higher-order focus. Continuing that focus-degeneracy branch while testing an outer-cycle fold is a more targeted next step than treating either point as a counterexample.
