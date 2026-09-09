# Original Q4: finite perturbations toward infinity

No five-cycle candidate emerged. Three exact rational perturbations along the previously computed three-anchor direction each give three numerical compact sign brackets and no further sampled crossing up to section radius 10⁸.

Starting at (X,Y) = (0,−R), I integrated forward and backward to the oppositely oriented crossing of X = 0 with 0 < Y < 1. There were 91 sampled radii for each of δ = 1/1000, 1/10000 and −1/1000. All 582 half passages, including tighter endpoint checks, completed. The three common brackets are approximately

- R ∈ (0.278256, 0.359381);
- R ∈ (0.599484, 0.774264);
- R ∈ (1.291550, 1.668101).

Tighter endpoint checks preserved every selected sign, with changes below 1.30×10⁻¹⁴. SciPy clamped the requested relative tolerance 2×10⁻¹⁴ to its double-precision floor, approximately 2.22×10⁻¹⁴. These remain nonvalidated numerical calculations.

There is also an independent consistency check on the boundary separation. The previously computed boundary integral divided by the energy gradient at the positive section predicts the derivative −0.004268103997. At radius 10⁸, the measured separation divided by δ is −0.004263665891, −0.004267659896 and −0.004272547578 for the three perturbations respectively. Their behavior is consistent with the predicted first-order derivative.

This does not prove the infinite-radius limit, control the perturbation remainder, or exclude additional endpoint cycles. It supplies reproducible finite-field support for three compact crossings and the separation calculation. A simultaneous endpoint return expansion remains the key missing step.
