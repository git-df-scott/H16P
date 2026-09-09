# Resonant reversible family: logarithmic first-order test

No five-cycle candidate was found. This is a numerical finite-sample result, not an upper bound on cyclicity.

The inherited moment scan explicitly omitted a = −1. For that case, with t = |y| and σ = sign(y), the base field has first integral

H = x²/t + bt + σ(1−b) log(t) + (2−b)/(4t).

Its derivative along the vector field vanishes identically, checked symbolically. This supplies a direct quadrature for the previously omitted logarithmic shapes. The perturbation is the same on both annuli: ε₀ = −τ/2, ε₁ = −cτ, ε₂ = −mτ. Its normalized first-order integral is σ(w−mv)−c, where v = J₋₁/J₋₂ and w = J₋₃/J₋₂ are weighted area ratios.

For b = 0.2, 0.5, 0.8, 1, 1.2, 1.5, 1.8, I computed 61 energy levels per annulus and swept the line arrangements of the sampled moment curves using the previously inspected repository routine. The largest sampled simultaneous crossing count was three, except at b = 1 where it was two. No sampled five-crossing direction appeared.

Gauss orders 256 and 512 differed by less than 2.62×10⁻¹⁵ in normalized moment values. At b = 1, the independent exact identities v = 1/2 and w = h+1 agreed within 3.2×10⁻¹⁵. These comparisons are sanity checks, not rigorous error bounds. The finite arrangement routine discards narrow numerical cells, and neither its grid nor its parameter samples cover all possibilities.

Marín and Villadelprat’s *The cyclicity of hyperbolic hemicycles*, arXiv:2501.16924v1 (2025 submission), Theorems B–C, leave the resonant upper bound open; Lemma 3.1 supplies the local five-parameter normalization. I checked those statements, not the entire proof. [Primary manuscript](https://arxiv.org/html/2501.16924v1), accessed September 7, 2026, Edmonton time.

The most useful remaining question is whether a joint near-boundary, higher-order unfolding can yield additional cycles while preserving compact ones. The present first-order sampled test does not answer it. Code, full profiles, convergence comparisons and the inherited arrangement source hash are preserved beside this report.
