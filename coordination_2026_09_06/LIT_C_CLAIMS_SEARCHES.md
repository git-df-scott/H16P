# LIT_C — Every claim, search, and certification effort for ≥4 limit cycles in planar quadratic systems, 1979–2026

Compiled 2026-09-06. Scope: slice C of the H(2) literature sweep — four-cycle constructions, all ≥5 claims and their refutations, computer-assisted certification, explicit restricted upper bounds, and Melnikov/Abelian-integral cyclicity of the four quadratic center families.

PDFs are in `/Users/scottg/Claude_all/papers/`. Text extraction used `/Users/scottg/Claude_all/h16p_env/.venv/bin/python` + pypdf 6.17.0 (helper `papers/x.py`).

**Notation.** A *nest* is the set of limit cycles surrounding one focus. `(p,q)` = p cycles around focus 1, q around focus 2. A quadratic system's limit cycles surround exactly one singular point, which is a focus, and a quadratic system has at most two foci — so at most two nests.

---

## 0. Executive summary of the headline question

- **Nobody has ever exhibited four limit cycles in one nest.** Every single one of the ~dozen independent four-cycle constructions from 1979 to 2026 is **(3,1)** (or its mirror (1,3)). Not one is (4,0) or (4,1).
- The distribution is *not* what blocks it: Zhang Pingguang (2002), verified and gap-filled by Zegeling (2024), proves the only possible distributions are **(n,0)** and **(n,1)** — so (4,1) and (4,0) are permitted by the distribution theorem. The obstruction is purely that **no construction has ever produced a fourth cycle in a single nest**.
- Bautin (1952) caps *small-amplitude* cycles at one focus at 3. So a 4-in-one-nest example needs at least one normal-size cycle in that nest, produced by a global mechanism (separatrix/graphic bifurcation, semistable-cycle bifurcation, or Poincaré bifurcation from a period annulus).
- H(2) ≥ 4 is the state of the art. **H(2) < ∞ is not known** (Gasull–Santana, arXiv:2407.13465, Oct 2024, verbatim: *"Even for the quadratic case, it is not known if H(2) < ∞."*).
- Every claim of ≥5 has failed or is unverified; every claim of "H(2)=4 proved" (Gaiko, Pedregal, da Silva et al.) is unaccepted, and two carry published/reviewed objections.

---

# PART 1 — The four-cycle constructions (all of them)

## 1.1 Chen Lansun & Wang Mingshu 1979 — the first example

**Citation.** Chen Lansun, Wang Mingshu, *The relative position and number of limit cycles of quadratic differential equations*, Acta Math. Sinica **22** (1979) 751–758. **Paywalled / Chinese-language original not obtained.**

**Explicit system** (verbatim, via Yu & Zeng arXiv:2002.09987 eq. (12), `papers/yu_zeng_2021_visualization_four_limit_cycles_near_integrable_2002.09987.pdf`):

```
ẋ = −δ₂ x − y − 3 x² + (1 − δ₁) x y + y²
ẏ =  x + (2/9) x² − 3 x y
```
with `0 < δ₁, δ₂ ≪ 1`, **values not specified in the original**.

Verbatim from Yu–Zeng p.4: *"In [3] Chen and Wang constructed two trapping regions (see Figure 1(b)) and used Poincaré-Bendixson theory to prove the existence of the big limit cycle around (0,1) and a small one around (0,0). Then they showed that further perturbing the 2nd-order fine focus (0,0) using the parameters δ₁ and δ₂ to obtain two more small limit cycles, but did not specify the values of δ₁ and δ₂."*

**Distribution: (3,1).** Big **stable** cycle around unstable focus (0,1); three small cycles around stable focus (0,0), outer unstable / middle stable / inner unstable. Origin is a **2nd-order fine focus** when δ₁=δ₂=0.

**Sizes.** The three small cycles are infinitesimal; the big one is normal-size. Yu–Zeng supply working values `δ₁ = 0.01, δ₂ = 0.00002`, giving `v₀ = −10⁻⁵, v₁ = 0.0025, v₂ = −172948799/2332800000 ≈ −0.074138`, hence amplitudes `r₁ = 0.068102, r₂ = 0.170538` — the first published numerically visualizable version of Chen–Wang.

## 1.2 Shi Songling 1980 — the Songling system

**Citation.** Shi Songling, *A concrete example of the existence of four limit cycles for plane quadratic systems*, Sci. Sinica **23** (1980) 153–158. (Chinese version: Sci. Sinica **11** (1979) 1051–1056.) **Paywalled.** Full statement recovered verbatim from Galias–Tucker 2022 and Perko 1984, both held.

**Explicit system** (verbatim from Galias–Tucker, `papers/galias_tucker_2022_songling_exactly_four_limit_cycles_AMC415.pdf` eq. (1)):

```
ẋ = λx − y − 10x² + (5 + δ)xy + y²
ẏ = x + x² + (−25 + 8ε − 9δ)xy
```

**Coefficient vector — note the literature disagrees on λ:**

| source | δ | ε | λ |
|---|---|---|---|
| Galias–Tucker 2022 (AMC 415, eq. 1) | −10⁻¹³ | −10⁻⁵² | **−10⁻²⁰⁰** |
| Perko 1984 (RMJM 14, eq. 5) | −10⁻¹³ | −10⁻⁵² | −10⁻²⁰⁰ (OCR "IO-20»", i.e. 10⁻²⁰⁰) |
| Yu–Zeng 2020 (arXiv 2002.09987, eq. 11) | −10⁻¹³ | −10⁻⁵² | **−10⁻²⁵⁰** |

Galias–Tucker's is the value actually used in the rigorous proof; treat **λ = −10⁻²⁰⁰** as canonical.

**Distribution: (3,1).** Two equilibria (0,0) and (0,1). Verbatim Yu–Zeng: *"both systems have a big stable limit cycle around the unstable focus (0,1), and three small limit cycles around the stable focus (0,0)."* Origin is a **3rd-order fine focus** when λ=ε=δ=0; Shi perturbs it with the three parameters λ, ε, δ.

**Method.** Verbatim Yu–Zeng p.4: *"In [8] Shi explicitly constructed four trapping regions (see Figure 1(a)) and applied Poincaré-Bendixson theory to prove the existence of four limit cycles, one of them around (0,1) and three of them around (0,0)."*

**Scales.** Perko 1984 p.635 (verbatim): *"the size of L₂ and L₁ cannot be determined by this method. All that is known at this time is that 10⁻⁶¹ < |y₂| < 10⁻¹⁹ and that 0 < |y₁| < 10⁻⁶⁰; this was established by Songling in [10], pp. 155-156."* Galias–Tucker later pinned all four exactly (§3.1 below).

## 1.3 Perko 1984 — the first NORMAL-SIZE four cycles

**Citation.** L. M. Perko, *Limit cycles of quadratic systems in the plane*, Rocky Mountain J. Math. **14**(3) (1984) 619–645. **PDF: `papers/perko1984_RMJM_limit_cycles_quadratic.pdf`.**

Perko lists six known configurations (a)–(f); **(f) is the (3,1) configuration**, and is the only one with four cycles. Verbatim p.620: *"And Shi Songling [10] gave an example with at least four limit cycles in the configuration (f), thereby disproving the assertion of Petrovskii and Landis [9] that a quadratic system can have at most three limit cycles. These examples represent all of the known limit cycle configurations for quadratic systems in the plane."*

**§5, "A Modification of Songling's Example", verbatim:** *"In order to obtain a quadratic system with four limit cycles of 'normal size' in the configuration (f), Songling's example was modified as follows. The system (5) was first embedded in a uniformly rotated vector field ... ẋ = P cos α − Q sin α, ẏ = P sin α + Q cos α (6)."*

**Explicit normal-size coefficient vector (verbatim p.639):** *"Figures 20 and 21 then show the four limit cycles of the system (6) with δ = −.5, ε = −.01, λ = −.005, and α = −.0023 in the configuration (f) of Figure 1."*

Poincaré-section intersections with the negative y-axis, four-figure accuracy: **y₁ = −.0425, y₂ = −.2160, y₃ = −1.3838** (plus the fourth cycle L₄ around the other focus).

**Rigour caveat (important).** Perko's §5 is a *numerical* study resting on an assumption. Verbatim p.637: *"The dashed curve indicates the qualitative behavior of the limit cycles L₁(α), L₂(α), and L₃(α) around the origin of (6) which follows from Songling's results [10] and the theory of rotated vector fields [5,7] **under the assumption that there are at most three limit cycles around the origin**."* Perko also records p.631: *"It is still an open problem to determine whether or not this system has exactly four limit cycles in the configuration (f)."* — settled only in 2022 by Galias–Tucker.

**Chin (Qin Yuanxun) et al.'s system**, recorded by Perko as (7), for λ = 10⁻⁸²⁰, δ = −10⁻³⁶⁶ (OCR uncertain), ε = −10⁻⁷⁸:
```
ẋ = λx − y + (2 + δ)xy − y²
ẏ = x + λy + x² − (5 + ε)xy − ...
```
This is a **configuration (c) = (3,0)** example (three cycles, one nest, all small-amplitude), not four. Perko: *"the smallest limit cycle in Chin's example (7) is very nearly a circle with a diameter of O(10⁻⁸²⁰)."*

## 1.4 Cherkas–Artés–Llibre 2003 — exact counts, normal size, via Dulac functions

**Citation.** L. A. Cherkas, J. C. Artés, J. Llibre, *Quadratic systems with limit cycles of normal size*, Bul. Acad. Ştiinţe Repub. Mold. Mat. **1(41)** (2003) 31–46. **PDF: `papers/cherkas2003.pdf` (= `papers/cherkas_artes_llibre_2003_normal_size_BASM41.pdf`). OPEN ACCESS.**

Verbatim abstract: *"In the class of planar autonomous quadratic polynomial differential systems we provide 6 different phase portraits having exactly 3 limit cycles surrounding a focus, 5 of them have a unique focus. We also provide 2 different phase portraits having exactly 3 limit cycles surrounding one focus and 1 limit cycle surrounding another focus. The existence of the exact given number of limit cycles is proved using the Dulac function. All limit cycles of the given systems can be detected through numerical methods; i.e. the limit cycles have 'a normal size' using Perko's terminology."*

**Normal form (2), verbatim:**
```
dx/dt = 1 + xy
dy/dt = a₀₀ + a₁₀x + a₂₀x² + a₀₁y + a₁₁xy + a y²,      a₀₀ = a₀₁ + a₁₁ − a₁₀ − a₂₀ − a
```

**Table 1 — the complete explicit coefficient vectors (verbatim):**

| No | a | a₂₀ | a₁₁ | a₀₁ | a₁₀ | Singular points | Cycle distr. |
|---|---|---|---|---|---|---|---|
| 1 | 3 | −12 | −1.398 | 8.4 | 15.28 | 1F + 1N + 2S∞ + 1N∞ | 3 |
| 2 | 1.5 | −15 | 0.79993 | 3.2 | 9.17 | 2F + 2S∞ + 1N∞ | (3,0) |
| 3 | −2 | 12 | 10.999 | −14 | −26.1 | 1F + 3S + 3N∞ | 3 |
| 4 | −2 | −1 | 9.49965 | −12.5 | 6.955 | 1F + 1S + 2N∞ + 1S∞ | 3 |
| 5 | −4 | −1 | 13.9987 | −21 | 12.4 | 1F + 1N + 2S + 2N∞ + 1S∞ | 3 |
| 6 | 5 | −50 | −5.49995 | 16.5 | 76.45 | 1F + 2N + 1S + 1N∞ + 2S∞ | 3 |
| **7** | **8/11** | **−12** | **2.1502** | **67/220** | **−26.5** | **2F + 1S∞** | **(3,1)** |
| **8** | **1.04** | **−120** | **1.51997** | **1.56** | **−79.6** | **2F + 2S∞ + 1A∞** | **(3,1)** |

(F=focus, N=node, S=saddle, subscript ∞ = at infinity, A∞ = antisaddle at infinity.)

**Rows 7 and 8 are the two four-cycle systems, both (3,1).** All six others are three-cycle. **No row is (4,0).**

**Critical context, verbatim from the paper's introduction (p.31):** *"Recently, Zhang Pingguang [20,21] has proved that if nᵢ > 0 for i = 1,2, then either n₁ = 1, or n₂ = 1. The following distributions of limit cycles for quadratic systems (1) are known: (a) 1 and (1,0); (b) 2 and (2,0); (c) 3 and (3,0); (d) (1,1); (e) (2,1); (f) (3,1)."* — the known-distribution list tops out at (3,1) and (3,0).

**Method:** Perko's rotated-parameter method systematized, plus Dulac functions Ψ(x,y,C) built as an optimization over a uniform net (their eq. (11)), plus a "method of reduction to the global uniqueness of a limit cycle" using the Andronov–Hopf function AH(x) = a₁₁. The exactness (upper bound) is the genuinely new content.

Foundation: Artés & Llibre, *Quadratic vector fields with a weak focus of third order*, Publ. Mat. **41** (1997) 7–39 (`papers/pubmat41.pdf`) — 20 global phase portraits from 16 local ones, of which only three have any limit cycle. Verbatim from its intro: *"A quadratic vector field has at most 4 limit cycles, and when..."* (the standing conjecture).

## 1.5 Leonov / Kuznetsov–Kuznetsova–Leonov — the "hidden oscillations" school, 2010–2013

All **paywalled**; abstracts recovered verbatim via OpenAlex.

**(a) G. A. Leonov, *A criterion for the existence of four limit cycles in quadratic systems*, J. Appl. Math. Mech. (PMM) **74**(4) (2010) 421–435, doi:10.1016/j.jappmathmech.2010.05.002. Paywalled, no abstract in OpenAlex.**
Reported content: uses asymptotic integration of trajectories of the Liénard equation to prove two *large* limit cycles in a quadratic system with a weak focus, then standard parameter perturbation for two *small* ones. *"The criterion obtained for the existence of four limit cycles generalizes the well known Shi theorem."* The criterion covers both **3 small + 1 large** and **2 small + 2 large**. Note carefully: "2 small + 2 large" is still **(3,1)** overall (2 small + 1 large in the first nest, 1 large in the second) — it is not a (2,2) distribution, which is provably impossible.

**(b) G. A. Leonov, *Four limit cycles in quadratic two-dimensional systems with a perturbed first-order weak focus*, Doklady Mathematics **81**(2) (2010) 287–289, doi:10.1134/S1064562410020237. Paywalled.**

**(c) G. A. Leonov, *Four normal size limit cycles in two-dimensional quadratic systems*, Int. J. Bifurcation and Chaos **21**(2) (2011) 425–429, doi:10.1142/S0218127411028532. Paywalled. Abstract verbatim (OpenAlex):**
> *"The existence criterion of three normal size limit cycles in quadratic systems with a weak focus of first order is obtained. Further, giving a finite disturbance for weak focus, the fourth normal size limit cycle is obtained. Bifurcation of appearance of two limit cycles via semistable cycle is given."*

This is the strongest existing statement: **three normal-size cycles around a first-order weak focus** — i.e. a nest of 3 that is *not* small-amplitude — plus a fourth from a finite disturbance. Still (3,1).

**(d) N. V. Kuznetsov, O. A. Kuznetsova, G. A. Leonov, *Visualization of Four Normal Size Limit Cycles in Two-Dimensional Polynomial Quadratic System*, Differ. Equ. Dyn. Syst. **21**(1–2) (2013) 29–34, doi:10.1007/s12591-012-0118-6. Paywalled (Springer 303-redirects to auth); no abstract in OpenAlex.**

**(e) Kuznetsov / Leonov, *Visualization of four limit cycles of two-dimensional quadratic systems in the parameter space*, Differential Equations **49**(13) (2013), doi:10.1134/S0012266113130028. Paywalled.**

**Explicit system and coefficients — recovered from Yu & Zeng arXiv:2002.09987 §3 (open access), which reproduces Kuznetsov et al.'s example verbatim as their system (8):**

```
ẋ = y + x² + x y
ẏ = a₂ x² + b₂ x y + c₂ y² + α₂ x + β₂ y
```

Existence conditions (Yu–Zeng eq. (16), quoting [7] = Kuznetsov et al.):
```
b₂ ∈ (1,3),  c₂ ∈ (1/3, 1),  4a₂(c₂ − 1) > (b₂ − 1)²,  b₂c₂ > 1,
α₂ ∈ ( a₂(b₂+2)/(b₂c₂−1),  a₂(b₂+2)/(b₂c₂−1) + δ ),   β₂ ∈ (0, ε),   0 < ε ≪ δ ≪ 1.
```

Chosen values (Yu–Zeng eq. (17), *"were chosen in [7] to obtain four big size limit cycles"*):
```
a₂ = −10,  b₂ = 2.2,  c₂ = 0.7,  α₂ = −72.7778,  β = 0.0015
```
Fixed points: **E₀ = (0,0)** and **E₁ = (−6.2596, 7.4498)**.

**Distribution: (3,1).** Verbatim Yu–Zeng p.8–9: *"the simulated three limit cycles around E₀ ... and one big unstable limit cycle around E₁."* Sizes: all four **normal ("big") size** — this is the most easily visualizable four-cycle quadratic system in the literature (the outer one has radius ~4000 in the figures).

## 1.6 Yu & Han 2010/2012 — first four cycles in a NEAR-INTEGRABLE quadratic system

**Citation.** Pei Yu, Maoan Han, *Four limit cycles from perturbing quadratic integrable systems by quadratic polynomials*, arXiv:1002.1055 (4 Feb 2010); Int. J. Bifurcation Chaos **22**(10) (2012) 1250254. **PDF: `papers/yu_han_2012_four_limit_cycles_perturbing_quadratic_integrable_1002.1055.pdf`. OPEN ACCESS.**

Verbatim abstract: *"In this paper, we give a positive answer to the open question: Can there exist 4 limit cycles in quadratic near-integrable polynomial systems? It is shown that when a quadratic integrable system has two centers and is perturbed by quadratic polynomials, it can generate at least 4 limit cycles with (3,1) distribution. The method of Melnikov function is used."*

**Perturbed system** (their eq. (9), a perturbation of the reversible class Q₃^R):
```
ẋ = y(1 + a₁x) + ε a₁₀ x
ẏ = −x + x² + a₄ y² + ε (b₀₁ y + b₁₁ x y)
```
Unperturbed (ε=0) centers at **(0,0)** and **(1,0)**.

Center classification used (their Theorem 1.1, in these coordinates):
```
Q₃^R  (reversible)      : a₃ = a₂ = 0
Q₃^H  (Hamiltonian)     : a₃ = a₁ + 2a₄ = 0
Q₃^LV (Lotka–Volterra)  : a₂ = 1 + a₄ = 0
Q₄    (codim-4)         : a₃ − 5a₂ = a₁ − (5 + 3a₄) = a₄ + 2(1 + a₂²) = 0
```

**Theorem 2.1 (verbatim):** *"When a₁ < −1, the quadratic near-integrable system can have small limit cycles bifurcating from the two centers (0,0) and (1,0) with distributions: (3,0), (0,3), (2,0), (0,2) and (1,1). (2,1)- or (1,2)-distribution does not exist."*

**Theorem 4.1 (verbatim):** *"For the case of bifurcation of small limit cycles from the two centers (0,0) and (1,0) with (3,0)-distribution (respectively, (0,3)-distribution) there exists at least one large limit cycle ... For the case of limit cycles with (2,0)-distribution (respectively, (0,2)-distribution) there exist at least two large limit cycles..."*

Yielding four cycles by four routes — **all totalling (3,1) or (1,3), never (4,0)**:

| case | small cycles | large cycles | total | parameter conditions |
|---|---|---|---|---|
| (A) | (3,0) | (0,1) | **(3,1)** | a₁=−30/7, a₄=(a₁−5)/3−ε₁, b₁₁=[(a₁+2a₄)(1+a₄−a₁)/(1+a₄)]a₁₀−ε₂, b₀₁=−a₁₀−ε₃ |
| (B) | (0,3) | (1,0) | **(1,3)** | a₁=−70/51, a₄=(6a₁+5)/3−ε₁, b₁₁=[(a₁+2a₄)(2a₁−a₄+1)/((1+a₁)²(a₁−a₄+1))]a₁₀−ε₂, b₀₁=−[(b₁₁+2a₄−1)/(1+a₁)]a₁₀−ε₃ |
| (C) | (2,0) | (1,1) | **(3,1)** | a₁=−4, a₄=−18/5−ε₁, b₁₁=(392/65)a₁₀−ε₂, b₀₁=−a₁₀−ε₃ |
| (D) | (0,2) | (1,1) | **(1,3)** | a₁=−4/3, a₄=−6/5−ε₁, b₁₁=(1176/65)a₁₀−ε₂, b₀₁=−(513/65)a₁₀−ε₃ |

**The decisive theoretical paragraph, verbatim, Yu–Han p.4:**
> *"It should be mentioned that Zhang [28] has proved that the possible cycle distributions in general quadratic systems with two foci must be (0,1)-distribution or (1,i)-distribution, i = 0,1,2,3,···. **So far, no results have been obtained for i ≥ 4.** This result also rules out the possibility of (2,2)-distribution. **It is conjectured that at most 3 limit cycles can exist around one focus point.** The problem of bifurcation of 3 limit cycles near an isolated homoclinic loop is still open."*

Also note the related Yu–Han IJBC 2013, *Eight limit cycles around a center in quadratic Hamiltonian system with third-order perturbation* (`papers/yu_han_eight_limit_cycles_around_center_quadratic_IJBC2013.pdf`) — **cubic** perturbation of a quadratic Hamiltonian, so it does not bear on H(2); it does show what a single nest can do once you leave degree 2.

## 1.7 Yu & Zeng 2020/2021 — visualization; four NORMAL-SIZE near-integrable cycles

**Citation.** Pei Yu, Yanni Zeng, *Visualization of Four Limit Cycles in Near-Integrable Quadratic Polynomial Systems*, arXiv:2002.09987 (23 Feb 2020). **PDF: `papers/yu_zeng_2021_visualization_four_limit_cycles_near_integrable_2002.09987.pdf`. OPEN ACCESS.** (Also duplicated in the dir as `kuznetsov_visualization_four_lc_near_integrable_2002.09987.pdf`.)

This is the single most useful open-access source for explicit four-cycle coefficient vectors: it reproduces and simulates **Shi (1980)**, **Chen–Wang (1979)**, **Kuznetsov et al. (2013)**, and produces a new one.

**Their new fully-explicit near-integrable system, eq. (20)** — realizes Yu–Han case (A) with `a₁₀ = 0.005, ε₁ = 0.1, ε₂ = 3×10⁻⁵, ε₃ = 10⁻⁸, ε = 1/100`, i.e.
`a₁ = −30/7, a₄ = −671/210, a₁₀ = 1/200, b₀₁ = −500001/100000000, b₁₁ = 49182857/968100000`:

```
ẋ = y(1 − (30/7) x) + (1/100)·(1/200) x
ẏ = −x + x² − (671/210) y² + (1/100)·( −(500001/100000000) y + (49182857/968100000) x y )
```

Focus values at (0,0):
```
v₀ = −1/200000000
v₁ = 461/56000000
v₂ = −36481804571/14817600000000
v₃ = 11326448548182069181/25092716544000000000
```
Approximate small-cycle amplitudes: **r₁ ≈ 0.026385, r₂ ≈ 0.079134, r₃ ≈ 0.140966**.

**Distribution: (3,1), normal size.** Verbatim from the abstract: *"finally provide a concrete near-integral quadratic polynomial system to show four normal size limit cycles."*

Their reported simulation difficulty is itself evidence about the geometry: even for these normal-size cycles, *"the simulation for this near-integrable system is extremely difficult since the convergence speed is too slow"* — R-K 4th order, with per-cycle time steps, iteration counts, and CPU times tabulated in their Table 1.

## 1.8 Gaiko's four-cycle construction (constructive part only; the ≤4 claim is §2.2)

**Citation.** V. A. Gaiko, *Geometry of Planar Quadratic Systems*, arXiv:math/0611142. **PDF: `papers/gaiko_geometry_planar_quadratic_math0611142.pdf`. OPEN ACCESS.**

**Theorem 3.1 (verbatim): "A quadratic system can have at least four limit cycles in the (3:1)-distribution."**

Explicit canonical system used (his eq. (3.1), the case a = 1/2, c = −1):
```
ẋ = −y(1 + x + α y)
ẏ = x + (λ + β + γ) y + (1/2) x² + (α + β + γ) x y − γ y²
```
Construction: start with α=β=γ=λ=0, giving two centers symmetric about the x-axis at **(0,0)** and **(−2,0)**. Then successively switch on four field-rotation parameters:
1. `0 < γ ≪ 1` → centers become foci: (0,0) unstable, (−2,0) stable.
2. `−1 ≪ λ < −γ < 0` → (0,0) changes stability, generating an **unstable** cycle.
3. `γ + λ ≪ α < 0` → destroys the line x = 1, generating **two** cycles from separatrix cycles: a **stable** one around (0,0) and an **unstable** one around (−2,0).
4. `0 < −γ − λ < β ≪ 1` → (0,0) changes stability again, generating a **stable** cycle.

Net: **three around (0,0), one around (−2,0) — (3,1).** No small-amplitude/normal-size distinction is drawn; all four come from field rotation, so all are normal size.

---

# PART 2 — Every claim of ≥5 (and every claim to have *proved* H(2)=4), with objections

## 2.1 Shi Songling's alleged 1978 five-cycle preprint

**NOT FOUND.** Repeated targeted searching turned up **no evidence that a Shi Songling preprint claiming five limit cycles ever existed**, and no rebuttal by Qin Yuanxun to such a claim. What exists:
- Shi Songling, Sci. Sinica **23** (1980) 153–158 — the **four**-cycle example. (Chinese version: Kexue Tongbao/Sci. Sinica 1979.)
- Qin Yuanxun, Shi Songling, Cai Suilin, *On limit cycles of planar quadratic systems*, Sci. Sinica **25** (1982) 41–50 — a **collaboration**, not a rebuttal.
- Qin Yuanxun ("Chin Yuan-shun") et al.'s own example is a **(3,0)** system (three cycles, one nest), recorded verbatim in Perko 1984 §6.

**Assessment: this appears to be a conflation.** The genuinely retracted five-cycle-adjacent claim in this history is **Petrovskii–Landis (1955/1957)**, who claimed to *prove* H(2)=3 and later acknowledged an error. Verbatim from Koditschek–Narendra 1984 (`papers/koditschek_narendra_1984_limit_cycles_planar_quadratic_JDE54.pdf`, p.181): *"Three years later, a paper by Petrovskii and Landis [3] purported to [show] that a quadratic system could support no more than three cycles on the plane. Although this result was called into question by several [authors] (and the authors later acknowledged an error in the proof [4]) it inspired a number of attempts..."* — an erroneous *upper*-bound claim, refuted by Chen–Wang and Shi's four-cycle examples. If a five-cycle Shi preprint exists it is not visible in any indexed source reachable from here.

## 2.2 Gaiko — "H(2) = 4 and only (3:1)"

**Citations.**
- V. A. Gaiko, *Global Bifurcation Theory and Hilbert's Sixteenth Problem*, Mathematics and Its Applications **562**, Kluwer/Springer, 2003, 208 pp. **Book, paywalled.**
- V. A. Gaiko, *Geometry of Planar Quadratic Systems*, arXiv:math/0611142. **PDF held.**
- Related, all held and open access: `gaiko_lienard_smale_math0611143.pdf`, `gaiko_quadratic_two_parallel_isoclines_0803.3055.pdf`, `gaiko_1104.3019.pdf`, `gaiko_general_lienard_1202.3540.pdf`, `gaiko_1504.03353.pdf`, `gaiko_1611.08113.pdf`, `gaiko_field_rotation_ihes.pdf`, `gaiko_hilbert16_quadratic_lienard_transformation_math0701869.pdf`.

**The claim, verbatim from the math/0611142 abstract:**
> *"Using geometric properties of four field rotation parameters of a new canonical system which is constructed in this paper, we present a proof of our earlier conjecture that the maximum number of limit cycles in a quadratic system is equal to four and the only possible their distribution is (3:1) [10]. Besides, applying the Wintner–Perko termination principle for multiple limit cycles to our canonical system, we prove in a different way that a quadratic system has at most three limit cycles around a singular point (focus) and give another proof of the same conjecture."*

**Theorem 3.2 (verbatim): "A quadratic system has at most four limit cycles and only in the (3:1)-distribution."**

**Theorem 4.3 (verbatim, the direct four-in-one-nest claim):**
> *"There exists no quadratic system having a swallow-tail bifurcation surface of multiplicity-four limit cycles in its parameter space. In other words, **a quadratic system cannot have neither a multiplicity-four limit cycle nor four limit cycles around a singular point (focus)**, and the maximum multiplicity [is three]."*

**Argument sketch.** Two routes. (i) Geometric: reduce to a canonical form with four field-rotation parameters λ, α, β, γ; introduce them successively; argue from the geometry of the spirals filling the interior/exterior of the cycles that no fifth cycle can be produced and that the (2:2) case is impossible. (ii) Wintner–Perko: *"Any one-parameter family of limit cycles belongs to a maximal one-parameter family which is either open or cyclic. If the family is open, then it terminates as the parameter or the orbits become unbounded, or it terminates at a critical point or on a (compound) separatrix cycle."* Gaiko argues a maximal one-parameter family of multiplicity-four cycles must terminate somewhere, and every possible termination contradicts Bautin's cyclicity-3 result.

**THE PUBLISHED OBJECTION — Roussarie's MathSciNet review MR2023976 (2005d:37102), quoted verbatim by Gaiko himself in §5 of math/0611142:**
> *"I just mention the hazardous claim made in Theorem 4.12, page 137, that there exists no quadratic system having a swallow-tail bifurcation surface of multiplicity-four limit cycles. Looking at the proof, it seems that the author unfortunately confuses two different notions: paths of limit cycles, as defined in Definition 4.7, page 112, and lines of multiple limit cycles, as defined by Perko (and recalled in Definition 4.13, page 127). In fact, there is nothing forbidding that a path begin at a parameter value with a multiplicity-four limit cycle and end at a focus point."*

**Gaiko's counter-rebuttal (verbatim, same section):** he calls the review *"awkward"*, notes an apparent definition-number mismatch, and replies that Roussarie *"unfortunately does not see (or does not want to see) Bautin's result [1] (Theorem 2.1, page 45) on the cyclicity of a singular point of the focus or center type, which is an obstacle on such a path. Or, maybe, Bautin's result is also 'questionable'?"* He also attacks the competing Roussarie/Dumortier/Rousseau programme as unable to determine cyclicity of non-monodromic separatrix cycles.

**Status.** Not accepted by the community. **The exact point Roussarie disputes is precisely the four-in-one-nest question.** Note also Gaiko's own admission at the end of §4: *"to complete the solution of Hilbert's Sixteenth Problem for quadratic systems (1.1), it is sufficient to prove impossibility of the (2:2)-distribution of limit cycles only in the case of two finite foci and a saddle at infinity"* — but (2,2) impossibility was independently established by Zhang Pingguang, and it does **not** yield H(2)=4; it only reduces H(2) to the one-nest question (see §4.1). So the reduction Gaiko relies on does not close the problem even if granted.

## 2.3 Pedregal — "H(2) = 4" via variational global analysis

**Citation.** Pablo Pedregal, *A variational approach to Hilbert's 16th problem within the framework of global analysis*, arXiv:2103.07193 (v1 12 Mar 2021; **v3 30 Aug 2024**), 69 pp. **PDF: `papers/pedregal_2021_H2_equals_4_claim_2103.07193.pdf`. OPEN ACCESS. Not published in a peer-reviewed journal as of this sweep.**

**The claim, verbatim abstract:**
> *"We focus on the second part of Hilbert's 16th problem and provide an upper bound on the number of limit cycles that a polynomial, differential, planar system may have, depending exclusively on the degree n of the system. Such a bound turns out to be a polynomial of degree 4 in n. More specifically, if H(n) indicates the maximum number of limit cycles among planar, differential, polynomial systems of degree n, then*
> *H(n) ≤ (5/2)n⁴ − (23/2)n³ + (43/2)n² − (37/2)n + 7 if n is even, and*
> *H(n) ≤ (5/2)n⁴ − (23/2)n³ + (41/2)n² − (33/2)n + 6 if n is odd.*
> ***For quadratic systems, we find H(2) = 4.*** *Our proof is entirely variational and utilizes in a fundamental way tools and facts from global analysis to the point that no particular expertise in dynamical systems is necessary or required."*

**Central theorem (verbatim, Theorem 1.1).** Under three hypotheses — (1) P, Q have no common non-trivial factor; (2) all connected components of the algebraic curve `Pₓ + Q_y = 0` are homeomorphic to a line or an oval (no singular points), M of them; (3) the system `P(Pₓₓ + Q_yₓ) + Q(Pₓ_y + Q_yy) = 0, Pₓ + Q_y = 0` has only N simple solutions — then
```
H(n) ≤ 1 + (n − 1)²(M + N).
```
Bezout and Harnack then bound M and N in terms of n, giving the quartic. For n=2: `(n−1)² = 1`, so `H(2) ≤ 1 + M + N`, and the divergence curve `Pₓ + Q_y = 0` is a **line**, so M = 1; the contact system has N ≤ 2 simple contacts; 1 + 1 + 2 = 4.

**Objections.**
1. **It is not the paper refuted by Buzzi–Novaes.** arXiv:2411.09594 rebuts a *different* H(2)=4 claim (§2.4). No paper specifically rebutting Pedregal was found.
2. **Prior-version history.** Pedregal's earlier attempts on the same programme — arXiv:1411.6814 (*Hilbert's 16th problem. When variational principles meet differential systems*, with Llibre) and arXiv:1904.01300 (*Hilbert's 16th problem. II. Pfaffian equations and variational methods*) — were **withdrawn after an error was found in the counting method for limit cycles**. That is the same step (§5.9 "Counting method" in the current version) on which the whole argument rests.
3. **Structural.** The hypotheses (2) and (3) are *genericity* assumptions; the paper carries a §5.10 "Non-generic situation". A uniform bound must survive the non-generic stratum, and the reduction of the non-generic case to the generic one is where a variational multiplicity argument is most fragile.
4. **The quartic bound does not contradict the n² log n growth**, so the Buzzi–Novaes asymptotic objection does not apply to Pedregal — this claim is at least *not obviously* false, unlike §2.4.
5. **Community status:** unrefereed, uncited as established, and Gasull–Santana (Oct 2024, §2.7) — writing after Pedregal v3 — still state flatly that *"it is not known if H(2) < ∞"*, which is a de facto non-acceptance.

## 2.4 da Silva, Vieira & Leonel — "H(n) = 2(n−1)(4(n−1)−2)", hence H(2)=4 — REFUTED

**Citation.** V. B. da Silva, J. P. Vieira, E. D. Leonel, *Exploring limit cycles of differential equations through information geometry unveils the solution to Hilbert's 16th problem*, Entropy **26**(9) (2024) 745. Open access at MDPI / PMC11431191 (host blocked download here; abstract and content recovered via the rebuttal, which quotes it in full).

**Claim (verbatim as quoted in the rebuttal, their Theorem 4):** `H(n) = 2(n − 1)(4(n − 1) − 2)` for n ≥ 2. For n = 2 this gives `2·1·(4−2) = 4`, i.e. **H(2) = 4**.

Method: a scalar curvature R of a Fisher information metric,
```
R = (1/√G)[ ∂/∂x( (1/√G) ∂G₂₂/∂x ) + ∂/∂y( (1/√G) ∂G₁₁/∂y ) ],
G₁₁ = 2[(∂P/∂x)² + (∂Q/∂x)²],  G₂₂ = 2[(∂P/∂y)² + (∂Q/∂y)²],  G = G₁₁G₂₂,
```
with limit cycles *redefined* as states where R is positive near equilibria and |R| is singular; the count of distinct divergences of |R| is asserted to be the maximum number of limit cycles.

**REFUTATION — Claudio A. Buzzi & Douglas D. Novaes, *A note on a recent attempt to solve the second part of Hilbert's 16th Problem*, arXiv:2411.09594 (14 Nov 2024). PDF: `papers/note_on_recent_attempt_hilbert16_rebuttal_2411.09594.pdf`. OPEN ACCESS.**

Two independent kill shots, both verbatim:
1. **Asymptotic contradiction.** *"H(n) grows asymptotically as fast as n² log n. A direct consequence of this growth estimation is that H(n) cannot be bounded from above by any quadratic polynomial function of n. ... Since this expression is quadratic in n, it contradicts the established asymptotic behavior and, therefore, cannot hold."* Concretely, against Christopher–Lloyd / Li–Chan–Chung's `H(2k−1) ≥ S_k = 4^(k−1)(k − 13/6) + (2k − 1)/3`, the claimed formula gives `H(2k−1) = 4(2^k − 2)(2^(k+1) − 5)`, *"which contradicts (4) for k ≥ 35."*
2. **The definition is neither necessary nor sufficient.** Four explicit counterexamples:
   - `ẋ = −y + x(x²+y²−1), ẏ = x + y(x²+y²−1)` — unique limit cycle (the unit circle), but **R(0,0) = −1 < 0**, so the criterion misses it.
   - `ẋ = −y + x(x²+y²−1)(x²+y²−4), ẏ = x + y(...)` — two limit cycles, **R(0,0) = −80/289 < 0**.
   - The first system after `(x,y) = (u, u + v/2)`: **R(0,0) = 6/5 > 0** but R₂ never vanishes, so **|R| has no singularities** despite a limit cycle existing.
   - `ẋ = −y + x², ẏ = x + xy` (isochronous quadratic centre S₂): **no limit cycles at all**, yet `R = 1/[(x²+1)²(4x²+(y+1)²)]` has R(0,0)=1 > 0 and |R| singular at (0,−1) — the criterion fires falsely.

**Status: definitively refuted.** Note the irony for our purposes: this refuted paper's *value* for H(2) happens to be 4.

## 2.5 Voronin & Lebedev 2022 — SIX limit cycles in a (3:3) distribution

**Citation.** A. V. Voronin, S. S. Lebedev, *Quadratic system of two differential equations with six limit cycles* / *Example of a quadratic system of two differential equations with six limit cycles*, ~June 2022. Hosted at `https://repository.hneu.edu.ua/bitstream/123456789/29182/1/Voronin_A.V._Quadratic_system_ResearchQate.pdf` (Simon Kuznets Kharkiv National University of Economics repository).

**PDF: INACCESSIBLE.** The host (212.111.204.10) refused all connections during this sweep (`ECONNREFUSED` direct, `TimeoutError` via proxy, curl silent failure over both http and https). Almost certainly wartime infrastructure. Not on arXiv (`au:Voronin AND all:"limit cycles"` returns nothing; no OpenAlex record).

**Claim, as recovered from the indexed abstract text:** two foci and **six limit cycles in a (3:3) arrangement**; a follow-up paper verifies the earlier result by a second method, *"both the method known from literature and the method proposed by the authors give the same result: in such nonlinear dynamic systems there are two foci and six limit cycles in the 3:3 arrangement."*

**OBJECTION — this is excluded by a theorem, twice over.**
- **Zhang Pingguang (2002), Theorem 1** (`papers/zhang_pingguang_2002_QTDS_abstract.pdf`, verbatim): *"A quadratic system having two foci has at most one limit cycle surrounding one of its two foci."* Hence *"the limit cycles of the quadratic system with two foci must be (0,1)–distribution or (1,i)-distribution (i = 0,1,2,...)."* A (3,3) distribution requires 3 cycles around **each** focus and is flatly impossible.
- **Zegeling (2024), Theorem 1.2** (independent verification, gaps in Zhang's proof filled): *"A quadratic system (1.1) can only have an (n,1) or (n,0) distribution of limit cycles where n ∈ ℕ ∪ {0}."*
- Zegeling's historical remark reads almost as if written about this claim, verbatim: *"Several attempts were made in the following years to prove that a (2,2) distribution is not possible or to find a (3,2) or even a (3,3) distribution. **These attempts remained futile.**"*

**Assessment.** Unless Voronin–Lebedev's example somehow evades the hypotheses of Zhang/Zegeling (both of which cover the generic two-foci case exhaustively), the six-cycle (3:3) claim is **false**. The most likely failure mode is numerically mistaking non-isolated or spurious closed curves — or trajectories in a slowly-converging annulus — for limit cycles, which the extreme scale separation in quadratic systems makes very easy (this is exactly the failure Galias–Tucker's rigorous method is designed to preclude). Verification is blocked pending the PDF becoming reachable.

## 2.6 "Malyarets et al. 2026, six cycles" — NOT FOUND

Exhaustive search returns **nothing**: `au:Malyarets` on the arXiv API returns zero results; `all:Malyarets AND all:"limit cycles"` zero; web search returns no such paper. Note that **L. M. Malyarets is a colleague of A. V. Voronin at Simon Kuznets Kharkiv National University of Economics** — this is very likely a garbled reference to the Voronin–Lebedev six-cycle claim in §2.5 (same institution, same claimed count of six). Treat §2.5 as the real item.

## 2.7 "Hernandez Rosales 2026, H(2)=7" — NOT FOUND

**No such paper exists in any reachable index.** arXiv `all:"Hernandez Rosales"` returns only Maribel Hernández Rosales's computational-biology papers (phylogenetics, best match graphs, Mexico City mobility networks) — a completely different researcher and field. No 2026 H(2)=7 claim on arXiv, in web search, or in OpenAlex. **Report as non-existent / fabricated reference.**

**What actually exists in 2026 on arXiv touching Hilbert 16:**
- **arXiv:2604.12883** (14 Apr 2026), O. Eshkobilov, S. Kadyrov, K. Mamayusupov, *Limit-Cycle Replication via Chebyshev Pullbacks and a Quadratic Ceiling for Separable Schemes*. **PDF held.** Proves `H(nm + m − 1) ≥ m² H(n)` for m ≥ 2 via the separable Chebyshev covering `Φ(u,v) = (T_m(u), T_m(v))`, and a **negative** result: pure separable replication can only give quadratic-in-degree growth, so superquadratic lower bounds need other mechanisms. New explicit lower bounds `H(14) ≥ 252, H(29) ≥ 1080, H(31) ≥ 1380, H(39) ≥ 2012`. **Says nothing new about H(2)** — it *consumes* H(2) ≥ 4 as a seed.
- **arXiv:2602.22558** (Feb 2026), Yovani Villanueva & Warwick Tucker, *Darboux-type center conditions for families of planar polynomial vector fields*. **PDF held.** **This is NOT a limit-cycle certification paper** — see §3.2.

## 2.8 Other H(2)-adjacent claims encountered

- **Gaiko, arXiv:math/0701869**, *Hilbert's 16th Problem for Quadratic Systems. New Methods Based on a Transformation to the Lienard Equation*. **PDF held.** Same programme as §2.2.
- Gaiko's Liénard/Smale companion `math/0611143` and the later Kukles-cubic and general-Liénard papers all rest on the same Wintner–Perko machinery and inherit the Roussarie objection in spirit.

---

# PART 3 — Computer-assisted and rigorous-numerical certification

## 3.1 Galias & Tucker 2022 — the gold standard: EXACTLY four, rigorously

**Citation.** Z. Galias (AGH Kraków), W. Tucker (Monash), *The Songling system has exactly four limit cycles*, Applied Mathematics and Computation **415** (2022) 126691. **OPEN ACCESS (CC-BY). PDF: `papers/galias_tucker_2022_songling_exactly_four_limit_cycles_AMC415.pdf`.**

**Theorem 1 (verbatim): "The Songling system (1) has exactly four limit cycles."**

Tools: **CAPD library** for rigorous Taylor integration (order 100), **MPFR** for multiple precision (up to **1024 bits** for existence, **2048 bits** for the interval-Newton uniqueness step). Runtime 4 s–14 s per fixed point for existence; **45 s to 40 minutes** per fixed point for the Newton step.

**Setup.** On Σ = {x = 0}: `ẋ = −y + y², ẏ = 0`, so ẋ < 0 for y ∈ I = (0,1). Every periodic orbit must cross `(0, I)`. Define the return map P : I → I.

**Non-rigorous scan (§3.1).** Sampling `y_k = 10^(−k/3)`, k=1..300, and plotting `f(y) = y − P(y)`: four sign changes in `y ∈ [10⁻¹⁰⁰, 10^(−1/3)]`, in
`[4.64·10⁻⁷⁵, 10⁻⁷⁴]`, `[2.15·10⁻²¹, 4.65·10⁻²¹]`, `[4.64·10⁻⁸, 10⁻⁷]`, `[2.15·10⁻², 4.64·10⁻²]`.
**This is the single best published picture of just how brutal the scale separation is: the four cycles live at radii spanning 73 orders of magnitude.**

**Lemma 2 (existence, verbatim intervals containing fixed points of P):**
```
[y₁l, y₁r] = 0.042689603882085(75..85)
[y₂l, y₂r] = 6.666660148(1..2) · 10⁻⁸
[y₃l, y₃r] = 2.24780594(7..8) · 10⁻²¹
[y₄l, y₄r] = 7.07106781186547524(4..5) · 10⁻⁷⁵
```
with verified inequalities e.g. `P(y₁l) − y₁l < −5.06·10⁻¹⁵`, `P(y₄r) − y₄r < −6.23·10⁻²⁹³`.

**Lemma 3 (uniqueness + stability, interval Newton).** Enclosures to ~90+ digits, with multipliers:
```
y₁ = 0.0426896038820800619842959753055429655870014298074 9(6..8)      STABLE, P'(y₁) ⊂ 9.11(8..9)·10⁻⁵
y₂ = 6.66666014815265057395069851799647931628168529042733641712740 6(88..91)·10⁻⁸   UNSTABLE, P'(y₂) ⊂ 1.0…0049(1..2)
y₃ = 2.2478059477961305860583886189574201301744379417437915330332524455975299877707946920093780(1..3)·10⁻²¹   STABLE, P'(y₃) ⊂ 0.99…9936(5..6)
y₄ = 7.0710678118654752440084436210484903928483593768847403658833986899536623915772018609186193383048745303417359328082824292486479850624247482966280728718679061243 5(4..6)·10⁻⁷⁵   UNSTABLE, P'(y₄) ⊂ 1.0…062(8..9)
```
Verbatim: *"One can see that for fixed points in y₂,₃,₄ the derivative of P is very close to one. As a consequence a very accurate integration method has to be used."* — this is the quantitative reason four-cycle quadratic systems resist naive numerics.

**Distribution, stated verbatim in the Fig. 2 caption:** *"Note that three of the four limit cycles surround the equilibrium at the origin."* → **(3,1)**, confirmed rigorously.

**The uniqueness half (§4).** Not self-contained: verbatim, *"From [12] it follows that we only have to prove that no limit cycles of (1) intersect the line segment (x,y) ∈ {0} × [0.004, 0.04]."* Reference [12] is the normal-form argument establishing exactly 3 fixed points of P in (0, 0.004) and exactly 1 in (0.04, 1) — the latter because the line `1 + x − 25y = 0` is transversal to the field and *"the Songling system ... has at most one limit cycle around one of the equilibria."* Galias–Tucker's contribution is closing the gap segment [0.004, 0.04] by interval methods. **So the "exactly three around the origin" part is inherited from normal-form theory, not from interval arithmetic** — worth knowing if one wants a fully independent certification.

**Techniques (reusable for a five-cycle hunt).** Four exclusion lemmas:
- **Lemma 1** (topological): if `g(x) ⊂ x` or `x ⊂ g(x)` then g has a fixed point in x.
- **Lemmas 4 & 5** ("iteration based method"): if `ȳ < P(ȳ)` there are **no** fixed points in `[ȳ, P(ȳ)]` (and the analogue for P⁻¹). Chaining these sweeps out fixed-point-free intervals, *including across points where P is undefined* — they use this to jump the escape-to-infinity point `y* ∈ [0.03689093, 0.03689096]` by showing `P(0.0364) > 0.3766`, hence no fixed points in `[0.0364, 0.3766] ∋ y*`.
- **Lemma 6** ("derivative based method"): if `1 ∉ P'(y)` then `f = y − P(y)` is strictly monotone on y, so sign(f) at the endpoints decides between zero and exactly one fixed point.
- **Lemma 7** (Lyapunov function method): `P⁻¹(y) > y` for `y ∈ [0.98, 1)`, proved after shifting (0,1) to the origin via `z = y − 1`.

Verbatim on generality: *"The techniques presented here are applicable to the much wider class of real-analytic planar differential equations."* **This is the toolchain a five-cycle search should be built on.**

Verbatim on why the problem is hard: *"One of the main difficulties is that the limit cycles can reside within areas of vastly different scales. This makes numerical explorations very hard to perform, requiring high precision computations, where the necessary precision is not known in advance. Using rigorous computations, we can dynamically determine the required precision, and localize all limit cycles of a given system."*

## 3.2 Villanueva & Tucker 2026 — arXiv:2602.22558 — MISCHARACTERIZED IN THE BRIEF

**Citation.** Yovani Villanueva, Warwick Tucker, *Darboux-type center conditions for families of planar polynomial vector fields*, arXiv:2602.22558, math.DS, Feb 2026, 19 pp. **PDF: `papers/villanueva_tucker_2026_2602.22558.pdf`. OPEN ACCESS.**

**Verbatim abstract:**
> *"We study the center-focus problem for planar polynomial vector fields, which can be viewed as a local version of Hilbert's 16th problem. Based on a Lyapunov function approach, we establish novel results regarding the center-focus conditions for two families of differential systems. More precisely, we find an enclosure of the Bautin ideal generated by the Lyapunov constants of these systems. Our results hold for any degree n ≥ 2."*

**This is NOT a limit-cycle certification paper and contains no four-cycle result.** It is about the **center-focus problem** and **Bautin ideal enclosures** via Lyapunov functions `V = Σ V_k`, solving linear systems for the coefficients `v_{i,j}` and the Lyapunov constants `L_k`. Relevance to a five-cycle hunt is **indirect but real**: Bautin-ideal structure is exactly what caps small-amplitude cyclicity at a focus (Bautin's 3 for quadratic), and any (4,·) example needs a mechanism outside that cap.

## 3.3 The Poincaré–Bendixson / transversal-curve toolchain

- **A. Gasull, H. Giacomini, M. Grau, *Proving the existence of numerically detected planar limit cycles*, arXiv:1602.00113 (30 Jan 2016). PDF held. OPEN ACCESS.** Verbatim: *"We provide a method to construct Poincaré–Bendixson regions by using transversal curves, that enables us to prove the existence of a limit cycle that has been numerically detected."* Applied to the Brusselator, Liénard systems, and (with other tools) to sharp saddle-node bifurcation values for the Rychkov system. **Key limitation for our problem, noted by Galias–Tucker verbatim about the analogous [7]: "the method presented in [7] cannot be used to prove the uniqueness of a limit cycle in a specified region."** — it gives existence (lower bounds), never upper bounds. That is precisely the right tool for a five-cycle *discovery* effort.
- **Follow-ups held:** `gasull_transversal_conics_existence_limit_cycles_1410.4480.pdf` (*Transversal conics and the existence of limit cycles*); `garcia_saldana_gasull_giacomini_new_approach_limit_cycles_1910.08098.pdf` (JDE 2020 — star-like limit cycles as solutions of an associated non-autonomous planar system or heteroclinics of a 3-D polynomial system); `gasull_extended_bendixson_dulac_1305.3402.pdf` (*Some Applications of the Extended Bendixson–Dulac Theorem*); `number_limit_cycles_planar_invariant_algebraic_curves_2210.15803.pdf` (QTDS 2023, Gasull–Giacomini, non-existence of periodic orbits off a given algebraic curve, applied to quadratic systems); `smooth_transformations_ruling_out_closed_orbits_2309.02513.pdf`. Also Gasull–Giacomini, *Effectiveness of the Bendixson–Dulac theorem*, JDE **305** (2021) — paywalled, not obtained.
- **Cherkas's Dulac-function optimization** (as used in Cherkas–Artés–Llibre 2003, §1.4) is the one method in the literature that has actually delivered **exact** counts (upper *and* lower) for normal-size four-cycle quadratic systems by hand. See also L. A. Cherkas, *Bendixson–Dulac criterion and reduction to global uniqueness in the problem of estimating the number of limit cycles*, Differ. Equ. **46**(1) (2010) — paywalled.

## 3.4 Large-scale numerical exploration of the 5-parameter space

**Finding: no systematic exhaustive sweep of the quadratic parameter space for ≥5 cycles has ever been published.** What exists is structurally adjacent but not this:

- **Artés, Llibre, Schlomiuk, Vulpe, *Geometric Configurations of Singularities of Planar Polynomial Differential Systems: A Global Classification in the Quadratic Case*, Birkhäuser/Springer, 2021, 699 pp.** Book, paywalled. A **complete bifurcation diagram in the 12-parameter space** of global geometric configurations of *singularities* of quadratic systems, in invariant form; **1765 distinct global geometric configurations of singularities**. This is the definitive atlas of *where in parameter space to look*, but it classifies singularities, **not limit cycles**. A five-cycle hunt should use it as the stratification.
  - Related open-access material held: `schlomiuk_vulpe_geometry_neighbourhood_infinity_math0405026.pdf`; `llibre_schlomiuk_2004_CJM56_weak_focus_third_order.pdf` (Canad. J. Math. 56, the geometry of quadratic systems with a **third-order weak focus** — the exact stratum every (3,·) example sits in).
- **Leonov's "hidden oscillations" group (Kuznetsov, Kuznetsova, Leonov)** did map *"domains of parameters corresponding to existence of different configurations of large limit cycles"* and *"the domain of parameters of quadratic systems for which four limit cycles can be obtained"* (DEDS 2013; Differential Equations 49 (2013)). This is the closest thing to a parameter-space search, but it is confined to the (3,1)-producing region of a 5-parameter subfamily and, being an *existence* search, could not have detected a fifth cycle it did not construct.
- **Galias–Tucker's Fig. 1 scan** is a rigorous one-dimensional sweep of the *return map* of a **single** system — not a parameter sweep. Extending it to a parameter sweep is technically straightforward and, as far as this survey can determine, **has not been done**.

**This is the clearest gap in the literature and the most obvious place for new work.**

---

# PART 4 — Explicit upper bounds: restricted and by subclass

## 4.1 The distribution theorems — the key structural constraints

**Zhang Pingguang (2002).** *On the Distribution and Number of Limit Cycles for Quadratic Systems with two Foci*, Qual. Theory Dyn. Syst. **3** (2002) 437–463. **PDF: `papers/zhang_pingguang_2002_QTDS_abstract.pdf`.**

Verbatim abstract: *"In this paper, we study the distribution and number of limit cycles for quadratic systems with two foci. It is proved that **a quadratic system with two foci has at most one limit cycles around one of the two foci**, and hence the limit cycles of the quadratic system with two foci must be **(0,1)–distribution or (1,i)-distribution (i = 0,1,2,...)**."*

Verbatim Theorem 1: *"A quadratic system having two foci has at most one limit cycle surrounding one of its two foci."* And verbatim: *"This theorem has been proved by the author in several papers published in Chinese. Here, by first time, we present its complete proof in English."*

Normal form used: `ẋ = −y + δx + lx² + mxy + y², ẏ = x(1 + ax + by)` with `1 + b < 0`, `a ≥ 0`, `|δ| < 2`, `(m+δ)² + 4(b+l) < 0`.

**Consequences:** (2,2) impossible. (3,2), (3,3) impossible. **(4,1) and (4,0) NOT excluded.**

**Zegeling (2024)** — the independent verification. A. Zegeling, *Nests of limit cycles in quadratic systems*, Advances in Nonlinear Analysis **13** (2024) 20240012, doi:10.1515/anona-2024-0012. **OPEN ACCESS (CC-BY). PDF: `papers/zegeling2024.pdf`.**

Verbatim abstract: *"We give a proof of the distribution property of limit cycles in so-called quadratic systems. We prove that the possible limit cycle distributions are either (n,0) or (n,1) (where n ∈ ℕ ∪ {0}). The aim of this article is to simplify and fill gaps in the original proof by Zhang... **A consequence of the distribution property is that it reduces the study of H(2) to the study of the maximum number of limit cycles surrounding one singularity.**"*

**Theorem 1.2 (verbatim): "A quadratic system (1.1) can only have an (n,1) or (n,0) distribution of limit cycles where n ∈ ℕ ∪ {0}."** Hence, verbatim, *"H(2) = max{p₁ + 1, 1 + p₂}"* — i.e. **H(2) is exactly the one-nest maximum plus one.**

**Why this paper matters enormously.** Verbatim: *"It was mentioned in Section 2.2 of the book [1] (= Artés–Llibre–Schlomiuk–Vulpe) that '**this result of Zhang Pingguang has so far not been checked by the mathematical community.**' Our aim is to show that the conclusion of Zhang [22] for the distribution of limit cycles is essentially correct."* And on the specific gap: *"On page 456 in the study by Zhang [22], the check of the monotonicity of a critical function F(x)f(x)/g(x) was omitted with the justification that it was easy to see. It was not mentioned in the original Chinese papers either, and **in our opinion, the suggested way of proving the monotonicity is not the correct way.** In this article, we will show the details of this critical part of the proof."*

**The single most relevant historical sentence in the entire literature for our question, verbatim (Zegeling 2024 §1):**
> *"The first famous examples of a nest distribution of limit cycles were given in the 70's in China where p ≥ 3, q ≥ 1 [2,11]. Several attempts were made in the following years to prove that a (2,2) distribution is not possible or to find a (3,2) or even a (3,3) distribution. **These attempts remained futile.**"*

**Zegeling's method** (relevant to any new bound): reduce to Liénard form and apply two classical theorems — a **non-existence** theorem (used by Li Chengzhi for third-order weak foci; requires non-existence of solutions to a 2×2 algebraic system) and **Coppel's uniqueness theorem** (requires that same algebraic system to have ≤1 solution, plus fixed sign of `d/dx [F(x)f(x)/g(x)]` on an interval). Zegeling reformulates both so the conditions reduce to a **function of one variable only**. Cases covered: §4 two finite real strong foci + one complex pair + ≥1 elementary saddle at infinity, no invariant line (*"the most difficult one and contains the examples of four limit cycles in a (3,1)-distribution"*); §5 four real singularities (simplifying Zegeling–Kooij [17] and Zhang [22]); §6 the non-generic remainder via rotated vector fields.

**Zhang Pingguang & Cai Suilin (1991)**, *Quadratic systems with a weak focus*, Bull. Austral. Math. Soc. **44** (1991). **PDF: `papers/zhang_pingguang_cai_suilin_1991_quadratic_systems_weak_focus_BAMS44.pdf`.** Proves **at most one limit cycle not surrounding the weak focus**, under any of: (i) ≥2 finite saddles; (ii) >2 finite singular points and >1 singular point at infinity; (iii) exactly 2 finite singular points, >1 at infinity, and the weak focus surrounded by ≥1 limit cycle. This is the precursor result feeding into the 2002 distribution theorem.

**Zegeling & Kooij (1999)**, *The distribution of limit cycles in quadratic systems with four finite singularities*, J. Differential Equations **151** (1999) 373–385. Paywalled. Subsumed by Zegeling 2024 §5.

## 4.2 Ilyashenko & Llibre 2009 — the only explicit uniform-type bound

**Citation.** Yu. Ilyashenko, J. Llibre, *A restricted version of the Hilbert's 16th problem for quadratic vector fields*, arXiv:0910.3443 (19 Oct 2009); Mosc. Math. J. **10**(2) (2010) 317–335. **PDF: `papers/ilyashenko_llibre_restricted_version_H2_0910.3443.pdf`. OPEN ACCESS.**

Normalization: `ż = μz + Az² + Bzz̄ + Cz̄²` with `μ = λ₁ + i`, `λ₁ ≥ 0`, and one of `A=1, |B|≤2, |C|≤1` / `B=2, |A|≤1, |C|≤1` / `C=1, |A|≤1, |B|≤2`. Parameter space Λ ≅ three glued copies of `ℝ₊ × D² × D²`.

Three restriction parameters:
- **δ-tame:** the cycle lies in `B(λ,δ)` = the disc `|z| ≤ δ⁻¹` minus the open δ-neighbourhoods of all singular points (real and complex) except 0. *(δ-distant from singular points and infinity.)*
- **σ-distant from centers:** distance σ from the center variety (center conditions in Żołądek's complex form).
- **κ-distant from singular quadratic fields:** a *singular* quadratic field is `ż = μz l(z)` with a line of singularities; after normalization `ż = μz + z² + (μ/μ̄)zz̄ =: v_s(z)`. Writing `v = v_s + u`, `u = bzz̄ + cz̄²`, the condition is `‖r⁻²u‖₂ = √(|b|² + |c|²) > κ`.

**Theorem 5 (Main Theorem), verbatim:**
> *"For any {δ, σ, κ} ⊂ (0, 0.1), the number of δ-tame limit cycles of a normalized quadratic vector field which is σ-distant from centers and κ-distant from singular quadratic vector fields is no greater than*
> ***H(2, δ, σ, κ) = |log σ| · exp(exp(10²⁵ δ⁻³¹ κ⁻²)).***
> *This estimate is irrealistic but this is the only known estimate of this kind."*

Machinery: Growth-and-Zeros / Bernstein-index theorem `#{z ∈ K | f(z) = 0} ≤ B_{K,U}(f) exp(...)` applied to the displacement function expanded by **Bautin's algorithm**, plus explicit lower bounds for trigonometric polynomials.

Sequel promised in the paper, verbatim: *"In a subsequent paper we prove that for κ sufficiently small: κ ≤ κ₀(δ,σ), the vector field (8) has only **one** δ-tame limit cycle."* — i.e. **near the singular stratum, exactly one cycle**, so no five-cycle example can live there.

**Practical reading for a five-cycle hunt.** The bound is astronomically far from 4, so it excludes nothing. Its real content is *structural*: any hypothetical five-cycle example is not obstructed by anything in the "tame" region, and the three degeneracy directions (near-center, near-singular, cycles near singular points or infinity) are the only places where finiteness is currently uncontrolled. Combined with Galias–Tucker's observation that the four known cycles span 73 orders of magnitude, this says a fifth cycle, if it exists, most likely lives **very close to a singular point or to infinity — precisely the region δ-tameness excludes.**

## 4.3 Subclass upper bounds — the table

| Subclass | Bound | Source | Status |
|---|---|---|---|
| Quadratic with an **invariant straight line** | **≤ 1** | classical (Coppel; Bautin theory) | Sharp; the line is a Dulac/transversality obstruction |
| **Bounded** quadratic systems | conjectured **≤ 2** | Coppel programme; Zegeling & Kooij | Conjecture, not theorem |
| Quadratic with a **weak focus of 3rd order** | **≤ 3 in that nest**, and a full topological classification | Li Chengzhi (non-existence thm); Artés–Llibre, Publ. Mat. 41 (1997) (`pubmat41.pdf`); Llibre–Schlomiuk, CJM 56 (2004) (`llibre_schlomiuk_2004_CJM56_weak_focus_third_order.pdf`) | Established; 20 global phase portraits from 16 local, only 3 with any limit cycle |
| **Two foci**, one nest | **≤ 1** in the smaller nest | Zhang Pingguang 2002, Thm 1; Zegeling 2024, Thm 1.2 | Established (Zegeling verified Zhang) |
| **(2,2), (3,2), (3,3)** distributions | **impossible** | same | Established |
| **Four finite singularities** | distribution property holds | Zegeling & Kooij, JDE 151 (1999); Zegeling 2024 §5 | Established |
| Near-**singular** quadratic fields (κ small) | **exactly 1** δ-tame cycle | Ilyashenko–Llibre (announced sequel) | Announced |
| **Small-amplitude at one focus** (any quadratic) | **≤ 3** | **Bautin 1952** | Sharp; the fundamental cap |
| **Any nest** | conjectured **≤ 3** | folklore; stated as conjecture in Yu–Han 2010 | **OPEN — this is the whole problem** |
| General quadratic, δ/σ/κ-restricted | `\|log σ\| exp(exp(10²⁵δ⁻³¹κ⁻²))` | Ilyashenko–Llibre 2009 | Only known explicit bound of this type |
| **General quadratic, unrestricted** | **UNKNOWN — H(2) < ∞ not proved** | Gasull–Santana, arXiv:2407.13465, Oct 2024 | Open |

Verbatim, Gasull & Santana, *A note on Hilbert 16th Problem*, arXiv:2407.13465v2 (1 Oct 2024), `papers/note_on_hilbert_16th_problem_2407.13465.pdf`:
> *"Even for the quadratic case, it is not known if H(2) < ∞. However, advances has been made and lower bounds for H(n) have been found. For small values of n, the best lower bounds so far are **H(2) ≥ 4** [3,18], H(3) ≥ 13 [10] and H(4) ≥ 28 [15]."*
> Theorem 1: *"Given n ∈ ℕ, it holds H(n+1) ≥ H(n) + 1."* And: H(n) is realizable by structurally stable vector fields with only hyperbolic limit cycles, and is strictly increasing whenever finite.

**Note the date: October 2024, after Pedregal v3 (Aug 2024) and after Gaiko's book (2003). Neither claim is treated as settling H(2).**

---

# PART 5 — Melnikov / Abelian-integral cyclicity of the four quadratic center families

**The four families** (Żołądek's classification; verbatim from Zhao arXiv:1011.2253, in the complex coordinate `z = x + iy`):
```
ż = −iz − z² + 2|z|² + (b + ic)z̄²                        Hamiltonian            (Q₃^H)
ż = −iz + az² + 2|z|² + b z̄²                              reversible             (Q₃^R)
ż = −iz + 4z² + 2|z|² + (b + ic)z̄²,  |b + ic| = 2         codimension four       (Q₄)
ż = −iz + z² + (b + ic)z̄²                                 generalized Lotka–Volterra (Q₃^LV)
ż = −iz + z̄²                                              Hamiltonian triangle
```

The *cyclicity of the period annulus* = the number of zeros of the first non-vanishing Melnikov function `M_k(h)` in the expansion of the displacement function `d(h,ε) = εM₁(h) + ε²M₂(h) + ···`. Melnikov functions for all quadratic centers were determined by **Iliev**, *Perturbations of quadratic centers*, Bull. Sci. Math. **122** (1998) 107–161 (**paywalled**).

## The table

| Family | Best known **upper** bound (cyclicity of open period annulus, quadratic perturbations) | Best known **lower** bound / realized | Conjectured sharp | Open? |
|---|---|---|---|---|
| **Q₃^H** (Hamiltonian, generic) | **2** — Gavrilov, *The infinitesimal 16th Hilbert problem in the quadratic case*, Invent. Math. **143** (2001) 449–497: the number of isolated zeros of any Abelian integral of a real quadratic 1-form over closed level curves of a real cubic Hamiltonian is **at most 2**. Earlier: Horozov–Iliev. | 2 | **2** | **SOLVED (sharp)** |
| **Q₃^LV** (generalized Lotka–Volterra, generic) | **2** — Żołądek, *Quadratic systems with center and their perturbations*, J. Diff. Eqns. **109** (1994) 223–273 | 2 | **2** | Essentially solved |
| **Q₃^R** (reversible, generic) | **OPEN.** Verbatim, Gavrilov–Iliev arXiv:0811.4602: *"Almost nothing is known about the generic reversible case (Q₃^R)."* Partial results only: specific subfamilies give **2** or **3** (Li Chengzhi; Zhao Yulin; Gautier–Mañosas–Villadelprat; Chen–Li–Llibre–Zhang). The subcase **Q₃^R ∩ Q₄** is bounded by **3** ("up to three limit cycles can emerge from the period annulus surrounding the centre"). | 3 (in known subfamilies) | **3** | **OPEN — the main gap** |
| **Q₄** (codimension four, generic) | **5** — Yulin Zhao, arXiv:1011.2253, Theorem 1 (verbatim): *"the perturbed quadratic system (2) has at most **five** limit cycles which emerge from the period annulus around the center."* Improving Gavrilov–Iliev (arXiv:0811.4602, Theorem 1, verbatim): *"The cyclicity of the open period annulus surrounding the center of any generic codimension-four plane quadratic system is less than or equal to **eight**."* | **3** — Zhao, Theorem 1 (verbatim): *"there exists the quadratic polynomials X₂ and Y₂ such that system (2) has at least **three** limit cycles produced by the period annulus of system (1)."* | **3** — Gavrilov–Iliev verbatim: *"The conjectural exact upper bound, as it is well known, is **three** [11,7]."* (Żołądek, Iliev) | **OPEN (gap 3 vs 5)** |
| **Hamiltonian triangle** | solved | — | — | Solved (several authors) |
| **Center itself** (inner boundary, any quadratic) | **3** — **Bautin 1952** | 3 | **3** | **SOLVED (sharp)** |

**Held PDFs for this section:**
- `gavrilov_iliev_quadratic_perturbations_codim_four_centers_0811.4602.pdf` — Gavrilov & Iliev, *Quadratic perturbations of quadratic codimension-four centers*, ≤ 8. Method: orbits are **affine elliptic curves**, `I(h)` is a complete elliptic integral (*"This remarkable fact (it seems to have gone unnoticed by the specialists) is the starting point of the paper"*); Picard–Fuchs equation `M₂ ∘ L₂ ∘ L₁(I) = 0` with `L₁ = h d/dh − 1`; Petrov's argument-principle method; Chebyshev-space property of ker L₂.
- `gavrilov_iliev_quadratic_perturbations_codim_four_centers_1011.2253.pdf` — Zhao, ≤ 5 and ≥ 3. Generating function `I(h) = μ₁hI₀,₀ + μ₂I₁,₀ + μ₃I₀,₁ + μ₄(2I₋₁,₀ + 3κhI₋₁,₁)` with `I_{i,j}(h) = ∬_{H<h} x^i y^j dx dy`, `h ∈ (−2/3, −2/(3√κ))`, `H = (2/3)(κ−1)x³ − (κ−1)x²y + (κ/3)y³ − y`, κ > 1; plus a 5-equation Picard–Fuchs system.
- `gavrilov_iliev_perturbations_quadratic_centers_genus_one_0705.1609.pdf` — the genus-one programme (*"the authors of the paper [9] propose a program for finding the cyclicity of period annuli of quadratic centers of genus one"*).
- `gavrilov_iliev_quadratic_hamiltonian_two_saddle_cycles_1306.2340.pdf` — outer-boundary (two-saddle polycycle) cyclicity for Q₃^H.
- `francoise_gavrilov_xiao_h16_period_annulus_nash_1610.07582.pdf` — Françoise–Gavrilov–Xiao, H16 on a period annulus and the Nash space of arcs.
- `binyamini_novikov_yakovenko_zeros_abelian_integrals_0808.2952.pdf` — the constructive solution of the *infinitesimal* Hilbert 16th problem (explicit double-exponential bound), for context on what "solved" means at the infinitesimal level.
- `yakovenko_tangential_hilbert16_lectures_math0104140.pdf` — lecture notes on the tangential H16.

**Boundary-of-annulus (outer) cyclicity — the newest result:**
**D. Marín, J. Villadelprat, *The cyclicity of hyperbolic hemicycles*, arXiv:2501.16924 (28 Jan 2025), 48 pp. PDF: `papers/marin_villadelprat_2025_2501.16924.pdf`. OPEN ACCESS.** Verbatim abstract:
> *"We consider families of planar polynomial vector fields of degree n and study the cyclicity of a type of unbounded polycycle Γ called hemicycle. Compactified to the Poincaré disc, Γ consists of an affine straight line together with half of the line at infinity and has two singular points, which are hyperbolic saddles located at infinity. We prove four main results. **Theorem A** deals with the cyclicity of Γ when perturbed without breaking the saddle connections. For the other results we consider the case n = 2. More concretely they are addressed to the **quadratic integrable systems belonging to the class Q₃^R** and having two hemicycles, Γ_u and Γ_ℓ, surrounding each one a center. **Theorem B** gives the cyclicity of Γ_u and Γ_ℓ when perturbed inside the whole family of quadratic systems. In **Theorem C** we study the number of limit cycles bifurcating **simultaneously** from Γ_u and Γ_ℓ when perturbed as well inside the whole family of quadratic systems. Finally, in **Theorem D** we show that for three specific cases there exists a **simultaneous alien limit cycle bifurcation** from Γ_u and Γ_ℓ."*

**Why Theorems C and D matter for a five-cycle hunt.** They are about **simultaneous** bifurcation from two limit periodic sets — exactly the mechanism that would be needed to beat (3,1) — and "alien" limit cycles are cycles not detected by the first Melnikov function, i.e. produced by higher-order or non-generic mechanisms. **This is the most promising *theoretical* line in the recent literature for the four-in-one-nest question.** Companion: `marin_villadelprat_criticality_reversible_quadratic_2203.12966.pdf` (*The criticality of reversible quadratic centers at the outer boundary of its period annulus*, arXiv:2203.12966); and `mardesic_marin_villadelprat_perko_conjectures_BT_1303.2065.pdf` (on Perko's conjectures for the Bogdanov–Takens system).

---

# PART 6 — HAS FOUR-IN-ONE-NEST EVER BEEN SEEN, AND WHAT IS THE STRONGEST EVIDENCE FOR AND AGAINST IT EXISTING

## 6.1 The direct answer

**No. Never. Not once, in 47 years and every construction ever published.**

Complete audit of every four-cycle quadratic system in the literature:

| # | Source | Year | Distribution | Sizes | Rigour |
|---|---|---|---|---|---|
| 1 | Chen Lansun & Wang Mingshu | 1979 | **(3,1)** | 3 infinitesimal + 1 normal | Poincaré–Bendixson, existence |
| 2 | Shi Songling | 1980 | **(3,1)** | 3 infinitesimal (down to 10⁻⁷⁵) + 1 normal | Poincaré–Bendixson, 4 trapping regions, existence |
| 3 | Perko (rotated Shi) | 1984 | **(3,1)** | **all 4 normal** | numerical, *assumes* ≤3 in the origin's nest |
| 4–5 | Cherkas–Artés–Llibre, Table 1 rows 7 & 8 | 2003 | **(3,1)** ×2 | all normal | Dulac function — **exact** count |
| 6 | Gaiko, Thm 3.1 | 2006 | **(3,1)** | all normal (field rotation) | constructive, existence |
| 7 | Leonov, criterion | 2010–11 | **(3,1)** (as "3 small+1 large" or "2 small+2 large") | 3 normal + 1 | analytic criterion |
| 8 | Kuznetsov–Kuznetsova–Leonov | 2013 | **(3,1)** | **all 4 "big"** | analytic + visualization |
| 9–12 | Yu & Han, cases (A)(B)(C)(D) | 2010/12 | **(3,1), (1,3), (3,1), (1,3)** | small + large mix | Melnikov, existence |
| 13 | Yu & Zeng, eq. (20) | 2020 | **(3,1)** | **all normal** | Melnikov + simulation |
| 14 | Galias & Tucker | 2022 | **(3,1)** | 3 at radii ~10⁻⁷⁵/10⁻²¹/10⁻⁸ + 1 at ~0.043 | **rigorous, EXACTLY four** |

**Fourteen constructions. Zero (4,0). Zero (4,1). Zero anything but 3-in-a-nest.**

Every single one produces its nest of three the same way: **Bautin's cap of three small-amplitude cycles at a fine focus of order ≤3**, then one more cycle in the *other* nest from a global mechanism. The "normal size" variants (Perko, Cherkas et al., Leonov, Kuznetsov et al., Yu–Zeng) do not change the count — they inflate the same three cycles by rotating the field or choosing parameters far from the degenerate limit. Leonov's IJBC 2011 is the one construction whose *three* are genuinely normal-size from the start (a **first**-order weak focus, three normal-size cycles, then a fourth from a finite disturbance) — and even that stops at (3,1).

Explicit targeted searches for `"four limit cycles surrounding one focus" quadratic`, `"(4,1) distribution"`, `"four limit cycles in one nest"` return **nothing** — every hit resolves to a (3,1) paper.

## 6.2 Evidence AGAINST four-in-one-nest existing

**A1. Bautin's theorem (1952) — the hard local cap.** At most **three** limit cycles can bifurcate from a focus or center of a quadratic system. This is sharp and unassailable. Every known route to three cycles in one nest is exhausted by it. A fourth cycle in the same nest **cannot** come from the local Hopf/fine-focus mechanism; it must come from a distinct global source (separatrix cycle, semistable-cycle bifurcation, Poincaré bifurcation from an annulus), and then the two mechanisms must be made to coexist in the same nest without one destroying the other.

**A2. Forty-seven years of failure by the strongest people in the field, using every available method.** Poincaré–Bendixson trapping regions (Shi, Chen–Wang), rotated vector fields (Perko), Dulac functions with numerical optimization (Cherkas–Artés–Llibre), field-rotation parameters (Gaiko), Liénard asymptotic integration (Leonov), Melnikov functions to first and second order (Yu–Han, Iliev, Gavrilov), and rigorous interval arithmetic (Galias–Tucker). None produced a fourth cycle in one nest. Zegeling 2024's verbatim verdict on the systematic attempts of the 80s and 90s: ***"These attempts remained futile."***

**A3. The conjecture is universal.** Stated as the community expectation in essentially every survey. Verbatim, Yu–Han 2010: *"**It is conjectured that at most 3 limit cycles can exist around one focus point.**"* Verbatim, Artés–Llibre 1997 (`pubmat41.pdf`): *"the main open problem ... it seems that the answer must be the following: A quadratic system has at most 4 limit cycles."* Verbatim, Scholarpedia: *"It is widely conjectured that H(2)=4."* Given Zegeling's reduction `H(2) = one-nest max + 1`, "H(2)=4" **is** "at most 3 in one nest".

**A4. Gaiko claims a proof.** Theorem 4.3 of math/0611142 asserts *"a quadratic system cannot have neither a multiplicity-four limit cycle nor four limit cycles around a singular point (focus)"*, via the Wintner–Perko termination principle plus Bautin. **Discount heavily** — see B4.

**A5. The Melnikov ceilings in every center family are ≤ 3 (conjecturally).** Q₃^H: **2** (Gavrilov 2001, proved sharp). Q₃^LV: **2** (Żołądek). Q₄: conjectured **3**, realized **3**, proved ≤ **5**. Q₃^R: known subfamilies realize **3**. **No center family is even conjectured to produce four cycles from one period annulus.** Since perturbation of a center is the most systematic generator of nested cycles available, this is substantial circumstantial evidence.

**A6. The Cherkas–Artés–Llibre exactness results.** Their eight systems are proved to have **exactly** 3 (or exactly 3+1) — the Dulac function certifies an upper bound, not just existence. In the six one-focus cases with three cycles, the Dulac construction closes at three; the ovals of `Ψ(x,y,C*) = 0` number exactly three.

**A7. The known-configuration list has been closed since 1984.** Perko, verbatim: *"These examples represent all of the known limit cycle configurations for quadratic systems in the plane"* — configurations (a) 1 and (1,0), (b) 2 and (2,0), (c) 3 and (3,0), (d) (1,1), (e) (2,1), (f) (3,1). **Forty-two years later the list is unchanged.**

## 6.3 Evidence FOR four-in-one-nest possibly existing

**B1. No theorem forbids it.** This is the decisive point. Zhang Pingguang's and Zegeling's distribution theorems permit **(n,0)** and **(n,1)** for **every** n. Verbatim Yu–Han: *"the possible cycle distributions ... must be (0,1)-distribution or (1,i)-distribution, i = 0,1,2,3,···. **So far, no results have been obtained for i ≥ 4.**"* The literature's own framing is that i ≥ 4 is **untouched**, not excluded. Zegeling's own abstract frames this as the *remaining* problem: the distribution property *"reduces the study of H(2) to the study of the maximum number of limit cycles surrounding one singularity"* — a reduction he would not present as progress if the one-nest bound were known.

**B2. Gaiko's claimed proof is specifically and expertly rejected on exactly this point.** Roussarie's MathSciNet review MR2023976, verbatim: *"I just mention the **hazardous claim** made in Theorem 4.12 ... it seems that the author unfortunately **confuses two different notions**: paths of limit cycles ... and lines of multiple limit cycles ... In fact, **there is nothing forbidding that a path begin at a parameter value with a multiplicity-four limit cycle and end at a focus point.**"* Roussarie — the author of the standard monograph on this exact subject — is saying, in a peer review, that **nothing known forbids a multiplicity-four limit cycle in a quadratic system.** That is the strongest single piece of expert testimony *for* the possibility. Gaiko's counter (that Bautin's cyclicity result blocks such a path) is unpersuasive: Bautin bounds cycles bifurcating *from the focus*, not the multiplicity of a cycle elsewhere in the nest that later shrinks to it under a *different* parameter path.

**B3. H(2) < ∞ is not even known.** Gasull–Santana, Oct 2024, verbatim: *"Even for the quadratic case, it is not known if H(2) < ∞."* If finiteness itself is open, "at most 3 in one nest" is not close to being settled — it is a conjecture with **no** partial upper bound whatsoever. There is not even a proof that a quadratic system has at most 10 cycles in one nest, or 10⁶.

**B4. Every "H(2)=4" proof has failed.** Petrovskii–Landis 1955 claimed H(2)=3 and was **wrong** (author-acknowledged error) — refuted by an explicit counterexample precisely at the count the community then believed. **This is the exact historical shape of the situation we are in now.** Gaiko 2003: reviewed as "hazardous", unaccepted. Pedregal 2021–24: unrefereed, and its two predecessor versions were **withdrawn for an error in the limit-cycle counting method** — the same step the current argument rests on. da Silva et al. 2024: **explicitly refuted**. The track record of confident H(2)=4 proofs is 0 for 3, plus a 0-for-1 on H(2)=3.

**B5. Extreme scale separation is the actual obstruction — and it is a *detection* obstruction, not an existence obstruction.** Galias–Tucker's four cycles sit at return-map values `4.3×10⁻², 6.7×10⁻⁸, 2.2×10⁻²¹, 7.1×10⁻⁷⁵` — **73 orders of magnitude**, with multipliers `P'` equal to 1 to 90+ decimal places, requiring **2048-bit** arithmetic. Their verbatim framing: *"the limit cycles can reside within areas of vastly different scales. This makes numerical explorations very hard to perform, requiring high precision computations, **where the necessary precision is not known in advance.**"* A fifth cycle at radius ~10⁻²⁰⁰ would be **completely invisible** to every method used before 2022, and to essentially all standard double-precision exploration since. Perko's 1984 remark is the historical proof of this: *"the size of L₂ and L₁ cannot be determined by this method"* and *"A(c) could only be carried down to c ≈ .008 since |A(c)| < 10⁻¹¹."* **Absence of evidence here is unusually weak evidence of absence.**

**B6. Nobody has ever systematically searched.** No published exhaustive sweep of the 5-parameter quadratic space for ≥5 cycles exists (§3.4). Leonov's group mapped only the region they already knew produced (3,1). The complete singularity atlas (Artés–Llibre–Schlomiuk–Vulpe, 1765 configurations) exists but was never coupled to a limit-cycle count. Galias–Tucker's rigorous return-map scan was applied to **one** system. The combination — atlas-guided stratified parameter sweep + rigorous multi-precision return-map scan over 200+ decades of scale — has **never been run**.

**B7. The (3,1) monoculture is a methodological artifact, not an observation.** Every construction reaches three via Bautin's fine-focus mechanism and one via a global mechanism, because that is the only recipe anyone has. It is a recipe that *structurally cannot* yield four in one nest. This is not evidence that four in one nest is impossible; it is evidence that **no one has tried a mechanism that could produce it**. What would be needed is 3 small-amplitude cycles from a 3rd-order fine focus **plus** one normal-size cycle in the *same* nest from a separatrix/semistable bifurcation. Leonov's *"Bifurcation of appearance of two limit cycles via semistable cycle"* (IJBC 2011) is the closest anyone has come to the required second mechanism, and Marín–Villadelprat's **simultaneous** and **alien** limit cycle bifurcations (arXiv:2501.16924, Theorems C and D) are the closest current theory.

**B8. Higher-order Melnikov analysis is essentially untouched for Q₃^R.** Verbatim, Gavrilov–Iliev: *"Almost nothing is known about the generic reversible case (Q₃^R)."* Yu–Han, verbatim: *"a different method was used in [27] with Melnikov function up to second order, but no more limit cycles were found"* — one attempt, on one subcase. The Q₄ gap (proved ≤ 5, realized 3) is likewise a genuine three-unit gap in which nobody has ruled out 4.

## 6.4 Net assessment

**The conjecture "≤3 in one nest" is very likely true, but it rests on nothing but 47 years of failure to find a counterexample using a family of methods that structurally could not have found one.** Its only claimed proofs are (i) reviewed as confusing two distinct notions by the leading expert in the area, or (ii) unrefereed with a history of withdrawal for an error at the critical step, or (iii) explicitly refuted. The distribution theorems — the one piece of genuinely solid structural knowledge, now double-verified by Zegeling — are carefully **neutral** on i ≥ 4.

**Highest-value directions for a refutation attempt, in order:**

1. **Q₃^R, higher-order Melnikov, simultaneous/alien bifurcation.** The one center family where the cyclicity is genuinely unknown, combined with the one recent theoretical tool (Marín–Villadelprat 2025, Theorems C & D) designed for simultaneous multi-source bifurcation. Target: 3 from the annulus + 1 alien, in the same nest.
2. **Q₄, closing the 3-vs-5 gap from below.** Gavrilov–Iliev proved ≤ 8, Zhao improved to ≤ 5, but only **3** are realized and the conjectured sharp value is 3. If the true value is 4, that is 4 cycles from **one** period annulus around **one** center — a (4,0) or, with a companion cycle, a (4,1). Nobody has seriously attacked the lower bound here.
3. **The atlas-guided rigorous sweep.** Stratify by Artés–Llibre–Schlomiuk–Vulpe's 1765 singularity configurations; on each stratum, run Galias–Tucker's exact toolchain (CAPD + MPFR, `f(y) = y − P(y)` sign-change scan over `y = 10^(−k/3)` for k up to 600+, then Lemmas 4/5/6 to exclude and Lemma 1 / interval-Newton to certify). This is the search that has never been run and is entirely within reach of existing software.
4. **Combine two mechanisms in one nest deliberately.** Take a 3rd-order fine focus (three small cycles guaranteed by Bautin, e.g. the Artés–Llibre 1997 / Llibre–Schlomiuk 2004 classification of exactly this stratum) and *simultaneously* engineer a semistable-cycle or separatrix bifurcation in the **same** nest, à la Leonov's IJBC 2011 semistable mechanism. Every existing construction spends the global mechanism on the *other* focus. Nobody has spent it on the same one.
5. **Search near the boundaries Ilyashenko–Llibre's δ-tameness excludes** — near singular points, near infinity, near the center variety. That is where finiteness is least controlled and where a 10⁻²⁰⁰-scale fifth cycle would hide.

---

# Appendix — PDF inventory (`/Users/scottg/Claude_all/papers/`)

**Four-cycle constructions**
- `galias_tucker_2022_songling_exactly_four_limit_cycles_AMC415.pdf` — Galias & Tucker, AMC 415 (2022) 126691. Rigorous, exactly four. OPEN.
- `perko1984_RMJM_limit_cycles_quadratic.pdf` — Perko, RMJM 14(3) (1984) 619–645. Normal-size four; configuration atlas (a)–(f).
- `cherkas2003.pdf` / `cherkas_artes_llibre_2003_normal_size_BASM41.pdf` — Cherkas, Artés, Llibre, BASM 1(41) (2003) 31–46. Table 1 with 8 explicit coefficient vectors. OPEN.
- `yu_han_2012_four_limit_cycles_perturbing_quadratic_integrable_1002.1055.pdf` — arXiv:1002.1055. Near-integrable (3,1), four cases. OPEN.
- `yu_zeng_2021_visualization_four_limit_cycles_near_integrable_2002.09987.pdf` (dup: `kuznetsov_visualization_four_lc_near_integrable_2002.09987.pdf`) — arXiv:2002.09987. Reproduces Shi, Chen–Wang, Kuznetsov et al. explicitly. OPEN. **Best single source for coefficient vectors.**
- `gaiko_geometry_planar_quadratic_math0611142.pdf` — arXiv:math/0611142. Thm 3.1 construction + Thm 3.2/4.3 claims + Roussarie quote. OPEN.
- `pubmat41.pdf` — Artés & Llibre, Publ. Mat. 41 (1997) 7–39, third-order weak focus. OPEN.
- `llibre_schlomiuk_2004_CJM56_weak_focus_third_order.pdf` — Canad. J. Math. 56 (2004).
- `koditschek_narendra_1984_limit_cycles_planar_quadratic_JDE54.pdf` — JDE 54 (1984) 181–195; documents the Petrovskii–Landis retraction.

**Claims of ≥5 / claimed H(2) proofs and refutations**
- `pedregal_2021_H2_equals_4_claim_2103.07193.pdf` — arXiv:2103.07193v3, 69 pp. OPEN. Unrefereed.
- `note_on_recent_attempt_hilbert16_rebuttal_2411.09594.pdf` — Buzzi & Novaes, arXiv:2411.09594. Refutes da Silva et al. OPEN.
- `note_on_hilbert_16th_problem_2407.13465.pdf` — Gasull & Santana, arXiv:2407.13465. "H(2) < ∞ not known". OPEN.
- `gaiko_hilbert16_quadratic_lienard_transformation_math0701869.pdf`, `gaiko_lienard_smale_math0611143.pdf`, `gaiko_quadratic_two_parallel_isoclines_0803.3055.pdf`, `gaiko_1104.3019.pdf`, `gaiko_general_lienard_1202.3540.pdf`, `gaiko_1504.03353.pdf`, `gaiko_1611.08113.pdf`, `gaiko_field_rotation_ihes.pdf`.
- `chebyshev_pullbacks_quadratic_ceiling_2604.12883.pdf` — arXiv:2604.12883 (Apr 2026). Replication lower bounds; H(2) unaffected. OPEN.
- **NOT OBTAINED:** Voronin & Lebedev 2022 six-cycle (3:3) claim — host `repository.hneu.edu.ua` unreachable. **NOT FOUND / NON-EXISTENT:** Malyarets 2026; Hernandez Rosales 2026 "H(2)=7"; Shi 1978 five-cycle preprint.

**Certification / rigorous numerics**
- `villanueva_tucker_2026_2602.22558.pdf` — arXiv:2602.22558. **Center conditions, not cycle certification.** OPEN.
- `proving_existence_numerically_detected_planar_limit_cycles_1602.00113.pdf` — Gasull, Giacomini, Grau, arXiv:1602.00113. OPEN.
- `gasull_transversal_conics_existence_limit_cycles_1410.4480.pdf`, `garcia_saldana_gasull_giacomini_new_approach_limit_cycles_1910.08098.pdf`, `gasull_extended_bendixson_dulac_1305.3402.pdf`, `number_limit_cycles_planar_invariant_algebraic_curves_2210.15803.pdf`, `smooth_transformations_ruling_out_closed_orbits_2309.02513.pdf`.

**Distribution / upper bounds**
- `zegeling2024.pdf` — Zegeling, Adv. Nonlinear Anal. 13 (2024) 20240012. Thm 1.2: only (n,0)/(n,1). **OPEN (CC-BY). Critical.**
- `zhang_pingguang_2002_QTDS_abstract.pdf` — Zhang, QTDS 3 (2002) 437–463. Thm 1.
- `zhang_pingguang_cai_suilin_1991_quadratic_systems_weak_focus_BAMS44.pdf` — BAMS 44 (1991).
- `ilyashenko_llibre_restricted_version_H2_0910.3443.pdf` — arXiv:0910.3443, Thm 5. OPEN.
- `schlomiuk_vulpe_geometry_neighbourhood_infinity_math0405026.pdf`, `msp_pjm2006_existence_limit_cycles_real_quadratic.pdf`, `nlin0502049_xiao_zhang_uniqueness_theorem_quoted.pdf`, `sansone_legacy_lienard_uniqueness_1101.2761.pdf`.

**Melnikov / Abelian integrals**
- `gavrilov_iliev_quadratic_perturbations_codim_four_centers_0811.4602.pdf` — Q₄ ≤ 8. OPEN.
- `gavrilov_iliev_quadratic_perturbations_codim_four_centers_1011.2253.pdf` — Zhao, Q₄ ≤ 5, ≥ 3. OPEN.
- `gavrilov_iliev_perturbations_quadratic_centers_genus_one_0705.1609.pdf`, `gavrilov_iliev_quadratic_hamiltonian_two_saddle_cycles_1306.2340.pdf`, `francoise_gavrilov_xiao_h16_period_annulus_nash_1610.07582.pdf`.
- `marin_villadelprat_2025_2501.16924.pdf` — arXiv:2501.16924, hemicycle cyclicity in Q₃^R, simultaneous + alien bifurcation. OPEN. **Most promising recent theory.**
- `marin_villadelprat_criticality_reversible_quadratic_2203.12966.pdf`, `mardesic_marin_villadelprat_perko_conjectures_BT_1303.2065.pdf`.
- `binyamini_novikov_yakovenko_zeros_abelian_integrals_0808.2952.pdf`, `yakovenko_tangential_hilbert16_lectures_math0104140.pdf`.
- `yu_han_eight_limit_cycles_around_center_quadratic_IJBC2013.pdf` — cubic perturbation of quadratic Hamiltonian (not H(2)).

**Paywalled, not obtained:** Shi Songling Sci. Sinica 23 (1980); Chen & Wang Acta Math. Sinica 22 (1979); Leonov JAMM 74 (2010), Doklady Math. 81 (2010), IJBC 21 (2011); Kuznetsov–Kuznetsova–Leonov DEDS 21 (2013); Differential Equations 49 (2013); Gaiko book (Kluwer 2003); Artés–Llibre–Schlomiuk–Vulpe book (Birkhäuser 2021); Gavrilov Invent. Math. 143 (2001); Iliev Bull. Sci. Math. 122 (1998); Żołądek JDE 109 (1994); Zegeling & Kooij JDE 151 (1999); Cherkas Differ. Equ. 46 (2010); Gasull–Giacomini JDE 305 (2021).
