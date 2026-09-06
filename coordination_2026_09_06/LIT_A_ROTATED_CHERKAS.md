# LIT_A — Rotated vector fields, Andronov–Hopf functions, multiplicity, Dulac–Cherkas
### Literature slice A for the H(2) ≥ 5 hunt. Compiled 2026-09-06.

Scope: theory of rotated vector fields (Duff, Perko), the Wintner–Perko termination
principle, multiplicity of limit cycles, and Dulac–Cherkas methods — **specifically as
they bear on how many limit cycles can live in ONE nest of a quadratic system.**

All PDFs saved under `/Users/scottg/Claude_all/papers/`.

---

## 0. Executive summary of what this slice establishes

1. **No accepted proof exists that a quadratic system has at most 3 limit cycles in one
   nest.** Every claimed proof (Gaiko 2003/2006/2008) rests on a non-rigorous
   "the spirals untwist so the distance between coils increases" argument, or on a
   misuse of the Wintner–Perko termination principle that Roussarie explicitly
   rejected in MathSciNet review MR2023976 (2005d:37102).
2. **The number 3 in the literature is Bautin's LOCAL number (cyclicity of a focus),
   not a global bound on a nest.** Perko himself (1984, closing remarks) lists
   "Determine whether or not a quadratic system can have more than three limit cycles
   around a single critical point" as an *open problem*, and says only that
   "there is no known example."
3. The rotated-field machinery gives you a **monotonicity + termination** toolkit that
   is exactly what a 4-in-a-nest search needs: a field-rotation parameter turns the
   search for a 4th cycle into a search for an extra **semistable (multiplicity-2)**
   cycle appearing in a specific annulus, and the Andronov–Hopf function
   `α = φ(x₁)` must acquire **three extrema** instead of two.
4. Dulac–Cherkas functions give **upper** bounds on a nest, but the bound they give is
   an artifact of the ansatz degree (`n/2` for a degree-`n`-in-`y` ansatz), not an
   intrinsic bound. **No Dulac–Cherkas theorem forbids 4 cycles in a nest.**
5. Explicit coefficient vectors for quadratic fields with **3 normal-size cycles in one
   nest** are given in §2.5 (Perko 1984) and §3.2 (Cherkas–Artés–Llibre 2003, 8 systems).
   These are the natural launch points for a 4-in-a-nest search.

---

# PART 1 — ROTATED VECTOR FIELDS AND PERKO'S THEORY

## 1.1 Duff 1953

**Citation.** G. F. D. Duff, *Limit-cycles and rotated vector fields*, Annals of
Mathematics (2) **57** (1953), 15–31.
(Warning: Perko's own reference lists cite this as vol. **67**; Gaiko cites
"Ann. Math. 67 (1953), 15–31". The correct volume is 57, issue Jan. 1953,
JSTOR stable/i307260.)

**PDF: PAYWALLED (JSTOR).** Not obtainable without a subscription.

Duff's content is nevertheless recoverable **verbatim through Perko's quotations**,
which are reproduced below.

### Definition (family of rotated vector fields) — Perko 1984, p. 621, verbatim:

> "A vector field is said to belong to a family of rotated vector fields with
> parameter α ∈ (a, b) if each vector in the vector field rotates through a positive
> angle as α increases in (a, b). Any vector field, such as the vector field
> ẋ = P(x,y), ẏ = Q(x,y) defined by (1), can be embedded in a family of rotated
> vector fields with parameter α ∈ (−π, π) by setting
>
>   ẋ = P(x, y) cos α − Q(x, y) sin α
>   ẏ = P(x, y) sin α + Q(x, y) cos α."

The *infinitesimal* form of the same condition (Perko 1992, §2, verbatim):

> "Since, by definition, any one-parameter family of rotated vector fields f(x, λ)
> satisfies
>
>   f ∧ f_λ(x, λ) > 0"

where `x ∧ y = x₁y₂ − x₂y₁`.

Cherkas–Artés–Llibre 2003 use the weaker "rotating parameter" version (verbatim):

> "We say that the parameter α rotates the vector field f(x,α) associated to system
> (3), or that α is a rotating parameter, if one of the two inequalities
>
>   (f₁)′_α f₂ − f₁(f₂)′_α ≥ 0 (≤ 0),  x ∈ R², α ∈ R,
>
> holds, and the inequality never becomes an identity equal to zero on any limit cycle
> of L(α). ... Moreover, if α is a rotating parameter, then L(α₁) ∩ L(α₂) = ∅ if
> α₁ ≠ α₂."

Gaiko's determinant form (Gaiko, *Geometry of Planar Quadratic Systems*, Lemma 2.1
proof, verbatim):

> "Using the definition of a field rotation parameter we can calculate the following
> determinants: Δ_λ = P Q′_λ − Q P′_λ ... Since, by definition, the vector field is
> rotated in positive direction (counterclockwise) when the determinant is positive
> and in negative direction (clockwise) when the determinant is negative."

### Duff's Table I (p. 21) — reproduced verbatim in Perko 1984, p. 622:

| Orientation | Stability | Motion as α ↑ |
|---|---|---|
| + | − | Contracts |
| + | + | Expands |
| − | − | Expands |
| − | + | Contracts |

### Duff's monotonicity / generation theorem — Perko 1984 p. 622, verbatim:

> "In his theory for limit cycles of a family of rotated vector fields, [5], Duff
> showed that a unique limit cycle is generated at the origin of (1_α), with the
> determinant of the linear terms a₁₀b₀₁ − b₁₀a₀₁ ≠ 0, at that value of α where the
> trace of the linear part, τ(α) = (a₁₀ + b₀₁)cos α + (a₀₁ − b₁₀)sin α, vanishes
> provided that the origin is not a center at this value of α. Furthermore, he showed
> that this limit cycle **expands monotonically with monotonically varying α, covering
> a deleted neighborhood of the origin until it either (i) intersects one or more
> critical points of (1_α) and forms a separatrix cycle, or (ii) intersects a second
> limit cycle of (1_α) and forms a semi-stable limit cycle, or (iii) becomes
> unbounded.**"

### Duff's Theorem 8 (splitting of a semistable cycle) — as used by Perko 1984 p. 624:

> "Since a semistable limit cycle can be made to split into two limit cycles with a
> proper variation of α, cf. Theorem 8, p. 23 in [5], this is equivalent to assuming
> that there are at most three limit cycles around the critical point (1,0)."

**This is the load-bearing logical step for the whole "3-in-a-nest" folklore, and it is
an ASSUMPTION, not a theorem.** See §5.

Duff's other numbered results referenced by Perko: **Theorem 7** (monotone motion of a
simple limit cycle with α), **Theorem 9** (a limit cycle expands to infinity /
terminates), **Theorem 10** (generation of a unique limit cycle at a weak focus),
**eq. (3.17), p. 23** (Duff's estimate for the growth rate of a limit cycle with α;
Perko 1992 reproduces its analogue as `∂s/∂λ(0,τ) = −d_λ(0,0,τ)/d_s(0,0,τ)`).

---

## 1.2 Perko 1975 (JDE)

**Citation.** L. M. Perko, *Rotated vector fields and the global behavior of limit
cycles for a class of quadratic systems in the plane*, J. Differential Equations
**18** (1975), 63–86. (Crossref-confirmed. NOTE: the task brief called this
"...for a class of nonlinear systems in the plane"; the actual title says
**quadratic** systems.)

**PDF: PAYWALLED (Elsevier; ScienceDirect bot-blocks all requests).**

Its content survives in Perko 1984's citations. The results Perko labels
**Theorems D, G, H** in [7] = Perko 1975:

- **Theorem D** — generalizes Duff's Theorem 7: monotone expansion/contraction of a
  simple limit cycle under a field rotation parameter.
- **Corollary to Theorem G** (Perko 1984, p. 627, verbatim):
  > "as the parameter α varies monotonically, any limit cycle L(α) of (4) covers an
  > annular neighborhood of its initial position and the inner and outer boundaries of
  > this region consist of either a single critical point, a separatrix cycle on the
  > Poincaré Sphere, or a semistable limit cycle."
- **Corollary to Theorem H** (Perko 1984, p. 628 & p. 637, verbatim):
  > "a unique limit cycle, L₃, is generated at the critical point (−2,0) at
  > α = −Tan⁻¹(1/3)" ... "it follows from the corollary to Theorem H in [7] that the
  > limit cycle L₁ is generated at the origin at that value of α for which the trace of
  > the linear part of (6) is zero; i.e., at α = tan⁻¹(λ/1)."

Companion paper: L. M. Perko and Shu Shih-lung, *Existence, uniqueness, and
non-existence of limit cycles for a class of quadratic systems in the plane*,
J. Differential Equations **53** (1984), 1–26. (Paywalled.)

---

## 1.3 Perko 1990 (TAMS) — the planar termination principle

**Citation.** L. M. Perko, *Global families of limit cycles of planar analytic
systems*, Trans. Amer. Math. Soc. **322** (1990), 627–656.

**PDF: `/Users/scottg/Claude_all/papers/perko1990_TAMS_global_families.pdf`**
(AMS open archive; **only the first 10 pages, pp. 627–636, are served** — the
appendix containing the explicit `d_λ` formula is not in the free portion, but the
same formula appears complete in Perko 1992, §1.5 below.)

### Abstract, verbatim:

> "The global behavior of any one-parameter family of limit cycles of a planar analytic
> system ẋ = f(x, λ) depending on a parameter λ ∈ R is determined. It is shown that
> any one-parameter family of limit cycles belongs to a maximal one-parameter family
> which is either open or cyclic. If the family is open, then it terminates as the
> parameter or the orbits become unbounded, or it terminates at a critical point or on
> a (compound) separatrix cycle of the system. This implies that the periods in a
> one-parameter family of limit cycles can become unbounded only if the orbits become
> unbounded or if they approach a degenerate critical point or (compound) separatrix
> cycle of the system. This is a more specific result for planar analytic systems than
> Wintner's principle of natural termination for n-dimensional systems where the
> periods can become unbounded in strange ways. This work generalizes Duff's results
> for one-parameter families of limit cycles generated by a one-parameter family of
> rotated vector fields. In particular, it is shown that the behavior at a nonsingular,
> multiple limit cycle of any one-parameter family of limit cycles is exactly the same
> as the behavior at a multiple limit cycle of a one-parameter family of limit cycles
> generated by a one-parameter family of rotated vector fields."

### **Planar termination principle** (Theorem 3 of §2) — verbatim, p. 629:

> "**Planar termination principle.** Any maximal one-parameter family of limit cycles of
> a planar analytic system (1_λ) is either open or cyclic. If it is open, then it
> terminates as either the parameter or the orbits become unbounded; or the family
> terminates either at a critical point or on a (compound) separatrix cycle of (1_λ)."

Contrast, verbatim p. 630:

> "**Wintner's principle of natural termination.** Any maximal one-parameter family of
> periodic orbits of analytic system (1_λ) with x ∈ Rⁿ is cyclic or terminates as
> either the periods, the parameter, or the orbits become unbounded, or the family
> terminates at an equilibrium point or at a period-doubling bifurcation orbit of
> (1_λ)."

### **Theorem (Behavior at a Nonsingular Multiple Limit Cycle)** — verbatim, p. 630:

> "**Theorem (Behavior at a Nonsingular Multiple Limit Cycle).** If L₀ is a nonsingular,
> multiple limit cycle of (1_{λ₀}), then L₀ belongs to a unique one-parameter family of
> limit cycles of (1_λ); furthermore,
> (1) if the multiplicity of L₀ is odd, then the family either expands or contracts
> monotonically as λ increases through λ₀, and
> (2) if the multiplicity of L₀ is even, then L₀ bifurcates into a stable and an
> unstable limit cycle as λ varies from λ₀ in one sense and L₀ disappears as λ varies
> from λ₀ in the opposite sense; i.e., there is a saddle-node bifurcation at λ₀."

### Multiplicity apparatus — verbatim, pp. 635–636:

> "**Definition 1.1.** The function d(n, λ) = h(n, λ) − n is called the displacement
> function along the transversal l."

> "**Lemma 1.3.** Let d(n, λ) be the displacement function along a straight line normal
> to the cycle L₀. Then
> (3)  d_n(0,0) = exp ∫₀^{T₀} [P_x(φ₀(t), ψ₀(t), 0) + Q_y(φ₀(t), ψ₀(t), 0)] dt − 1."

> "**Definition 1.2.** A cycle L₀ of (1₀) is called a simple limit cycle if
> ∫₀^{T₀}[P_x(φ₀(t), ψ₀(t), 0) + Q_y(φ₀(t), ψ₀(t), 0)] dt ≠ 0. If L₀ is a limit cycle and
> the above condition is not satisfied, then L₀ is called a multiple limit cycle."

> "**Definition 1.3.** If L₀ is a (limit) cycle of (1₀) and if
> d_n(0,0) = ··· = d_n^{(k−1)}(0,0) = 0, d_n^{(k)}(0,0) ≠ 0 then L₀ is called a limit
> cycle of multiplicity k."

> "**Remark 1.3.** Theorem 42 on p. 277 in [6] states that if L₀ is a limit cycle of
> multiplicity k ≥ 2 of the analytic system (1₀), then there exists an analytic
> perturbation of (1₀) which has exactly k limit cycles in a neighborhood of L₀."

**Consequence directly relevant to the hunt:** a **multiplicity-4** limit cycle of a
quadratic system, if one exists, splits into **4 limit cycles in one nest** under
perturbation — but Remark 1.3 only guarantees an *analytic* perturbation, not a
*quadratic* one. That gap is exactly where the H(2) question lives.

---

## 1.4 Perko 1992 (Proc. AMS) — geometric theory, and the `d_λ` formula

**Citation.** L. M. Perko, *Bifurcation of limit cycles: geometric theory*,
Proc. Amer. Math. Soc. **114** (1992), 225–236.

**PDF: `/Users/scottg/Claude_all/papers/perko1992_ProcAMS_bifurcation_geometric_theory.pdf`**
(AMS open archive; **complete, all 12 pages.**)

### THE FORMULA FOR THE PARAMETER DERIVATIVE OF THE DISPLACEMENT FUNCTION

**Lemma 2** — verbatim, p. 227 (attributed to Andronov et al. [1, p. 384], and
highlighted by Chicone & Jacobs):

> "**Lemma 2.** For δ > 0, |s| < δ, |λ| < δ, and τ ∈ R, let d(s, λ, τ) denote the
> displacement function for the system (1_λ) along the normal line l_τ and let ω₀
> denote the orientation of Γ₀. Then
>
> (3)  d_λ(0,0,0) = −(ω₀ / |f₀(0)|) ∫₀^{T₀} e^{−∫₀^t ∇·f₀(u) du} f ∧ f_λ(γ₀(t), 0) dt."

**Lemma 3** (the τ-dependent version) — verbatim, p. 227:

> "**Lemma 3.** Under the hypotheses of Lemma 2,
> (4)  d_τ(0,0,τ) = −(ω₀ e^{∫₀^τ ∇·f₀(t)dt} / |f₀(τ)|) ∫_τ^{τ+T₀}
>      e^{−∫₀^t ∇·f₀(u)du} f ∧ f_λ(γ₀(t), 0) dt."
> *(the OCR of the printed paper is degraded at this equation; the structure is
> λ-derivative along the normal at the point γ₀(τ), with lower limit τ.)*

**Melnikov connection** — verbatim, Remark 2, p. 228:

> "The integral containing the wedge product in equation (4) is related to the Melnikov
> function, which plays such an important role in the theory of perturbed dynamical
> systems. We define the function
> M(τ) = ∫_τ^{T₀} e^{−∫₀^t ∇·f₀(u)du} f ∧ f_λ(γ₀(t), 0) dt for later use."

### Multiplicity and singular/nonsingular multiple cycles — verbatim:

> "**Definition 1.** The limit cycle Γ₀ is a multiple limit cycle of (1₀) if
> ∫₀^{T₀} ∇·f₀(t) dt = 0."

> "**Lemma 1.** For δ > 0, |s| < δ, and |λ| < δ, let d(s, λ) denote the displacement
> function for the system (1_λ). Then the derivative d_s(0,0) is independent of the
> point γ₀(τ) on Γ₀ and (2) d_s(0,0) = e^{∫₀^{T₀} ∇·f₀(t) dt} − 1."

> "**Definition 2.** If d(0,0) = d_s(0,0) = d_s^{(2)}(0,0) = ··· = d_s^{(m−1)}(0,0) = 0
> and d_s^{(m)}(0,0) ≠ 0 then, for m > 1, Γ₀ is called a multiple limit cycle of
> multiplicity m. If m = 1, i.e., if d(0,0) = 0 and d_s(0,0) ≠ 0, then Γ₀ is called a
> simple limit cycle or hyperbolic limit cycle."

> "**Remark 1.** It follows from Theorem 42 [1, p. 277] that the multiplicity m of a
> multiple limit cycle Γ₀ is independent of the point γ₀(τ) on Γ₀; in fact, according to
> Theorem 42 in [1], **the multiplicity m of Γ₀ is equal to the maximum number of limit
> cycles that can bifurcate from Γ₀ under a perturbation of (1₀)**."

> "**Definition 3.** The limit cycle Γ₀ is a singular, multiple limit cycle of (1₀) if
> d_s(0,0) = 0 and d_λ(0,0) = 0. If d_s(0,0) = 0 and d_λ(0,0) ≠ 0, then Γ₀ is a
> nonsingular, multiple limit cycle of (1₀)."

> "**Remark 3.** Note that Γ₀ is a singular, multiple limit cycle of (1₀) if and only if
> ∫₀^{T₀} ∇·f₀(t) dt = 0 and ∫₀^{T₀} e^{−∫₀^t ∇·f₀(u)du} f ∧ f_λ(γ₀(t),0) dt = 0."

> "**Definition 4.** ρ(τ) = sgn M(τ) = sgn ∫_τ^{τ+T₀} e^{−∫₀^t ∇·f₀(u)du}
> f ∧ f_λ(γ₀(t),0) dt."

### **THE MONOTONICITY THEOREM** — verbatim, Theorem 3, p. 229:

> "**Theorem 3.** If Γ₀ is a nonsingular, multiple limit cycle of (1₀), then Γ₀ belongs
> to a unique, one-parameter family of limit cycles of (1_λ) and
> (1) if the multiplicity of Γ₀ is odd, then the family either expands or contracts
> monotonically as λ increases through zero as determined by Table 1 while
> (2) if the multiplicity of Γ₀ is even, then Γ₀ bifurcates into a simple stable limit
> cycle and a simple unstable limit cycle as λ varies in a certain sense, determined by
> Table 1, and Γ₀ disappears as λ varies in the opposite sense."

**Table 1** (verbatim caption): "The change in λ, Δλ, that causes the expansion of a
nonsingular, multiple limit cycle Γ₀ of odd multiplicity or the bifurcation of a
nonsingular, multiple limit cycle Γ₀ of even multiplicity." The table is indexed by
`ω₀` (orientation), `σ₀` (±1 as Γ₀ is unstable/stable on its exterior) and `ρ₀ = ρ(0)`;
the operative sign is `sgn[b₀(τ)] = ω₀σ₀ρ₀` (Remark 4).

> "**Theorem 2.** If f(x,λ) defines a one-parameter family of rotated vector fields in
> R² with parameter λ ∈ R, then **any multiple limit cycle of (1_λ) is a nonsingular,
> multiple limit cycle** of (1_λ)."

> "**Theorem 4.** If Γ₀ is a simple limit cycle of (1₀), then Γ₀ belongs to a unique,
> one-parameter family of limit cycles Γ_λ of (1_λ) and at any point γ₀(τ) on Γ₀,
> increasing the parameter λ causes the limit cycle Γ_λ to expand or contract along the
> normal line l_τ if and only if ω₀σ₀ρ(τ) = ±1 respectively."

### Singular multiple limit cycles (Puiseux series) — verbatim, Theorem 5, p. 232:

> "**Theorem 5.** Suppose that Γ₀ is a singular, multiple limit cycle of (1₀) that
> belongs to a one-parameter family of limit cycles Γ_λ of (1_λ), corresponding to a
> branch s(λ, τ) of d(s, λ, τ) = 0, of reduced multiplicity m. Then for all but possibly
> a finite number of τ ∈ [0, T₀), there is a δ > 0 such that s(λ, τ) can be expanded in a
> Puiseux series, (5) with a₀(τ) ≠ 0, which converges for 0 < σλ < σδ where σ = ±1;
> furthermore,
> (1) if m is even, then Γ₀ bifurcates into a simple stable limit cycle and a simple
> unstable limit cycle belonging to the family Γ_λ as σλ increases and Γ₀ disappears as
> σλ decreases;
> (2) if m is odd and κ is odd, then the limit cycles of the family Γ_λ expand or
> contract along the normal line l_τ to Γ₀ as σλ increases according to whether a₀(τ) is
> positive or negative respectively; and
> (3) if m is odd and κ is even, then the limit cycles of the family Γ_λ expand or
> contract along the normal line l_τ to Γ₀ as λ increases in (0, δ) or as λ decreases in
> (−δ, 0) according to whether a₀(τ) is positive or negative respectively."

### Conclusions, verbatim, p. 235:

> "To summarize, any bifurcation at a periodic orbit of a planar analytic system (1_λ)
> depending on a parameter λ ∈ R occurs at a multiple limit cycle of (1_λ) or at a cycle
> belonging to a continuous band of cycles of (1_λ). **The only bifurcation that occurs
> at a nonsingular, multiple limit cycle is the saddle-node bifurcation** with the stable
> and unstable bifurcating limit cycles expanding and contracting monotonically. On the
> other hand, **as many as m one-parameter families of limit cycles can bifurcate from a
> singular, multiple limit cycle of multiplicity m**; however, any bifurcating
> one-parameter family whose reduced multiplicity is even corresponds to a saddle-node
> type of bifurcation..."

**Search-relevant corollary.** Under a *field rotation* parameter (Theorem 2), every
multiple cycle is nonsingular, so **only fold (saddle-node) bifurcations occur**. So in
a rotated family the number of cycles in a nest can only change by ±2 at a semistable
cycle, or by ±1 at a Hopf/graphic. That is the exact combinatorial constraint governing
a 4-in-a-nest search.

---

## 1.5 Perko 1984 (Rocky Mountain J. Math.) — the quadratic paper

**Citation.** L. M. Perko, *Limit cycles of quadratic systems in the plane*,
Rocky Mountain J. Math. **14** (1984), no. 3, 619–645.

**PDF: `/Users/scottg/Claude_all/papers/perko1984_RMJM_limit_cycles_quadratic.pdf`**
(Project Euclid, **free, complete**.)

### Perko's own statements about the maximum number of cycles in a nest

**(a)** p. 624, verbatim (the key epistemic statement):

> "Since a semistable limit cycle can be made to split into two limit cycles with a
> proper variation of α, cf. Theorem 8, p. 23 in [5], this is equivalent to assuming
> that there are at most three limit cycles around the critical point (1, 0).
> **This seems like a reasonable assumption in view of Bautin's result [2] that at most
> three limit cycles can disappear into a critical point of focus or center type for a
> quadratic system and since there is no known example of a quadratic system with more
> than three limit cycles around a single critical point.**"

**(b)** §7 Closing Remarks, p. 644, verbatim:

> "Two important mathematical problems are suggested by this study.
> 1. Determine the exact number of limit cycles in the examples of Songling [10] and
> Tung Chin-chu [11], a problem previously suggested by Chicone and Jinghuang [3].
> **2. Determine whether or not a quadratic system can have more than three limit cycles
> around a single critical point.**
> The second problem would have a direct bearing on Hilbert's problem 16 for quadratic
> systems. **In fact, if a quadratic system can have at most three limit cycles around
> any one critical point, it would then follow from Theorem 2.9(c) in [3] that the
> maximum number of limit cycles possible for a quadratic system is less than or equal
> to six.** It follows from Songling's example that this number is greater than or equal
> to four."

> **Reading this carefully: Perko does NOT assert H(2)=4. He asserts that "≤3 per nest"
> would only give H(2) ≤ 6.** The gap between 4 and 6 is closed only later by
> Zhang Pingguang's result (a two-nest configuration must have min(n₁,n₂)=1); with that,
> "≤3 per nest" ⟹ H(2) ≤ 4.

**(c)** Perko's theorem for Tung Chin-chu's (2:1) example, verbatim, pp. 625–626 —
note the *hypothesis*:

> "**THEOREM.** The quadratic system (4) forms a complete family of rotated vector fields
> with parameter α ∈ [−π, π). A unique, negatively-oriented limit cycle, L₁, is
> generated at the critical point (1, 0) at α = 0. **Under the assumption that there are
> at most three limit cycles around the critical point (1, 0)**, it follows that there
> exists a unique α* ∈ (−π/2, 0) such that for all α ∈ (α*, 0), there is a second,
> negatively-oriented limit cycle, L₂, around (1, 0). The limit cycle L₁ is unstable and
> expands monotonically as α decreases in the interval (α*, 0), and the limit cycle L₂
> is stable and contracts monotonically as α decreases in the interval (α*, 0). The two
> limit cycles intersect in a negatively-oriented semistable limit cycle, stable on its
> exterior and unstable on its interior, at α = α*. This semistable limit cycle
> disappears as α is decreased below the critical value α*."

α* ≈ **−0.0093** (numerically determined, Fig. 11).

### EXPLICIT COEFFICIENT VECTORS FROM PERKO 1984

**(P1) Tung Chin-chu's system: 3 limit cycles, distribution (2 : 1)** — Perko's eq. (4):

```
ẋ = P(x,y) cos α − Q(x,y) sin α
ẏ = P(x,y) sin α + Q(x,y) cos α
P(x,y) = x y
Q(x,y) = −(1/3)(x − 1)(x + 2) + (1/2)y² + (1/3)x y − (1/3) y
```
Critical points (1,0) and (−2,0). For α ∈ (α*, 0) with α* ≈ −0.0093 there are **two
cycles around (1,0)** and **one around (−2,0)**. Perko's plots use **α = −0.009**.
(This is a *complete* family of rotated vector fields with vectors of constant length.)

**(P2) Shi Songling's system — 4 limit cycles, distribution (3 : 1)** — Perko's eq. (5):

```
ẋ = P(x,y) = λ x − y − 10 x² + (5 + δ) x y + y²
ẏ = Q(x,y) = x + x² + (8ε − 25 − 9δ) x y
```
Shi's original values: **λ = −10⁻²⁰⁰ (Perko prints −10⁻²⁰⁰; Yu–Zeng print −10⁻²⁵⁰),
ε = −10⁻⁵², δ = −10⁻¹³.** Three cycles around (0,0), one large around (0,1).
Sizes: 10⁻⁶¹ < |y₂| < 10⁻¹⁹ and 0 < |y₁| < 10⁻⁶⁰. **Not numerically visible.**

**(P3) ★ PERKO'S NORMAL-SIZE (3 : 1) SYSTEM ★** — Perko's eq. (6) with modified
coefficients (p. 620 and §5, pp. 636–639):

```
ẋ = P(x,y) cos α − Q(x,y) sin α
ẏ = P(x,y) sin α + Q(x,y) cos α

P(x,y) = λ x − y − 10 x² + (5 + δ) x y + y²
Q(x,y) = x + x² + (8ε − 25 − 9δ) x y

α = −0.0023,  λ = −0.005,  ε = −0.01,  δ = −0.5
```
i.e. numerically
```
P(x,y) = −0.005 x − y − 10 x² + 4.5 x y + y²
Q(x,y) =  x + x² − 20.58 x y
```
Perko, p. 620, verbatim: *"the numerical results in §5 of this paper strongly indicate
that the quadratic system ... has exactly four limit cycles in the configuration (f)."*

**Cycle locations on the negative y-axis (four-figure accuracy, Perko p. 638):**
`y₁ = −0.0425, y₂ = −0.2160, y₃ = −1.3838` (three cycles around the origin), plus the
large cycle L₄ around (1,0).

> **THIS IS THE SINGLE MOST USEFUL STARTING POINT IN THIS SLICE.** It is an explicit,
> normal-size quadratic vector field with 3 hyperbolic cycles in one nest, embedded in a
> complete rotated family with 4 free parameters (α, λ, ε, δ). A 4th cycle in the nest
> would have to appear as a semistable cycle in the annulus between two of
> `y₁, y₂, y₃`, or outside `y₃`, or inside `y₁`. Perko's Figure 19 plots `Δ(c)/c` on the
> negative y-axis and shows exactly three sign changes — the displacement-function
> curve you would want to re-compute at higher precision over a 4-parameter sweep.

**(P4) Chin Yuan-Shun et al.'s system — 3 cycles in ONE nest (configuration (c))** —
Perko's eq. (7):

```
ẋ = λ x − y + (2 + δ) x y − y²
ẏ = x + λ y + x² − (5 + ε) x y − y²     [OCR of the ẏ line is partly damaged; see below]
```
with **λ = 10⁻³²⁰ (Perko prints 10⁻⁸²⁰), δ = −10⁻³⁰⁰ (Perko prints −10⁻³ˣˣ),
ε = −10⁻⁷⁸.**
Verbatim, p. 641: *"Chin et al. [4] showed that for λ = 10⁻⁸²⁰, δ = −10⁻³ˣˣ and
ε = −10⁻⁷⁸, the quadratic system (7) has three limit cycles around the origin in the
configuration (c) of Figure 1."*
Perko, p. 642, verbatim: *"it follows ... that the smallest limit cycle in Chin's
example (7) is very nearly a circle with a diameter of O(10⁻⁸²⁰)."*
Perko, p. 643, verbatim: **"A modification of Chin's example with three normal-size
limit cycles in the configuration (c) of Figure 1 has not yet been obtained."**

**(P5) Mieussens' normal-size modification of Chin — 2 cycles in one nest** —
system (7) with **ε = 0.1, δ = −0.06, λ = −0.0001**; cycles at
`y₁ ≈ −0.146, y₂ ≈ −0.403` (Perko p. 642, refining Mieussens' x₂ = 0.213).

**(P6) One-cycle examples (Perko & Shu 1984), useful as controls:**
```
(2)  ẋ = y + y²,      ẏ = −x + a y − x y + (1 + a) y²        0 < a < 1: exactly 1 LC
(3)  ẋ = y + y²,      ẏ = −0.5x + a y − x y + (0.8 + a) y²   0 < a < 0.8: exactly (1,1)
```

### The Andronov–Hopf function in Perko 1984

Perko calls it "the function α(x)". Verbatim, p. 629:

> "The inverse of each of the relationships shown in Figures 10 and 11 describes a
> continuous function, α(x), defined on [0, ∞). This is the function described in
> Theorem 3, p. 340 in [15] for a general class of quadratic systems."

([15] = M. Mieussens, *Sur les cycles limites des systèmes quadratiques*, C. R. Acad.
Sci. Paris **291** Série A (1980), 337–340. **PAYWALLED / not located.**)

**Operational reading:** for a rotated family, the map
`α ↦ {cycle radii}` inverts to a single-valued function `α = φ(r)` — the
**Andronov–Hopf function**. Number of cycles in the nest at parameter α₀ = number of
solutions of `φ(r) = α₀`. **3 cycles ⟺ φ has 2 interior extrema; 4 cycles ⟺ φ has 3
interior extrema** (an additional fold). Perko's Figure 18 is the graph with 2 extrema.
**Hunting a 4th cycle = hunting a THIRD extremum of the Andronov–Hopf function.**
Cherkas–Artés–Llibre make exactly this point (§3.2 below).

---

## 1.6 Perko 1995 (JDE) — multiple limit cycle bifurcation surfaces

**Citation.** L. M. Perko, *Multiple limit cycle bifurcation surfaces and global
families of multiple limit cycles*, J. Differential Equations **122** (1995), 89–113.
DOI 10.1006/jdeq.1995.1140.

**PDF: PAYWALLED (Elsevier).** No abstract available via Crossref.

The paper's content **is restated verbatim** by Gaiko in the open-access preprint
`gaiko_1104.3019.pdf` (arXiv:1104.3019), §4, "Bifurcation surfaces of multiple limit
cycles". Gaiko's preamble, verbatim:

> "In this section, we restate Perko's theorems on the local existence of
> (n−m+1)-dimensional surfaces, C_m, of multiplicity-m limit cycles for the polynomial
> system (3.1) with µ ∈ Rⁿ and n ≥ m ≥ 2. These results describe the topological
> structure of the codimension (m−1) bifurcation surfaces C_m. For m = 2, 3, 4, C₂, C₃,
> and C₄ are the familiar fold, cusp, and swallow-tail bifurcation surfaces; for m ≥ 5,
> the topological structure of the surfaces C_m is more complex. For instance, C₅ and C₆
> are the butterfly and wigwam bifurcation surfaces, respectively."

### Displacement-function derivative, as restated by Gaiko (Perko's eq.), verbatim:

> "d_s(0, µ₀) = exp ∫₀^{T₀} ∇·f(φ₀(t), µ₀) dt − 1   (3.2)
>
> d_{µ_j}(0, µ₀) = −ω₀/‖f(φ₀(0), µ₀)‖ ∫₀^{T₀} exp[ −∫₀^t ∇·f(φ₀(τ), µ₀) dτ ]
>                   f ∧ f_{µ_j}(φ₀(t), µ₀) dt   (3.3)
>
> for j = 1, …, n, where ω₀ = ±1 according to whether L₀ is positively or negatively
> oriented, respectively, and where the wedge product of two vectors x = (x₁,x₂) and
> y = (y₁,y₂) in R² is defined as x ∧ y = x₁y₂ − x₂y₁."

> "**Definition 3.1.** A limit cycle L₀ of (3.1) is a multiple limit cycle iff
> d(0, µ₀) = d_r(0, µ₀) = 0 and it is a simple limit cycle (or hyperbolic limit cycle) if
> it is not a multiple limit cycle; furthermore, L₀ is a limit cycle of multiplicity m
> iff d(0,µ₀) = d_r(0,µ₀) = … = d_r^{(m−1)}(0,µ₀) = 0, d_r^{(m)}(0,µ₀) ≠ 0."

### **Definition 4.1 (fold surface C₂ of multiplicity-two limit cycles)** — verbatim:

> "An (n−1)-dimensional analytic surface C₂ ⊂ Rⁿ is an (n−1)-dimensional fold
> bifurcation surface of multiplicity-two limit cycles of (3.1) through a point µ₀ ∈ Rⁿ,
> if for all ε > 0 there exists a δ > 0 such that for each µ ∈ C₂ with ‖µ − µ₀‖ < δ, the
> system (3.1) has a unique multiplicity-two limit cycle L_µ in an ε-neighborhood of L₀
> and the system (3.1) undergoes a fold bifurcation at L_µ; i.e., for ‖µ − µ₀‖ < δ, L_µ
> splits into a simple stable and a simple unstable limit cycles in an ε-neighborhood of
> L₀ for µ on one side of C₂ and L_µ vanishes for µ on the other side of C₂."

> "**Theorem 4.1.** Suppose that n ≥ 2, that for µ = µ₀ ∈ Rⁿ the system (3.1) has a
> multiplicity-two limit cycle L₀, and that d_{µ₁}(0, µ₀) ≠ 0. Then given ε > 0, there is
> a δ > 0 and a unique function g(µ₂,…,µₙ) with g(µ₂⁽⁰⁾,…,µₙ⁽⁰⁾) = µ₁⁽⁰⁾, defined and
> analytic for |µ_j − µ_j⁽⁰⁾| < δ, such that C₂ : µ₁ = g(µ₂,…,µₙ) is an (n−1)-dimensional,
> analytic fold bifurcation surface of multiplicity-two limit cycles of (3.1) through the
> point µ₀."

### **Definition 4.2 (cusp surface C₃ of multiplicity-three limit cycles)** — verbatim:

> "An analytic surface C₃ ⊂ Rⁿ is an (n−2)-dimensional cusp bifurcation surface of
> multiplicity-three limit cycles of (3.1) through a point µ₀ ∈ Rⁿ, if for all ε > 0
> there exists a δ > 0 such that for each µ ∈ C₃ with ‖µ − µ₀‖ < δ, the system (3.1) has
> a unique multiplicity-three limit cycle L_µ in an ε-neighborhood of L₀ and the system
> (3.1) undergoes a cusp bifurcation at L_µ; i.e., C₃ is the intersection of two
> (n−1)-dimensional fold bifurcation surfaces of multiplicity-two limit cycles of (3.1),
> C₂^±, which intersect in a cusp along C₃; **for ‖µ − µ₀‖ < δ and for µ in the cuspidal
> region between C₂⁺ and C₂⁻, the system (3.1) has three simple limit cycles in an
> ε-neighborhood of L₀**; and for ‖µ − µ₀‖ < δ and µ outside the cuspidal region, the
> system (3.1) has one simple limit cycle in an ε-neighborhood of L₀."

> "**Theorem 4.2.** Suppose that n ≥ 3, that for µ = µ₀ ∈ Rⁿ the system (3.1) has a
> multiplicity-three limit cycle L₀, that d_{µ₁}(0,µ₀) ≠ 0, d_{rµ₁}(0,µ₀) ≠ 0 and for
> j = 2,…,n, Δ_j ≡ ∂(d, d_r)/∂(µ₁, µ_j)(0, µ₀) ≠ 0. Then … C₃ … is an (n−2)-dimensional,
> analytic, cusp bifurcation surface of multiplicity-three limit cycles of (3.1) through
> the point µ₀ and C₂^± : µ₁ = g^±(µ₂,…,µₙ) are two (n−1)-dimensional, analytic, fold
> bifurcation surfaces of multiplicity-two limit cycles of (3.1) which intersect in a
> cusp along C₃."

### **Definition 4.3 (swallow-tail surface C₄ of multiplicity-four limit cycles)** — verbatim:

> "An analytic surface C₄ ⊂ Rⁿ is an (n−3)-dimensional swallow-tail bifurcation surface
> of multiplicity-four limit cycles of (3.1) through a point µ₀ ∈ Rⁿ, if for all ε > 0
> there exists a δ > 0 such that for each µ ∈ C₄ with ‖µ − µ₀‖ < δ, the system (3.1) has
> a unique multiplicity-four limit cycle L_µ in an ε-neighborhood of L₀ and the system
> (3.1) undergoes a swallow-tail bifurcation at L_µ; i.e., C₄ is the intersection of two
> (n−2)-dimensional cusp bifurcation surfaces of multiplicity-three limit cycles C₃^±
> which intersect in a cusp along C₄; furthermore, there are three (n−1)-dimensional fold
> bifurcation surfaces of multiplicity-two limit cycles of (3.1), C₂^{(i)}, i = 0,1,2,
> such that C₂^{(0)} and C₂^{(1)} intersect in a cusp along C₃⁺, C₂^{(0)} and C₂^{(2)}
> intersect in a cusp along C₃⁻, and C₂^{(1)} and C₂^{(2)} intersect along an
> (n−2)-dimensional surface on which (3.1) has two multiplicity-two limit cycles;
> finally, **for ‖µ − µ₀‖ < δ and for µ in the swallow-tail region, the system (3.1) has
> four simple limit cycles in an ε-neighborhood of L₀**; for ‖µ − µ₀‖ < δ and µ above the
> surfaces C₂^{(i)}, i = 0,1,2, the system (3.1) has two simple limit cycles in an
> ε-neighborhood of L₀; and for ‖µ − µ₀‖ < δ and µ below the surfaces C₂^{(i)},
> i = 0,1,2, the system (3.1) has no limit cycles in an ε-neighborhood of L₀."

> "**Theorem 4.3.** Suppose that n ≥ 4, that for µ = µ₀ ∈ Rⁿ the system (3.1) has a
> multiplicity-four limit cycle L₀, that d_{µ₁}(0,µ₀) ≠ 0, d_{rµ₁}(0,µ₀) ≠ 0,
> d_{rrµ₁}(0,µ₀) ≠ 0, and that for j = 2,…,n,
> ∂(d,d_r)/∂(µ₁,µ_j)(0,µ₀) ≠ 0, ∂(d,d_rr)/∂(µ₁,µ_j)(0,µ₀) ≠ 0,
> ∂(d_r,d_rr)/∂(µ₁,µ_j)(0,µ₀) ≠ 0. Then … [C₄ exists as an (n−3)-dimensional analytic
> swallow-tail bifurcation surface of multiplicity-four limit cycles]."

**★ THIS IS THE PIVOT OF THE WHOLE QUESTION.** Perko's Theorem 4.3 says: *if* a
quadratic system had a multiplicity-4 limit cycle at some µ₀ (with the nondegeneracy
conditions), *then* the quadratic parameter space (n = 12, or 5 after normalization —
in either case n ≥ 4) contains a swallow-tail region where the system has **four simple
limit cycles in one nest**. Gaiko's entire "H(2)=4" program is an attempt to prove that
no such µ₀ exists for quadratics. See §4.

---

# PART 2 — THE WINTNER–PERKO TERMINATION PRINCIPLE AND GAIKO'S CLAIM

## 2.1 The principle, verbatim

Sources (all open-access, all downloaded):
- `gaiko_geometry_planar_quadratic_math0611142.pdf` (arXiv:math/0611142), Theorem 4.1
- `gaiko_field_rotation_ihes.pdf` (IHES preprint), Theorem 5.1
- `gaiko_1611.08113.pdf` (Kukles), Theorem 4.1
- Also L. Perko, *Differential Equations and Dynamical Systems*, 3rd ed., Springer 2001,
  **§4.6 "Multiple Limit Cycles"** and §4.7 — book, not obtainable here.

> "**Theorem 4.1 (Wintner–Perko termination principle).** Any one-parameter family of
> multiplicity-m limit cycles of relatively prime polynomial system (4.1_µ) can be
> extended in a unique way to a maximal one-parameter family of multiplicity-m limit
> cycles of (4.1_µ) which is either open or cyclic. If it is open, then it terminates
> either as the parameter or the limit cycles become unbounded; or, the family terminates
> either at a singular point of (4.1_µ), which is typically a fine focus of multiplicity
> m, or on a (compound) separatrix cycle of (4.1_µ), which is also typically of
> multiplicity m."

> "**Theorem 4.2.** If L₀ is a nonsingular multiple limit cycle of (4.1₀), then L₀
> belongs to a one-parameter family of limit cycles of (4.1_λ); furthermore:
> 1) if the multiplicity of L₀ is odd, then the family either expands or contracts
> monotonically as λ increases through λ₀;
> 2) if the multiplicity of L₀ is even, then L₀ bifurcates into a stable and an unstable
> limit cycle as λ varies from λ₀ in one sense and L₀ disappears as λ varies from λ₀ in
> the opposite sense; i.e., there is a fold bifurcation at λ₀."

Note the word **"typically"** in Theorem 4.1. This is precisely the hedge that makes the
principle non-usable as a hard bound — see §2.3.

## 2.2 Gaiko's canonical quadratic system with FOUR field-rotation parameters

**Citation.** V. A. Gaiko, *Geometry of planar quadratic systems*, arXiv:math/0611142
(2006). **PDF: `/Users/scottg/Claude_all/papers/gaiko_geometry_planar_quadratic_math0611142.pdf`**

> "**Theorem 2.1.** A quadratic system with limit cycles can be reduced to the canonical
> form
>
>   ẋ = −y(1 + x + α y) ≡ P,
>   ẏ = x + (λ + β + γ)y + a x² + (α + β + γ)xy + c γ y² ≡ Q        (2.4)
>
> or
>
>   ẋ = −y(1 + ν y), ν = 0; 1,
>   ẏ = x + (λ + β + γ)y + a x² + (β + γ)xy + c γ y².                (2.5)"

Derived from Erugin's two-isocline reduction to
`ẋ = −y + m x y + n y², ẏ = x + λ y + a x² + b x y + c y²`, m = −1 or 0.

> "**Lemma 2.1.** Each of the parameters λ, β, γ, and α rotates the vector field of (2.4)
> in the domains of existence of its limit cycles, under the fixed other parameters of
> this system, namely: when the parameter λ, β, γ, or α increases (decreases), the field
> is rotated in positive (negative) direction, i.e., counterclockwise (clockwise), in the
> domains, respectively:
>   1 + x + α y < 0 (> 0);
>   (1 + x)(1 + x + α y) < 0 (> 0);
>   (1 + x + c y)(1 + x + α y) < 0 (> 0);
>   (λ + β + γ) y + (a − 1) x² + (β + γ) x y + c γ y² < 0 (> 0)."

with determinants
```
Δ_λ = P Q′_λ − Q P′_λ = −y²(1 + x + α y)
Δ_β = P Q′_β − Q P′_β = −y²(1 + x)(1 + x + α y)
Δ_γ = P Q′_γ − Q P′_γ = −y²(1 + x + c y)(1 + x + α y)
Δ_α = P Q′_α − Q P′_α =  y²((λ + β + γ)y + (a − 1)x² + (β + γ)xy + c γ y²)
```

**★ HIGHLY USEFUL FOR THE HUNT.** This is an explicit 4-parameter rotated embedding of
the full quadratic family (plus the two "shape" parameters a, c). Any 4-in-a-nest search
can be organized as a sweep over (a, c) with (λ, α, β, γ) as rotation knobs.

Gaiko's explicit 4-cycle (3:1) construction (Theorem 3.1 proof) uses
**a = 1/2, c = −1**, i.e.
```
ẋ = −y(1 + x + α y)
ẏ = x + (λ + β + γ)y + (1/2)x² + (α + β + γ)x y − γ y²        (3.1)
```
with the ordering `γ > 0` fixed, then `−1 ≪ λ < −γ < 0`, then `γ + λ ≪ α < 0`, then
`0 < −γ − λ < β ≪ 1`. Finite singularities at (0,0) and (−2,0); invariant line x = −1
when α = 0.

> "**Theorem 3.1.** A quadratic system can have at least four limit cycles in the
> (3 : 1)-distribution."

> "**Theorem 3.2.** A quadratic system has at most four limit cycles and only in the
> (3 : 1)-distribution."

## 2.3 Gaiko's claim of "at most 3 in a nest" — and why it is not accepted

> "**Theorem 4.3.** There exists no quadratic system having a swallow-tail bifurcation
> surface of multiplicity-four limit cycles in its parameter space. In other words, a
> quadratic system cannot have neither a multiplicity-four limit cycle nor four limit
> cycles around a singular point (focus), and the maximum multiplicity or the maximum
> number of limit cycles surrounding a focus is equal to three."

Gaiko's proof, verbatim (arXiv:math/0611142, pp. 10–11):

> "Suppose that system (2.4) with four field rotation parameters, λ, α, β, and γ, has
> four limit cycles around the origin … Thus, there is a domain bounded by three fold
> bifurcation surfaces forming a swallow-tail bifurcation surface of multiplicity-four
> limit cycles in the space of the field rotation parameters λ, α, β, and γ.
> The corresponding maximal one-parameter family of multiplicity-four limit cycles cannot
> be cyclic, otherwise there will be at least one point corresponding to the limit cycle
> of multiplicity five (or even higher) in the parameter space. Extending the bifurcation
> curve of multiplicity-five limit cycles through this point and parameterizing the
> corresponding maximal one-parameter family of multiplicity-five limit cycles by a
> field-rotation parameter, according to Theorem 4.2, we will obtain a monotonic curve
> which, by the Wintner–Perko termination principle (Theorem 4.1), terminates either at
> the origin or on some separatrix cycle surrounding the origin. Since we know absolutely
> precisely at least the cyclicity of the singular point (Bautin's result [1]) which is
> equal to three, we have got a contradiction with the termination principle stating that
> the multiplicity of limit cycles cannot be higher than the multiplicity (cyclicity) of
> the singular point in which they terminate.
> If the maximal one-parameter family of multiplicity-four limit cycles is not cyclic, on
> the same principle (Theorem 4.2), this again contradicts to Bautin's result not
> admitting the multiplicity of limit cycles higher than three."

### ★ ROUSSARIE'S OBJECTION — quoted verbatim by Gaiko himself (p. 12):

> "Thereupon, it makes sense to say some words on Roussarie's review MR2023976
> (2005d:37102) on [10]. The only concrete remark in this 'awkward' review is the
> following: **'I just mention the hazardous claim made in Theorem 4.12, page 137, that
> there exists no quadratic system having a swallow-tail bifurcation surface of
> multiplicity-four limit cycles. Looking at the proof, it seems that the author
> unfortunately confuses two different notions: paths of limit cycles, as defined in
> Definition 4.7, page 112, and lines of multiple limit cycles, as defined by Perko (and
> recalled in Definition 4.13, page 127). In fact, there is nothing forbidding that a
> path begin at a parameter value with a multiplicity-four limit cycle and end at a focus
> point'.**"

Gaiko's rebuttal is rhetorical, not mathematical. **Assessment: Roussarie is right.**
The termination principle says a family of multiplicity-m cycles terminates at a focus
"**typically** of multiplicity m". "Typically" is not "always". A path of limit cycles
that *begins* at a multiplicity-4 cycle can perfectly well *end* at an ordinary focus,
with the multiplicity dropping along the way. Bautin's cyclicity-3 result constrains only
what happens **in a neighbourhood of the singular point**, not the multiplicity of a
cycle of normal size elsewhere in the phase plane.

Gaiko's second (§3) "geometric" proof of Theorem 3.2 is even weaker: it repeatedly
argues *"a semi-stable limit cycle cannot appear in the domain D_i because of increasing
the distance between the spiral coils filling these domains"*. This is a heuristic about
spiral geometry with no proof; it is the same argument he applies verbatim to Liénard
and Kukles systems to "prove" their Hilbert numbers, and it has never been accepted.

**Bottom line for the hunt: the literature contains NO valid proof that a quadratic
system has at most 3 limit cycles in one nest.**

## 2.4 Related Gaiko items downloaded

- `gaiko_field_rotation_ihes.pdf` — IHES preprint, *Field Rotation Parameters and Limit
  Cycle Bifurcations*. Contains the general "k field-rotation parameters ⟹ at most k−1
  limit cycles around the origin" argument (Section 4), again via the untwisting-spirals
  heuristic. **Not rigorous, but the parameter-ordering bookkeeping is worth mining as a
  search strategy.**
- `gaiko_quadratic_two_parallel_isoclines_0803.3055.pdf` — *Limit Cycles of a Quadratic
  System with Two Parallel Straight Line-Isoclines* (with Zegeling; NWO-supported).
- `gaiko_lienard_smale_math0611143.pdf`, `gaiko_general_lienard_1202.3540.pdf`,
  `gaiko_1104.3019.pdf` (FitzHugh–Nagumo), `gaiko_1504.03353.pdf` (Holling),
  `gaiko_1611.08113.pdf` (Kukles cubic).
- V. A. Gaiko, *Wintner–Perko termination principle, parameters rotating a field, and
  limit-cycle problem*, J. Math. Sci. **126** (2005), 1259–1266. **PAYWALLED (Springer).**
- V. A. Gaiko, *Global Bifurcation Theory and Hilbert's Sixteenth Problem*, Kluwer,
  Boston, 2003 — the book Roussarie reviewed. **Chapter 4 "Multiple Limit Cycles and
  Wintner–Perko Termination Principle" is PAYWALLED (Springer).**

---

# PART 3 — CHERKAS AND THE DULAC–CHERKAS METHOD

## 3.1 Cherkas–Artés–Llibre 2003

**Citation.** L. A. Cherkas, J. C. Artés, J. Llibre, *Quadratic systems with limit
cycles of normal size*, Buletinul Academiei de Ştiinţe a Republicii Moldova. Matematica
**1(41)** (2003), 31–46. ISSN 1024-7696.

**PDF: `/Users/scottg/Claude_all/papers/cherkas2003.pdf`** (already present; complete.)

### Abstract, verbatim:

> "In the class of planar autonomous quadratic polynomial differential systems we provide
> 6 different phase portraits having exactly 3 limit cycles surrounding a focus, 5 of
> them have a unique focus. We also provide 2 different phase portraits having exactly 3
> limit cycles surrounding one focus and 1 limit cycle surrounding another focus. The
> existence of the exact given number of limit cycles is proved using the Dulac function.
> All limit cycles of the given systems can be detected through numerical methods; i.e.
> the limit cycles have 'a normal size' using Perko's terminology."

### Structural facts, verbatim (p. 31):

> "It is known (see, for instance [17]) that a quadratic system can have only limit
> cycles enclosing a unique singular point, which is a focus. As system (1) has no more
> than two foci [17], only the following distributions of limit cycles are allowed:
> n, (n₁, n₂), where n ∈ N, and n₁, n₂ ∈ N ∪ {0} with n₁ + n₂ > 0. … **Recently, Zhang
> Pingguang [20, 21] has proved that if nᵢ > 0 for i = 1, 2, then either n₁ = 1, or
> n₂ = 1.**"

([20] Zhang Pingguang, *On the distribution and number of limit cycles for quadratic
systems with two foci*, Acta Math. Sinica **44** (2001), 37–44 (Chinese);
[21] same title, Qualitative Theory of Dynamical Systems, 2003.)

> "The following distributions of limit cycles for quadratic systems (1) are known:
> (a) 1 and (1,0); (b) 2 and (2,0); (c) 3 and (3,0); (d) (1,1); (e) (2,1); (f) (3,1)."

> "**Remark 1.** Kooij and Zegeling proved in [18,19] that the distribution of limit
> cycles (3,1) is possible only for quadratic system of the type 2A + 1S∞,
> 2A + 2S∞ + 1S∞ which we have considered."

### **The normal form used** (p. 32, verbatim):

> "By means of an affine transformation of the phase variables and a change of the time
> scale, a quadratic system (1) generically can be written as
>
>   dx/dt = 1 + x y,
>   dy/dt = a₀₀ + a₁₀ x + a₂₀ x² + a₀₁ y + a₁₁ x y + a y²,          (2)
>
> where a₀₀ = a₀₁ + a₁₁ − a₁₀ − a₂₀ − a."

> "We remark that for the quadratic systems (2), **a₁₁ is a rotating parameter.**"

The weak focus / center is at **A = (1, −1)**, with conditions (verbatim):

> "System (2) has a weak focus or a center at the point A = (1, −1), if the conditions
>   L = 2a − a₀₁ − a₁₀ − 2a₂₀ > 0,  V₁ = a₁₁ + a₀₁ − 2a − 1 = 0,          (13)
> hold. The last condition says that the divergence of system (2) at A is zero."

Focal values (verbatim, eq. (15)):
```
V₃ = W₀ − a₁₀ W
V₅ = (4 − 2a − a₁₁) V / W
V₇ = −(a₁₁ + 2a + 1) U V / W
```
with
```
W₀ = a₁₁²(a+1) + a₁₁(2a² + a − 1) − a₂₀(a₁₁(2a−1) + (2a+1)(2a−3))
W  = −1 + 2a² + a(a₁₁ − 1)
V  = −a₁₁² a(a+1) + a₂₀(a−1)(2a+1)²
U  = (8 − 2a²)(a₁₁ + 2a + 1)² − 35(2a+1)(a₁₁ + 2a + 1) + 35(2a+1)²
```
Third-order weak focus condition (verbatim, eq. (16)):
```
a₁₁ = ã*₁₁ = 4 − 2a,  a ≠ 2
a₀₁ = ã*₀₁ = 2a + 1 − ã*₁₁
a₁₀ = ã*₁₀ = (6(a² − a − 2) + a₂₀(6a − 7)) / (1 − 3a)
(a − 3 − a₂₀)/(1 − 3a) < 0
```
> "In short, we note that we have a 2-parameter family of quadratic systems (2) with a
> weak focus of order three at A, the two parameters are a and a₂₀."

### **Cherkas's Liénard transformation of a quadratic system** — verbatim (p. 36):

> "Since the straight line x = 0 is transversal for the vector field associated to system
> (2), its limit cycles do not intersect x = 0. Therefore, in the half–planes x < 0 and
> x > 0 its limit cycles can be studied separately. In the half–plane x > 0 the
> transformation **x = 1/ξ, y = (ỹ − F(ξ))ξ^{−a} − ξ** writes system (2) into the Liénard
> system
>
>   dξ/dt = ỹ − F(ξ),  dỹ/dt = −g(ξ),        (12)
>
> where
>   f(ξ) = (a₁₁ + a₀₁ ξ − (2a+1)ξ²) ξ^{a−2},
>   g(ξ) = (a₀₀ + a₁₀ ξ + (a₀₀ − a₁₁)ξ² − a₀₁ ξ³ + a ξ⁴) ξ^{2a−3},
>   F(ξ) = ∫₁^ξ f(t) dt = P̃₂(ξ) ξ^{a−1} − P̃₂(1)."

> "Clearly, for a = 2, 3, … system (12) is a Liénard polynomial differential system.
> Moreover, for a = −2, −3, … system (12), under the transformation ξ = 1/x, ỹ = −y,
> goes over to
>   dx/dt = y + P̂₂(x)x^{−a−1} − P̂₂(1),  dy/dt = P̂₄(x)x^{−2a−3}.   (14)"

**★ This is the operational core of the whole Dulac–Cherkas program: it converts the
quadratic nest-counting problem into a Liénard nest-counting problem, where the
Dulac-function ansatz `Ψ = Σᵢ Ψᵢ(x) y^{n−i}` is tractable.**

### **THE UPPER-BOUND THEOREMS** — verbatim (pp. 34–35):

> "**Theorem 1.** Assume that system (1) is structurally stable in a connected region
> Ω ⊂ R². Then, there exist a function Ψ(x,y) ∈ C¹(Ω) and a constant k < 0, such that the
> inequality
>
>   Φ = k Ψ div f + (∂Ψ/∂x) P + (∂Ψ/∂y) Q > 0,  f = (P, Q),          (4)
>
> is satisfied in the region Ω. Moreover, the limit cycles of system (1) do not intersect
> the set W = {(x,y) ∈ Ω : Ψ(x,y) = 0}, and **in every two–dimensional connected
> subregion of Ω where either Ψ(x,y) > 0 or Ψ(x,y) < 0, system (1) has at most one limit
> cycle** γ, and if exists, is hyperbolic and stable (respectively unstable) if kΨ|_γ < 0
> (respectively > 0)."
>
> "If the function Ψ(x,y) satisfies the condition (4), the function B(x,y) =
> |Ψ(x,y)|^{1/k} is a Dulac function in each subregion Ψ(x,y) > 0 or Ψ(x,y) < 0, and we
> have that div(Bf) = Φ|Ψ|^{1/k−1}(sign Ψ)/k."

> "**Theorem 2.** Let Ω be a simple connected region where system (1) is defined and has
> a unique singular point, the antisaddle A with div f(A) ≠ 0. Assume that there exist a
> function Ψ and a number k < 0 satisfying the assumptions of Theorem 1. Suppose that the
> equation Ψ(x,y) = 0 determines in the region Ω a nest of m of ovals surrounding the
> point A. Then, in each of the m − 1 annulus limited by two adjacent ovals, system (1)
> has exactly one limit cycle. **Moreover, system (1) has in the region Ω at most m limit
> cycles.**"
>
> "By Theorem 2 it follows that the ovals are transversal to the vector field associated
> to system (1), and that the annulus limited by two adjacent ovals satisfies the
> Bendixson principle … An additional m-th limit cycle can exist between the most
> external oval and the boundary of the region Ω."

> "**Theorem 3.** Suppose that the function g(x) of the Liénard system (5) satisfies that
> g(0) = 0, and that its two nearest zeros at 0 are x₁ and x₂ with x₁ < 0 < x₂. Assume
> that there exist the constants k < 0 and Cᵢ for i = 1,…,n, such that the function Φ
> given in (9) is positive for x ∈ (x₁, x₂). **Then, system (5) has at most n/2 limit
> cycles surrounding the singular point (0,0).**"

**★ CRITICAL READING OF THEOREM 3.** The bound `n/2` is a function of the DEGREE OF THE
ANSATZ `Ψ = Σ_{i=1}^{n} Ψᵢ(x) y^{n−i}`, chosen by the user. It says: *if you can find a
degree-n Dulac–Cherkas function with Φ > 0, then ≤ n/2 cycles*. Cherkas et al. use
n = 10, 11, 12 for the 3-cycle examples — giving bounds 5, 5.5, 6 — and then sharpen to
exactly 3 by counting the ovals of `Ψ = 0` (Theorem 2) plus "reduction to global
uniqueness". **There is NO theorem in this literature bounding a quadratic nest by 3 in
general.** The Dulac–Cherkas method is a *certification* tool applied case by case, and
its failure to certify ≤3 for a given system is exactly what one would see if that
system had 4.

### The Andronov–Hopf function and the "reduction to global uniqueness" method

> "By definition the **Andronov–Hopf function** F : ∪_{α∈R} L(α) → R associates to the
> points of L(α) the value α. Therefore, the surface of limit cycles is determined by the
> equation α = F(x) running α in R.
> If the limit cycles surrounding the singular point x = 0 … intersect the half–axis
> x₂ = 0, x₁ > 0 only in one point, instead of function F(x) it is more convenient to
> consider the function **α = φ(x₁) = F(x)|_{x₂=0, x₁>0}**, which provides a full
> information about the limit cycles of system (3) surrounding the point x = 0, and their
> bifurcations when the parameter α varies."

> "We consider the Andronov–Hopf function AH(x) = a₁₁ … We fix all the parameters and we
> move only the parameter a₁₁. The function AH(x) is considered on the interval
> I₁ = [x₀, x_max] where the endpoints satisfy x₀ < 1 and x_max > 1, and x_max corresponds
> to the bifurcation of a limit cycle from a loop of the saddle S. **If in a subinterval
> I₀ = [x₁, x₂] of I₁ the number of zeros of the function AH(x) = a⁰₁₁ is 2p, then the
> number of limit cycles of system (2) in the strip x₁ < x < x₂ is p.** Now, suppose that
> the equation AH(x) = a¹₁₁ with a¹₁₁ < a⁰₁₁ provides a unique limit cycle which is
> localized in the strip x₃ < x < x₄ with [x₃,x₄] ⊂ I₀ ⊂ I₁, then the function AH(x)
> cannot take the value a⁰₁₁ outside the interval I₁. Consequently, for the value a⁰₁₁
> system (2) has exactly p limit cycles. **This is the method of reduction to the global
> uniqueness of a limit cycle.**"

### The construction recipe (verbatim, p. 38 — how they got 3 in a nest):

> "We fix the parameters a and a₂₀ of a system (2) having a weak focus of third order at
> A, and change the parameters a₁₁, a₀₁ and a₁₀ in order to obtain a quadratic system
> with one small limit cycle surrounding A, being A a weak focus of second order. We must
> change the parameters a₁₁, a₀₁ and a₁₀ in such a way that V₁ = 0 and V₅V₇ < 0. …
> [then a₁₀ alone for a 2nd cycle, then a₁₁ alone for a 3rd] … **and the appropriate
> Andronov–Hopf function a₁₁ = AH(x) will have two extrema.** … The limit cycles (which
> originally bifurcated from A) are not necessarily small."

> "**Remark 2.** For constructing the examples of quadratic system with the maximum
> number of limit cycles it is not necessary to use the function ã₁₁(x). It is enough to
> know that the function a₁₀(x) has an extremum, then the function AH(x) will have two
> extrema and provides the existence of an interval for the function a₁₁ in which the
> system has three limit cycles."

**★ THE 4-CYCLE CRITERION FALLS OUT IMMEDIATELY:** in this framework, **4 limit cycles
in a nest ⟺ the Andronov–Hopf function AH(x) = a₁₁ has THREE extrema on I₁** (i.e. the
horizontal line a₁₁ = const meets the graph 4 times). This is a completely concrete,
computable target: compute AH(x) numerically over the (a, a₂₀) parameter plane of
third-order weak foci and look for a third extremum.

## 3.2 ★ EXPLICIT COEFFICIENT VECTORS: 8 QUADRATIC SYSTEMS WITH 3 CYCLES IN ONE NEST ★

All in normal form (2): `ẋ = 1 + x y`, `ẏ = a₀₀ + a₁₀x + a₂₀x² + a₀₁y + a₁₁xy + a y²`,
with `a₀₀ = a₀₁ + a₁₁ − a₁₀ − a₂₀ − a`. Focus A = (1, −1). Cherkas–Artés–Llibre Table 1:

| N° | a | a₂₀ | a₁₁ | a₀₁ | a₁₀ | Singular points | Distr. |
|---|---|---|---|---|---|---|---|
| 1 | 3 | −12 | −1.398 | 8.4 | 15.28 | 1F + 1N + 2S∞ + 1N∞ | **3** |
| 2 | 1.5 | −15 | 0.79993 | 3.2 | 9.17 | 2F + 2S∞ + 1N∞ | **(3,0)** |
| 3 | −2 | 12 | 10.999 | −14 | −26.1 | 1F + 3S + 3N∞ | **3** |
| 4 | −2 | −1 | 9.49965 | −12.5 | 6.955 | 1F + 1S + 2N∞ + 1S∞ | **3** |
| 5 | −4 | −1 | 13.9987 | −21 | 12.4 | 1F + 1N + 2S + 2N∞ + 1S∞ | **3** |
| 6 | 5 | −50 | −5.49995 | 16.5 | 76.45 | 1F + 2N + 1S + 1N∞ + 2S∞ | **3** |
| 7 | 8/11 | −12 | 2.1502 | 67/220 | −26.5 | 2F + 1S∞ | **(3,1)** |
| 8 | 1.04 | −120 | 1.51997 | 1.56 | −79.6 | 2F + 2S∞ + 1A∞ | **(3,1)** |

**Cycle locations** (all cycles cross the line y = −1 at the listed x-values, verbatim
from the paper):

| N° | x₁ | x₂ | x₃ | second nest |
|---|---|---|---|---|
| 1 | 1.26 | 1.98 | 3.95 | — |
| 2 | 1.4 | 1.9 | 3.1 | (0 around the other focus) |
| 3 | 0.32 | 0.66 | 0.8 | — |
| 4 | 0.56 | 0.75 | 0.87 | — |
| 5 | 0.63 | 0.80 | 0.88 | — |
| 6 | 1.05 | 1.16 | 1.5 | — |
| 7 | 1.28 | 1.15 (sic — likely 2.15) | 4.43 | 1 cycle around B = (−3.2, 1/3.2) |
| 8 | 1.29 | 2.22 | 4.63 | 1 cycle around B = (−1.79, 1/1.79) |

Certifying data (Dulac–Cherkas ansatz orders used, all with k = −1 unless noted):
- Ex. 1: n = 10, k = −1, optimization net [−0.8, 0.5], N = 320; C* listed in paper.
- Ex. 2: n = 11, k = −1, net [0.2, 1.7], N = 200.
- Ex. 3: n = 10, k = −1, net [0.2, 1.72], N = 650; Liénard data
  `F(x) = −1001/3000 − 3x + 7x² − (10999/3000)x³`,
  `g(x) = 2x − 14x² − 2.1x³ + 26.1x⁴ − 12x⁵`.
- Ex. 4: n = 12, k = −1, net [0.3, 1.4], N = 750; Liénard data
  `F(x) = −782/9375 − 3x + (25/4)x² − (118747/375)x³`,
  `g(x) = 2x − (25/2)x² + (3291/200)x³ + (1391/200)x⁴ + x⁵`.
  Plus reduction-to-uniqueness at `a₁₁ = 9.4993` with n = 5, k = −2/3.
  **`AH(x) ≈ 8.89863 + 4.39482x − 13.5991x² + 22.9703x³ − 22.4248x⁴ + 11.9886x⁵
  − 2.72941x⁶` on I₀ = [0.6, 0.9].** ← *the only explicitly published Andronov–Hopf
  function in this literature. Its degree-6 form has exactly two interior extrema.*
- Ex. 5: n = 11, k = −1, net [0.5, 1.33], N = 750; reduction at `a₁₁ = 13.998`, n = 7,
  k = −2/3. Liénard data
  `F(x) = −17539/150000 − (7/3)x³ + (21/4)x⁴ − (139987/50000)x⁵`,
  `g(x) = x⁵(4 − 21x + (142/5)x² − (62/5)x³ + x⁴)`.
- Ex. 6: n = 10, k = −1, net [0.6, 1.21], N = 450; reduction at `a₁₁ = −5.4997`, n = 5,
  k = −2/3. Liénard data
  `F(ξ) = −22003/240000 − (1099999/80000)ξ⁴ + (33/10)ξ⁵ − (11/6)ξ⁶`,
  `g(ξ) = ξ⁷(−50 + (1529/20)ξ − (299/20)ξ² − (33/2)ξ³ + 5ξ⁴)`.
- Ex. 7: n = 11, k = −1, net [0.001, 4], N = 790; reduction at `a₁₁ = 2.156`, n = 3,
  k = −2/3. Non-polynomial Liénard (a = 8/11):
  `F(ξ) = 10130461/5700000 − 118261/(75000 ξ^{3/11}) + (67/800)ξ^{8/11} − (7/95)ξ^{19/11}`,
  `g(ξ) = −(−12)/(25 ξ^{17/11}) − 53/(5 ξ^{6/11}) + (8377/5500)ξ^{5/11}
  − (67/5500)ξ^{6/11} + (8/275)ξ^{27/11}`.
- Ex. 8: n = 10, k = −1, net [0.1, 2.2], N = 400; reduction at `a₁₁ = 1.5198`, n = 7.
  Non-polynomial Liénard (a = 1.04):
  `F(ξ) = −7749847/2040000 + (151997/40000)ξ^{1/25} + (3/20)ξ^{26/25} − (77/510)ξ^{51/25}`,
  `g(ξ) = −(−6)/(5 ξ^{23/25}) − (119/250)ξ^{2/25} + (5003/2500)ξ^{27/25}
  − (39/2500)ξ^{52/25} + (13/1250)ξ^{77/25}`.

> **★★ The nearness of the a₁₁ values to the "reduction" values is striking:
> Ex. 4: 9.49965 vs 9.4993 (Δ ≈ 3.5·10⁻⁴). Ex. 5: 13.9987 vs 13.998. Ex. 6: −5.49995 vs
> −5.4997. Ex. 7: 2.1502 vs 2.156. Ex. 8: 1.51997 vs 1.5198.
> These are the widths of the a₁₁-intervals for which 3 cycles exist. The intervals are
> TINY — 10⁻⁴ to 10⁻³ in a parameter of size O(10). If a 4-cycle window exists it would
> plausibly be narrower still, which is exactly why random search fails and why the
> Andronov–Hopf-function approach (find a 3rd extremum) is the right instrument.**

## 3.3 Cherkas's other work and the Dulac–Cherkas chain

- **L. A. Cherkas, I. L. Shevtsov**, *Normal-size limit cycles of quadratic systems with
  a structurally unstable focus*, Differential Equations **40** (2004), no. 8, 1076–1084
  (Springer, DOI 10.1023/B:DIEQ.0000049831.36933.b7). **PAYWALLED.**
  From the publisher abstract/description: quadratic systems with limit cycles of normal
  size were previously obtained with distributions **3, (3,0), (3,1)**; this paper treats
  a **structurally unstable focus** and obtains distributions **2, (2,0), (2,1)**, aiming
  at the maximum number of normal-size cycles for all configurations of singular points.
  The method: an **algebraic construction of a Dulac–Cherkas function in a neighbourhood
  of the focus as a polynomial of degree 4**.
  **No claim of 4 in a nest.**
- **L. A. Cherkas, A. A. Grin'**, *Bendixson–Dulac criterion and reduction to global
  uniqueness in the problem of estimating the number of limit cycles*, Differential
  Equations **46** (2010), 61–69. **PAYWALLED (Springer, DOI 10.1134/S0012266110010076).**
  Publisher description: "a regular method for localizing and estimating the number of
  limit cycles surrounding the unique singular point … dividing the phase plane into
  annulus-shaped domains with transversal boundaries in each of which a Dulac function is
  constructed by solving an optimization problem … the principle of reduction to global
  uniqueness … in the case of existence of an Andronov–Hopf function of limit cycles to
  obtain a sharp global estimate of the number of limit cycles."
- **L. A. Cherkas, A. A. Grin', K. R. Schneider**, *Dulac–Cherkas functions for
  generalized Liénard systems*, Electron. J. Qual. Theory Differ. Equ. **2011**, no. 35,
  1–23. Open access in principle; **the EJQTDE site did not serve the PDF at the time of
  this sweep.** Content: constructs a class of Dulac–Cherkas functions for
  `ẋ = y, ẏ = Σ_{j=0}^{l} h_j(x) y^j`; for `1 ≤ l ≤ 3` linear ODEs suffice, for `l ≥ 4`
  one must solve an overdetermined linear differential-algebraic system.
- **L. A. Cherkas**, *Dulac function for polynomial autonomous systems on a plane*,
  Differential Equations **33** (1997), 692–701. **PAYWALLED.**
- **L. A. Cherkas**, *Methods for counting the number of limit cycles of autonomous
  systems*, Differential Equations **13** (1977), 779–801. **PAYWALLED.** (This is
  reference [4] of Cherkas 2003 and the origin of Theorem 1 above.)
- **A. A. Grin, L. A. Cherkas**, *Dulac function for Liénard systems*, Proc. Inst. Math.
  Belarus NAS **4** (2000), 29–38 (Russian). **Not located.**
- **A. A. Grin', K. R. Schneider**, *Dulac–Cherkas function in a neighborhood of a
  structurally unstable focus of an autonomous polynomial system on the plane*,
  Differential Equations **50** (2014). **PAYWALLED.**

### The general Bendixson–Dulac upper-bound theorem (open-access substitute)

**Citation.** A. Gasull, H. Giacomini, *Effectiveness of the Bendixson–Dulac theorem*,
arXiv:2101.03874. **PDF: `/Users/scottg/Claude_all/papers/gasull_effectiveness_bendixson_dulac_2101.03874.pdf`**
(Companion: `gasull_extended_bendixson_dulac_1305.3402.pdf`, *Some Applications of the
Extended Bendixson–Dulac Theorem*.)

> "**Definition 1.1.** Given a function V : R² → R of class C¹ we will say that it is
> admissible if: (i) The vector ∇V vanishes on {V(x,y) = 0} at finitely many points.
> (ii) The set {V(x,y) = 0} has finitely many connected components. (iii) The set
> R² \ {V(x,y) = 0} has j connected components, U_i, i = 1,…,j, and for all of them
> ℓ(U_i) < ∞. Associated to V, we define the non negative integer number
> L(V) := Σ_{i=1}^{j} ℓ(U_i)."

> "**Theorem 1.2 (Bendixson–Dulac theorem).** Consider a C¹ planar differential system
> ẋ = P(x,y), ẏ = Q(x,y), and denote by X = (P,Q) its associated vector field. Let
> V : R² → R be an admissible function such that there exists s ∈ R⁺ for which the
> function
>
>   M_s := (∂V/∂x)P + (∂V/∂y)Q − s(∂P/∂x + ∂Q/∂y) V             (2)
>
> does not change sign and vanishes only on a null measure set. Define
> L_X(V) := N + L(V), where N is the number of periodic orbits of X contained in the set
> V = {V(x,y) = 0}. **Then, the differential system (1) has at most L_X(V) periodic
> orbits, which are limit cycles.** Moreover, each limit cycle not contained in V is
> hyperbolic, it is contained in one of the connected components U_i of R² \ V and, for
> each i = 1,2,…,j, **there are at most ℓ(U_i) limit cycles in the component U_i**. The
> stability of each of these limit cycles is given by the sign of −V M_s on the region
> U_i."

(Here ℓ(U) = number of "holes" of U. `V` is precisely a Dulac–Cherkas function with
`s = −1/k`.)

> "Observe also that, somehow, this version of the Bendixson–Dulac theorem relates the
> second part of the Hilbert's 16th problem, which deals with the number of limit cycles,
> with the first part, that deals with the number and distribution of the ovals of a
> planar algebraic curve."

**★ ASSESSMENT OF THE METHOD'S POWER.** The bound is `Σ ℓ(U_i)`, i.e. the number of ovals
of the Dulac–Cherkas curve. To *disprove* a 4-in-a-nest example you would need a
Dulac–Cherkas function whose zero set has only 3 nested ovals around the focus — which
is exactly what fails to exist if the system really has 4 cycles. **The method can
certify ≤3 but can never prove that ≤3 always holds.** It is therefore a *falsifier*, not
a *prover*: for a candidate 4-in-a-nest system, the failure of every degree-n
Dulac–Cherkas optimization to certify ≤3 is weak positive evidence.

---

# PART 4 — RELATED EXPLICIT SYSTEMS FOUND IN THIS SWEEP

**Citation.** P. Yu, Y. Zeng, *Visualization of Four Limit Cycles in Near-Integrable
Quadratic Polynomial Systems*, arXiv:2002.09987.
**PDF: `/Users/scottg/Claude_all/papers/kuznetsov_visualization_four_lc_near_integrable_2002.09987.pdf`**

Standard classification normal form, verbatim:
```
ẋ = −y + l x² + m x y + n y²
ẏ = x(1 + a x + b y)                                      (5)
```
(with n ≠ 0, rescale to n = 1). Equivalent (under (x,y)→(y,x)) to
```
ẋ = y(1 + a₁ x + a₂ y)
ẏ = −x + x² + a₃ x y + a₄ y²                              (7)
```
Center conditions for (7) at the origin (verbatim):
> "Q₃^R − Reversible system: a₃ = a₂ = 0;
> Q₃^H − Hamiltonian system: a₃ = a₁ + 2a₄ = 0;
> Q₃^{LV} − Lotka–Volterra system: a₂ = 1 + a₄ = 0; and
> Q₄ − Codimension-4 system: a₃ − 5a₂ = a₁ − (5 + 3a₄) = a₄ + 2(1 + a₂²) = 0."

**(Y1) Shi Songling's system with LARGE parameters that still gives 3 in a nest:**
```
ẋ = λ x − y − 10 x² + (5 + δ) x y + y²
ẏ = x + x² + (−25 + 8ε − 9δ) x y

λ = −2/10⁸ = −2×10⁻⁸,  ε = −1/1000,  δ = −1/10
```
Focus values at these values (verbatim):
`v₀ = −10⁻⁸, v₁ = 10⁻³, v₂ = −2079402109/112500000, v₃ = 59143866813736153313/1.62×10¹⁶`
Amplitudes of the three cycles around the origin: **r₁ ≈ 0.003636, r₂ ≈ 0.006431,
r₃ ≈ 0.070769.** (Plus the large cycle around (0,1).)
When λ = ε = δ = 0, `v₀ = v₁ = v₂ = 0, v₃ = 35625/8` (third-order fine focus).

**(Y2) Chen & Wang's system — 4 cycles, (3:1):**
```
ẋ = −δ₂ x − y − 3 x² + (1 − δ₁) x y + y²
ẏ = x + (2/9) x² − 3 x y                                  0 < δ₁, δ₂ ≪ 1
```
Origin is a **2nd-order** fine focus when δ₁ = δ₂ = 0 (note: only 2nd order, yet three
cycles are produced around (0,0) by perturbing — one from a trapping region plus two
from the fine focus).

**(Y3) Kuznetsov et al.'s form used for "four big size limit cycles":**
```
ẋ = y + x² + x y
ẏ = a₂ x² + b₂ x y + c₂ y² + α₂ x + β₂ y                  (8)
```

**(Y4) Near-integrable perturbation of the reversible center (Yu–Zeng eq. (9)):**
```
ẋ = y(1 + a₁ x) + ε a₁₀ x
ẏ = −x + x² + a₄ y² + ε(b₀₁ y + b₁₁ x y)
```
— yields four **normal-size** limit cycles.

Also present in `papers/` from parallel sweeps: `yu_han_2012_four_limit_cycles_perturbing_quadratic_integrable_1002.1055.pdf`,
`yu_han_eight_limit_cycles_around_center_quadratic_IJBC2013.pdf`,
`gasull_open_problems_low_dim_dynamical_systems_2012.02524.pdf`,
`gasull_santana_note_hilbert16_2407.13465.pdf`, `pubmat41.pdf` (Artés–Llibre, *Quadratic
vector fields with a weak focus of third order*, Publ. Mat. 41 (1997), 7–39 — the
parameter-space partition Cherkas et al. build on), `zegeling2024.pdf`.

Related paywalled items worth noting:
- **"Four normal size limit cycles in two-dimensional quadratic systems"**, Int. J.
  Bifurcation and Chaos **21** (2011), DOI 10.1142/S0218127411028532 (Leonov/Kuznetsov
  school). **PAYWALLED.** Abstract line: *"the existence criterion of three normal size
  limit cycles in quadratic systems with a weak focus of first order, and … a fourth
  normal size limit cycle … by giving a finite disturbance for weak focus, with
  bifurcation of appearance of two limit cycles via semistable cycle."*
  ★ This is a **(3:1)** result, not 4 in a nest, but the "**bifurcation of appearance of
  two limit cycles via semistable cycle**" is precisely the mechanism a 4th nest cycle
  would need.
- **"Visualization of Four Normal Size Limit Cycles in Two-Dimensional Polynomial
  Quadratic System"**, Differ. Equ. Dyn. Syst. (Springer, 2012), DOI 10.1007/s12591-012-0118-6.
  **PAYWALLED.**
- **"Limit cycles of quadratic systems with a perturbed weak focus of order 3 and a
  saddle equilibrium at infinity"**, Doklady Math. (2010), DOI 10.1134/S1064562410050042.
  **PAYWALLED.**

---

# PART 5 — WHAT THIS SLICE SAYS ABOUT WHETHER A NEST CAN HOLD 4 CYCLES

## 5.1 The state of the proof landscape

**There is no theorem, anywhere in this slice of the literature, that forbids 4 limit
cycles around a single focus of a quadratic system.** What exists is:

| Claim | Status |
|---|---|
| Bautin (1952): cyclicity of a focus/center of a quadratic system is exactly 3 | **THEOREM.** But it is *local*: it bounds the number of cycles that can bifurcate *out of the singular point*, in an arbitrarily small neighbourhood, under arbitrarily small quadratic perturbation. It says nothing about a nest of normal-size cycles. |
| Perko (1984): "there is no known example of a quadratic system with more than three limit cycles around a single critical point" | **EMPIRICAL OBSERVATION**, explicitly flagged by Perko as an open problem. |
| Perko (1984, Tung theorem): the (2:1) analysis | **PROVED ONLY UNDER THE HYPOTHESIS** "at most three limit cycles around (1,0)". |
| Gaiko (2003–2008): "a quadratic system cannot have four limit cycles around a focus" | **NOT ACCEPTED.** Roussarie's MathSciNet review MR2023976 identifies the fatal conflation of *paths of limit cycles* with *lines of multiple limit cycles*: "there is nothing forbidding that a path begin at a parameter value with a multiplicity-four limit cycle and end at a focus point." The alternative "geometric" proof rests on an unproved claim about the spacing of spiral coils. |
| Dulac–Cherkas certifications (Cherkas–Artés–Llibre 2003, Cherkas–Shevtsov 2004, Cherkas–Grin 2010) | **THEOREMS**, but they are *per-system certificates*. Each certifies ≤3 for a specific coefficient vector by exhibiting a Ψ with 3 ovals. They give no universal bound. |
| Zhang Pingguang: if both nests are non-empty, min(n₁,n₂) = 1 | **THEOREM.** This is what turns "≤3 per nest" into "H(2) ≤ 4". Conversely it means: **if you find 4 in a nest, H(2) ≥ 5 follows immediately if you can also keep 1 cycle in the second nest — and if the system has only one focus, 4 in a nest already gives H(2) ≥ 4 with a new configuration, while 4 in a nest + 1 elsewhere gives H(2) ≥ 5.** |

## 5.2 What Perko's own machinery says about the possibility

Perko 1995 Theorem 4.3 (Definition 4.3) is an **existence** theorem, not a
non-existence theorem: *if* a quadratic vector field somewhere in parameter space carries
a **multiplicity-4** limit cycle satisfying the three Jacobian nondegeneracy conditions,
*then* there is a nearby open swallow-tail region in which the system has **four simple
limit cycles in one nest**. Since the quadratic family has n = 5 essential parameters
after normalization (e.g. Cherkas's `(a, a₂₀, a₁₁, a₀₁, a₁₀)`, or Gaiko's
`(a, c, λ, α, β, γ)` with 6), the codimension-3 condition `n ≥ 4` for C₄ is comfortably
satisfied. **Nothing in the codimension count obstructs a swallow-tail in the quadratic
family.**

The only obstruction anybody has proposed is Bautin's number 3 — and Bautin's number
constrains the **singular point**, not a **multiplicity-4 cycle of normal amplitude**.
Roussarie's objection is precisely that the termination principle's word "typically"
allows a multiplicity-4 family to terminate at an ordinary (multiplicity ≤ 3) focus.

## 5.3 The rotated-field constraints a 4-cycle nest must satisfy

From Perko 1992 Theorems 2 & 3 combined with Duff's Table I, in a rotated family
`ẋ = P cos α − Q sin α, ẏ = P sin α + Q cos α`:

1. Every multiple cycle is **nonsingular** (Perko 1992 Thm 2), hence
2. The **only** bifurcation of cycles as α varies is a **fold at a semistable cycle**
   (Perko 1992 Thm 3(2), Conclusions §4), plus Hopf birth/death at the focus and
   birth/death on a graphic.
3. Cycles in a nest alternate in stability, and adjacent cycles move in *opposite*
   directions as α increases (Duff's table: an unstable and a stable cycle in the same
   nest have the same orientation but opposite stability, hence one expands and one
   contracts).
4. Therefore **in a rotated family, a 4th cycle in a nest can arise in exactly three
   ways**:
   (a) a **new semistable (multiplicity-2) cycle** is born in one of the annuli
       `(0, L₁), (L₁, L₂), (L₂, L₃), (L₃, ∂)` and splits — this is a *fold surface C₂*;
   (b) a **multiplicity-4 cycle** exists at an isolated parameter and unfolds
       (swallow-tail, Perko 1995 Thm 4.3);
   (c) a cycle enters from the Hopf bifurcation at the focus or from a separatrix
       graphic while the other three persist.
   Route (c) is blocked at the focus by Bautin (only 3 can come out of the focus at once),
   so **route (c) requires the 4th to come from a graphic**, and routes (a)/(b) require a
   semistable/quadruple cycle of *normal* amplitude.

**★ This gives a concrete search target: in a rotated quadratic family that already has
3 cycles in a nest, look for a value of the rotation parameter at which the displacement
function `Δ(c)` on the transversal develops a FOURTH sign change — equivalently, at which
`Δ(c)` acquires a double root (semistable cycle) in one of the annuli between the three
existing cycles or outside the outermost one.** Perko's Figure 19 for system (P3) plots
exactly this quantity `Δ(c)/c`; it shows three sign changes at
`y = −0.0425, −0.2160, −1.3838`. **Recomputing that curve at high precision across the
(α, λ, ε, δ) box is the single most direct experiment this literature suggests.**

## 5.4 The Andronov–Hopf-function formulation (Cherkas's framing) — the cleanest target

Cherkas–Artés–Llibre reduce the whole question to a single scalar function. With `a₁₁`
as the rotating parameter and `x` the abscissa where a cycle meets `y = −1`:

- number of limit cycles in the nest at `a₁₁ = c` = **number of solutions of AH(x) = c**;
- 1 cycle ⟺ AH monotone; 3 cycles ⟺ AH has **2** interior extrema;
- **4 cycles ⟺ AH has 3 interior extrema** (an N-shaped-plus-one graph).

For Example 4 they publish `AH(x)` explicitly as a degree-6 polynomial fit on [0.6, 0.9].
A degree-6 polynomial *can* have 5 critical points — **so nothing in the shape of AH
obstructs three extrema a priori.** Cherkas et al. reach 2 extrema by moving three
parameters in sequence (`a₁₁` → 1 cycle, `a₁₀` → 2 cycles, `a₁₁` again → 3 cycles); the
natural next move — **use the fourth free parameter (`a₀₁`, or `a₂₀`, or `a`) to induce a
third extremum** — is *not attempted anywhere in this literature*. They stop at 3 because
Bautin's 3 is the target they are aiming at, not because anything obstructed a 4th step.

> **This is, in my assessment, the most promising concrete gap in the literature.**
> The construction is a greedy sequential unfolding of a third-order weak focus: V₁, V₃,
> V₅, V₇ with `V₇ ≠ 0` gives at most 3 small cycles. But the cycles *"are not necessarily
> small"* (their words), and the 3 cycles they get are of normal size. A 4th cycle would
> not come from the focus (Bautin forbids it) but from a **fold in the annulus** — i.e.
> from a third extremum of AH. There is no theorem in this slice saying AH can have only
> two extrema.

## 5.5 My assessment

**Can a nest hold 4 cycles? The honest answer from this slice is: nobody knows, and the
belief that it cannot rests on (i) Bautin's local theorem misapplied globally, (ii) forty
years of failed search, and (iii) two published arguments that experts reject.**

Specifically:

- **The strongest formal barrier is absent.** There is no analogue of Bautin's theorem for
  cycles of normal amplitude, no bound on the multiplicity of a limit cycle of a quadratic
  system, and no Dulac-type universal bound. Perko's own framing (1984, closing remarks)
  and Perko 1995's swallow-tail theorem both treat 4-in-a-nest as *a priori possible*.
- **The rotated-field theory does NOT forbid it — it tells you how it would look.**
  Perko 1992 Thm 3 + Duff's table imply that in a rotated family the 4th cycle must
  appear via a fold at a semistable cycle in a specific annulus. That is a codimension-1
  event in a 5-parameter family: **it should be generic, not exotic, IF the fold surface
  C₂ reaches into the region where 3 cycles already exist.**
- **The evidence against is purely negative.** Every certified example tops out at 3, and
  every Dulac–Cherkas certificate found so far has exactly 3 ovals. But the parameter
  windows in which 3 cycles exist are of width ~10⁻⁴ in a coefficient of size ~10, and
  Shi's original 3-cycle window required coefficients of size 10⁻⁸²⁰. **A 4-cycle window,
  if it exists, could plausibly be far below the resolution of every search performed to
  date.** The historical pattern is instructive: Petrovskii–Landis "proved" H(2) ≤ 3;
  Shi Songling found 4 in 1979 with coefficients of order 10⁻²⁵⁰.
- **The two most actionable leads from this slice:**
  1. Take Perko's normal-size (3:1) system (P3) and Cherkas Examples 1–8, and compute
     the **full displacement function / Andronov–Hopf function** to high precision over
     the surrounding parameter box, looking for a fourth sign change of `Δ(c)` or a third
     extremum of `AH(x)`. Both are cheap, well-posed numerical problems and neither has
     been published.
  2. Search directly for a **multiplicity-4 limit cycle** in Gaiko's 4-rotation-parameter
     canonical form (2.4), i.e. solve
     `d = d_r = d_rr = d_rrr = 0` in `(λ, α, β, γ)` at fixed `(a, c)` — four equations,
     four unknowns, so **generically a discrete set of solutions**, and a solution would
     immediately give 4 cycles in a nest by Perko 1995 Theorem 4.3. If no real solution
     exists for any `(a, c)`, that would be the first *evidence-based* argument for
     "≤3 in a nest" — and would be a genuinely new result either way.

**Confidence:** I would put maybe 25–35% on "a quadratic system with 4 cycles in one nest
exists". The strongest reason to doubt it is sociological (many strong people have looked)
rather than mathematical; the strongest reason to believe it possible is that the only
formal barrier ever proposed (Bautin's 3) is provably about the wrong object, and the
swallow-tail machinery that would deliver it is fully available in the quadratic
parameter space's dimension count.

---

## Appendix: PDFs saved by this slice

```
/Users/scottg/Claude_all/papers/
  perko1984_RMJM_limit_cycles_quadratic.pdf            (complete, free)
  perko1990_TAMS_global_families.pdf                   (pp. 627-636 only)
  perko1992_ProcAMS_bifurcation_geometric_theory.pdf   (complete, free)
  cherkas2003.pdf                                      (complete)
  gaiko_geometry_planar_quadratic_math0611142.pdf
  gaiko_field_rotation_ihes.pdf
  gaiko_lienard_smale_math0611143.pdf
  gaiko_quadratic_two_parallel_isoclines_0803.3055.pdf
  gaiko_general_lienard_1202.3540.pdf
  gaiko_1104.3019.pdf   (restates Perko 1995 Defs 4.1-4.3, Thms 4.1-4.3 verbatim)
  gaiko_1504.03353.pdf
  gaiko_1611.08113.pdf
  gasull_effectiveness_bendixson_dulac_2101.03874.pdf
  gasull_extended_bendixson_dulac_1305.3402.pdf
  kuznetsov_visualization_four_lc_near_integrable_2002.09987.pdf
```

**Confirmed paywalled / unobtainable:** Duff 1953 (JSTOR); Perko 1975 & 1995 (Elsevier,
ScienceDirect bot-blocked); Perko & Shu 1984 (Elsevier); Mieussens 1980 (CRAS);
Cherkas–Shevtsov 2004, Cherkas–Grin 2010, Cherkas 1977/1997, Grin–Schneider 2014
(Springer); Gaiko 2003 Kluwer book & 2005 J. Math. Sci.; IJBC 2011 "Four normal size
limit cycles"; Perko's textbook §4.6.
