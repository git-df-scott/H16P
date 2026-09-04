# Rigorous certification protocol

## Objective

Given exact rational coefficients, or coefficient intervals intended to prove
a whole box, certify **at least five distinct limit cycles**. Exact global
counting is optional.

The preferred proof object is a validated one-dimensional Poincare return map.
CAPD::DynSys supports interval ODE integration, event-defined Poincare maps,
derivatives, monodromy matrices, and arbitrary-precision MPFR intervals.

## Per-cycle proof

For cycle `i`, provide an affine transversal

\[
\Sigma_i=\{p_i+s v_i:s\in I_i\}
\]

with rational `p_i`, `v_i`, and rational interval `I_i`.

1. **Regularity and orientation.** Interval-evaluate the section normal against
   the vector field on the section flow box and prove `n_i dot F` has a fixed
   nonzero sign. Prove the box contains no equilibrium.
2. **Validated return.** Use an interval ODE solver with directed rounding to
   enclose the first return time and `P_i(J)` for every subinterval `J` used in
   the proof. Rule out an earlier section crossing.
3. **Existence.** Either enclose the displacement `D_i=P_i-id` at rational
   endpoints with opposite strict signs, or use interval Newton/Krawczyk:

   \[
   N(m,I_i)=m-\frac{D_i(m)}{D_i'(I_i)}\subset\operatorname{int}(I_i).
   \]

4. **Isolation and hyperbolicity.** Prove `0 notin D_i'(I_i)`, equivalently
   `1 notin P_i'(I_i)`. This gives a unique fixed point in `I_i` and hence an
   isolated periodic orbit. If this derivative test fails, subdivision or a
   higher-multiplicity analytic argument is mandatory.
5. **Floquet enclosure.** Record `P_i'(I_i)`. Independently, when convenient,
   enclose

   \[
   \mu_i=\exp\left(\int_{\gamma_i}\operatorname{div}F\,dt\right).
   \]

   The two enclosures must overlap and must exclude `1` for the hyperbolic
   gate.

## Distinctness

Five fixed intervals on one section are enough only if they correspond to
different first-return orbits. The certificate must additionally provide one
of:

- five pairwise disjoint validated flow tubes;
- five pairwise disjoint isolating annuli;
- a nesting proof using disjoint section intervals and Jordan-curve ordering,
  together with proof that each map is the first return;
- different enclosed equilibria/nests plus disjointness within each nest.

No trajectory sampling is accepted as a distinctness argument.

## Optional isolating-annulus route

For a rational polygonal annulus, interval-evaluate the outward normal
component of `F` on every edge. Opposite trapping orientations on the two
boundaries, plus absence of equilibria, allow Poincare--Bendixson to produce a
periodic orbit. To certify a **limit cycle**, add a return-map derivative or a
strict contraction/expansion argument; a trapping picture alone need not rule
out a period annulus.

## Global issues not required for a counterexample

The verifier need not exclude a sixth cycle. It need not cover the entire
Poincare sphere. Galias--Tucker needed global exclusion because their theorem
was “exactly four”; a counterexample needs only five local existence and
distinctness certificates. This makes verification materially easier than
global discovery.

## Replayable certificate format

Directory layout:

```text
certificate/
  manifest.json
  coefficients.json
  build-lock.txt
  proof_driver.cpp
  cycles/
    01.json ... 05.json
    01.log  ... 05.log
  hashes.sha256
```

Minimum `manifest.json` fields:

```json
{
  "claim": "at least five distinct limit cycles",
  "coefficient_encoding": "reduced rationals or hexadecimal interval endpoints",
  "vector_field_sha256": "...",
  "toolchain": {
    "capd_commit": "...",
    "compiler": "...",
    "mpfr_version": "...",
    "rounding_mode_test": "pass"
  },
  "precision_bits": 256,
  "cycles": ["cycles/01.json", "cycles/02.json", "cycles/03.json", "cycles/04.json", "cycles/05.json"]
}
```

Each cycle record must contain the exact section, initial interval, oriented
crossing bound, return-time interval, `P(I)`, `P'(I)`, displacement endpoint
signs or interval-Newton image, equilibrium exclusions, and a reference to the
distinctness object.

## Hostile checks

- Recompile with a second compiler and replay at higher precision.
- Perturb every printed decimal outward by one ulp; exact rational parsing must
  be invariant.
- Re-run with two set representations (e.g. doubleton and tripleton) to detect
  hidden wrapping assumptions.
- Split every accepted initial interval and verify compatible enclosures.
- Independently compute a Floquet enclosure from the divergence integral.
- Ensure no event is accepted at the initial time and every map is first
  return, not an iterate.

## Expected cost

For a well-separated hyperbolic candidate, five existence proofs should cost
minutes to hours, not a campaign. Cost can explode near a multiplier of one,
a saddle connection, or scale separation like the Songling `10^-200`
hierarchy. Adaptive precision and coordinate blow-ups are then essential.
Failure of interval enclosure is `UNRESOLVED`, never evidence that a numerical
cycle is false.
