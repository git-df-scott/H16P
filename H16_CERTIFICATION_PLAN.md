# Executable certification plan for at least five quadratic limit cycles

## Mathematical contract

Input: one exact rational real planar polynomial vector field of degree <=2, five proposed transverse section intervals and their full-return itineraries. Output: five independently replayable proofs of isolated periodic orbits. An exact global cycle count is unnecessary. Existence of at least five suffices, even if a sixth is present.

For Ji=[ai,bi], set D=P−identity. The basic gate is

    validated P exists on all of Ji,
    sup D(ai)<0<inf D(bi), or the reverse,
    1 not in P'(Ji).

Opposite endpoint signs give existence. Derivative separation supplies uniqueness and hyperbolicity in the interval. Alternatively, for mi in Ji, verify

    Ni = mi − (P(mi)−mi)/(P'(Ji)−1) is contained in interior(Ji).

Use outward-rounded interval operations throughout. Excluding zero from the denominator is required. An interval Newton or Krawczyk solver does not repair an invalid or discontinuous return map.

## What Galias–Tucker actually established

Galias–Tucker, *Applied Mathematics and Computation* 415 (2022), 126691, accepted September 27, 2021, prove their exact Songling field has exactly four limit cycles. [Open author paper](https://www.zet.agh.edu.pl/~galias/ps/amc2022.pdf), [DOI](https://doi.org/10.1016/j.amc.2021.126691).

Their Lemma 2 uses validated endpoint returns and existence of P over the complete section interval; MPFR-backed CAPD computations reach 1024 bits, taking 4–14 seconds per fixed point on their stated hardware. Lemma 3 uses interval Newton and derivative bounds, up to 2048 bits and 45 seconds–40 minutes per cycle. Taylor order 100 is used. These timings are historical measurements, not estimates for a new candidate.

The expensive part of their exact-count proof is excluding all other cycles. They combine interval iterations, polar-coordinate derivative estimates and analytic normal forms/Lyapunov functions. A fifth-cycle lower-bound proof does not need that global exclusion stage. Their three tiny multipliers are extremely close to one, so double precision cannot resolve them reliably.

The audit found a likely sign-label typo in Lemma 2: the printed signs match f=y−P(y) used in Figure 1 and Lemma 3, but are labelled P(y)−y. No inverse-map explanation appears there. Sign reversal leaves the existence brackets intact. Use Lemma 3's multipliers for stability and define displacement orientation explicitly in every replay.

## Toolchain and implementation

Use CAPD::DynSys with MPFR/MpInterval and a validated Taylor/Lohner solver, its Poincaré map support and C1 variational propagation. [Library paper](https://arxiv.org/abs/2010.07097), [official map documentation](https://capd.sourceforge.net/capdDynSys/docs/html/a05226.html). Treat decimal constants as named interval parameters initialized from exact rationals; do not silently parse them as binary double constants. In particular, setting 5−10^-13 and −25−8·10^-52+9·10^-13 from double precision loses the relevant focus data.

A clean implementation should expose the following executable operations:

```text
certify-field field.json --check-degree --check-equilibria
certify-return field.json sections.json --precision 256 --order 40
certify-return field.json sections.json --precision 512 --order 60
certify-distinctness field.json returns.json
verify-certificate certificate-directory/
```

These names specify the required command contract, not already installed executables. **The current repository contains numerical replay scripts, not a finished CAPD verifier.** Implement and pin this small verifier before accepting a new numerical candidate as a certificate. CAPD APIs, build options and the current release must be checked against its pinned source. No maintained one-command source archive for the exact Galias–Tucker proof was located.

Required algorithm:

1. Parse rational coefficient numerators/denominators. Construct their outward-rounded enclosures at the chosen precision. Verify all monomials have total degree <=2.
2. Isolate equilibria in the relevant domains. Enclose F·n on each section and prove a strict sign. Validate entry/exit times so the initial t=0 crossing is not counted as a return.
3. Propagate each entire section interval using a validated flow set with wrapping control. Use a sequence of local transverse sections or compactified charts near infinity. Prove no premature/alternate return, missed event or blow-up occurs.
4. Propagate the variational equation. For section s(z)=0 use the projected derivative `(Id − F n^T/(n^T F)) Dphi`; reduce it to the scalar section coordinate. On y=0 the scalar formula is `P'=Q(initial)/Q(return)*exp(integral div F dt)`. The divergence exponential alone is a multiplier only at a fixed point.
5. Evaluate endpoint displacement signs and an interval Newton inclusion. Record strict interval margins. Increase precision/order or subdivide the initial section if a bound overlaps zero; do not replace an inconclusive enclosure by a floating-point sign.
6. Prove the five tubes are distinct by common-section ordering, one crossing per full lap, two half-planes and nesting, or pairwise disjoint tubes. Merely listing five roots on unrelated sections is insufficient.
7. Replay all gates from a clean process with fixed software hashes. Have an independent reviewer verify the exact coefficients, displacement sign convention and orbit identities.

## Precision and chart strategy

Start the KKL controls at 128–256 interval bits, increasing to 512 when multipliers approach one. Original GT is a stress control requiring 1024–2048 bits in the published implementation. Scale tiny cycles as (x,y)=r(X,Y); use logarithmic or blown-up section coordinates. Near a focus, prove angular monotonicity before integrating dr/dtheta. Near saddles, combine a validated local Dulac map with regular transition maps to control long flight times. A numerical asymptotic expansion at infinity needs an explicit remainder before certification.

Taylor models help retain parameter dependence in persistence boxes; interval Newton/Krawczyk helps solve fixed-point or multi-shooting equations. An isolating annulus with verified inward/outward flow and no equilibria is an alternative Poincaré–Bendixson existence certificate. To establish isolated hyperbolic cycles robustly, prefer return-derivative bounds. Guckenheimer's geometric boundary approach can reduce long integration costs when a suitable annulus is available. [Primary paper](https://pi.math.cornell.edu/~gucken/PDF/planar_proofs.pdf).

## Required manifest

For each orbit save:

| Field | Meaning |
|---|---|
| exact_field | Rational coefficients and polynomial basis in a declared order |
| section | Equation, orientation, coordinate interval and parameterization |
| itinerary | Ordered charts/sections, complete-return and no-escape bounds |
| return_time | Positive interval enclosing the first full return |
| D_left, D_right | Strict outward-rounded endpoint sign intervals |
| Pprime | Derivative enclosure excluding 1, or Newton/Krawczyk proof data |
| flow_tube | Validated enclosures and local transition certificates |
| equilibrium_nest | Which unique focus lies inside; distinctness proof |
| arithmetic | Rounding mode, MPFR precision, Taylor order and remainder policy |
| provenance | Input/source/build hashes, versions, commands and complete logs |

For the Hopf strike additionally save a common parameter box proving persistence of J1–J4, the exact K>0 and nonzero trace-derivative bounds, and the final explicit beta<0. Either bound the analytic Hopf remainder or directly certify the new J0. Existence for “sufficiently small beta” cannot substitute for exact final coefficients.

## Controls and acceptance levels

1. Replay the ordinary KKL four-cycle numerical control and confirm S/U/S plus remote U.
2. Replay the 900-bit GT small-cycle sign checks as a numerical stress test, then replicate the published **interval** brackets before claiming the new verifier handles extremely weak multipliers.
3. Test an escape trajectory near GT section y≈0.03689094: the verifier must reject continuity/return existence rather than infer a cycle from a jump.
4. Test a quadratic center: periodic ovals must not pass the derivative-isolation gate as distinct limit cycles.
5. Optionally reproduce the genuine quadratic Hamiltonian alien-two control in the graphics appendix. It exercises two-saddle itineraries; its nearby global bound prevents treating it as a five-cycle seed.

Current evidence levels: published GT theorem **COMPUTER-ASSISTED RIGOROUS / ACCEPTED**; new repository replays **NUMERICAL ONLY**; no new five-cycle certificate. The certificate plan is executable in mathematical and software-interface detail, but implementing/running the interval verifier remains a required gate for Strike #5, not work falsely reported as completed here.
