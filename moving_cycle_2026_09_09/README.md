# Executed moving-cycle / outer-fold strike

**No counterexample found.** Neither bounded branch produced an additional outer stationary point or a five-cycle field. The strongest newly replayed fields have four numerical cycle brackets. A prototype interval verifier passed existence and isolation gates for **one known cycle of the original seed**, with a higher-precision replay; it does not certify four or five cycles.

Source: main commit `498b58fc1378486601f348f67e0ab3a5362e7933`. This directory executes and refines `H16P_NEXT_STRIKE.md` and the accompanying `PLAN.md`. Only this directory is added; inherited results remain unchanged. The run stopped at its declared branch caps. No background research or automation is started by this package.

## What actually ran

1. Reproduced the six saved seed witnesses with the Cartesian field and a logarithmic field. Refined four numerical roots, checked complete physical-time first returns, and computed divergence/Floquet diagnostics.
2. Replayed all six witnesses using an independent 60-digit Cartesian Taylor method, with an additional tighter-tail replay of a weak witness. This method has numerical tail estimates, not interval remainder bounds.
3. Implemented analytic state and parameter variations, including initialization and section-event timing. Checked the vector-field Jacobians by complex-step differentiation and selected complete-return derivatives by independent central differences.
4. Sampled the upper outer displacement/derivative at 15 seed heights. No outer stationary point was detected.
5. Continued two explicit moving-root slices. At each accepted step, retained the three upper and one lower roots, local opposite-sign witnesses, nonzero numerical derivatives, and complete half-return gates. Examined the upper outer profile for a derivative sign change.
6. Diagnosed two failed far-outer returns using an additional positive time rescaling and a switch back to Cartesian coordinates. Localized the returning/nonreturning transitions with 12 bisections per branch and replayed four points on each returning side.
7. Continued inside the valid part of each section domain until **500 charged paired evaluations per branch**, including the intervening boundary probes. An accepted-step cap of 50 was also imposed but was not reached.
8. Independently replayed the four brackets in the two intermediate fields and the two final fields using exact decimal coefficient encodings and 60-digit Cartesian Taylor calculations. Full Cartesian first returns were also checked at the numerical roots.
9. Built an experimental interval-Taylor verifier and tested one original seed cycle. Fixed a failed Picard-tube enlargement strategy without weakening its inclusion gate. The revised implementation passed at 80 interval digits/order 32 and in a fresh process at 100 interval digits/order 36.

## Fields and results

All fields use

\[
\dot x=(b-2)/4+(1-b)y+ax^2+by^2+e_1x+e_2xy,
\qquad \dot y=e_0-2xy.
\]

The original exact seed has

\[
(a,b,e_0,e_1,e_2)=
(-7/4,1/3,-1/27500,-31379/250000000,-7517/50000000).
\]

Its numerical section heights are approximately 2.66510266, 4.01464333 and 7.73075810 in the upper nest, and −9355.42474 in the lower nest. The first height was subsequently refined with arbitrary precision for the interval smoke test.

The final replay fields are defined by the following **decimal strings interpreted exactly as rationals**, including the displayed finite approximation to b=1/3. Every bracket in a column was replayed in that one fixed field.

| Coefficient | Shape branch, final field | Amplitude branch, final field |
|---|---:|---:|
| a | −1.285811438560486 | −1.75 |
| b | 0.3333333333333333 | 0.3333333333333333 |
| e0 | −0.000036363636363636364 | −0.12283156115176493 |
| e1 | 0.000016233149616318536 | −0.6263511200791534 |
| e2 | −0.000364948697341617 | −0.46491694211759743 |

| Outcome | Shape branch | Amplitude branch |
|---|---:|---:|
| Accepted continuation steps | 18 | 17 |
| Charged discovery/boundary pair attempts | 500 | 500 |
| Upper section heights, approximately | 0.551004, 0.830019, 1.598317 | 5.488640, 8.267949, 15.921094 |
| Lower section height, approximately | −265.456800 | −1,238,704.297482 |
| Numerical brackets retained | 4 | 4 |
| Additional outer stationary points detected | 0 | 0 |
| End-of-run issue | Increasingly weak inner-cycle derivative; step reductions | Corrector/return failures beyond the last accepted point |

The amplitude parameter reached exp(8.125) ≈ 3377.867932 times its seed value. The other perturbation ratios were corrected to retain the cycles; this was not uniform multiplication of every perturbing coefficient by the same factor.

In the amplitude field, the approximate log multipliers of the four tracked cycles are 0.00194441, −0.00137868, 0.00414861 and −1.46302. These are numerical stability diagnostics. The remote cycle is still present, despite its large section height.

The final high-precision upper witnesses have alternating signs across each of the three local brackets. The lower witnesses also have opposite signs. See `shape_resumed_independent_replay.json` and `amplitude_resumed_independent_replay.json` for complete coefficient strings and values. No cycles from different fields have been combined.

## What the continuation covered

Write q=(a,b,e0/10^-4,e1/10^-4,e2/10^-4), and let s=log|y| on x=0. The three upper root coordinates have the form

\[
s_i=s_i^{\rm seed}+t\quad (i=1,2,3),
\]

where the common translation t is an unknown. Their two log-coordinate gaps are fixed; the lower root moves independently.

- **Shape slice:** increase a from −7/4 while holding b and e0 fixed; solve the three upper return equations for e1, e2 and t, then correct the lower root.
- **Amplitude slice:** hold a and b fixed, prescribe e0=e0(seed) exp(lambda), and solve for the two scaled perturbation ratios and t, then correct the lower root.

These choices supply definite local paths through the larger family. They do not cover all five parameter directions, varying upper gap ratios, additional normal-form charts, or disconnected components. Fixed upper gaps remain a restriction even though all four root positions move.

For each field, the outer profile used a small set of logarithmically separated section heights, beginning beyond the outermost tracked root. Derivative sign changes would have triggered stationary-point isolation. None occurred. Sparse negative samples do **not** establish derivative negativity throughout the intervening intervals. Consequently, this result is a bounded search failure, not an exclusion theorem for either slice or family.

## The outer-domain failures and a conditional endpoint lemma

At the first stopping points, shape lambda=0.3953125 and amplitude lambda=6.5, the smooth logarithmic chart failed at the farthest sample. A further positive time rescaling removed the artificial finite-time blow-up of that chart. Switching to Cartesian coordinates before y became too small showed the forward trajectories crossing y=0 before returning to x=0 in the upper half-plane.

For e0<0, the exact identity Q(x,0)=e0 means every crossing of y=0 is downward. A periodic orbit cannot cross this line, since closing the orbit would require an upward crossing. The numerical crossing diagnostics therefore identify a meaningful obstruction to these particular starting trajectories, rather than a new cycle. The location of the transition has not been interval certified.

The last returning/first failing log-height intervals were:

| Intermediate field | Numerical transition interval in s |
|---|---|
| Shape | [20.325946002150765, 20.328387408400765] |
| Amplitude | [21.001049118793674, 21.003490525043674] |

On the returning side, D becomes sharply negative and D_s becomes sharply negative. The data are consistent with a regular logarithmic endpoint singularity, not an approaching outer fold.

**Conditional local lemma.** Suppose an upper starting point (0,exp(s*)) reaches (0,0) in finite forward time, with the intended first half-return itinerary, and the backward half-return remains at a strictly positive height B(s*). Assume the nearby half-maps are defined up to these transversal endpoints. If 0<b<2 and exp(s*)>1/2, then the forward endpoint height A satisfies A(s*)=0 and A'(s*)<0. Indeed, the derivative of a planar half-map has the sign of the ratio of its normal velocities; here P(0,exp(s*))>0 and P(0,0)=(b−2)/4<0, while the exponential divergence factor is positive. With delta=s*−s,

\[
D(s)=\log A(s)-\log B(s)
=\log\delta+\log\!\left(\frac{-A'(s_*)}{B(s_*)}\right)+O(\delta),
\qquad
D_s(s)=-\frac1\delta+O(1).
\]

Thus no stationary point lies sufficiently close to this regular endpoint from the returning side. This is an analytic statement under its stated itinerary and regularity hypotheses; it does not provide a uniform neighborhood over the searched parameter curves. Our data support, but do not rigorously certify, those hypotheses at the two numerical transitions. In particular, the derivative-based estimates s−1/D_s converge toward approximately 20.32712628 and 21.00228811. This does not settle singular hemicycle or parameter-dependent endpoint regimes of the broader H16P campaign.

After identifying this boundary mechanism, the continuation resumed with inward-adjusted outer probes. Failed exterior probes remain in the records; their disappearance from an active profile was not counted as absence of cycles.

## Interval-verifier smoke test

The prototype uses directed interval arithmetic, a Picard inclusion to bound each flow segment, Taylor coefficients derived directly from the quadratic ODE, and a Lagrange remainder coefficient enclosed over that segment's entire tube. It rejects section-ambiguous steps, encloses both half-return events, and checks the full derivative by the normal-velocity ratio times the exponential of the divergence integral.

For the original exact seed, it uses a rational section center

```
2.6651026557169748663668514641384580595273000830967678977738577064
```

and radius 10^-20. The left and right full-return displacements have strict opposite signs, approximately −5.7814993×10^-27 and +5.7814993×10^-27. Over the whole initial interval, the derivative enclosure is contained in

```
[1.000000578149921, 1.000000578149940]
```

and therefore excludes one. Both the 80-digit/order-32 run and the fresh 100-digit/order-36 run passed all prototype gates. Each used 1,287 validated Taylor steps and six event enclosures across the two endpoint tests and the whole-interval test.

This is a local existence/isolation verification for one already-known seed cycle within an experimental implementation. The numerical Cartesian and logarithmic calculations provide independent cross-checks, but the interval code has not received an independent implementation audit. We do not label the entire four-cycle seed, either new four-cycle field, or any five-cycle field as certified by this smoke test.

The first prototype failed its Picard inclusion check because its enclosure iteration enlarged an already-sufficient coordinate while the other coordinate was still catching up. The repair keeps sufficient coordinates fixed and inflates only insufficient coordinates. The inclusion criterion was retained. Failed-attempt JSON files and logs remain archived. Higher-precision replay logs use explicitly outward-widened decimal serialization; the proof decisions themselves use binary interval endpoints.

## Accounting, limitations and reproducibility

- The two branch caps cover **1,000 charged paired discovery/boundary attempts**. Seed calibration, independent Cartesian/Taylor replays and interval verification are separate verification work, not hidden additional continuation steps. Per-process counters are retained; no account-wide usage percentage is inferred from them.
- There is one provenance gap: the original shape ledger retains evaluations 1–272, while its result charges 273 attempts. Its final failed outer attempt lacks its detailed original JSONL row. The location and failure appear in the checkpoint and were replayed explicitly by `shape_boundary.json`. The replacement replay is not represented as the original missing record.
- Limits include fixed upper log gaps, sparse outer sampling, numerical branch correctors, unvalidated parameter-domain boundaries, and partial corrector attempts stopped by the evaluation cap. No global H16P conclusion follows.
- The apparent approach to a weak inner-cycle derivative on the shape slice is not a proved Hopf collision. The amplitude corrector failures do not establish nonexistence beyond its last accepted point.
- All four numerical brackets of both final fields survived exact-decimal, arbitrary-precision replay. The extra pair required by this particular construction was not detected.

Scripts write their named outputs beside themselves. **To replay without replacing this evidence, copy the experiment directory and retain the repository's `reversible_reseed/data/verified_control.json` at its relative location.** Use the pinned source commit and the versions in `requirements.txt`. Typical commands in that copy are:

```bash
python seed_replay.py
python mp_replay.py
python outer_seed_profile.py
python continuation.py shape
python continuation.py amplitude
python boundary_replay.py shape
python boundary_replay.py amplitude
python final_field_replay.py
python continuation.py shape --resume
python continuation.py amplitude --resume
python final_field_replay.py --resumed
python interval_smoke.py
python interval_smoke.py 100 36
```

The directory manifest records artifact hashes; hashes establish integrity, not mathematical correctness. The useful next research change would release an upper gap ratio or identify a different endpoint organizer with a new compatibility calculation. Simply increasing the two exhausted branch budgets is not supported by the present evidence.
