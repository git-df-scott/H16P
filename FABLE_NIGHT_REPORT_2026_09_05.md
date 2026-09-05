# Fable night report, 2026-09-05 (compiled 08:10 UTC)

No field with five limit cycles was found. No candidate exists to verify. Every lane opened during the
night reached a decisive negative, with the exception of two questions that are now theory rather than search.
All numbers below are NUMERICAL unless marked; scripts and data are in `audit/fable_engine/` and the
per-lane verdicts with figures are in `FABLE_LANES_2026_09_05.md`.

## Closed tonight

| Lane | Question | Verdict |
|---|---|---|
| F3/F5 | Random and Sobol sweeps of the Shi chart (Cartesian engine) | At most one cycle at trace zero, three otherwise; lower bounds only (first engine settings missed small cycles) |
| F6 | First-order Melnikov of the reversible center at the Shi loop point | Span three, at most two zeros |
| F11 | Order-two weak focus with a neutral homoclinic loop | Does not exist as a focus: trace is proportional to -eta_2 along the loop branch, the neutral point is a center |
| F12 | KKL double center (beta = K = J = 0) | Span three, at most two zeros (one three-zero direction in 400k at the noise floor) |
| F13 | Yu-Han reversible family on its degenerate curve, both annuli, a1 from -1.1 to -15 | The triple-zero element never has an interior zero: (3,1) is the ceiling of that mechanism |
| F14/F4/F20 | Descents from twelve (3,1) seeds, KKL incumbent and KKL c* seeds, Cartesian and compactified | Maximum four; the displacement between the outer roots is a single hump, no forming fold |
| F15 | Neutral-hemicycle unfoldings: two-center chart (20k) and KKL near c* (10k), compactified counting | Maximum three and four respectively; every four is (3,1) with the remote cycle at radius 1e6 to 1e9 |
| F16 | Dulac coefficients at the neutral hemicycle (transverse section) | First order: D = pi e0 - pi (2e1+e2) y + O(y^2); one cycle, two with the ratio parameter; exactly Marin-Villadelprat |
| F17 | Second-order Melnikov along the null direction at the holomorphic point (a,b) = (-1,1) | The direction is center-preserving; no second-order route |
| F18/F18b/F18c | Unfolding the hemicycle at its codimension-two point, neutral and at the proven alien points, including the physical scaling e0 ~ delta^2, da ~ delta (11k fields) | Two cycles from the graphic and never a third |
| F19 | Q4 alien test at its neutral resonant graphic: fifteen constructed three-zero directions, real cycles counted to radius e^40 | Exactly three real cycles every time, including zeros within 2e-4 of the boundary; no alien; no four-zero direction in 20k targeted trials |
| F21 | Full Shi chart with compactified counting, 30k fields | Maximum three |

## What the night established

1. Rigidity: every degenerate boundary object placed next to a degenerate focus collapsed onto an
   integrable stratum (order-three loop, order-two neutral loop, order-two neutral infinity graphic). A fifth
   cycle cannot come from stacking degeneracies in the field.
2. The two boundaries that stay degenerate for free, the equator connection and saddle resonance at
   infinity, emit exactly what first order predicts: two cycles at the neutral hemicycle, three interior
   zeros at Q4, no aliens in any tested direction or scaling.
3. Four-cycle fields are not rare: about 0.02 to 0.3 percent of random fields in natural five-parameter
   families are (3,1). Sixty generations of descent from more than forty such seeds never produced a fifth.
4. The instrument problem is solved: the compactified return map counts cycles at radius 1e12 in
   milliseconds and reproduces every published four-cycle field and Astra's KKL table.

## Open lanes for the day (two)

- **D1. Joint KKL fold work with Astra, cross-verified.** Astra's numerical fold at K = 1/512, c = 0.96889
  is the only pair-creating event anyone has seen. Continue it in the compactified counter through the
  neutrality line and down to K = 0, with the Hopf and remote gates, and hostile-verify whatever Astra
  reports. Expected outcome: (3,1) at best, but this is the one place a fold exists.
- **D2. Turn F11 into a theorem.** The observed law "saddle trace = -k eta_2 along a homoclinic loop
  around a second-order weak focus" is the order-two analogue of Li-Cherkas. Proving it (or the weaker
  sign law) closes K1's neighbour rigorously and explains why four is the ceiling in one nest. Start from
  the Zhang 1999 Liénard chart and the loop's Melnikov integral of the divergence.

Everything else is closed at the numerical level. A counterexample, if it exists, is not near any
integrable stratum, any degenerate graphic, or any published four-cycle field, and it is not reachable by
descent from them. The honest prior after tonight is well under one in twenty.
