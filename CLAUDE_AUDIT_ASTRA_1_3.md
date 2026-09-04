# Claude adversarial audit of Astra Strikes 1–3

Date: 2026-09-04. Audited base: main commit `8eb89c66c43b9570233774fb8bbd853a007634af`.
Role: adversarial verification half of the FASTRA pipeline. No Astra artifact
was modified or deleted. Everything below was reconstructed from the formulas
and scripts, not from the prose summaries, and every numerical statement here
was produced by the independent scripts in [`audit/`](audit/).

## 0. Verdict in one paragraph

All three strikes survive hostile review. No false theorem, no sign error,
no silent quantifier gap, and no numerical-to-analytic promotion was found.
The one place where the campaign is exposed is not an error but a strategy:
the proposed Strike #4 (search the `Y=Y'=0` locus) targets a set that a new
necessary condition, derived and measured here, shows to be empty everywhere
tested and empty in every asymptotic regime of the late-root corner. The
recommendation is therefore to convert Strike #4 into an exclusion theorem.
Details are in section 7 and in [FASTRA_H16_HANDOFF.md](FASTRA_H16_HANDOFF.md).

```text
ASTRA STRIKE 1 SOUND: YES
ASTRA STRIKE 2 SOUND: YES
ASTRA STRIKE 3 SOUND: YES
r > 5/11 CERTIFICATE: VERIFIED
LARGE-KAPPA OBSTRUCTION: VERIFIED
Q4 STILL LIVE: YES
Y=Y'=0 FOURTH STRIKE JUSTIFIED: NO
CANON CORRECTED: NO
```

## 1. What was replayed and what was built

Replayed from the repository (all exit 0):

| Command | Result |
|---|---|
| `python3 q4/test_q4.py` | 7 tests OK |
| `python3 q4/test_q4_second.py` | 4 tests OK |
| `python3 q4/test_q4_third.py` | 4 tests OK |
| `python3 q4/q4_structure_checks.py` | all exact checks pass |
| `python3 q4/q4_green_endpoint_third.py` | all exact checks pass |
| `python3 q4/q4_reconstruction.py` | bounded diagnostic passes, exact exclusions replayed |

Note for future replays: the repository test modules install a ten-second
per-process CPU fuse at import time, so they must be run as separate
processes, not collected into one `pytest` session (which is killed by the
fuse and reports nothing).

Built by Claude (all in `audit/`, all numerical unless stated; the exact
rational re-derivation in the threshold checker is exact):

| Script | What it attacks |
|---|---|
| `claude_check_strike1.py` | Stieltjes representation of `M`, Wronskian signs, `R(t)` monotone with `R(0)=54/31`, `R(1-)=1`, strip mapping, endpoint identity `H(1)`, closed moments `J_n`, `K_j` |
| `claude_check_reconstruction.py` | (R1)–(R3) plus coefficient transport on random `(A,B,eta)` at `kappa=1.15,1.7,3,9`, against direct area quadrature |
| `claude_check_threshold.py` | independent quadrature and an independently written exact rational series with its own tail bound for the frozen late-root point; equality and both-side controls at `r=5/11` |
| `claude_check_strike2.py` | `sup C_a`, `r(a)<-1/2`, `y` bounds, the `601/136136` ratio at box corners, the `8/9` integral, brute-force `Z(p_1)` with `tau_1<=5/11` |
| `claude_check_boundaries.py` | cusp construction and cusp exclusion mechanism, anchor-collision boundary, P2 cubic identity, `w''(kappa)`, orientation |
| `claude_check_endpoint_identities.py` | coefficient limit `(94/77,-17/77,1)`, (G1), the beta moments (G2), `Rcal` closed form, variation-of-parameters identity `Y=Phi*y+P*y_2` against the repository ODE |
| `claude_check_large_kappa.py` | scaled `P(1-eps)` along `gamma(1-eps)` at fixed `lambda`, versus the claimed limit (G5) |
| `claude_corner_map.py`, `claude_corner_derivatives.py`, `claude_shots_phi.py`, `claude_lobe_scan_phi.py` | the Strike-4 question (section 7) |
| `test_claude_hostile.py` | wrapper running the six checkers as separate processes |

## 2. Hostile review of Strike 1

Every claim was recomputed from the formulas.

- **Period pair and companion.** `K=(1-t)(F+6tF')` with `K={}_2F_1(-1/6,1/6;1;t)` holds to 35 digits at `t=0.3, 0.77`.
- **Stieltjes representation.** `M(t)=\int_0^1 rho(u)/(1-tu)du` with `rho(u)=3/(2\pi^2|F(1/u+i0)|^2)` was evaluated by direct contour-limit quadrature and agrees with `1-6(1-t)F'/F` to `10^{-26}` at `t=0.1,0.5,0.9`. This is the load-bearing analytic input for everything downstream; it is correct.
- **Wronskians.** For the ordered basis `(1,t,M,tM)`: `W_3=M''>0`, `W_4=3M''^2-2M'M'''<0` on a 20-point grid; the first-strike summary's signs `+,+,+,-` are consistent with Q4_STRUCTURE's `(1,t,w,v)` signs `W_3<0`, `W_4>0` by the constant basis change of determinant `-1`. No inconsistency.
- **Inflection function.** `R(t)=t+2M'/M''` is strictly decreasing, `R(0)=54/31` exactly from the rational moments, `R(0.999)=1.00134`, `R(0.99999)=1.0000126`. The strip `1<eta<54/31` and its original form `(54-23kappa)/31<beta_0<1` follow; the mapping was checked exactly. The withdrawal of the `kappa<85/23` cutoff is correct.
- **Cusp exclusion mechanism.** At `b=1.3` the triple contact is at `t_*=0.67932`, `q_0'''(t_*)=-0.492<0`, `H_{q_0}(t_*)=2.36\cdot10^{-4}>0`, and after unfolding with `lambda=10^{-3}` the primitive has exactly one sign change on `(0,1)`. The mechanism works as stated. Its scope (open neighborhoods of each fixed interior cusp, size not uniform toward the endpoints) is honestly stated.
- **Weighted-lobe inequalities.** The strictness question (does a double primitive zero, counted with multiplicity, evade the strict form?) was examined. It does not: five distinct `I` zeros force, by the anchored Rolle steps through `G`, the weighted `L_2` factorization and `Z(g)<=3`, three *distinct simple* zeros of `H`, which is exactly the strict form (L). The argument in Q4_ZERO_GEOMETRY is complete and I found no gap.
- **P2 cubic filter.** The identity `(g/(s-beta))''=2P_2(beta)/(s-beta)^3+w''(s)` was verified numerically at `kappa=3`, and `w''(kappa)=-25/(216(kappa-1)^2)`, `w'''>0` in `s` hold. The corrected filter is sound; the old linear one is unsupported.

Verdict: **VERIFIED**, no correction.

## 3. Hostile review of Strike 2

- **Lobe region nonemptiness and the certificate.** The exact rational series certificate is a correct proof of three primitive zeros of the frozen point (the tail majorant was re-derived independently and is valid: `f_n` decreasing, `d_n/f_n=6n/(6n-1)` decreasing, geometric tail). Combined with `Z(H)<=Z(q)<=3` (anchored Rolle plus the auxiliary ECT property) the four alternating signs give exactly three simple primitive zeros. The anchor-map theorem (unique normalized coefficient vector for three prescribed primitive roots) is a clean consequence of the multiplicity bound; surjectivity and openness follow as stated. I also checked that the anchor matrix determinant goes to zero linearly as two anchors collide (`delta=10^{-2},10^{-4},10^{-6}` give `det=-8.8\cdot10^{-7},-8.6\cdot10^{-9},-8.6\cdot10^{-11}`) with `H(midpoint)\to0`: the lobe-equality boundary behaves exactly as described.
- **Reconstruction (R1)–(R3), transport and the sign correction.** This is the single most important thing to verify, and it was verified in the strongest available way: random `(A,B,eta)` at four values of `kappa` (`1.15, 1.7, 3, 9`), three values of `t` each, comparing the scalar reconstruction with direct 40-digit area quadrature of the original four-term integral. Worst relative discrepancy `4.6\cdot10^{-11}` (at a value of size `10^{-6}`); typically `10^{-14}`. Since the four coefficient directions are exercised simultaneously, this confirms the linear transport map, the center data (R2), the forcing sign in (R1) and the multiplier in (R3) at every tested `kappa`. Astra's claim that Zhao's printed sign in equation (24) is wrong is correct in the repository's conventions.
- **Variation of parameters.** The decomposition `Y=Phi\,y+P\,y_2` with `y_2=y\,Rcal`, `Phi=Z-P\,Rcal`, was verified against the repository ODE solution to `10^{-12}`; `Rcal`'s hyperbolic closed form matches its integral to 25 digits.
- **`P_0<=0` sign chain, (S1)–(S3).** Re-derived. The monotone-lobe arguments are correct, including the endpoint behaviours `P\to+\infty`, `Z\to+\infty` (which need `H(1)<0`, true on the lobe region) and finiteness of `X(1)` for `a<1`.
- **The `5/11` theorem (E0).** Every ingredient was checked: `sup_a C_a=2.2114<9/4` (the supremum is interior, near `a\approx0.9`; the limit `a\to1` is `13/6`); `max_a r(a)=-0.50007<-1/2`; `(1-t)^{5/6}<y<(1-t)^{1/2}` on a grid; `\int_0^{5/11}(1-t)^{-13/6}dt=0.88133<8/9`; and the ratio bound `|Y_0|/eta<601/136136` is attained only in the limit at the box corner (`0.0044147029411` versus `0.0044147029441`), so it is valid but with no slack. The integer inequality `11\cdot3^{11}<2^7\cdot5^6` is right. The quantifiers are correct: the theorem says that for every `kappa>1` and every lobe point with first primitive root at most `5/11`, no five distinct original zeros exist; it uses only `P<=P_0` on `[0,p_1]`, which is forced by `H>0` there.
- **The `kappa>21636/19043` theorem.** Same ingredients plus `r(a)<-1/2`, checked; the arithmetic `a<=2593/21636` was replayed exactly by the repository and by hand.
- **Attempted counterexample.** A brute-force scan of lobe points with `tau_1<=5/11` and `a` on a grid found no configuration with a `P` root before `tau_1` at all (the `S1` band is `~10^{-7}` thin in `P_0`), consistent with, though weaker than, the theorem. No counterexample.

Verdict: **VERIFIED**, no correction.

## 4. Hostile review of Strike 3

- **Definition of `r`.** `r` is the first zero of the universal primitive `H`, identically equal to the first anchor along `gamma(r)=T(r,(1+r)/2,(3+r)/4)`. Verified: the coefficients from `coefficients_from_r` reproduce `H(r)=0` to `10^{-52}` (repository test) and my independent bisection finds the first root at `0.45454445`, `0.45454545`, `0.45454645` for `r=5/11\mp10^{-6}`, with `H(5/11)` changing sign through `7\cdot10^{-43}` at `r=5/11` exactly. The crossing is transverse and exact.
- **Rational late-root certificate.** Independently re-derived. Direct quadrature values at `23/32, 13/16, 29/32, 31/32` lie inside the JSON enclosures; my own rational series with my own tail bound reproduces the four signs with margins `1.7\cdot10^{-5}, 1.5\cdot10^{-5}, 1.0\cdot10^{-5}, 9.3\cdot10^{-5}` (all `>10^{-5}` as claimed); `H(5/11)=1.8667\cdot10^{-4}>0`; the box perturbation bound `122047/3072\cdot10^{-8}` was checked against the actual integral (`1.49\le39.7`, conservative). Because `Z(H)<=3`, the four signs prove exactly three simple primitive roots with the first strictly beyond `23/32>5/11` and none earlier. **VERIFIED.**
- **Endpoint asymptotics and (G1)–(G3).** Coefficients along the path converge to `(94/77,-17/77,1)` (at `eps=10^{-12}` the agreement is `10^{-9}`); `H_*` closed form (G1) holds to the precision of the evaluator; the beta moments (G2) are `1` and `25` to `10^{-12}` and `4\cdot10^{-4}` respectively (the second integrand is `(1-t)^{-5/6}\log`-singular; my quadrature is the limitation, the identity is exact and also proved in the repository by positive series).
- **The scoped fixed-`lambda` obstruction (G5).** The exact quantifiers, as stated by Astra and as I read the proof: for fixed strict anchor-ratio triple, fixed `lambda=1/(kappa\,eps)\in(0,\infty)`, and `eps\to0`, the scaled quantity `2304\pi P(1-eps c)/(eps^{5/6}L)` converges to `\mathcal P_lambda(c)`, locally uniformly in `lambda` on compact subsets; and `\mathcal P_lambda(1)>0`. The sign lemma (G4) was re-derived and is correct for all strict ratios (it only uses `D>0`, which is universal, and Rolle on `h(c)/c`). The matching constant is fixed by the two exact moments and the finite-part primitive; the repository's exact script confirms the primitive has no constant term at infinity and I confirmed `\mathcal B_lambda'=omega_lambda c^2` numerically. Direct quadrature of `P(1-eps)` (no ODE, no series) along the path gives:

| `lambda` | claimed limit `\mathcal P_lambda(1)` | `eps=10^{-2}` | `10^{-3}` | `10^{-4}` | `10^{-5}` | `10^{-6}` | `10^{-7}` |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.5 | 23.755 | 7.33 | 9.78 | 11.78 | 13.34 | 14.54 | 15.50 |
| 2 | 36.446 | 8.72 | 13.41 | 16.88 | 19.47 | — | — |

  The scaled values are positive at every `eps`, increase monotonically, and are consistent with a relative correction of order `8/L` (`L=\log(432/eps)`), which is the expected size of the neglected `O(eps L^2)/ (eps^{5/6}L)` and `1/L` terms. The sign conclusion (`P(tau_1)>0`, so (S1) fails) holds at every tested `eps`, not only asymptotically. **VERIFIED.** A scope remark that Astra did not make explicit: the convergence is logarithmic, so "sufficiently small `eps`" is not quantified by the proof; in practice the sign is already right at `eps=10^{-2}`.
- **What the theorem does not cover** is exactly what Astra listed. Section 7 addresses those regimes numerically.

Verdict: **VERIFIED**, no correction. The eight tuned shots and four confluent shots were also reproduced (my quadrature gives `Z(p_1)` to eight digits with `P(p_1)\approx10^{-16}`).

## 5. Global multiplicity logic

The certificate architecture "six rigorously alternating signs at interior rational points imply five distinct simple zeros and no others" is correct given Zhao's bound `Z(I;(1,kappa))\le5` counting multiplicity. Endpoint zeros do not count (open interval), `I` is analytic on the open annulus (no poles), and the six points are interior. The chain `Z(I)\le Z(G)\le Z(\mathcal F)+2\le Z(g)+2\le5` was re-derived from the repository's own objects: `I/h\to0` at the center gives the anchored Rolle step; the homogeneous `L_s` equation is disconjugate on `(1,\infty)` because `O(x)>0` for `x>0` is a positive solution, which gives the `+2`; `\mathcal F(kappa)=0` gives `Z(\mathcal F)\le Z(g)`; the auxiliary ECT property gives `Z(g)\le3`. Consistent with the repository's own numerical reconstruction. **VERIFIED.**

## 6. Claim ledger

| Claim | Strike | Status | Independently reproduced? | Proof type | Failure mode tested | Verdict |
|---|---|---|---|---|---|---|
| Universal chart `t=(kappa-s)/(kappa-1)`, `F={}_2F_1(1/6,5/6;1;t)`, `J_2/C=(1-t)(F+6tF')` | 1 | proved | yes (35 digits) | exact | wrong companion, normalization | VERIFIED |
| Positive Stieltjes representation of `M` | 1 | proved | yes (26 digits, direct contour quadrature) | analytic | wrong density, branch | VERIFIED |
| Wronskian signs `+,+,+,-` for `(1,t,M,tM)`; ECT | 1 | proved | yes (grid) | analytic | sign, basis-order confusion | VERIFIED |
| `R` decreasing, `R(0)=54/31`, `R(1-)=1` | 1 | proved | yes | analytic + exact moments | monotonicity, endpoint | VERIFIED |
| Corrected strip `(54-23kappa)/31<beta_0<1`; no `kappa<85/23` | 1 | proved | yes (exact mapping) | exact | reversed sign | VERIFIED |
| Cubic P2 filter | 1 | proved | yes (identity at `kappa=3`) | exact | missing factor 2, wrong power | VERIFIED |
| Cusp-neighborhood exclusion (`Z(I)\le3` near each interior cusp, all `kappa`) | 1 | proved | mechanism reproduced at `b=1.3` | analytic | `H(t_*)` sign, unfolding | VERIFIED |
| Strict lobe inequalities (L) necessary for five distinct zeros | 1 | proved | logic re-derived | analytic | multiplicity/strictness gap | VERIFIED |
| Lobe region = bounded analytic cell, anchor map diffeomorphism | 2 | proved | boundary behaviour reproduced | analytic | anchor collision | VERIFIED |
| Rational lobe box certificate | 2 | rigorous computation | repository replay + logic | exact rational | tail bound | VERIFIED |
| (R1)–(R3), (R2) center data, forcing sign correction | 2 | proved | yes, random coefficients, 4 `kappa` | exact + numerical | any sign/transport error | VERIFIED |
| Green kernel positive; `P`, `Z` factorization | 2 | proved | identity re-derived; `Y=Phi y+P y_2` checked | analytic | integrating factor | VERIFIED |
| `Y_0<0` on lobe region | 2 | proved | algebra re-derived | exact | box bounds | VERIFIED |
| `P_0\le0` excludes five zeros; (S1)–(S3) necessary and sufficient | 2 | proved | re-derived | analytic | endpoint limits | VERIFIED |
| `kappa>21636/19043` necessary | 2 | proved | ingredients checked | analytic + exact | `r(a)` bound | VERIFIED |
| First primitive root `>5/11` necessary (E0) | 2 | proved | ingredients checked, `C_a` sup measured | analytic + exact | `C_a`, `y` bounds, ratio bound | VERIFIED |
| Certified box excluded for all `kappa` | 2 | proved | follows from E0; margin replayed | exact | — | VERIFIED |
| Three-simple-plus-double reduced target | 2 | conditional | logic checked | analytic | — | CONDITIONAL (no point exists) |
| Path `gamma(r)` crosses `5/11` transversally | 3 | proved | yes | analytic + numerical | equality case | VERIFIED |
| Rational late-root box certificate | 3 | rigorous computation | yes, independent series and quadrature | exact rational | tail, margins, box | VERIFIED |
| Closed moments `J_n`, `K_j`, endpoint `H(1)` identity | 3 | proved | yes (30+ digits) | exact | — | VERIFIED |
| Coefficient limit `(94/77,-17/77,1)`, (G1), (G2), (G3) | 3 | proved | yes | exact | — | VERIFIED |
| Sign lemma (G4) | 3 | proved | re-derived | analytic | degenerate ratios | VERIFIED |
| Fixed-`lambda` obstruction (G5) with stated scope | 3 | proved | sign and trend at `lambda=0.5,2` | matched asymptotics | order of limits, constant | VERIFIED (convergence is logarithmic; sign robust at all tested `eps`) |
| Eight tuned shots fail (S2); confluent shots; reverse-tangency lines | 3 | numerical | `Z(p_1)` reproduced to 8 digits | numerical | — | VERIFIED as diagnostics |
| Six alternating signs certify five simple zeros | 1–3 | theorem-level | logic checked | analytic | endpoint zeros, poles | VERIFIED |
| "Q4 still live" | 1–3 | status | — | — | — | VERIFIED as of the audited proofs; see section 7 for what remains |

No claim received UNSUPPORTED, FALSE, or UNKNOWN. Two prose-level remarks:
the `601/136136` bound has no slack (it is sharp at the box corner), and the
fixed-`lambda` theorem's "sufficiently small `eps`" is not quantified because
the corrections are of relative size `~8/L`. Neither affects any conclusion.

## 7. Strike #4 (`Y=Y'=0`): killed as a search, replaced by a target theorem

### 7.1 What `Y` is and why `Y=Y'=0` is the boundary event

`Y=G/C` where `G=hI_h-I`, `C=\pi/\sqrt{kappa-1}`. With the variation-of-parameters decomposition
`Y=Phi\,y+P\,y_2` (`y` the positive homogeneous solution vanishing at the loop, `y_2=y\,Rcal` the one vanishing at the center),

\[
Phi(t)=Y_0+\int_0^t Rcal\,\Omega\,H\,du,\qquad P(t)=P_0-\int_0^t\Omega H\,du,
\]

and since `W(y,y_2)=1/p\ne0`, `Y(t_*)=Y'(t_*)=0` is exactly `Phi(t_*)=P(t_*)=0`: a critical point of `Z=Y/y` at height zero, i.e. the boundary of the first (S2) inequality `Z(p_1)>0`. So Astra's proposed event is the right boundary for the missing positive first maximum. This part of the proposal is correct.

### 7.2 A tuning-independent necessary condition

`Phi` does not depend on `P_0`, hence not on the delicate (S1) tuning of `kappa`; it depends on the universal point and on `a` only through `Rcal\,\Omega`. On `(0,tau_1)` we have `H>0`, so `Phi` is strictly increasing there, and `p_1<tau_1`. Therefore, for every five-zero point,

\[
\boxed{\;\Phi(\tau_1)=Y_0+\int_0^{\tau_1}Rcal(u)\,\Omega(u)\,H(u)\,du>0\;}\qquad\text{(N1)}
\]

This is a single scalar inequality on `(A,B,eta,a)`, with no shooting. It is weaker than (S2) but far easier to evaluate and to bound. Claude derived it; it is not in the Astra artifacts (they have the equivalent moment identity at `p_1`, which still requires the tuned `p_1`).

### 7.3 (N1) fails by 97–99.99% at every S1-tuned point Astra found

Recomputed by direct quadrature from the frozen JSON (my `Z(p_1)` agrees with Astra's to eight digits):

| path | `r` | `kappa` | `Phi(tau_1)/|Y_0|` |
|---|---:|---:|---:|
| affine | 0.5 | 2.177 | `-0.99998` |
| affine | 0.75 | 2.585 | `-0.99979` |
| affine | 0.9 | 3.273 | `-0.99868` |
| affine | 0.99 | 5.564 | `-0.98823` |
| affine | 0.9999 | 8.107 | `-0.96865` |
| power | 0.6 | 2.412 | `-0.99993` |
| power | 0.75 | 2.975 | `-0.99958` |
| power | 0.9 | 4.231 | `-0.99692` |

The gain `\int_0^{tau_1}Rcal\,\Omega H` is between `2\cdot10^{-5}` and `3\cdot10^{-2}` of `|Y_0|`. This is not a near miss; it is a different order of magnitude, and it explains the uniform `Z(p_1)\approx Y_0` that all of Astra's tables show.

### 7.4 The late-root corner: the only place where (N1) approaches equality, and there (S1) dies first

At the corner point `(94/77,-17/77,1)` and lift parameter `kappa`, by direct quadrature:

| `kappa` | `P_*(1;kappa)` | `Phi_*(1;kappa)/|Y_{0,*}|` |
|---:|---:|---:|
| 3 | `-5.76e-4` | `-0.983` |
| 5 | `-1.92e-4` | `-0.976` |
| 8 | `-8.86e-6` | `-0.968` |
| 8.5 | `+6.85e-6` | `-0.966` |
| 9 | `+2.03e-5` | `-0.965` |
| 16 | `+1.01e-4` | `-0.951` |
| 30 | `+1.19e-4` | `-0.933` |
| 100 | `+8.55e-5` | `-0.887` |
| 1000 | `+2.40e-5` | `-0.770` |
| 10000 | `+5.14e-6` | `-0.638` |

Facts established by these numbers together with Astra's exact results:

1. `P_*(1;kappa)` changes sign once, between `8` and `8.5`, and is positive for all larger `kappa` tested; asymptotically it behaves like `kappa^{-5/6}L`, which is exactly the finite part `D\,\mathcal B_lambda(0)>0` of Astra's matched expansion. So for fixed `kappa>kappa_c\approx8.3` and `eps\to0`, `P(tau_1)\to P_*(1;kappa)>0` and (S1) fails; the only fixed-`kappa` corner compatible with (S1) is `kappa=kappa_c`, where `Phi_*/|Y_0|=-0.966`.
2. `Phi_*(1;kappa)\to0^-` as `kappa\to\infty` (the two exact cancellations (G2)–(G3)), but empirically only like `kappa^{-1/6}L^{1.2}`: at `kappa=10^4` the deficit is still 64%.
3. In the joint regime `kappa=eps^{-theta}`: for `theta<1` the anchors sit inside the loop layer and `P(tau_1)=P_*(1;kappa)+O(eps L^2)+O(eps^{2-7theta/6}L)`, dominated by the positive `P_*(1;kappa)\sim eps^{5theta/6}L`, so (S1) fails; for `theta=1` Astra's theorem applies; for `theta>1` the anchors sit outside the loop layer and the leading inner term is `\mathcal P_0(1)=(6/5)D-\int_1^\infty v^{-13/6}(e-Vv+Qv\log v)dv>0` by (G4), so (S1) fails. Degenerating anchor ratios change `e,V,Q` but not the sign of any of these leading terms (the sign lemma survives every limit of the ratios, and `D=30/77` is universal).
4. The linearization at `kappa_c` along the universal deviation direction `(dA,dB,deta)=v(-643/462,1105/462,1)`, `v=eta-1>0`, which every late-root path follows to leading order, gives `dP(1)/dv=+0.0081>0` and `dPhi(1)/dv=-0.0165<0`: leaving the corner makes both (S1) and (N1) worse. (Further `kappa` values are appended to `audit/corner_derivatives.log` as they complete.)

Conclusion: the set `\{(S1)\}\cap\{Phi(p_1)=0\}` that Strike #4 wants to intersect with the lobe region is empty in every asymptotic regime of the late-root corner, and at every finite point examined (Astra's eight shots, four confluent shots, four reverse-tangency lines, and Claude's coarse landscape scan below) the necessary condition (N1) fails by more than 60%.

Coarse landscape scan (ten anchor triples, all with `y_1>5/11`, seven values of `kappa` from `1.2` to `1000`; `audit/claude_lobe_scan_phi.py`):

| anchors `(y_1,y_2,y_3)` | `eta` | worst-case `Phi(tau_1)/|Y_0|` over `kappa` |
|---|---:|---:|
| (0.46, 0.47, 0.48) | 1.5008 | `-0.999999` |
| (0.46, 0.7, 0.9) | 1.3291 | `-0.999986` |
| (0.46, 0.99, 0.999) | 1.1112 | `-0.999908` |
| (0.6, 0.8, 0.95) | 1.2540 | `-0.999923` |
| (0.7, 0.75, 0.8) | 1.3140 | `-0.999938` |
| (0.9, 0.95, 0.99) | 1.1069 | `-0.997143` |
| (0.5, 0.999, 0.9999) | 1.0746 | `-0.999835` |
| (0.99, 0.999, 0.9999) | 1.0139 | `-0.948445` |
| (0.999, 0.9995, 0.9999) | 1.0046 | `-0.879450` (at `kappa=1000`) |
| (0.95, 0.999, 0.99999) | 1.0275 | `-0.983856` |

The ratio increases monotonically toward the corner (all three roots late, `kappa` large) and nowhere else; the best value found is `-0.879`, i.e. a 12% deficit even at `kappa=1000` with roots within `10^{-3}` of the loop, where the (S1) band has already moved to `kappa\approx kappa_c\approx8.3` (deficit 97%).

### 7.5 Answers to the eight Strike-4 questions

1. `Y=G/C`, the center-anchored first Rolle derivative of the original integral, reconstructed as the solution of (R1) with center data (R2).
2. `Y=Y'=0` is `Phi=P=0`: the boundary of the first (S2) inequality. Correct boundary event.
3. Codimension: for fixed `(a,t_*)` it is an affine line in `(A,B,eta)`; the union over `(a,t_*)` is generically three-dimensional inside the four-dimensional `(a,A,B,eta)`, so there is enough freedom *in principle*. In practice the additional requirements (`t_*=p_1<tau_1`, (S1)) cut it down to the set discussed above.
4. It would allow the missing sign topology only if it met the lobe region with `t_*=p_1`; it does not, by (N1).
5. Yes: (N1). It already rules the locus out to the precision of everything computed, and analytically in the corner.
6. The locus can be searched in a bounded way, but there is nothing to find. The bounded object worth computing is instead `sup Phi(tau_1)/|Y_0|` over `L\times(0,1)`, which is a linear functional over an explicitly parametrized cell and is interval-certifiable.

### 7.6 What Strike #4 should become

Prove

> **Theorem N (target).** For every `kappa>1` and every point of the strict weighted-lobe region, `Phi(tau_1)<0`.

By (N1) this implies that no nonzero Q4 integral has five distinct zeros in the open annulus, i.e. at most four distinct zeros: a rigorous improvement of Zhao's bound (five with multiplicity) that closes Attack 1 with a publishable negative result. The proof architecture, launch conditions and the exact objects are in [FASTRA_H16_HANDOFF.md](FASTRA_H16_HANDOFF.md).

## 8. Silent-gap sweep

| Transition | Where it could hide | Finding |
|---|---|---|
| numerical → theorem | E0, `kappa` threshold, box exclusion | all analytic; numerics only motivate |
| transformed → original integral | (R1)–(R3), transport | random-coefficient test at four `kappa` |
| necessary → sufficient | (L), (S1)–(S3), thresholds | every document states "necessary, not sufficient"; the (S) chain is proved both ways |
| local → global | cusp exclusion | scoped as local; boundedness of the cell proved separately |
| pointwise → uniform asymptotic | (G5) | scope stated (compact `lambda`, fixed ratios); numerically consistent; logarithmic rate noted here |
| strict → non-strict | (L), thresholds | equality cases excluded explicitly and tested at `r=5/11` |
| auxiliary → original zeros | Strike 1 | never conflated; three-zero auxiliary mechanisms are labelled as such |

## 9. Files added by this audit

- `CLAUDE_AUDIT_ASTRA_1_3.md` (this file), `FASTRA_H16_HANDOFF.md`.
- `audit/*.py` checkers and `audit/test_claude_hostile.py`; logs of the long runs in `audit/*.log`.
- README and ASTRA_HANDOFF gained a pointer to this audit; no Astra artifact was altered.
