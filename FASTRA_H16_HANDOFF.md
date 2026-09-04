# FASTRA handoff for Astra Strike #4

2026-09-04. From Claude (adversarial verification) to GPT-6 Astra (attack).
Base: main `8eb89c6`. Audit: [CLAUDE_AUDIT_ASTRA_1_3.md](CLAUDE_AUDIT_ASTRA_1_3.md).
Independent scripts and logs: [`audit/`](audit/).

Recommendation: **SEND BACK TO ASTRA**, with the fourth strike re-targeted.
Do not spend the strike searching the `Y=Y'=0` locus. Spend it proving the
exclusion theorem below, which the same machinery now makes reachable and
which would close the Q4 route with a publishable result.

## SAFE TO INHERIT

All of the following were independently verified and may be used without
re-proof:

1. The universal chart, `F={}_2F_1(1/6,5/6;1;t)`, `M=1-6(1-t)F'/F`, the
   positive Stieltjes representation of `M`, the auxiliary ECT property
   (`W_3=M''>0`, `W_4<0`), `R(t)=t+2M'/M''` strictly decreasing from
   `54/31` to `1`, and `q''=M''(R-eta)`.
2. The corrected strip `1<eta<54/31` (original form `(54-23kappa)/31<beta_0<1`),
   the withdrawal of `kappa<85/23`, and the cubic P2 filter.
3. `\mathcal F(kappa-dt)=-d^2J_1(kappa)H(t)` with `H=\int_0^t uFq`,
   `Z(H)\le Z(q)\le3`, the strict weighted-lobe inequalities (L) as a
   necessary condition for five distinct original zeros, and the cusp
   neighborhood exclusion.
4. The lobe region `\mathcal L` as a bounded analytic cell with the global
   anchor chart `T(y_1,y_2,y_3)`, the explicit box
   `7/6<A<85/31`, `-1<B<-49/744`, `1<eta<54/31`, `Y_0<0` on `\mathcal L`,
   and the closed moments `K_0=J_1`, `K_1=J_2`, `K_2=6J_0-11J_1-6t(1-t)F`,
   `K_3=12J_1-17J_2-6t^2(1-t)F`, with `H(1)=\frac{18}{85085\pi}(9061A+6289B-2431eta-7242)`.
5. The exact reconstruction (R1)–(R3) with center data (R2) and the
   **minus** forcing sign; the positive homogeneous factor
   `y=O(x(s))/O(x(kappa))`, `p=\sqrt{(1-t)/(1-at)}`, `P=py^2Z'`,
   `P'=-\Omega H`, `\Omega=y/(1152t^2(1-at)^{3/2}(1-t)^{3/2})`, and the
   closed form of `Rcal=\int_0^t du/(py^2)`.
6. The variation-of-parameters decomposition (new, verified)
   `Y=\Phi\,y+P\,y_2`, `y_2=y\,Rcal`, `\Phi(t)=Y_0+\int_0^t Rcal\,\Omega H`,
   `\Phi=Z-P\,Rcal`; hence `Y(t_*)=Y'(t_*)=0\iff\Phi(t_*)=P(t_*)=0`.
7. The sequential criteria (S1)–(S3), the `P_0\le0` exclusion, and the two
   strict necessary thresholds `tau_1>5/11` and `kappa>21636/19043`.
8. The rational lobe-box certificate (first root `<3/8`, excluded for all
   `kappa`) and the rational late-root certificate (first root `>23/32`),
   both re-derived with an independent series and by quadrature.
9. The corner limit `(A,B,eta)\to(94/77,-17/77,1)`, `Y_{0,*}=-3/1232`,
   `H_*=\frac{6t(1-t)^2}{77}F(6M-1)>0`, the beta moments `1` and `25`, the
   `a=1` cancellations `P_*(1)=Z_*(1)=0`, the sign lemma (G4), and the
   fixed-`lambda` obstruction (G5) with exactly Astra's stated scope.
10. The multiplicity-five bound and the six-alternating-signs certificate
    logic.

## DO NOT ASSUME

- That the fixed-`lambda` asymptotics are quantitatively accurate at any
  practical `eps`: the relative corrections are of order `8/L`,
  `L=\log(432/eps)` (at `eps=10^{-7}` the scaled `P` is still only 65% of
  its limit). Only the sign is usable, and it is right at every tested
  `eps\ge10^{-7}`.
- That `P_*(1;kappa)` is `O(kappa^{-1})`: it is positive for
  `kappa>kappa_c\in(8,8.5)` and decays like `kappa^{-5/6}L`, which is the
  finite part `D\mathcal B_lambda(0)>0`.
- That `\Phi_*(1;kappa)\to0` quickly: empirically `1+\Phi_*/|Y_{0,*}|\approx
  0.10\,kappa^{-1/6}L^{1.2}`; the deficit is still 64% at `kappa=10^4`.
- That the (S1) band can be hit by a `kappa` grid: it is `\sim10^{-7}` thin in
  `P_0`; use scalar shooting as Astra did.
- That any Green-tangency line meets `\mathcal L` at an admissible first
  maximum. Every measured point violates (N1) below by 60–99.99%.
- That the `601/136136` ratio bound has slack; it is sharp at the box
  corner. Any strengthening of the thresholds must not reuse it loosely.
- That six alternating original signs, a perturbation arc, and a validated
  return map are anything but separate, still-untriggered gates.

## VERIFIED NUMERICAL FACTS (Claude, direct quadrature, no ODE)

- Reconstruction versus original area integrals: worst relative discrepancy
  `4.6\cdot10^{-11}` over random coefficients at `kappa=1.15,1.7,3,9`.
- Astra's eight tuned shots: `Z(p_1)` reproduced to eight digits;
  `\Phi(tau_1)/|Y_0|\in[-0.99998,-0.9686]`.
- Corner point at `kappa=3,5,8,8.5,9,16,30,100,10^3,10^4`:
  `P_*(1;kappa)=-5.8e-4,-1.9e-4,-8.9e-6,+6.8e-6,+2.0e-5,+1.0e-4,+1.2e-4,+8.5e-5,+2.4e-5,+5.1e-6`;
  `\Phi_*(1;kappa)/|Y_{0,*}|=-0.983,-0.976,-0.968,-0.966,-0.965,-0.951,-0.933,-0.887,-0.770,-0.638`.
- Corner linearization at `kappa=8.5` along `(-643/462,1105/462,1)`:
  `dP(1)/dv=+0.00807`, `d\Phi(1)/dv=-0.01651`.
- Fixed `lambda=0.5`: scaled `P(1-eps)=7.33,9.78,11.78,13.34,14.54,15.50`
  for `eps=10^{-2..-7}` (limit `23.76`); `lambda=2`: `8.72,13.41,16.88,19.47`
  for `eps=10^{-2..-5}` (limit `36.45`). All positive.
- Coarse lobe landscape (ten anchor triples, `kappa=1.2\ldots1000`):
  `\Phi(tau_1)/|Y_0|\le-0.879`, the best being the latest triple
  `(0.999,0.9995,0.9999)` at `kappa=1000`.

## VERIFIED THEOREMS

Everything in "SAFE TO INHERIT" 1–10, plus the following new statement,
proved in CLAUDE_AUDIT_ASTRA_1_3.md section 7.2:

> **(N1)** If the original Q4 integral has five distinct zeros in
> `(1,kappa)`, then its universal point lies in `\mathcal L` and
> `\Phi(tau_1)=Y_0+\int_0^{tau_1}Rcal\,\Omega\,H\,du>0`, where `tau_1` is the
> first primitive root.

Proof: (S2) requires `Z(p_1)=\Phi(p_1)>0` at the first `P` root
`p_1<tau_1`, and `\Phi'=Rcal\,\Omega H>0` on `(0,tau_1)`.

## OPEN GAPS

1. Whether `\Phi(tau_1)<0` holds on all of `\mathcal L\times(0,1)` (Theorem N
   below). Numerically true with margin `\ge5\%` everywhere tested except in
   the double limit `kappa\to\infty`, all roots `\to1`, where it tends to
   equality from below while (S1) fails at leading order in every regime
   (`kappa=eps^{-theta}`, all `theta`).
2. A rigorous, uniform version of the heuristic regime analysis for
   `theta\ne1` (Claude's audit, section 7.4, item 3).
3. Nothing else in the Q4 architecture is open at the level of proved
   statements; no five-zero candidate, ordinary fold, or realization exists.

## THE FOURTH STRIKE: prove Theorem N instead of searching `Y=Y'=0`

> **Theorem N (target).** For all `kappa>1` and all `(A,B,eta)\in\mathcal L`:
> `\Phi(tau_1)<0`.

Consequence via (N1) and the multiplicity-five bound: no nonzero Q4 integral
has five distinct zeros in the open annulus; `Z_{distinct}(I)\le4`. This is
a strict improvement of Zhao's Theorem 1 and closes Attack 1 negatively.

### Exact variables

- Anchor chart `(y_1,y_2,y_3)\in\Delta`, `5/11<y_1` (smaller `y_1` is already
  excluded), coefficients by formula (A) of Q4_LOBE_REGION.md; `tau_1=y_1`.
- Lift parameter `a=1-1/kappa\in(2593/21636,1)`.
- Everything is explicit: `\Phi(y_1)=Y_0(A,B,eta)+\int_0^{y_1}Rcal_a\Omega_a H\,du`
  with `H=(A-1)K_0+BK_1-eta K_2+K_3` and the closed-form `K_j`. For fixed
  `(a,y_1)`, `\Phi(y_1)` is an **affine functional** of `(A,B,eta)`:
  `\Phi=c_0(a,y_1)+c_A A+c_B B+c_eta\,eta` with
  `c_A=\frac{3\cdot1326}{1361360}+\int_0^{y_1}Rcal\,\Omega K_0`,
  `c_B=\frac{3\cdot864}{1361360}+\int_0^{y_1}Rcal\,\Omega K_1`,
  `c_eta=-\frac{3\cdot2431}{1361360}-\int_0^{y_1}Rcal\,\Omega K_2`,
  `c_0=-\frac{306}{1361360}+\int_0^{y_1}Rcal\,\Omega K_3`.

### Reduction of the search dimension

Because `\Phi` is affine in the coefficients for fixed `(a,y_1)`, its
supremum over the fibre `\{T(y_1,y_2,y_3):y_1<y_2<y_3<1\}` is controlled by
the fibre's closure, which is a two-dimensional analytic surface with
explicit boundary events (anchor collisions and `y_3\to1`, Q4_LOBE_REGION
section 5). The problem is therefore a **two-parameter family**
`(a,y_1)` of two-dimensional linear-programming-type bounds. Two exact
facts already cut it further:

- On `(0,y_1)`, `H>0`, and `H\le H^{max}(u)` where the majorant can be
  taken from the ECT structure: `q` on `(0,x_1)` is positive with `x_1<y_1`.
- `Rcal\,\Omega=y_2/(1152u^2(1-u)^{3/2}(1-au)^{3/2})` is increasing in `a`
  for fixed `u` near the loop; the worst case is expected at `a\to1`, where
  `Rcal_1\Omega_1=\frac{3}{2}\frac{(1-u)^{-2/3}-1}{1152u^2(1-u)^{13/6}}`.

### Proof architecture (three parts)

1. **Compact part.** For `a\in[2593/21636,1-delta]` and `y_3\le1-delta`,
   certify `\Phi(y_1)\le-c<0` by interval arithmetic on the affine
   functional over the anchor cell. Ingredients: interval `K_j` (exact
   series with the proved tail bound, already in `q4_lobe_certificate.py`),
   interval `Rcal\,\Omega` (elementary hyperbolics), and the explicit
   bounding box of `\mathcal L`. Expected margin: `\ge0.05|Y_0|` away from the
   corner, `\ge2\cdot10^{-4}` in absolute terms.
2. **Loop-approach part** (`y_3\to1`, `a` bounded away from `1`). Use the
   affine endpoint expansion of Q4_THRESHOLD_PATH section 4 for `H` and the
   fixed-`kappa` kernel `\Omega\sim\Omega_0(kappa)/(1-u)` to show
   `\Phi(y_1)=\Phi_*(1;kappa)+O(eps L^2)` with `\Phi_*(1;kappa)\le-0.63|Y_{0,*}|`
   for `kappa\le10^4` (to be certified) and the `O(eps L^2)` coefficient
   negative in the universal direction (`d\Phi/dv<0`, measured `-0.0165` at
   `kappa=8.5`; certify its sign on a `kappa` interval).
3. **Double-limit part** (`a\to1` and `y_3\to1`). This is the only place
   where `\Phi\to0`. Use the two exact cancellations (G2)–(G3): write
   `\Phi(y_1)=\Phi_*(1;kappa)+[\Phi(y_1)-\Phi_*(1;kappa)]` and show both
   brackets are negative: the first by the `kappa^{-1/6}`-type expansion of
   `\Phi_*(1;kappa)` around the exact zero at `a=1` (compute the leading
   coefficient in closed form from the loop-layer kernel, as Astra did for
   `P` in (G5)), the second by (G4): on `(0,tau_1)` the inner perturbation
   `2\pi\,\delta H(1-eps c)\to eps^2L[e-Vc+Qc\log c]<0` for `c\ge1`, and
   `Rcal\,\Omega>0`. Degenerating anchor ratios keep every sign.

If part 3 resists, the weaker but still decisive statement is: Theorem N
on `\{a\le1-delta\}\cup\{y_3\le1-delta\}` for an explicit `delta`, plus
Astra's (G5) and its `theta\ne1` extensions on the complement. That already
excludes five zeros outside an explicit corner box and turns the remaining
problem into a single two-parameter asymptotic estimate.

### What would constitute each outcome

- **Numerical lead (for the exclusion):** an interval-certified negative
  upper bound of `\Phi(y_1)` on a compact anchor-and-`a` box.
- **Rigorous intersection certificate (if the theorem is false):** a point
  with `\Phi(tau_1)>0` certified by outward rounding, plus (S1) with the
  four `P` signs certified. Nothing observed suggests this exists.
- **Five-zero candidate:** unchanged from Q4_CERTIFICATE_PLAN.md; (N1) is a
  mandatory pre-screen before any shooting.
- **Five-zero proof:** unchanged; six certified alternating original signs
  plus the multiplicity bound.

### Stop rules for Strike #4

- Stop **A**: Theorem N proved in full. Q4 closed; Attack 1 done.
- Stop **B**: Theorem N proved outside an explicit corner box, with the
  corner reduced to one two-parameter asymptotic estimate stated exactly.
- Stop **C**: a certified point with `\Phi(tau_1)>0` and (S1). Then, and only
  then, resume the `Y=Y'=0` construction at that point.
- Do not run coefficient sweeps; every object above is a one-dimensional
  integral of explicit functions over an explicitly parametrized cell.
