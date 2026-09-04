# Astra Strike 4: Theorem N proved globally

2026-09-04. Inherits the verified first three strikes and the FASTRA
handoff at `084fae3`. **Stop A is reached for the stated scalar theorem.**

Theorem N is proved by three global comparisons:

\[
Y_0<Y_{0,*}=-3/1232,\qquad
0<H(t)<H_*(t)\quad(0<t<\tau_1),\qquad
0<\mathcal R_a\Omega_a<\mathcal R_1\Omega_1.
\]

The exact inherited corner moment then gives

\[
\boxed{\Phi_a(\tau_1)
<-\int_{\tau_1}^1\mathcal R_1\Omega_1H_*\,dt<0.}
\]

The complete proof is [Q4_THEOREM_N.md](Q4_THEOREM_N.md). It covers every
finite `kappa>1`, every strict anchor triple, all partial loop approaches,
and every coupling or degeneration in the double limit. No asymptotic
uniformity is inferred from finite numerical evidence.

## New proof ingredients

1. The weighted moment curve `x=K1/K0`, `m=K2/K0` is strictly convex,
   with exact maximal slope `1105/462`. The transported center functional
   is positive on every variation with two roots and positive initial
   sign; its exact residual is `25/231` after normalization. Cardinal
   interpolation therefore makes both `Y0` and the primitive before its
   first root strictly increase with each anchor. Sending all anchors
   to one along a fixed-ratio comparison path proves the first two bounds.
2. Set `v_a=y_a Rcal_a/(1-at)^(3/2)`. Its homogeneous equation has
   limiting solution `v_1=(3/2)[(1-t)^(-4/3)-(1-t)^(-2/3)]`, whose residual
   for the finite-lift operator is exactly
   `(1-a)[22-7(1-t)^(2/3)]/[6(1-t)^(7/3)]>0`.
   The positive causal Green kernel gives the third bound globally.
3. The corner integral is exactly `3/1232`. Its integrable positive tail
   supplies strictness even in the limiting equation at `a=1`.

The requested compact part has an explicit analytic certificate. With
`delta=1/64`, throughout
`a in [2593/21636,1-delta]`, `y3<=1-delta`,

\[
\Phi_a(\tau_1)<-395/3784704.
\]

This replaces the proposed interval cover. The global comparison also
replaces both endpoint expansions, including the proposed two-parameter
corner estimate.

## The original-zero consequence needs its stated scope

Five distinct original zeros are impossible globally by the audited
necessary condition (N1). Thus the global distinct interior bound improves
from five to **four**. On the strict lobe region, the failed first Green
maximum and the complete monotone-lobe argument prove at most **three
distinct** original zeros, including nonsimple zeros.

The handoff's appended claim of a global bound of three omits the
coefficient directions outside the lobe region. The canon only places
five original zeros in that region. In particular, when `H` has two
interior zeros, the inherited `Z(I)<=Z(H)+2` bound can still allow four,
and Theorem N does not apply. This is a remaining proof obligation,
not a constructed four-zero example. The global three-zero conjecture
therefore remains open within the results established here.

The five-zero target in Attack 1 is excluded. A claim that all interior
questions are closed would go beyond this proof. No endpoint-born cycle
or Attack 2 work was undertaken.

## Verification and handoff correction

Independent derivations are in
[q4/notes_N_compact.md](q4/notes_N_compact.md),
[q4/notes_N_loop.md](q4/notes_N_loop.md), and
[q4/notes_N_double.md](q4/notes_N_double.md).
The final combined proof was checked for the kernel transformation,
convexity, cardinal signs, strictness, singular integrability, compact
constant, and distinct-zero conclusion.

Two small exact replays pass:

```
python q4/check_N_kernel.py
python q4/q4_N_loop_checks.py
```

Both use reduced priority and a ten-second CPU ceiling. They verify only
the new symbolic and rational identities. The proof of each global sign
is analytic. No old shots, coefficient sweeps, or tangency searches ran.

The new handoff's affine constant is missing `-K0`: the correct term is
`c0=-306/1361360+integral W_a(K3-K0)`. Existing scripts use the complete
primitive correctly. None of the three audited strikes or Claude's
source files was changed.

## Requested report

```
THEOREM N PROVED: YES
COMPACT PART CERTIFIED: YES (delta = 1/64; analytic certificate)
LOOP-APPROACH PART PROVED: YES
DOUBLE-LIMIT PART PROVED: YES
Phi(tau1) > 0 POINT FOUND: NO
Q4 INTERIOR BOUND: 4
Q4 STILL LIVE (interior): YES
```

The last two entries describe the global original family. The bound is
three on the strict lobe region; globally five is now excluded, and the
remaining four-zero possibility lies outside that region.
