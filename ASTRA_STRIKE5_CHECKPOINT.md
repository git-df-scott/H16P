# Astra Strike 5 — paused checkpoint

Saved 2026-09-04 at the user's request before an app refresh.
**Work is paused. This is a progress checkpoint, not a completed proof or
a counterexample certificate.**

Historical checkpoint: work subsequently resumed. The resulting proved
reductions and remaining open case are recorded in
[ASTRA_FIFTH_STRIKE.md](ASTRA_FIFTH_STRIKE.md).

## Repository and accepted state

- Repository: `/Users/scottg/Documents/Codex/2026-09-04/he/work/H16P`.
- Canonical base: `7def4d5`, fetched and fast-forwarded at the start of
  Strike 5. The working tree was clean before this checkpoint.
- Read in the requested order: `CLAUDE_AUDIT_ASTRA_4.md`,
  `FASTRA_H16_HANDOFF_5.md`, `CLAUDE_THOUGHT_SESSION.md`,
  `Q4_THEOREM_N.md`.
- Theorem N and its lemmas are accepted without re-proof. Global distinct
  original bound: four. Bound on strict lobe region: three. Five distinct
  original zeros are excluded. No four-zero example has been found.
- Strike 4's completed proof and exact checks are committed and were
  pushed; its artifacts also remain in `outputs/strike4`.

## Current user assignment and limits

The user's final clarification makes this a **single task**: prove or
refute the remaining four-distinct-original-zero case outside the strict
lobe region, concentrating on primitives with two interior zeros. The
earlier ambiguous request to do both tasks was superseded by the repeated
single-task instruction. Do not undertake the conditional secondary
loop/focus task.

Derive the two-anchor chart, determine the center sign, write the exact
two-root counterparts of (S1)–(S3), and attack the first Green maximum.
An actual counterexample needs five alternating exact original-integral
signs and the inherited bound; an auxiliary sign is not a certificate.

No coefficient sweeps, old lobe shots, tuned shots, searches for
`Y=Y'=0`, or corner asymptotics. Do not enter Attack 2, Attack 3,
reversible-center, or boundary-graphic work. Claude owns those lanes.
Do not modify Claude's files. Small exact replays only, one numerical
thread, reduced priority, and a ten-second CPU fuse. No numerical search
or exact replay has been run in Strike 5 so far.

On completion, commit proofs, notes, and exact checks to main and report:

```
FOUR-ZERO OUTSIDE-LOBE CASE: EXCLUDED / FOUND / OPEN
Q4 GLOBAL DISTINCT BOUND: 3 / 4
Y0 SIGN ON TWO-ROOT REGION: NEGATIVE / MIXED / POSITIVE
NEW THEOREM: YES / NO
```

Then give one paragraph on the strongest next task. Do not present this
checkpoint as completion of the research task.

## Work division at pause

All three agents are interrupted/paused, with no numerical jobs launched.
They had reported deductions in messages but had not yet written their
assigned notes. Resume them with a follow-up task when the user resumes.

- `zero_geometry`: two-anchor chart and center sign;
  planned file `q4/notes_fifth_two_anchor.md`.
- `pf_structure`: full sign chain, endpoint and multiplicity cases;
  planned file `q4/notes_fifth_sign_chain.md`.
- `audit_filters`: first-Green-maximum obstruction on the anchor fibres;
  planned file `q4/notes_fifth_green.md`.
- Root: independent analytic work on fibre comparisons and first-peak
  reduction; integration and eventual verification.

## Deductions to preserve and finish checking

These are new working deductions. They are not yet a final independently
checked Strike 5 theorem. The agents should write their complete proofs
and review them before publication.

### Two-anchor chart and mixed center sign

Fix `0<r<s<1`. Write `x=K1/K0`, `m=K2/K0`. Let the secant line through
`(x(r),m(r))` and `(x(s),m(s))` be `b+S*x`, and define

```
V = K2 - b K0 - S K1.
```

It has the two roots `r,s`, coefficient of `K2` equal to one, and signs
`+,-,+`. Its transported center functional `ell(V)` is positive by (N3).
Let `B=H_{r,s,1}` denote the `K3`-normalized endpoint primitive with
roots `r,s,1`, obtained as the third anchor tends to one. It has positive
first and last interior lobes and `Y_B<0` by the closure comparison from N.

Every `K3=1` primitive with the two anchors is

```
H_lambda = B + lambda V,
lambda = H_lambda(1)/V(1),
Y0(lambda) = Y_B + lambda ell(V),
eta(lambda) = eta_B - lambda.
```

The removable quotient `R=B/V` is strictly decreasing from `R(0)>0`
to `R(1)=0`. Repeated levels, or a critical point, would produce four
primitive zeros counted with multiplicity, contrary to the inherited
ECT bound. The two-simple-root branches are:

- `lambda>=0`: signs `+,-,+` (including the endpoint-zero boundary at zero).
- `lambda<=-R(0)`: signs `-,+,-` (with the center degeneracy at equality).

The intermediate interval has a third root, except at the two ordinary
double-contact levels. Define `lambda_c=-Y_B/ell(V)>0`.
The negative branch has `Y0<0`; the positive branch is mixed, with its
single center sign transition at `lambda_c`. The only possible four-zero
branch left by orientation is

```
0 < lambda < lambda_c, H signs +,-,+, H(1)>0, Y0<0.
```

### Stronger center-functional chord bounds

Root and geometry independently used the following improvement on the
coarse (N3) estimate. For an increasing convex graph, any interior secant
line evaluated at an interior abscissa lies below the chord joining the
two full endpoints. Apply this at `X=144/221` to `m(x)` and to
`n(x)=K3/K0`, which is also strictly convex because `tM(t)` is strictly
convex. Geometry reports the exact bounds

```
ell(V) > 231/50312,
Y0(B + eta_B V) > 27/12578.
```

The second primitive has `eta=0`. Thus `lambda_c<eta_B`, and the
center-zero primitive `C=B+lambda_c V` has `eta_C>0`.
Consequently `P_C(0)=-eta_C/192<0`; since `H_C>0` on `(0,r)`,
`P_C(t)<0` and `Z_C(t)<0` on `(0,r]`, with `Z_C(0)=0`.
The rational constants still need their small exact replay saved.

### Complete crossing conditions reported by the sign-chain agent

For `H` with signs `sigma,-sigma,sigma`, `H(1)!=0`, four original
zeros require the center sign `sign Y0=-sigma`, then

```
sign P0 = sigma,
sign P(r) = -sigma,
sign P(s) = sigma.
```

These give three P roots `p1<p2<p3`. Four crossings of `Z` require its
three extremum signs `sigma,-sigma,sigma`; its endpoint signs are both
`-sigma`. If the resulting Y roots are `v1<v2<v3<v4`, four original
crossings require

```
sign X(v2) = sigma,
sign X(v3) = -sigma,
sign X(v4) = sigma,
sign X(1) = -sigma.
```

For `H(1)=0`, `P(1)` is finite and the necessary additional final sign is
`sigma P(1)<0`. Cases `Y0=0` and `beta1=0` are excluded by the sign
argument (for beta1=0, (N3) makes Y0 share the initial H sign).
Multiple H zeros reduce the number of H sign changes to at most one;
the sign-chain proof must record the resulting bound carefully. The
anchored argument bounding distinct X zeros by Y sign changes remains
available and handles nonsimple original zeros.

### Why direct extension of Theorem N does not work

On the positive two-anchor fibre, `Phi_lambda(r)` increases affinely.
At `lambda=lambda_c`, its initial term is zero and its first-lobe
integral is strictly positive. Therefore `Phi(r)>0` already occurs with
`Y0<0` for some fibre points. This does not certify any original zeros:
the required P crossing conditions can fail. Do not claim Phi negativity
on the entire dangerous fibre.

Likewise N gives `Phi_B(t)<0` before r, not `Z_B(t)<0` when `P_B(t)>0`.
A convex-combination argument using B and C is incomplete without that
distinction. The normalized ratio `Z_lambda/|Y0(lambda)|` decreases with
lambda before r, but this alone does not settle the case.

### Latest first-peak reduction from the Green agent

Use the baseline B and center-zero endpoint C above. Set

```
K(t) = P_B(t) Z_C(t) - Z_B(t) P_C(t).
```

The quotient `B/C` is positive before r. Exact differentiation gives

```
K' = Omega C [ Z_B - (B/C) Z_C ].
```

At any K zero with `P_B>0`, the bracket equals
`[P_B-(B/C)P_C]*(Z_C/P_C)>0`. Thus K can cross zero only upwards while
`P_B>0`. If B has a first P zero `p_B<=r`, N gives
`Z_B(p_B)=Phi_B(p_B)<0`, hence `K(p_B)<0`. It follows that K is negative
before p_B. Every intermediate-fibre first P zero occurs there, so its
first Z maximum is negative. This excludes the branch with

```
P_B(r) <= 0.
```

The remaining baseline case is `P_B(r)>0`. This is the immediate next
analytic target, together with completing the endpoint and degeneracy
arguments. No claim of full exclusion or counterexample has been made.

## Resume instructions

Wait for the user to resume. Then restore the three bounded agent tasks,
have them save and cross-check their deductions, and attack the remaining
`P_B(r)>0` case. Keep the global bound at four unless the missing argument
or a certified example is actually obtained. This checkpoint contains no
credentials; authentication is not needed merely to resume the proof.
