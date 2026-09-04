# Astra Strike 5: an exact reduction; no four-zero certificate

2026-09-04. Canonical base `7def4d5`; inherits the accepted Theorem N.
Work resumed from the local pause checkpoint `4ee2e84`.

**The requested four-zero exclusion or counterexample is not established.**
The new result is a proved two-anchor reduction theorem, with complete
center-sign classification and an exact criterion for the first positive
Green maximum. The global distinct bound remains four.

[Q4_TWO_ROOT_REDUCTION.md](Q4_TWO_ROOT_REDUCTION.md) contains the combined
theorem and proof. Supporting independent derivations are
[q4/notes_fifth_two_anchor.md](q4/notes_fifth_two_anchor.md),
[q4/notes_fifth_sign_chain.md](q4/notes_fifth_sign_chain.md), and
[q4/notes_fifth_green.md](q4/notes_fifth_green.md).

## The proved progress

For fixed primitive roots `r<s`, the full normalized fibre is exactly
`H_lambda=B+lambda V`, where `B=H_{r,s,1}` and V is the `K2`-normalized
two-root primitive in the lower three-dimensional space. The quotient
`B/V` is strictly decreasing, which classifies both two-root branches
and their center and double-root boundaries.

The center sign is **mixed**. With `v=Y0(V)>0`,

\[
Y_0(\lambda)=Y_B+\lambda v,\qquad
\lambda_c=-Y_B/v,
\]
\[
\frac{231}{50312}<v<\frac3{616},\qquad
0<\lambda_c<\eta_B.
\]

Four original zeros are excluded on the negative-orientation branch,
at and above the center-sign transition, in the entire `beta1=0` sector,
on `H(1)=0`, and in all multiple-primitive-root cases. The exact two-root
counterparts of all three shooting gates are written out, including
nonsimple original zeros and finite endpoint momentum.

The only possible normalized coefficients have `0<lambda<lambda_c`,
primitive signs `+,-,+`, `H(1)>0`, and `Y0<0`. Put

\[
C=B+\lambda_c V,\qquad
\mathcal K(r)=P_B(r)Z_C(r)-Z_B(r)P_C(r).
\]

Here `P_C,Z_C<0` before the first anchor. A complete additional exclusion
holds for every such fibre with `P_B(r)<=0`. If `P_B(r)>0`, a positive
first Green maximum exists somewhere on the coefficient interval
**if and only if** `K(r)>0`. Every fibre with `K(r)<=0` is therefore
excluded as well. The equality cases are proved.

The determinant criterion follows from the exact identity
`K'=Omega(C Z_B-B Z_C)`: every zero while `P_B>0` crosses upwards.
The first maximum of a mixture has the same sign as K at its first
momentum zero. This removes the free coefficient from the first-maximum
question, leaving the three variables `a,r,s`.

## The precise remaining gap

No uniform sign of K has been proved on `P_B(r)>0`, and no positive K
point with all later gates has been certified. Even a positive K would
only satisfy the first Green-height test. It would not provide four
original zeros without the later Green extrema and the four X signs.

Directly extending Theorem N as `Phi(r)<0` would be incorrect: Phi is
strictly positive at the center-zero endpoint C, and hence also at some
nearby fibre points with `Y0<0`. These points can fail the momentum gate.
The proof preserves that distinction.

Two exact tools for the remaining original gates are recorded separately:
[the endpoint note](q4/notes_fifth_loop_gate.md) derives the original loop
functional as an elementary rational integral, and
[the original-height note](q4/notes_fifth_original_height.md) derives a
positive-kernel functional equal to X at a Y zero. Neither has been
assigned an unproved sign on the remaining interval.

Small exact replays check the rational center bounds, determinant
identities, and original endpoint algebra. The analytic proofs were
independently cross-checked for signs, endpoint limits, and multiplicity
cases. No coefficient sweeps, tuned shots, searches for `Y=Y'=0`, corner
asymptotics, or work in Claude's other lanes was performed. No Claude
file was modified.

The three exact replays all passed:

```
python q4/check_fifth_two_anchor.py
python q4/q4_fifth_green_checks.py
python q4/check_fifth_loop_gate.py
```

Their frozen output is
[q4/data/fifth_exact_checks.txt](q4/data/fifth_exact_checks.txt).
The determinant replay also checks the original-height derivative and
positive-kernel algebra. These identities do not substitute for the
analytic sign proofs or certify any proposed parameter point.

## Requested report

```
FOUR-ZERO OUTSIDE-LOBE CASE: OPEN
Q4 GLOBAL DISTINCT BOUND: 4
Y0 SIGN ON TWO-ROOT REGION: MIXED
NEW THEOREM: YES
```

The strongest next task is to determine the sign of
`K(r)=P_B(r)Phi_C(r)-Phi_B(r)P_C(r)` on the remaining baselines
`P_B(r)>0`. A global nonpositive sign would complete the sharp bound
three. A certified positive sign would identify the exact coefficient
interval on which to test the later Green heights and original X gates;
only five alternating exact signs of the original integral would certify
the requested four-zero counterexample.
