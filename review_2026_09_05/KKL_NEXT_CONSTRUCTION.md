# KKL: the next construction test after the 206-call checkpoint

Independent Astra review, 2026-09-05. Baseline: main `6048ed8`,
`STRIKE5_PRECURSOR.md`, council §9, and the committed KKL derivations and
return records. This review uses those records and elementary return-map
arguments. It launches no return integrations, parameter searches, or
symbolic experiments. It changes no earlier result.

**Recommendation.** Within the present assignment, spend at most one
64-call pilot on a newly seeded stationary-return curve with the remote
cycle's section coordinate held fixed. Do not spend the remaining budget
merely driving the old stable multiplier toward one. If the construction
scope is reconsidered, the most useful reframe is to keep the known four
incumbent cycles in the same fixed-coefficient KKL family and search for
an additional finite pair. Six cycles satisfy the objective too. That
second proposal requires changing the current beta/K restrictions; it is
not authorized continuation of the present precursor search.

Neither recommendation comes with an existence argument for the missing
pair. The current evidence excludes neither the assigned box nor the
larger KKL family.

## 1. What the completed continuation actually tested

The fixed family is

\[
 \dot x=(1+x)y+x^2,\qquad
 \dot y=-10x^2+\frac{11}{5}xy+cy^2+\alpha x+\beta y,
 \qquad K=-\alpha(11c/5-1)-42.
\]

The assigned precursor has beta zero, `1/2 <= c <= 3/2`,
`-200 <= alpha <= -10`, and `K >= 1/64`. At one field it needs three
origin cycles S/U/S and a remote U surrounding an attracting focus.
The final small negative-beta Hopf step is conditional on that precursor.
The section bounds and all return-validity gates remain in force.

There are 206 charged calls, leaving 3890 of the 4096 allowance. The main
path continued one stable origin cycle and the remote U. Its origin
multiplier rose from approximately 0.810 to 0.976. The last admissible
field was `c=9301/10000`, `K=1/64`; its remote section was approximately
`-1048286.51`, very close to the prescribed `-2^20` cutoff. The estimated
crossing of that cutoff is an experimental boundary, not a proof that
the remote cycle disappears or reaches infinity.

At that field the saved stationary-return calculation found

\[
 r_*\simeq28.17411716,\qquad D(r_*)\simeq0.24269325,
 \qquad D_r(r_*)=0,\qquad D_{rr}(r_*)\simeq-0.00112207.
\]

This is a positive maximum. It is not a cycle or a saddle-node of cycles.
The later four inward remote returns at the single `c=1001/1000` control
did not bracket a remote U. They do not establish its absence, even at
that one field.

## 2. The missing return geometry

Write `D=R-r` for the downward full return. At beta zero and positive K,
the proved small-return expansion is

\[
 D(r)=\frac{\pi K}{4(-\alpha)^{3/2}}r^3+O(r^4)>0
 \quad\hbox{for sufficiently small positive }r.
\]

If three simple origin roots are ordered S/U/S, then D crosses down,
up, down. Between the center and those roots there must be at least
three distinct finite stationary points: a positive maximum before the
first S, a negative minimum between S and U, and a positive maximum
between U and the last S. This follows directly from the signs and
Rolle's theorem on their common return domain. It does not assume that
these are the only stationary points or roots.

The pilot has identified one of the required types, a positive maximum.
It has not seeded the additional negative minimum. Two possible
finite-amplitude births would accomplish the desired count:

* An interior minimum crosses from positive to negative while the outer
  S remains: an additional S/U pair appears inside it.
* An exterior maximum crosses from negative to positive outside a
  retained S: an additional U/S pair appears outside it.

An ordinary fold requires `D=D_r=0`, `D_rr != 0`, and a nonzero
transverse parameter derivative. A creation of *stationary points*
instead has `D_r=D_rr=0`; unless D also vanishes, that is not a fold of
cycles. Tracking the latter event can seed the previously missing
minimum/maximum branches. Their heights must then be followed to the
actual cycle-producing event.

Following the old root toward multiplier one was a legitimate pilot,
but its first fold could remove that root. Conversely, the name or
historical identity of the “old S” is not a mathematical constraint:
root sheets can turn, reconnect, or exchange which cycle is innermost.
The real requirement is simultaneous S/U/S plus remote U at the final
precursor. No argument here proves that continuation from the known
root cannot reach it, or that every useful fold is connected to that root.

The saved second derivatives also explain why continuing directly from
the last corner is unattractive. On a simple stationary branch, its
height h satisfies

\[
 h_c\big|_K=R_c-R_\alpha\frac{11\alpha}{11c-5},
 \qquad h_K\big|_c=-\frac{5R_\alpha}{11c-5}.
\]

Using the frozen derivatives gives approximately `h_c|K = -4.71` and
`h_K|c = 1.74`. At the limiting corner where the remote cap is active,
the saved remote-radius derivatives impose approximately
`59.73 dc + 1.875 dK <= 0`, while the K margin imposes `dK >= 0`.
These first-order inequalities give no inward direction lowering this
maximum's height. The actual last field has a small amount of cap slack;
this is numerical directional evidence, not an exact local optimum or
a global obstruction.

## 3. A concrete next pilot within the assigned scope

Use a remote section coordinate as a continuation constraint instead of
allowing its radius to consume the remaining margin:

\[
 r_R=-2^{15},\qquad c_0=33/40.
\]

At this same c, the saved endpoint fields `K=1/64` and `K=6/5` have
remote U coordinates about `-26821.86` and `-36339.12`. They therefore
motivate a bracket for a field satisfying

\[
 D_R(-2^{15};c_0,K)=0.
\]

The bracket and remote hyperbolicity must be checked with valid returns;
the two saved coordinates alone are not a parameter-space certificate.
Both endpoints and the intended bracket lie within the assigned box.
No negative-K bridge or coefficient sweep is required.

At that field seed a positive origin stationary maximum before the
known stable origin cycle. Existence of at least one such maximum follows
from positive D near zero and the stable crossing; the numerical bracket
must still be obtained. Continue the coupled equations

\[
 D_R(-2^{15};c,K)=0,\qquad D_r(r;c,K)=0
\]

with c as the local parameter and `(K,r)` as unknowns. The required
Jacobian uses derivatives already available in the reviewed return
evaluator. Choose the first local direction from the calculated
stationary-height derivative; then use ordinary continuation corrections,
not a coefficient optimizer. A singular parameterization is an event to
record or resolve by pseudo-arclength, not permission for unbounded
retries.

The pilot's new output is the stationary height and curvature along a
curve with a retained remote orbit. Check the stable origin control at
accepted fields. Record any newly bracketed stationary minimum, loss of
curvature, or fold, and test the common-field root signs immediately if
one appears. A maximum becoming flatter or the multiplier approaching
one is progress information, not the success trigger.

Hard pilot stop: 64 additional charged evaluations or 16 accepted steps,
whichever comes first. Failed returns, derivative returns, corrections,
and control replays all count. Preserve the one-thread/ten-CPU-second
fuses. Stop at the K margin, section or itinerary failure, loss of the
required remote focus/cycle, or an unresolved singular event; do not
cross the next infinity stratum without its existing gate. If no new
stationary branch is seeded, report that this particular constrained
curve did not supply one. Do not spend the entire remaining allowance
extending the same negative result.

This curve can miss disconnected stationary sheets. Its advantage is
specific: it tests previously unrecorded critical-point geometry while
removing the particular remote-amplitude obstruction of the old path.
It does not assert that fixed remote amplitude is the optimal or a
sufficient family for finding five cycles.

## 4. What the incumbent negative control can and cannot seed

The exact incumbent is

\[
 (c,\alpha,\beta)=(7/10,-363889/5000,3/2000),
 \qquad K=-674997/250000<0.
\]

It supplies three detected origin cycles S/U/S and a remote U. The same
shape at beta zero supplies two detected origin roots U/S. These are
authorized controls; neither is an admissible positive-K precursor.

On the beta-zero negative control, a negative minimum before the U and
a positive maximum between U and S follow from the return signs, on
the common domain. Thus a fixed-field stationary calculation there is
a meaningful control for a critical-point evaluator. It is not a seed
already in the required box.

At fixed c with the proved quintic coefficient Delta positive, the
small-return unfolding gives, as K tends to zero from below,

\[
 r_U^2\sim-30K/\Delta,\qquad
 r_{\mathrm{min}}^2\sim-18K/\Delta.
\]

The second relation follows by differentiating the cubic-plus-quintic
expansion. Both the small unstable root and its negative minimum
collapse at K zero. It is not justified to identify the finite control
root with this local branch without following it. But **if** one follows
this local branch, it cannot be carried into positive K as the missing
finite minimum.

The special simultaneous zero `K=J(c)=0` is a proved double center,
not a third-order focus. Its local beta-zero unfolding has at most one
small isolated cycle and no small fold. This closes the proposed local
shortcut, not an entire finite-amplitude fold sheet. A finite fold could
remain wholly in positive K, be disconnected from the controls, or meet
some other boundary. No transport through the center, and no blanket
finite-annulus exclusion, follows from the local result.

New continuation through negative-K *parameter values* would require a
scope amendment. Reading/reproducing a fixed negative control, or using
its data to initialize an admissible predictor, does not grant that
amendment.

## 5. The strongest scope amendment: retain the incumbent four

The objective is **at least five**, not exactly five. If a finite
saddle-node pair is created around the origin while all four incumbent
cycles persist, the configuration becomes `(5,1)`: six limit cycles.
This is already a counterexample to `H(2)=4`. Its distribution is not
excluded by the `(n,0)`/`(n,1)` conclusion of
[Zegeling, Theorem 1.2](https://doi.org/10.1515/anona-2024-0012).
The theorem supplies no upper bound on n. This is a compatibility
statement, not evidence that six cycles occur in KKL.

The strongest bounded reframe by proximity to an existing control is:
keep the coefficients `-10` and `11/5` fixed; start at the exact finite-beta
incumbent; continue its three origin roots, remote root, and intervening
stationary branches in `(c,alpha,beta)`; look for an *additional* finite
pair while retaining the four roots at the candidate. This begins with
the richer three-stationary-point geometry, instead of first replacing
the four-cycle field by a two-cycle positive-K field.

This proposal must explicitly allow the incumbent's beta-positive,
K-negative neighborhood. The c/alpha box and section bounds can initially
stay as they are. It does **not** require freeing `11/5` or `-10` as a
first step. Any later widening of those coefficients would be another
decision, unsupported by the present 206 calls.

A sensible first block is again 64 calls: reconstruct stationary
brackets at the incumbent, calculate their transverse sensitivities,
and select one continuation direction compatible with the four
hyperbolic controls. Continuing an existing S/U annihilation boundary
is useful geometric information but does not by itself create a fifth
or sixth cycle. If no additional branch is identified, this is still
only a seeded-sheet result.

Why consider this reframe? It avoids treating the unproved beta-zero
four-finite-cycle precursor as a necessary condition for the original
problem. Why remain cautious? Six cycles is a stronger outcome than
five; having a closer known cycle count is not a probability estimate
for a new pair. The incumbent's two inner multipliers are also close to
one, so certified persistence and fold calculations may be harder than
the already explored positive-K controls.

There is a second, narrower amendment: start with the known positive-K
S plus remote U, take a sufficiently small negative beta to create the
inner U, and then search for a finite pair while those three persist.
That targets five directly and can retain the existing K margin and
final negative-beta range. It changes the prescribed **order** of the
construction, since beta is presently reserved for use after the
four-cycle precursor. The new small cycle must satisfy the lower section
bound and common persistence estimates; an arbitrary negative beta is
not automatically a valid control.

Both amendments address a logical restriction: a five-cycle field at
finite negative beta need not admit a deformation to beta zero that
preserves four finite cycles. The conditional Hopf theorem does require
that persistence for its own construction, but there is no established
theorem making it necessary for every possible counterexample.

## 6. Remote-U pruning and certificate cost

The proved remote Hopf sign prevents a *small* remote U from being born
at positive-K remote Hopf points with `J<=0`. A global exclusion of finite
remote U cycles in that attracting-focus region would prune part of the
box. It would not supply the missing origin pair or dispose of the
unexplored `J>0` finite-amplitude geometry. The four inward c>1 controls
are not such a proof. A short exact Dulac/return-comparison argument is
a useful parallel task if an actual sign mechanism is available;
otherwise it is a detour from construction.

There is also a concrete observation gap before interpreting remote
negative controls: y=0 does not pass through the remote focus, and no
proved result forces every candidate remote U to cross that line. The
independent [coverage review, §4](COVERAGE_AUDIT.md) derives a complete
beta-zero cycle section on the velocity nullcline
`(r,-r^2/(1+r))`, using its left branch `r<x*` for the remote nest.
Changing the event evaluator and translating the physical radius cap
are necessary before using it. This is an analytic repair, not a new
observed cycle or permission to evade the old radius restriction. It
strengthens the reason not to turn four inward y=0 returns into a
global remote-U exclusion.

For the original precursor route, a numerical S/U/S plus remote-U
candidate still needs a common negative-beta choice and five disjoint
original-return certificates. For the incumbent-plus-pair route, a
candidate has six gates; proving any five isolated cycles at the same
exact field suffices for the objective, although retaining all six
certificates would document the proposed mechanism. No cycle count can
be combined across different fields. The new pair near its fold is
poorly conditioned, so certification should use a transverse parameter
on the cycle-present side with resolved separation, not the double
cycle itself.

The completed pilot was useful: it identified a genuine constrained
boundary and eliminated the local double-center shortcut. It did not
establish that KKL lacks the required finite pair. The next useful
experimental object is an additional stationary-return branch; the
next useful scope decision is whether the beta-zero precursor remains
the exclusive construction mechanism.

## Reviewed sources

* [Current precursor report](../STRIKE5_PRECURSOR.md) and
  [council §9](../FASTRA_COUNCIL_2026_09_04.md).
* [Stationary control](../kkl/data/stationary_control.json),
  [boundary sensitivities](../kkl/data/fold_direction_at_boundary.json),
  and [single-field later control](../kkl/NEXT_BOUNDED_TEST.md).
* [Liénard/energy identities](../kkl/notes_lienard.md),
  [double-center and local unfolding](../kkl/notes_local_unfolding.md),
  and [remote/infinity gates](../kkl/notes_other_strata.md).
* Primary distribution theorem: André Zegeling, *Nests of limit cycles
  in quadratic systems*, Theorem 1.2, DOI above. The archived article's
  full text was checked for this statement. Its Corollary 6.8 bounds
  cycles around the **other** focus in the stated weak-focus cases;
  it is not used here as a bound on the weak focus's own nest.

Independent bounded review by the return-derivative auditor found no
actionable defect in the sign arguments, derivative conversion, local
scales, or scope distinctions. The stationary pilot still begins on the
old positive maximum; discovering an additional branch remains its
essential unresolved task.
