# KKL coverage audit: what the recorded experiment did and missed

Independent adverse review, 2026-09-05. Scope: the recovered handoff,
current KKL source, 206-call ledger, path records, checkpoint and council.
This review ran no ODE, numerical experiment, parameter search or new
literature search. Small read-only JSON summaries were used to inspect
the existing ledger. No existing file was modified.

## Finding

The recorded work establishes a functioning numerical return evaluator,
one continued pair of already-known cycles, and several useful exact
local restrictions. It does **not** constitute a substantial search of
the missing fold pair. The main algorithm follows the old stable origin
cycle and its remote companion; no additional stationary-return branch
or fold curve was continued. Reaching the selected radius and K limits
does not repair that mismatch with the handoff's central task.

The current checkpoint correctly says that the full box is open and the
allowance is unexhausted. Those qualifications are necessary. The adverse
conclusion is about how little of the actual construction problem has
been tested, not an accusation that the checkpoint claims a global proof.

## Sources and identification

The controlling recovered packet is
`/Users/scottg/Documents/Codex/2026-09-04/your-auditing-now-for-h16-x20/work/audit-only-publication/ASTRA_FIFTH_STRIKE_HANDOFF.md`,
at audit commit `de39ea78d56208a2a3267b594ce5c117b6b14c1e`.
Its numbered **Bounded numerical task** is used below.

Current-repository sources:

- `STRIKE5_PRECURSOR.md` and `FASTRA_COUNCIL_2026_09_04.md`, especially
  the corrective section 9, state the experiment and its claimed scope.
- `kkl/data/returns.jsonl` is the call record. Call numbers below are its
  `evaluation` values and also its current line numbers; they are
  consecutive from 1 through 206.
- `kkl/data/continuation_events.jsonl` aggregates the individual path
  files and preserves their names and line numbers.
- `kkl/return_map.py`, `kkl/continue_path.py`, `kkl/follow_segment.py`,
  `kkl/pilot.py`, and `kkl/README.md` establish what algorithms exist.
- `kkl/NEXT_BOUNDED_TEST.md` and
  `kkl/data/c_gt_1_remote_controls.json` establish the four-call c>1 test.
- The exact limitations are proved in `kkl/notes_lienard.md`,
  `kkl/notes_local_unfolding.md`, and `kkl/notes_other_strata.md`.
  They must be kept separate from numerical coverage.

## 1. Requested work versus execution

| Handoff obligation | Recorded execution | Coverage judgment |
|---|---|---|
| 1. Reproduce four incumbent returns and stability | Calls 3–6 evaluate the four supplied section coordinates at beta=3/2000 | Numerical control completed; no new interval root certificates |
| 2. Continue beta from the incumbent down to zero | Calls 7–8 evaluate two supplied beta-zero coordinates | **Not performed as continuation**; endpoints were substituted for the requested path |
| 3. Start K>0 from one origin S and remote U | Calls 1–2 reproduce those controls | Completed numerically |
| 3. Continue the root sheet in two shape variables and its fold curves | Calls 12–63 and 66–186 follow prescribed segments of the same two roots | One broken path in the two-parameter plane; **no executed fold-curve continuation** |
| 3–4. Track compactified escape boundaries and continuous return domains | Two section crossings, winding diagnostics, Cartesian guards, and exact infinity types are recorded | No compactified flow/transition implementation or return-domain atlas |
| 5. Establish S/U/S plus remote U at one parameter | Never obtained | Central construction gate remains unmet |
| 6. Validate the precursor, then select beta<0 and certify five | Not attempted | Correctly deferred; inventing the fifth now would double-count an absent precursor |

Only calls 3–6 have nonzero beta, and all four use the same value
3/2000. All other calls have beta=0. There is no intermediate positive-beta
trace in this ledger. Consequently it is valid to say that the static
controls agree with the published Hopf-loss mechanism; it is not valid
to say that this campaign followed that loss or located every event on
the incumbent's beta path. This gap is real but lower priority than the
missing fold experiment: another clean negative control is not itself a
route to the extra pair.

## 2. The largest missed step: a different object was optimized

The handoff seeks an additional S/U pair while the old S and remote U
persist. The actual continuation chooses directions that make the
**old origin cycle's multiplier** approach one. Calls 65 and 199 supply
the corresponding derivative diagnostics. That is a useful local
measurement, but it is not a search for a new pair coexisting with that
cycle. An ordinary fold of the old S would require another cycle and may
annihilate the old one; it does not automatically produce the desired
three-cycle nest.

This distinction appears directly in the code:

- `continue_path.py:correct` is an ordinary root corrector in log radius.
  It explicitly raises `near fold: switch to augmented fold equations`
  when its slope is small. No such augmented solver is implemented there.
- Both path drivers accept the tracked origin only with `R_r<1`, and
  preserve one selected remote root with `R_r>1`.
- `follow_segment.py` follows a specified straight segment in (c,K),
  using the previous roots as predictors. It does not explore the entire
  root surface, branch-switch, or continue through a fold.
- The second derivatives exist in the evaluator, but availability of
  derivatives is not execution of a fold algorithm.

Calls 187–196 make one ten-point origin section profile, **only at the
last field**. Calls 200–202 locate one stationary point of its displacement:
r approximately 28.1741, D approximately +0.242693, and D_rr<0. This lies
inside the known stable root at r approximately 48.6948. It is the positive
hump expected in a simple one-cycle displacement picture. It is not an
additional interior minimum or an exterior maximum demonstrating a new
pair, and it is not D=D_r=0. The report correctly calls it not a fold.

The central implementation gap is therefore still present: **no additional
stationary-return branch has been identified and followed, and no
ordinary cycle-fold curve has been solved and oriented.** Successive
algebraic reviews of the derivative formulas did not establish adequacy
of the discovery method. The exact local work rules out shortcuts; it
does not substitute for this missing finite-amplitude investigation.

## 3. Parameter coverage was one thin path, not a sheet survey

The accepted common-parameter records comprise:

1. Seven points on K=6/5 from c=7/10 to 33/40:
   `constant_K_path.jsonl`, calls 12–63.
2. Four points reducing K to 1/64 at c=33/40:
   `reduce_K.jsonl`, calls 66–90.
3. Points increasing c with K fixed at 1/64:
   `K_margin_to_c09.jsonl`, the repaired c=9/10 point, and
   `K_margin_toward_boundary.jsonl`, calls 91–173.
4. One accepted point just inside the remote section cap:
   `last_inside_boundary.jsonl`, calls 180–186.

There are 24 attempted common-parameter records and 22 accepted ones.
No parameter cell has a certified return-map cover. The available curves
do not investigate the rest of the rectangle, disconnected cycle sheets,
isolas, or roots that are not reached by the old S predictor. Nor do
sampled accepted points certify persistence over the intervals between
them. The distinction is explicit in the checkpoint but must also govern
the strategic interpretation of its negative result.

The selected direction stops at **two imposed conditioning boundaries**:
K=1/64 and remote |r|=2^20. At the last point both reported cycles remain
numerically hyperbolic; their multipliers are approximately 0.97634 and
1.90281. The stop is not an observed fold or destruction of either cycle.
Calls 174–179 locate the remote section-cap crossing numerically. The
original attempt at c=149/160 stopped before an admissible remote
corrector could be evaluated; its origin result remains in the log.

Keeping these limits was compliant with the handoff. Treating them as a
natural frontier of the mathematical mechanism would not be. Within the
existing experiment, another seeded direction must retain room inside
both limits. Any later decision to reduce the K margin or exceed the
physical radius cap should be an explicit new bounded protocol with
appropriate charts; renaming the coordinate does not enlarge the allowed
physical region or prove an escape theorem.

## 4. The remote section has an unproved completeness assumption

The origin section y=0 passes through its focus. The remote focus instead
has y*=s^2/(s-1)>0. A closed curve surrounding that remote focus can remain
entirely above y=0. The known very large remote cycles cross y=0, but no
recorded theorem forces **every possible admissible remote U cycle** to
cross this fixed line.

This is more serious than the finite radius cutoff. The present search
can miss a remote cycle even if its size is moderate and its coefficients
lie comfortably inside the box. This review does not claim that such a
missed U cycle exists. The new Hopf sign theorem blocks its simplest
positive-K local birth, but does not prove that every finite U cycle
intersects y=0.

A coverage repair should use a transverse ray/section through the remote
focus, with the moving equilibrium and all event projections accounted
for, or prove that the target U cycles must hit the existing section.
A section through the focus meets every surrounding Jordan curve, but
its transverse return domains still need to be partitioned; merely
changing the line is not a global return theorem. The known U orbit
provides a useful identity control for transporting between sections.

### An exact, more specific repair: the velocity nullcline

The parent proposed the rational curve

\[
 \sigma(r)=\left(r,-\frac{r^2}{1+r}\right),\qquad P(\sigma(r))=0.
\]

At beta=0, direct substitution gives

\[
 Q(\sigma(r))=-\frac{rW(r)}{(1+r)^2}
             =\frac{rT(r)}{(1+r)^2},\qquad
 \dot P\big|_{P=0}=(1+r)Q(\sigma(r)).
\]

Here W=-T is the cubic in `kkl/notes_lienard.md` (L3). Its proved
positivity on x>-1 makes every right origin turning point (r>0) a
strict maximum of x. On the remote side, write r=-s<-1. The identity
W(-s)=(s-1)^2[m-F(s,c)] and the proved F_s>0 give
Q=s[m-F(s,c)]<0 for r<x*. Thus the branch r<x* consists of strict
minima of x. The opposite remote branch consists of strict maxima.

Every cycle surrounding the corresponding focus has an x extremum on
the selected branch. These intersections are transverse: Q is nonzero
away from the focus and r=-1 is absent. More strongly, in the smooth
coordinates (x,v=P) on either half-plane, the selected branch becomes
a ray from the equilibrium on v=0. All intersections with it have the
same orientation. Transverse intersections of a Jordan curve with a
ray alternate entering and leaving, so there is exactly one such
intersection per surrounding cycle. This supplies a complete cycle
section for this beta-zero geometry, without assuming y=0 intersects
the cycle. It does not prove that every off-cycle starting point returns.

For a transverse full return on this curved section, its scalar
derivative still satisfies

\[
 R_r=\frac{\det(F_0,\sigma'(r))}
             {\det(F_T,\sigma'(R))}
          \exp\!\int_0^T\operatorname{div}F\,dt
     =\frac{Q_0}{Q_T}
          \exp\!\int_0^T\operatorname{div}F\,dt,
\]

because P=0 and the first component of sigma' is one, giving
det(F,sigma')=-Q at both endpoints. Initial tangent vectors and event
projection must nevertheless be changed in the evaluator; the present
y=0 implementation cannot be reused unchanged. This is an analytic
repair proposal only: no new evaluator or return was run. Its x-extremum
coordinate is not the old y=0 section coordinate, so the experimental
radius limits must be translated or explicitly revised, not silently
relaxed by the section change. Near the focus Q tends to zero, requiring
the usual conditioning and local-chart treatment.

For the origin section, the immediate omissions are the finite bounds,
unsampled radial gaps, and possibly disconnected return domains. The
last profile begins at 1/64, whereas the allowed lower radius is 2^-12.
The local no-collapse results have unquantified neighborhoods; they do
not certify that this omitted numerical interval contains no additional
root at the last field. Nor do ten section samples exclude a narrow
negative dip or positive bump between them.

## 5. Nonreturns and itinerary boundaries were not mapped

All 206 logged evaluator results have NUMERICAL_ONLY status. There are
no failed ODE-return records in this selected sample. That does **not**
establish a connected return domain between the samples.

`return_map.py` uses Cartesian integration, horizon 10, an absolute
coordinate guard of 10^7, and the remote x=-1 guard. It records the first
opposite crossing and then the desired downward crossing. These are
reasonable pointwise diagnostics. The source has no compactified
integrator, saddle transition map, boundary continuation or interval
return-domain proof. Exact classification of infinity directions does
not provide those missing orbit maps.

The c=9/10 rejection was a different event: the derivative discrepancy
at call 130 was numerical cancellation, not a nonreturn. Call 131 tightened
tolerance; calls 132–133 checked the determinant reformulation; calls
134–141 successfully replayed the common point. That repair was sound
and the rejection was not improperly promoted to a mathematical obstacle.
The remaining task is to distinguish and retain return-domain boundaries
when they actually arise, rather than interpreting a future timeout,
guard hit or sign jump as a fold or an absence proof.

## 6. The c>1 test was an algebraically motivated guess, not transported geometry

Calls 203–206 test the single field c=1001/1000, m=196/5 at remote
r=-8,-512,-32768,-1048576. All four returns point inward. No origin
controls run at that field because the protocol makes them conditional
on obtaining a remote U bracket.

The changed sign pattern of the quartic multiplier density N is a valid
necessary-amplitude observation. It is not evidence of a stationary
return branch, much less of a remote U orbit. The field was not reached
by continuation of a known periodic orbit through either infinity
transition. Thus it is an algebraically motivated fixed-field seeding
attempt, with no inherited cycle geometry. Calling it a second explored
cycle component would overstate the work.

Four radial values separated by factors as large as 64 leave large gaps,
the near-focus section omission above, and every other c>1 shape open.
Stopping after four of the allowed eight initial remote calls was a valid
bounded decision, but supplies especially weak adverse evidence about
this field. It is not grounds to close c>1 or to conclude that another
unseeded field is the best next use of the budget.

The stronger next obligation is to supply a reason a remote U branch is
present and accessible: transport a known branch with controlled itinerary,
establish an appropriate return/annulus sign bracket, or derive a theorem
that excludes a stated region. The remote Hopf calculation is relevant
adverse evidence: its positive-K small cycle is S on the wrong focus side,
so that familiar local mechanism does not provide the missing U seed.

## 7. Numerical and certification debt

The projected and determinant derivative equations were independently
checked and no algebraic error is alleged here. Numerical validation is
narrower than that statement:

- Calls 9–11 provide an off-root radial finite-difference check; call 64
  and the later control support the radial second derivative.
- No c/alpha finite-difference control of R_c, R_alpha, R_rc or R_ralpha
  appears in the ledger, although those quantities select parameter
  directions and feed predictor steps. Before relying on a delicate
  gradient cancellation, compare those derivatives at a well-conditioned
  common field and at an appropriately challenging remote control.
- The determinant radial equation and the logarithmic-divergence formula
  express the same scalar variational law. Their agreement is useful
  consistency evidence, not an independent validation of the orbit,
  itinerary or all parameter sensitivities.
- The code contains no MPFR interval ODE return verifier. Exact Fraction
  equilibrium gates do not validate periodic orbits. No common parameter
  neighborhood preserving even the starting two cycles was certified.

The lack of five certificates is not itself a failure once no precursor
exists to certify. Nevertheless a robust one-cycle interval benchmark
would expose certification costs before a discovery is urgently awaiting
proof. Final source hashes are available; contemporaneous source hashes
for calls 1–202 are absent and cannot be reconstructed by relabeling the
current engine. The newer calls record provenance fields. Current
documentation discloses this limitation.

## 8. Budget and status: no exhaustion claim is available

The ledger is consecutive and the stated counts agree. The pilot used
64 calls and seven accepted continuation steps, within its 64/16 limits.
The whole current KKL ledger has **206 of 4096 calls**, leaving **3890**;
22 common continuation points are far below the 256-step ceiling.
There is no evidence of a numerical resource ceiling being reached:
recorded evaluator CPU totals are 6.86162 seconds, and subprocess wall
time totals are about 67.07 seconds. These exclude analysis, symbolic
checks, tooling and writing; they are not the cost of the entire project.

No budget violation was found. The failure is one of coverage and method,
not excessive call consumption. This phase consumed about 5% of its call
allowance and then paused after a constrained path and one weak new-field
test. That is a legitimate checkpoint. It is not completion of the full
bounded construction/exclusion task. A whole-box negative would need
analytic bounds or a validated cover of all relevant return domains,
including disconnected sheets; the current machinery does neither.

## 9. Ranked fixes and a credible counterexample workflow

1. **Implement the missing finite-amplitude discovery object.** Work with
   stationary displacement branches D_r=0 and their heights, then solve
   D=D_r=0 with nonzero D_rr and a transverse shape derivative. Add
   pseudo-arclength and explicit branch identity. Preserve the old S and
   remote U while determining the two-cycle side. Do not keep using the
   old S multiplier approaching one as a proxy for the missing pair.
2. **Repair radial and remote-section coverage at a known common field.**
   Start where both old cycles have substantial numerical margins and
   room inside the imposed bounds, for example the original K>0 control,
   rather than at the nearly exhausted remote-radius point. Use bounded
   adaptive C1 section investigation to look for an interior minimum or
   an exterior maximum relevant to an extra pair. Partition continuous
   return domains and transport the remote cycle to the complete
   nullcline section proved in section 4. This is targeted seeding, not an
   unbounded coefficient grid.
3. **Use the two shape controls on that object, not only on the old root.**
   Continue a genuine additional stationary branch and its height within
   a stated small parameter region with the common remote ledger. A
   stationary-branch turning/degeneration event also needs a recorded
   local model; a single Newton extrapolation across a failed geometry
   gate is not a new seed. Unvisited isolas stay explicitly unresolved.
4. **Validate the directional derivatives and prepare one interval
   benchmark.** These are small, concrete implementation checks, not
   another long local-theory detour. Preserve source versions before
   future runs and label controls, candidates and validated results
   separately.
5. **Complete or explicitly retire the missing beta-path control.** A
   short genuine incumbent beta continuation would close handoff item 2
   and improve branch-identity handling. It should not displace the
   missing-pair task or be presented as new support for a fifth cycle.
6. **Re-enter other infinity strata only with a seed or an explicit
   transport problem.** The exact K_H and local-center results constrain
   that transport. Another arbitrary field chosen from N's sign chart
   without a remote U argument is low-value evidence.

The first decisive constructive result is one beta-zero rational shape
with **three distinct nonzero origin roots S/U/S and a remote U root**,
all on valid transverse full returns at the same parameter. A fold pair
must be checked on its two-cycle side; the old stable root cannot silently
become one of the pair and then be counted again. Refine the four gates
and validate hyperbolicity/persistence on one common neighborhood. Then
choose an explicit sufficiently small negative rational beta, establish
the additional unstable Hopf cycle inside the first old origin cycle,
and revalidate all four old cycles at that same beta. Five disjoint
isolated return gates suffice; a global exact cycle count is unnecessary.

None of those next steps guarantees a counterexample. They do address the
object the handoff actually asked for. At present the serious open question
is the finite-amplitude coexistence geometry. The recorded local theorems
remove several tempting shortcuts, while the numerical campaign has not
yet tested that remaining question with the required method.
