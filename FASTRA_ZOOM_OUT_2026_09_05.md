# FASTRA H16: transcript, campaign audit, and construction reset

Snapshot: 2026-09-05 UTC, still 2026-09-04 in Edmonton. Research baseline:
main `6048ed8b09e12094b0875fdc68c5a1fb6341775d`. This is a review of the
completed work, not a new numerical attack. No ODE evaluations were added.

## 1. The answer the record supports

**We have no five-cycle candidate and no five-cycle certificate. We also
have not performed an adequate negative test of the selected KKL family.**
The campaign proved useful restrictions, reproduced four-cycle controls,
and followed a two-cycle branch to its prescribed section limit. The
central finite-amplitude construction task—finding a new return-map
minimum/maximum pair and bringing its height through zero—was not carried
out. That gap matters more immediately than finding another distant
family to sample.

There is a second logical restriction to remove from our thinking. The
goal is **at least five cycles**, not precisely five, and the beta-zero
Hopf precursor is sufficient for one construction but is not necessary
for every possible counterexample. Preserving the incumbent four while
creating an additional pair would give six. Conversely, six is a stronger
outcome than five; proximity to a four-cycle control is not a probability
estimate. Both observations belong in the strategy.

The accepted lower bound remains four. A recent primary paper still
describes the quadratic maximum as unknown; nothing retrieved in this
review establishes a new accepted example. Uniform finiteness of the
maximum is also not known. [Artés–Cairó–Llibre, August 2026](https://link.springer.com/article/10.1007/s12346-026-01563-4),
[Gasull–Santana, 2025](https://ddd.uab.cat/pub/artpub/2025/309367/GasSan24-Postprint.pdf).

## 2. What “full transcript” contains

The companion local archive contains **78 visible messages: 14 user and
64 assistant messages**, through the user's request to zoom out. It also
contains four distinct recorded task objectives, explicitly separated
from conversation messages. It preserves the original visible wording,
including statements subsequently corrected. The GitHub credential is
redacted. A machine-readable event export and a SHA-256 manifest accompany
the readable Markdown.

This is the full available visible record of this local task, not a
fabricated transcript of the older Astra/Claude/Fable conversations.
Earlier strikes are reconstructed below from committed reports. Private
reasoning, internal instructions, tools and agent exchanges are not
conversation entries. The extractor reads only public message events and
the separate user-provided objective records.

Local artifacts are named `CONVERSATION_TRANSCRIPT.md`,
`CONVERSATION_TRANSCRIPT.jsonl` and `TRANSCRIPT_MANIFEST.json`, in the
delivered `outputs/review_2026_09_05/` directory. The raw conversation is
not included in this public research commit. The reproducible extraction
logic is [export_visible_transcript.py](review_2026_09_05/export_visible_transcript.py).

The repository review includes the user's named canon, the associated
proof/check scripts, the latest endpoint note, the KKL evaluator and all
206 saved return records. The independent post-Q4 packet is on a separate
commit, not main; it was read explicitly rather than assumed merged:
[frontier](https://github.com/git-df-scott/H16P/blob/de39ea78d56208a2a3267b594ce5c117b6b14c1e/H16_POST_Q4_FRONTIER_2026_09_04.md),
[historical claims](https://github.com/git-df-scott/H16P/blob/de39ea78d56208a2a3267b594ce5c117b6b14c1e/HISTORICAL_FIVE_CYCLE_CLAIMS.md),
[mechanisms](https://github.com/git-df-scott/H16P/blob/de39ea78d56208a2a3267b594ce5c117b6b14c1e/FIVE_CYCLE_MECHANISMS.md),
[four-cycle seeds](https://github.com/git-df-scott/H16P/blob/de39ea78d56208a2a3267b594ce5c117b6b14c1e/FOUR_CYCLE_SEED_LEDGER.md),
[KKL assignment](https://github.com/git-df-scott/H16P/blob/de39ea78d56208a2a3267b594ce5c117b6b14c1e/ASTRA_FIFTH_STRIKE_HANDOFF.md).

## 3. Campaign chronology, with corrections kept visible

Commit subjects and historical handoffs are records of what was believed
then; the last column states what survives now.

| Phase | Anchor commits | Work done | Current meaning |
|---|---|---|---|
| Q4 setup and Strike 1 | `9db4cb3`, `16f9406` | Universal chart, Stieltjes/ECT structure, corrected strip and cusp exclusion | Inherited exact restrictions; no five-zero example |
| Q4 Strike 2 | `6165190` | Lobe cell, reconstruction, Green/first-root constraints | Inherited proof machinery |
| Q4 Strike 3 | `8eb89c6`, audit `fe99a20` | Rational threshold certificates and scoped large-kappa obstruction | Verified within their stated scope; no extrapolation of slow asymptotics |
| Premature global claims | `084fae3`, `48ff5b1` | Handoff said global three; endpoint discussion treated the elliptic saddle as the original loop | Both claims were corrected later; do not inherit them from their subjects |
| Theorem N | `46bca95`, audit `ed5817b` | Negative first Green maximum on the strict lobe | Three distinct original interior zeros there; four globally. Not a global sharp-three theorem |
| Remaining two-root case | `4ee2e84`, `6c96d70` | Two-anchor reduction, mixed center sign, remaining Green determinant condition | No certified four-zero original integral, and no global exclusion of four |
| Other topology gates | `b85125c`, external audit | Tested maximal-focus graphic and quoted reversible seeds fail | Closures are tied to their strata/seed geometry, not all quadratic systems |
| Council | `5a2d653`, `2136896`, correction `4d5b53e` | KKL direct construction and original-coordinate Q4 endpoint allocation | Council §9 is the governing correction; two independent lanes |
| Fable endpoint checkpoint | `d168b7b` | First-order compatibility identities and numerical samples | Useful negative samples; composed original return and global compatibility remain unproved |
| KKL pilot | `bf04ad9` | 202 evaluations, 22 accepted common-cycle points, derivative repair, positive stationary maximum | One selected path; no precursor and no fold certificate |
| KKL exact gates and alternate control | `6048ed8` | Double-center/local-unfolding and remote-Hopf restrictions; four c>1 returns | 206 calls total; no additional cycle and no stratum-wide exclusion |
| Present review | this report | Full public transcript, three independent audits, complete curved-section proof | Identifies missing coverage and a concrete revised construction protocol; no new ODEs |

The [two-root reduction](Q4_TWO_ROOT_REDUCTION.md) deserves preservation,
but it is not the direct five-cycle lane by itself. Its center functional
is mixed on the two-root region. It reduces the remaining first Green
maximum question to an explicit determinant on the surviving fiber; it
does not supply four original integral zeros. Completing it remains a
mathematical result distinct from the ultimate construction target.

## 4. The KKL experiment we actually ran

The exact fixed family is

\[
 \dot x=y+x^2+xy,\qquad
 \dot y=-10x^2+\frac{11}{5}xy+cy^2+\alpha x+\beta y,
 \qquad K=-\alpha(11c/5-1)-42.
\]

S and U below mean numerically stable and unstable simple return roots.
All locations and multipliers in this section are ordinary numerical
results, not interval certificates.

| Common exact field | Detected origin controls | Detected remote control | Meaning |
|---|---|---|---|
| `c=7/10, alpha=-363889/5000, beta=3/2000` | S/U/S at roughly 0.68321, 2.18370, 15.96278 | U at roughly -3711.56 | Incumbent four-cycle numerical reproduction |
| Same c and alpha, `beta=0` | U/S at roughly 3.06885, 15.06407 | Not an additional beta-path record | Static negative-K control, not a demonstrated beta continuation |
| `c=7/10, alpha=-80, beta=0` | S at 64.55543, multiplier 0.80969 | U at -5391.14, multiplier 12.1680 | Positive-K starting pair; K=6/5 |
| `c=9301/10000, alpha=-8403125/209244, beta=0` | S at 48.69484, multiplier 0.97634 | U at -1048286.51, multiplier 1.90281 | Last accepted common field, K=1/64, close to the remote section cap |

The path increased c at K=6/5, reduced K at c=33/40, then increased c at
K=1/64. It reached the remote `-2^20` section-coordinate allowance before
detecting an origin fold. This is a bound on the chosen experiment. It
does not show the orbit disappears or that another path cannot retain it.

At the last field the stationary calculation found

\[
 r\simeq28.17411716,\quad D(r)\simeq+0.24269325,\quad
 D_r(r)=0,\quad D_{rr}(r)\simeq-0.00112207.
\]

This is the expected positive maximum inside the known S root. A fold
of cycles instead needs **both** D=0 and D_r=0. At the later exact
`c=1001/1000, alpha=-196/5, beta=0` checkpoint, four widely separated
remote starting coordinates all gave positive displacement. No remote U
was bracketed; absence of such a cycle was not proved.

The accounting is **206/4096 evaluations, 3890 unused**, with 22 accepted
points from 24 continuation events. All 206 return records completed;
rejected continuation points are not failed ODE integrations. Recorded
evaluator CPU time was about 6.86 seconds and subprocess wall time about
67.1 seconds. Those figures exclude analysis, derivations, reading and
writing. They cannot be used to promise how long an unsolved construction
will take. [Frozen summary](kkl/data/strike_summary.json),
[full ledger](kkl/data/returns.jsonl),
[coverage audit](review_2026_09_05/COVERAGE_AUDIT.md).

## 5. What we missed or overstated

### 5.1 We did not build the fold continuation promised by the council

The drivers correct a known simple root, require its stable multiplier,
and stop near a singular correction. There is no implemented augmented
fold solver or pseudo-arclength fold continuation in that pilot. A
multiplier approaching one is not discovery of a cycle-creating fold.
The fold of the root already being tracked could annihilate it.

For beta zero and K>0, D is positive sufficiently near the origin. S/U/S
requires at least three finite critical points: a positive maximum, a
negative minimum, and a further positive maximum. We located one positive
maximum. We did not locate the additional minimum or maximum. The useful
new object is a new stationary branch, not more samples of the old root.

### 5.2 A selected path and a sparse radial profile are not coverage

The ten-point final origin profile does not exclude zeros between its
samples. It also starts above the specified lower section cutoff. A
qualitative local no-collapse theorem supplies no numerical radius
covering that gap. The c>1 control was a single algebraically selected
field, not a transported four-cycle seed. No whole-box negative statement
follows from any of these records.

The code has Cartesian time/coordinate guards, but no implemented
compactified return atlas. Classifying infinity equilibria symbolically
does not validate passages near infinity or extend a stopped return.

### 5.3 The old remote section lacks a completeness argument

The remote focus lies above y=0. Four inward returns from that horizontal
line cannot exclude a surrounding cycle that never meets it. We have no
evidence that an actual U was missed; the issue is unproved observation
coverage.

This review supplies an analytic repair. Under the stated two-focus and
single-root gates, the rational nullcline

\[
 \sigma(r)=\left(r,-\frac{r^2}{1+r}\right)
\]

meets each origin cycle once on r>0 and each remote cycle once on r<x_*.
The proof uses the barrier x=-1, the diffeomorphism `(x,y)->(x,x')`, and
the common orientation of all crossings of a ray from the enclosed
focus. It applies also to nonzero beta when its explicit cubic gate
holds. The off-root derivative is

\[
 R'(r)=\frac{Q(\sigma(r))}{Q(\sigma(R(r)))}
       \exp\int_0^T\operatorname{div}F\,dt.
\]

The proof was independently checked and the polynomial/rational
identities passed an exact single-thread replay with a ten-CPU-second
fuse. This is a new **section theorem**, not an orbit certificate.
Its numerical evaluator has not been implemented, its controls have not
been transported, and old section-coordinate caps cannot silently be
reinterpreted as caps in the new coordinate.
[Proof and implementation requirements](review_2026_09_05/KKL_SECTION_REPAIR.md),
[exact check](review_2026_09_05/nullcline_exact_check.txt).

### 5.4 The Hopf construction became an unnecessary exclusive assumption

The intended conditional step is sound: beta-zero K>0 with origin S/U/S
and remote U, followed by a sufficiently small common negative beta,
would produce the inner Hopf cycle while preserving four hyperbolic
cycles. But no theorem says every possible five-cycle KKL field can be
deformed to that precursor without losing finite cycles.

The known four-cycle control has K<0 and beta>0. Replacing it with the
positive-K beta-zero starting pair gave up two detected cycles before
searching for the missing pair. That was a route choice, not a necessity.
An alternative is to retain all four and seek an additional finite pair
in the larger nest. Another is to create the small Hopf cycle first and
then continue finite-amplitude geometry at nonzero beta. Neither has
been tested by the saved pilot.

At least-five configurations must obey the distribution theorems. For
exactly five, the possibilities are (5,0) and (4,1); an added pair while
four persist can give the permitted (5,1), namely six. Adding remote
cycles to a retained (3,1) would instead violate the allowed
distributions. The theorem does not cap the larger nest at four.
[Zegeling, Theorems 1.2 and 5.4](https://d-nb.info/1332906729/34).

### 5.5 Local exact obstructions remain local

The simultaneous beta-zero K=J=0 point, with
`J(c)=305+634c-11c^2-1000c^3`, is a double center, not the desired
higher-order focus. The proved local unfolding prevents the proposed
small-cycle pair there. The remote Hopf sign likewise rules out the
desired *small* U birth at its stated positive-K points. Neither result
excludes a finite-amplitude U or origin fold elsewhere.

Proving a global remote-U obstruction at the single c>1 control could be
useful if a definite Dulac/Liénard sign is available. Starting another
large analytic detour without such a sign would postpone the missing
construction test. [Exact local result](kkl/notes_local_unfolding.md),
[remote gate](kkl/notes_other_strata.md).

### 5.6 Numerical derivatives and certification still need work

The transverse-determinant derivative repair addresses the observed
cancellation at large remote amplitude. It does not turn the integrations
into rigorous ones. Direct parameter-derivative controls and moving-event
second-derivative checks should precede an augmented solver. Agreement
between related variational identities on the same trajectory is useful
but is not independent validation of that trajectory. Historical rows
1–202 lack per-run source hashes; this is recorded rather than recreated.

There is no implemented interval ODE certifier in the current lane.
Certifying one known control should be an early benchmark, so discovery
does not outrun the method needed to prove a candidate.

## 6. Fable's lane: keep the original return, remove premature closures

The original Q4 boundary has two infinity saddles. The prior single-loop
closure in the elliptic chart was invalid for that original graphic.
The corrected lane remains open, but the latest first-order samples do
not finish it either.

Sampled positive endpoint values do not prove positivity for every shape
parameter. The 546 sampled two-root points and other sampled families do
not establish a universal first-order closed count. Theorem N does not
by itself append a boundary zero to its distinct interior bound, nor
provide the nonlinear composed displacement. Hence the reduction to
“only two alien cycles remain” is premature.

The next bounded task is exactly the original joint calculation: two
Dulac maps, both regular transitions, a common admissible quadratic
perturbation, and controlled remainders. Pair it with an actual direction
having three simple zeros of the original integral, as opposed to three
zeros of the auxiliary primitive H. Two actual endpoint cycles need not
be two alien cycles. Four interior plus one endpoint remains conditional
on the unresolved original four-zero question.

This does not reopen an unrestricted multiplicity loophole. Multiple
interior zeros often split into additional distinct zeros under
admissible perturbations; whether that splitting is available must be
checked before proposing them as an escape from a distinct-zero theorem.
The specific uncompleted issue here is original endpoint compatibility.

Published alien-cycle examples show that first-order counting can miss
cycles, not that arbitrary two-saddle graphics contribute two extras.
The quadratic Hamiltonian comparison itself has a closed-annulus upper
bound under its hypotheses. [Gavrilov–Iliev, Theorems 3–4 and Appendix A](https://www.math.univ-toulouse.fr/~gavrilov/publications/50.pdf),
[Zhao, three-cycle result](https://arxiv.org/pdf/1011.2253).
The detailed adverse review is [FRONTIER_AUDIT.md](review_2026_09_05/FRONTIER_AUDIT.md).
No Claude/Fable file was changed and no duplicate Q4 attack was run.

## 7. Where to search next, in order of readiness

The ordering below is an experimental judgment, not a success probability
or a theorem that other routes cannot work.

| Route | Actual starting evidence | Next bounded discriminating task | First useful success signal | What a negative result would mean |
|---|---|---|---|---|
| KKL finite-amplitude construction, with complete sections | Four numerical controls and reviewed derivatives in an explicit rational family | Implement/benchmark section transport, then continue stationary branches with explicit retained-cycle controls | A new minimum/maximum branch with a direction toward zero while required cycles persist | Failure of that seeded branch only |
| KKL incumbent-preserving amendment | Same exact four-cycle control, at finite beta and negative K | A 64-call stationary-geometry block in a written three-parameter neighborhood | Additional finite pair separated from the three old origin roots; remote U retained | Does not exclude a different sheet or the full KKL family |
| Fable: Q4 original closed annulus | Verified interior machinery and identified two-saddle infinity graphic | Derive composed original return and compatible scaling | One perturbation direction with three actual interior cycles and two actual endpoint cycles | An obstruction to that scaling; broader closure needs a uniform theorem |
| Conditioned GT/Songling backup | Published certificate at the original field; recorded relaxed four-cycle numerical endpoint | Interval-certify four relaxed controls and a preserving neighborhood | A certified usable four-cycle neighborhood, then a new finite-pair direction in the larger nest | Failure of conditioning or selected direction, not all Shi geometry |
| Resonant infinity hemicycle | Explicit family and a theorem omitting an upper bound on the resonance | Derive the missing joint resonant return with shared parameters | Five compatible sign brackets on one original return configuration | A result only for the analyzed graphic/unfolding |
| Published reversible two-center four-cycle family | A legitimate (3,1) mechanism distinct from the rejected three-saddle seed | Recover an exact conditioned control and its finite-amplitude critical geometry | A retained four-cycle seed plus a new larger-nest branch | Not a revival of the invalid center-plus-three-saddles proposal |

The GT relaxed field already recorded in the independent audit is

\[
 \dot x=-10^{-24}x-y-10x^2+\frac{499}{100}xy+y^2,\qquad
 \dot y=x+x^2-\frac{311375001}{12500000}xy.
\]

It is a numerical four-cycle control requiring validation, not a new
candidate. The published exact-four theorem concerns its original
Songling instance, not every nearby field. [Galias–Tucker publication](https://research.monash.edu/en/publications/the-songling-system-has-exactly-four-limit-cycles/).
The distinct reversible two-center construction is
[Yu–Han](https://arxiv.org/abs/1002.1055).

The explicit resonance backup is

\[
 \dot x=(b-2)/4+\epsilon_1x+(1-b)y+ax^2+\epsilon_2xy+by^2,
 \qquad\dot y=\epsilon_0-2xy.
\]

On a=-1 the cited theorem gives lower bounds rather than the sharp
upper bound available in other regions. This is an uncompleted
compatibility calculation, not evidence of five.
[Marín–Villadelprat, Theorems B–D](https://arxiv.org/html/2501.16924v1).

Historical five/six claims already contradicted by their own parameter
conditions remain frozen. The recent seven-cycle abstract still provides
no retrievable explicit field in the checked record; it is unverified,
not refuted by lack of access. Unavailable manuscripts and old rebuttal
details remain source-retrieval tasks, not construction seeds. No new
accepted witness was located in this bounded literature check.

## 8. The next experiment, made concrete

**Immediate priority: repair observation and stationary continuation,
then run one small discriminating construction block.** Do not spend the
remaining allowance extending the old broken-line path.

The [construction review](review_2026_09_05/KKL_NEXT_CONSTRUCTION.md)
specifies a pilot that keeps the existing restrictions: at c=33/40,
bracket K so the remote U has old section coordinate -2^15, then continue
that remote constraint jointly with D_r=0 for an origin stationary
point. At most 64 additional evaluations or 16 accepted steps. The
existing remote coordinates at K=1/64 and 6/5 motivate the bracket; its
endpoint displacements must actually be computed. The pilot starts on
the known maximum, so it becomes construction progress only if it seeds
new stationary geometry or a genuine cycle-producing event. It is not a
cover of disconnected branches.

For the broader reset, my recommended **scope amendment** is to preserve
the incumbent four at finite beta and inspect its critical branches
before returning to an exclusive beta-zero precursor. Keep the fixed
coefficients -10 and 11/5. The first block should be a bounded local
experiment, with the changed beta/K restrictions and section coordinates
written down before any search. A suggested exact initial box is

\[
 c\in[69/100,71/100],\quad
 \alpha\in[-74,-72],\quad
 \beta\in[1/1000,1/500].
\]

This box merely defines the pilot domain; it does **not** assert four
cycles throughout it. At every accepted field, topology and all retained
roots must be checked. Use the exact incumbent in section 4, one
continuation direction chosen from measured critical-height derivatives,
no coefficient grid, and the same 64-call/16-step cap. A new stationary
pair, its continuation to a nondegenerate D=D_r=0 event, and subsequent
separated simple roots are three different milestones. Following an
annihilation of the old pair is not success.

The complete curved section should be benchmarked before interpreting
new remote negatives. Preserve the original phase-space coordinate guard
and ten-CPU-second, one-thread limits; record control transport between
old and new sections and explicitly translate or replace the old section
coordinate bounds. Reserve at most 24 charged evaluations for this
validation; if it fails, fix the evaluator and do not consume the
construction block. Neither this amendment nor either pilot has been
executed by the present review.

These are alternatives within Astra's lane, not instructions to launch
two numerical lanes. My preference after validation is the finite-beta
incumbent block because it tests the unnecessary precursor restriction
directly. If the original council scope is retained, use the constrained
remote-amplitude pilot. The full reason for considering both, including
the weaker conditioning near the incumbent's inner cycles, is in
[KKL_NEXT_CONSTRUCTION.md](review_2026_09_05/KKL_NEXT_CONSTRUCTION.md).

An honest next-hour deliverable is a validated control transport plus
either a new stationary branch with a retained-cycle direction or a
clearly bounded failure record. A five-cycle candidate is a possible
outcome, not a scheduling promise. Fable's independent hour should produce
the original composed-return expansion or a stated missing coefficient,
not another unqualified global claim from samples.

## 9. How an actual counterexample gets proved

The trigger is one **exact rational parameter vector** with at least five
numerically isolated transverse full-return roots, consistent itineraries
and correct nesting at those same coefficients. Stop discovery and hand
that vector, sections, brackets, orbit bounds and software provenance to
the other lane for independent hostile reproduction. A fold itself is
not this trigger: move to the cycle-present side where roots separate.

Then obtain a rigorous original-field certificate:

1. Fix the exact quadratic coefficients; re-evaluate after any rounding.
2. Validate existence of each full return over its entire section bracket,
   including the event and exclusion of escape or itinerary changes.
3. Prove strict opposite displacement signs at its endpoints and an
   isolating/uniqueness condition, such as interval Newton or derivative
   separation from one.
4. Prove the brackets represent distinct cycles. The complete-section
   theorem or validated disjoint orbit neighborhoods supplies this step;
   multiple crossings of one orbit cannot be counted twice.
5. Independently replay all five certificates at the same coefficients.

The incumbent-plus-pair mechanism would predict six. Proving any five
distinct isolated cycles already establishes H(2)>=5; documenting all six
would establish the stronger lower bound. No local Hopf count, auxiliary
integral zero count, or cycles at differing coefficients can substitute
for these common-field claims.

## 10. Evidence status at the end of this review

| Claim | Evidence class |
|---|---|
| Known quadratic lower bound at least four | Published result; published computer-assisted certificates exist for particular fields |
| Global maximum four or uniform finiteness | Open; not asserted |
| Q4 Theorem N and two-root reduction | Audited repository proofs, inherited within their exact scope |
| KKL double-center/local and remote-Hopf gates | Analytic derivations with exact symbolic checks; local scope |
| New complete curved section | Analytic proof, independent review, exact rational identities checked |
| 206 KKL returns and listed roots | Numerical only, not interval-certified or a parameter-space cover |
| Five-cycle candidate in this campaign | None |
| Five-cycle quadratic field certified | No |
| KKL extra finite pair, Q4 joint endpoint compatibility | Open |
| Preference for finite-beta incumbent pilot | Experimental judgment, not an existence theorem |

The durable output of this review is a corrected map of what was actually
tested and a repair that makes future negative observations meaningful.
The next construction must produce new critical geometry and preserve
the required cycles at one field. That is the missing step between the
current four-cycle controls and a counterexample.
