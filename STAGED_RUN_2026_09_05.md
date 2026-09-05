# H16 staged strike: a finite KKL fold signal, with gates still open

Baseline: main `3b94f34`. Date: 2026-09-05. All new computations and failed
requests are in `staged_2026_09_05/`. The original 206-call ledger is unchanged.

**The augmented search located a numerical finite-radius fold and independently
reproduced the associated two-cycle sign pattern. It did not produce M1 or a
five-cycle field.** The old radius cutoff was also passed: both tracked cycles
persist numerically to c=0.9683, where the remote x-extremum is approximately
−6.5372 billion. Neither proposed family-wide kill condition has been proved.

## 1. Results against the requested stages

| Stage | Executed result | Remaining gate |
|---|---|---|
| 1: theory first | Bounded source/proof attempt, corrected attribution, exact failed comparison hypothesis, new amplitude restriction | K1 unresolved; full 1999 proof unavailable, complete two-session extension not performed |
| 2: KKL | 400 charged returns; compactified evaluator, terminal continuation, six K seeds for augmented Newton, 162 sign profiles; finite fold signal independently checked | No certified fold curve/whole-region exclusion; infinity connection and global transition coefficient unresolved |
| 3: outside KKL | 150 charged returns on rational Shi and Chen–Wang trace paths; both lose inner cycle into focus; independent 40-digit checks | Two paths do not cover the four-parameter chart; no M1 or interval verifier |
| 4: remote cycle | Exact remote stability gate checked at the new pair field | M1 prerequisite absent; remote-Hopf construction not run |
| 5: certification | Exact rational pair coefficients saved; independent numerical formulations agree | No five-cycle candidate; interval five-cycle certification not triggered |

The independent roles were separate review agents and distinct numerical
formulations. This record does not claim the external Fable model executed.

## 2. The concrete new construction signal

The beta-zero family is

\[
 x'=y+x^2+xy,\qquad
 y'=-10x^2+\frac{11}{5}xy+cy^2+\alpha x,\qquad
 K=-\alpha(11c/5-1)-42.
\]

On the complete curved section \(\sigma(r)=(r,-r^2/(1+r))\), augmented
Newton gives the numerical candidate

\[
 K=1/512,\quad c\simeq0.9688884793906646,\quad
 r\simeq6.949087993605231.
\]

The equations actually solved were

\[
 F_1=\log(R/r)=0,\qquad F_2=rR'(r)/R-1=0.
\]

They are equivalent to \(D=D_r=0\) at any finite positive fixed point.
The final follow-up divides both equations by \((r/(1+r))^2\) to remove
their automatic quadratic small-radius vanishing. The final normalized
residuals are approximately \((-1.8,1.9)10^{-14}\). An independent
original-time Cartesian return gives \(\log(R/r)=-2.02\cdot10^{-11}\),
with divergence exponential \(0.999999999943\). These differences measure
numerical agreement; they are not rigorous error intervals.

The preceding finite-difference Jacobian of the normalized equations, in
variables \((\log r,c)\) at fixed K, was approximately

\[
 \begin{pmatrix}
 2.41\cdot10^{-9} & -0.213746\\
 2.8914\cdot10^{-5} & -0.204861
 \end{pmatrix}.
\]

Its nonzero determinant supports ordinary-fold nondegeneracy numerically.
In particular the return minimum has nonzero curvature and moves downward
as c increases. This is not an interval nondegeneracy certificate.

The pair-side field was evaluated at **exact rational coefficients**:

\[
 c=\frac{4844447396953323}{5000000000000000},\qquad
 \alpha=-\frac{1050048828125000000}{28288921366486553},\qquad
 \beta=0,
\]

with the unchanged quadratic coefficients and exactly \(K=1/512\).
Its independent Cartesian returns give:

| Section r | Numerical log(R/r) | Sign |
|---:|---:|:---:|
| 4 | 0.00000170306870 | + |
| 6.949087993605231 | −0.000000163369528 | − |
| 12 | 0.00000517773575 | + |

This is the numerical signature of an S/U pair, conditional on the usual
return-domain continuity between the sampled points. The central sign at
the same radius but c one millionth **below** the fold is positive,
\(1.6332878\cdot10^{-7}\). A further pair-side sample at r=20000 is
positive, \(\log(R/r)=0.00971712488\). No third origin-cycle bracket was
obtained at this field, and intervals between or beyond samples are not
excluded.

The exact remote trace gate also fails the former stable-remote-focus
precursor condition: \(K_H\simeq0.03835816>K=0.001953125\), so its remote
equilibrium has **positive trace**. Thus this finite pair is not M1 and
cannot be relabelled the desired common-field precursor. See
[fold scope review](staged_2026_09_05/fold_scope_review.md) and the machine
[coefficient record](staged_2026_09_05/fold_pair_coefficients.json).

The useful advance is specifically the finite fold signal that the earlier
pilot never implemented a solver to seek. It is not a fifth-cycle claim.

## 3. The old cutoff was not the end of the cycles

Both old horizontal and complete curved sections were benchmarked. The new
engine integrates log radius and angle with positive time scaling
\(dt/d\tau=1/(1+\sqrt{x^2+y^2})\). It also integrates physical time and
original divergence. The bounded section coordinate is
\(q=1/(1+|r|)\); the old Cartesian section cutoff is removed. Numerical
guards remain: log radius <32, physical time <10, desingularized time <20000.
Guard stops are unresolved returns, not escape proofs.

The following are curved-section coordinates, not the old y=0 coordinates.
All listed points have K=1/64 and beta=0, with alpha computed exactly from c.

| c | Origin r | Origin multiplier | Remote r | Remote multiplier | Remote physical period |
|---:|---:|---:|---:|---:|---:|
| 0.9301 | 50.8726 | 0.976343 | −1.09990e6 | 1.902808 | 0.544221 |
| 0.9500 | 56.5679 | 0.987580 | −4.46920e6 | 1.427140 | 0.555090 |
| 0.9640 | 82.7328 | 0.995898 | −2.09928e7 | 1.106510 | 0.563099 |
| 0.9680 | 847.213 | 0.998288 | −3.95252e8 | 1.017363 | 0.566216 |
| 0.9682 | 2538.12 | 0.998378 | −1.64001e9 | 1.013017 | 0.566756 |
| 0.9683 | 7036.91 | 0.998405 | −6.53720e9 | 1.010887 | 0.567210 |

The final common point was replayed in original Cartesian coordinates.
Its origin and remote log-return residuals were −1.18e−9 and +3.55e−9.
The remote multiplier differed between formulations by about 1.8e−7.
These are numerical consistency checks, not certification at billion-scale
amplitudes. The follow-up allowance ended during the c=0.9684 origin
corrector; that field was not accepted and no common pair is claimed there.

The exact infinity audit explains why the original diagnostic was too strong:

* Below c=241/250, the vertical graphic's interior connection is impossible
  because it would cross x=−1 right-to-left, while x'=1 on that line.
* For 241/250<c<1, the relevant two-saddle eigenvalue product is neutral on
  J(c)=305+634c−11c²−1000c³=0, giving c*=0.968620633553494…, independent of K.
* This is a **candidate eigenvalue-neutrality line**, not a proved graphic.
  Connection splitting and the global transition coefficient C−1 remain
  uncomputed. Nonhyperbolic boundaries need separate treatment.
* A physical period can stay bounded near infinity despite divergent
  desingularized residence time. A nonunit multiplier also does not exclude
  a graphic. Neither diagnostic alone is a valid terminal classification.

See [exact infinity analysis](STAGED_INFINITY_2026_09_05.md).

## 4. K coverage and what the sign table means

The finite profile design used c in {0.6,0.7,0.825,0.93,0.95,0.965,1,1.2,1.5},
K in {1/512,1/64,1/16,1/4,3/5,119/100}, and r in {2,20,20000}.
These are 162 return requests, spanning six K values in the requested open
interval. The profile is an exploratory finite design, not coverage of all
K>0 or the inherited parameter box. The upper c strata can fail the remote
precursor gate; their origin returns must not be counted as common precursors.

[The exact sign table](staged_2026_09_05/PROFILE_SIGNS.md) and
`profiles.json` preserve every result. At c=0.965 the outer sign changes
across the sampled K values; at c=1, small K gives the sequence −,−,+.
Six r=20000 returns at c=1.2 hit numerical guards. These are missing returns,
not negative displacements. The samples do not prove displacement monotonicity.

Six augmented-Newton seed attempts ran for four iterations each. The smallest
K attempt supplied a follow-up that converged to the finite fold candidate.
A separate follow-up used the accepted path cycle with smallest |P'−1| from
the first continuation block. It moved toward large amplitude and did not
converge in its allowance. No continuous fold curve across K was completed;
disconnected fold sheets remain untested.

## 5. Theory and the outside-KKL controls

[Stage 1](STAGED_K1_THEORY_2026_09_05.md) distinguishes Zhang–Cai's
weak/strong-focus distribution result from Zhang's 1999 general order-two
own-nest theorem. The latter's full proof could not be retrieved. The
accessible monotone-ratio proof fails exactly because its derivative changes
sign in the KKL region. This is a failed proof hypothesis, not evidence of
three cycles.

A new exact Bernstein-coefficient certificate gives a stronger amplitude
restriction throughout 0.9≤c≤1, 0<K≤1.2: the multiplier polynomial N is
nonnegative on −1≤x≤0 and has a unique positive root a(c,K). Every stable
or multiplier-one origin cycle must reach x>a. The new fold candidate passes
this necessary test; the test supplies no existence or uniqueness theorem.

[Stage 3](STAGED_SHI_2026_09_05.md) gives both explicit rational trace paths,
their root brackets, and independent 40-digit direct-polar checks. Both
initial fields have three numerical origin cycles; each endpoint has two
tracked cycles after the inner Hopf collapse. No remote cycle was recomputed
on these paths. The Shi control is a conditioned rational field, not the
published tiny-scale Galias–Tucker field.

Three simple cycles at trace zero would not require a double cycle at that
same field. A codimension-one fold surface is also not guaranteed to be easy
to detect by random sampling. These corrected premises prevent treating the
two successful trace diagnostics as a global K1 kill test.

## 6. Reproduction and honest stopping point

New evaluations: **400 KKL (394 completed, 6 unresolved) + 150 Shi/Chen–Wang
(139 completed, 11 failed) = 550**. With the original 206 calls, the conservative
shared 4096 accounting is **756 used, 3340 remaining**. Stage 2 spent exactly
its 400-call allocation. Recorded successful KKL evaluator CPU totals about
8.77 seconds; Shi records about 13.72 seconds. These exclude reading, symbolic
work, process startup and reporting, and are not total task runtimes.

Exact checks and static artifact validation consume no ODE calls:

```sh
python staged_2026_09_05/theory_exact_checks.py
python staged_2026_09_05/infinity_check.py
python staged_2026_09_05/validate_artifacts.py
```

Discovery replay sources are `run_kkl.py`, `refine_kkl.py`, `finalize_kkl.py`,
`shi_run_continuation.py`, and the low-level evaluators. The initial KKL
benchmark requests are fully recorded in the ledger. The KKL supervisor now
refuses additional calls against its full 400-call ledger. Reproduction must
retain the historical record and use a separately budgeted output location.

Every KKL row contains a request, source hash, result/status, time and purpose.
The exact original evaluator source is archived as `compact_return_v1.py`;
its hash matches rows 1–392. The current version adds explicit flux and remote
branch guards, as requested by the independent source review. Rows 393–400
contain the new guarded profile and independent Cartesian checks.
Shi's initial seven successful rows predate its diagnostic correction; they
are labelled historical. Per-evaluation Shi source hashes were not captured,
so the final manifest must not be mistaken for historical execution provenance.

No M1 point was certified. Consequently Stages 4 and 5 did not run their
conditional constructions. No full four-parameter fold search or interval
Poincaré certifier was completed. The first global graphic coefficient is
also still missing. These are explicit remaining tasks, not hidden completions.

The next mathematical question sharpened by this run is whether the located
finite pair can be continued to coexist with a third origin cycle, while
respecting the remote stability/distribution constraints. The present data
establish neither that coexistence nor an obstruction to every such path.
