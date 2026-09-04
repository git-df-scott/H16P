# Third-strike completion-scope audit

Audit stage: final mathematical, artifact, and regression audit, completed
before the external publication receipt. The inherited base is
`616519009cdec6617c2969254a44ea05f86bdb42`. All 59 inherited tracked files
were verified unchanged. The three required third-strike root artifacts
exist and have been inspected. This audit performs no new search.

Publication is deliberately a separate check: the final commit hash and
remote branch state must be verified in the external delivery receipt.
A committed audit cannot contain its own eventual commit hash without
changing that hash. This file therefore proves content completion and
records the publication requirement without asserting a publication result.

## The actual objective remains unachieved

No current artifact proves five distinct simple zeros of the original Q4
integral, exhibits a five-zero numerical candidate, or forces a positive
first Green maximum. The threshold certificate is a certificate of three
universal primitive zeros only. It must not be presented as a smaller
replacement for the requested original-five-zero result.

The documented bounded strike ends under **G: the bounded high-value
tasks are exhausted**. Conditions A, B, C, D, E, and F are not
established. In particular:

- E is disproved by the exact threshold path and rational late-root box.
- D requires both a forced positive maximum and resolution of the +2
  mechanism; neither has been obtained.
- F requires closure of the whole positive-Green-maximum architecture.
  A fixed-shape, fixed-positive-lambda endpoint obstruction is narrower.

A stop under G must be described as exhaustion of this documented bounded
strike, with the original construction still open. It is not mathematical
exhaustion of every remaining Q4 parameter or asymptotic regime.

## Requirement-by-requirement status

| Objective requirement | Current authoritative evidence | Audit result |
|---|---|---|
| Continue from the specified first/second-strike state | `git rev-parse HEAD` returns the exact specified base; old tracked artifacts unchanged | Satisfied at audit time |
| Preserve corrected normalization, beta strip, cubic filter, orientation, and reconstruction sign | Inherited exact modules are reused; two new original-area checks record the corrected filter as surviving | Satisfied within tested scope |
| Identify what the first-root and kappa thresholds mean and their strictness | `Q4_THRESHOLD_PATH.md`; `notes_green_third.md`; exact threshold test in `test_q4_third.py` | Proved/documented |
| Construct a continuous lobe-admissible path crossing 5/11 | Exact anchor path `T(r,(1+r)/2,(3+r)/4)` and global analytic-cell theorem | Proved |
| Produce explicit admissible coefficients on the greater-than side | `third_threshold_certificate.json`; exact rational point and radius 1e-8 box | Rigorously certified for the primitive |
| Verify positivity before 5/11 and simple late primitive roots | Four exact signs, auxiliary ECT plus anchored Rolle, and an additional positive threshold enclosure | Proved |
| Analyze the exact Green/PF lift and use kappa above its strict threshold | `q4_green_max_3.py`, original transport, inherited reconstruction, and 16 frozen initial shots | Completed bounded diagnostic; not a zero certificate |
| Begin near the kappa threshold and continue beyond it | Initial shots use kappa 1.137, 2, 4, 16; tuned shots select one-dimensional kappa balances | Completed bounded diagnostic |
| Track required Green extrema and original witness signs | Eight tuned S1 records store four P roots and their Z heights; all first heights are negative | Completed bounded diagnostic; positive maximum unachieved |
| Understand the actual +2 sign mechanism | The sequential S1/S2/S3 criteria are retained; shots realize S1 numerically and fail S2 | Mechanism tested, saturation unresolved |
| Use a minimal boundary configuration when naturally exposed | Four primitive-confluent shots and four reverse Green-tangency lines | Tested narrower boundary routes; no original three-simple-plus-double point |
| Prefer signs to poorly conditioned root claims | Frozen primitive signs; sampled original signs explicitly numerical; no candidate promotion | Satisfied |
| Deliberately analyze the late-root / kappa-infinity regime | Exact endpoint direction and cancellations; matching theorem and uniform remainder in `Q4_GREEN_MAX_3.md`, independently checked in `notes_audit_third.md` section 5 | Proved for fixed strict anchor ratios and fixed finite positive lambda, locally uniformly on compact lambda intervals |
| Certify any five-original-zero lead immediately | No lead exists in any frozen record | Conditional gate not triggered |
| Independently check a five-zero candidate in original normalization | No five-zero candidate; two representative new lifts nevertheless receive independent area checks | Conditional gate not triggered; extra control completed |
| Produce exact original five-zero parameters and six witnesses | No such evidence | Unachieved; not replaced by primitive witnesses |
| Attack the realization gate after five zeros are certified | No five-zero certificate | Conditional gate not triggered |
| Produce required third-strike deliverables | `ASTRA_THIRD_STRIKE.md`, `Q4_GREEN_MAX_3.md`, and `Q4_THRESHOLD_PATH.md` all exist and were inspected | Complete |
| Add meaningful regression tests | Four new tests cover exact thresholds, closed moments, path residuals, and the full rational certificate/hash/box replay | Passed in final sequential replay |
| Run all inherited checks and new checks | `third_verification.txt`: old 7 tests, old 4 tests, new 4 tests, exact structure and endpoint checks all exited zero | Complete; source hashes independently rechecked |
| Preserve failed informative trials and normalization checks | Frozen numerical records, two independent original comparisons, exact threshold certificate | Completed for current trials |
| Commit and push | The immutable commit and remote main hash must be checked after this content audit is saved | External delivery receipt required; no self-embedded publication claim |
| Give the exact requested nine-field report and one next task | `ASTRA_THIRD_STRIKE.md` contains all nine fields with the correct truth values and one reverse-tangency follow-up | Artifact complete; final delivery must preserve these truth values |

## Bounded construction work inspected directly

1. `third_initial_shoot.json`: 16 shots at four prescribed path points and
   four kappa values. Every record labels itself numerical; all candidate
   flags are false and all sampled original crossing counts are zero.
2. `third_tuned_shoot.json`: five affine-path points, including r=0.9999,
   tuned to four P crossings. All four Z extrema in every row are negative.
3. `third_shape_shoot.json`: three non-affine anchor shapes, again with four
   P crossings and negative Z extrema. These are the other three of the
   eight S1-tuned tests.
4. `third_confluent_shoot.json`: four primitive triple-contact boundary shots.
   Each has numerically vanishing P at the contact and Z below -0.002.
   They are not contacts of the excluded auxiliary triple-root cusp, and
   are not original-integral double roots.
5. `third_reverse_tangency.json`: four selected Green tangency lines from
   Y(t*)=Y'(t*)=0. The two late lines include a one-dimensional ratio
   determinant diagnostic. No three-root primitive was detected. The
   source and JSON explicitly retain the possibility of missed roots or
   narrow windows; there is no rigorous whole-line exclusion.
6. `third_independent_checks.json`: two original-normalization checks against
   the positive-area evaluator. Absolute differences are approximately
   1.08e-18 and 1.16e-16. They validate the tested implementations and do not
   establish a global root count.
7. `third_threshold_certificate.json`: the exact rational primitive point
   and coefficient box. Hashes and a full rational replay are checked by
   the new regression suite.

These cover the requested analytic path, threshold crossing, one-dimensional
kappa shooting, a non-affine variation, primitive-boundary continuation,
reverse tangency, exact endpoint structure, and independent original
normalization controls. No generic coefficient sweep or long CPU run is
needed to complete the remaining administrative verification.

## Final verification and publication boundary

The formerly open matching remainder has been independently proved in
`q4/notes_audit_third.md`, section 5. The final `Q4_GREEN_MAX_3.md` incorporates
that proof and labels the result with its fixed-strict-ratio,
fixed-finite-positive-lambda scope. It contains no stale pending-proof
qualification and no global delayed-lobe closure claim.

The required root artifacts have all been inspected:

- `ASTRA_THIRD_STRIKE.md`: outcome, exact nine-field report, failed bounded
  construction routes, scoped asymptotic theorem, and one next task.
- `Q4_THRESHOLD_PATH.md`: exact path, strict threshold meanings, rational
  certificate and box, closed moments, and endpoint coefficient asymptotics.
- `Q4_GREEN_MAX_3.md`: exact reconstruction identities, scoped matched
  theorem with remainder proof, numerical diagnostics, and remaining domain.

`q4/data/third_verification.txt` records the final sequential successful
replay. This audit independently compared its script hashes to the current
source files and found every one equal:

| Check | Result | Recorded process CPU |
|---|---|---:|
| Exact inherited structure checks | Exit 0 | 0.3650 s |
| Seven original unit tests | Exit 0 | 0.2947 s |
| Four second-strike unit tests | Exit 0 | 7.0908 s |
| Four third-strike unit tests | Exit 0 | 2.3189 s |
| Exact endpoint / finite-part checks | Exit 0 | 0.5649 s |

The metadata assembly's initial assumption that every JSON used the key
`script_sha256` was corrected to accept the reverse-tangency record's
`source_sha256`. This was a metadata issue, not a failed mathematical check.
Every frozen replay now matches its current script hash. The late-root
certificate was fully recomputed within the four new tests, including
rational tails, directed rounding, and the coefficient-box estimate.
All 59 inherited files remain unchanged, as confirmed by the preservation
check and an independent base-to-worktree diff.

No mathematical replay or artifact assembly remains outstanding in this
bounded strike. Publication must still be confirmed by the external delivery
receipt: record the actual committed hash, verify the pushed remote main
hash agrees, and report any delivery failure accurately. This audit does
not certify publication in advance.

The finite numerical data do not exclude the complete path, a complete
kappa interval, the whole lobe region, or the full Green architecture. The
bounded stop-G report preserves that limitation and the original objective.

The following mathematical escapes remain open after the documented work:

- finite, non-asymptotic points in the late-root lobe region outside the
  particular tested paths;
- anchor ratios that themselves degenerate with epsilon;
- joint limits with lambda=(1-a)/epsilon tending to zero or infinity,
  which are not covered by a theorem for fixed positive lambda;
- other reverse-Green-tangency lines in the two-variable (a,t*) family.

These are not evidence that a candidate exists. They delimit what the
current proofs do not establish. No straightforward unexecuted candidate
promotion or certification step is visible: there is no positive first-peak
lead to promote. Extending any of the listed escapes requires another
specified mathematical construction or proof, rather than more copies of
the failed bounded shots.

## Verified final report truth values

```
FIVE Q4 ZEROS CERTIFIED: NO
FIVE-ZERO NUMERICAL CANDIDATE: NO
PRIMITIVE ROOT r > 5/11: YES
KAPPA > 21636/19043 USED: YES
POSITIVE GREEN MAXIMUM FORCED: NO
RECONSTRUCTION SATURATES +2: UNKNOWN
THREE-SIMPLE+DOUBLE POINT: NO
Q4 STILL LIVE: YES
FIVE LIMIT CYCLES REALIZED: NO
```

The final report names one surviving construction problem: solve the exact
reverse Green-tangency line's intersection with the strict late-root lobe
region, using a certified ratio-determinant analysis, and then unfold an
ordinary first-maximum crossing. This is a proposed next task, not work
already performed or a claimed tangency point. Even a successful first
maximum would still need the remaining S2-S3 sign conditions to yield five
original zeros.
