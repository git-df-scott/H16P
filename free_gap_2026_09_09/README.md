# Executed free-spacing counterexample hunt

**No five-cycle counterexample found.** This follow-up released both upper-cycle spacing restrictions from the previous experiment. It then pursued a finite-field a=−2 boundary test and a search adjusted to avoid an explicitly identified centre degeneracy. Four retained numerical brackets in each final field survived independent 65-digit Cartesian replay. The extra outer sign needed for a fifth bracket remained absent.

This is a new bounded experiment based on local commit `c836b00b5c5c68ab9bc3408f19925182660c85b9`, which was itself based on repository main `498b58fc1378486601f348f67e0ab3a5362e7933`. The earlier GitHub publication remains pending explicit authorization; this follow-up does not retry that push. No global nonexistence conclusion follows from these searches.

## Search progression

The vector fields remain

\[
\dot x=ax^2+by^2+e_2xy+e_1x+(1-b)y+(b-2)/4,
\qquad \dot y=e_0-2xy.
\]

All cycle counts refer to one complete coefficient vector at a time.

| Experiment | Paired search attempts | Accepted corrected steps | Outcome |
|---|---:|---:|---|
| Direct free-spacing inequalities, strong seed | 28 | Not a corrected continuation | Line search failed; no new feasible field beat the seed |
| Direct free-spacing inequalities, medium seed | 553 | Not a corrected continuation | Improved outer signs in some trials, but old sign constraints failed |
| Corrected free-spacing ascent, strong seed | 391 | 12 | Four brackets retained; return map approached a centre degeneracy |
| Corrected free-spacing ascent, medium seed | 476 | 22 | Four brackets retained; reached the chosen parameter/gap bounds |
| Extended a=−2 boundary attempt | 225 | 3 | Four brackets retained up to a≈−1.998524; no accepted crossing of a=−2 |
| Centre-distance-normalized ascent | 498 | 21 | Four brackets retained; no positive outer witness |

There were also **28 paired attempts in two failed startup runs** before a NumPy Boolean serialization issue was repaired. Their logs and evaluations are in `startup_failure/`. The total ordinary search accounting is therefore **2,199 paired attempts**. Independent gradient checks added 56 ordinary paired evaluations. Four final-field replays added 36 arbitrary-precision paired evaluations. A pair consists of forward/backward half-return attempts; it is not itself a cycle or a certificate.

The initially declared pair caps were 700 for each direct search, 900 for each corrected search, and 500 for each targeted follow-up. Runs stopped earlier when the optimizer or return/isolation gates failed persistently. No cap was silently enlarged and no background process is needed to obtain the recorded results.

## What changed from the previous hunt

The first attempt optimized all five field parameters and seven section-witness positions. Four ordered upper witnesses had to retain +,−,+,− displacement signs, and the lower pair had to retain its opposite signs. A fifth upper witness, separated from the fourth, was pushed toward positive displacement. A successful sign pattern would have supplied four upper brackets and one lower bracket, subject to complete-return, distinctness and rigorous isolation verification.

The feasible set was very narrow. In the medium run, some trials made the outer displacement positive, but those same fields failed the old witness constraints. They were not treated as five-cycle candidates. Loss of a particular witness is also not a proof that its cycle ceased to exist elsewhere; it only invalidates this acceptance test.

The corrected method then used upper log-section coordinates

\[
s_1=t,\qquad s_2=t+g_1,\qquad s_3=t+g_1+g_2.
\]

Both g1 and g2 were free controls. With e0 fixed away from zero in each corrected run, the controls were (a,b,g1,g2,d), where d places the new witness beyond s3. Each trial solved the three upper return equalities for e1, e2 and t, then corrected the lower root and checked four actual opposite-sign brackets before acceptance.

The implicit derivative is obtained by differentiating those three root equations, including the previously checked event-time derivatives. The new centre-distance objective derivative was independently checked against corrected central differences: maximum absolute discrepancy across the four checks was approximately 2.83×10^-10. No derivative-sign census between the original roots was mistaken for an additional outer fold.

This is still a local search: it retains the seed ordering and derivative orientations, minimum log-gap sizes, coefficient bounds and numerical return gates. It does not cover all possible four-cycle fields or all continuations through nonhyperbolic cycles.

## Independent final-field results

The following are **high-precision numerical** log-height mismatches at the final outer witnesses. They are all negative; smaller magnitudes do not establish a double zero.

| Final field | Independently replayed outer D, approximately | Four original brackets replayed? |
|---|---:|---|
| Strong, corrected | −6.56852566461×10^-10 | Yes |
| Medium, corrected | −3.72040623147×10^-8 | Yes |
| a=−2 boundary attempt | −8.57747919745×10^-10 | Yes |
| Centre-distance-normalized | −6.96259760378×10^-8 | Yes |

Each replay parsed one full set of decimal coefficients as exact rationals and used it unchanged on both nests. The Cartesian Taylor method ran at 65 digits with an estimated tail tolerance of 10^-32. Estimated tails are not validated interval remainder bounds. The replay JSON files retain the coefficient strings, all eight old witness values, the extra outer witness, and the four individual sign checks.

The boundary attempt widened only the numerical engine's a guard; its ODE and positive time-rescaling formulas are unchanged. Trials approaching and crossing a=−2 failed the required old-cycle derivative, orientation or sign gates. The last accepted value was a=−1.998524392813805. This does not close the a=−2 regime or prove absence on its other side.

## Exact centre-degeneracy finding

The strong branch's nearly flat return map lay close to an explicit affine-reversible locus. The full derivation is in [CENTER_GEOMETRY.md](CENTER_GEOMETRY.md). For −2<a<0, a≠−1 and b>0, define

\[
k=\sqrt{\frac{-a(a+2)}{b(2-a)}},\quad
e_1=-\frac{(a-1)k(1-b)}{a+1},\quad
e_2=\frac{2(a-1)(a+2)}{(2-a)k},\quad
e_0=\frac{k(a+b)(b-a-2)}{4b(a+1)^2}.
\]

An explicit affine reflection reverses the field. Exact rational reflection tests support the implementation, and the accompanying algebra establishes the parameter-family identity. There is no novelty claim.

Interval calculations establish two nondegenerate elliptic equilibria on the reflection axis throughout the small box a∈[−1.908,−1.907], b∈[0.365,0.366]. Their determinant enclosures lie inside [0.6507,0.7910] and [4.3023,4.7080], respectively. Reversibility gives local period annuli at those equilibria: these are centres with nonisolated periodic ovals.

At the strong search's fixed b and e0, interval endpoint signs establish a reversible parameter point with

```
-1.907584173654 < a_center < -1.907584173652.
```

The corresponding changes from the final searched field are only a few millionths in a, e1 and e2. This explains a concrete way the outer mismatch can become very small without producing an isolated extra cycle. The actual searched field is near this locus; it was not relabeled as an exact centre.

The follow-up objective divides D by |e0| times a normalized distance to this explicit locus and imposes a minimum distance. It is a numerical distance normalization, not an analytic factorization theorem. Every acceptance decision still uses raw same-field return signs. The final outer sign remained negative.

## Evidence levels and limitations

- There is **no new interval certificate of five cycles**, or of all four cycles in any final field. The previous one-cycle interval prototype is not being counted as additional cycles here.
- The centre-locus algebra and interval geometry checks are separate from the numerical limit-cycle searches.
- A failed optimizer, a failed return, a very small residual or a centre's nonisolated ovals do not meet the counterexample target.
- Sparse section witnesses cannot exclude roots elsewhere. Fixed normal-form sectors, minimum gaps and old-cycle persistence gates may miss other configurations.
- No candidate simultaneously passed the five required same-field brackets. Therefore no candidate-level five-cycle certification was initiated.
- In the original divided-run log, the printed label `D/amp` denotes the normalized objective D/(|e0|r); the stored `outer.D` is always the raw mismatch and `center_distance` records r. The code's future print label is corrected to `objective`.
- The startup serialization failure was unrelated to the ODE calculation. Its 28 paired evaluations are included in the total rather than omitted.
- Four detailed direct-search ledger rows are missing: calls 27–28 in `strong_evaluations.jsonl` and calls 13–14 in `startup_failure/medium_evaluations.jsonl`. Run summaries charge those calls, so the accounting remains 2,199 attempts, while 2,195 detailed JSONL rows are retained. The cause of the missing tails was not established. No original rows were reconstructed. All four corrected-search ledgers and the independent final-field replays are complete.

## Reproduction

Use a copy of this directory because the scripts write their named result files. Retain the sibling `moving_cycle_2026_09_09/` dependencies: `engine.py`, `mp_replay.py`, `amplitude_result.json`, and `amplitude_resumed_result.json`. The delivered replay ZIP includes these four dependencies. Install `requirements.txt` in a Python 3.12-compatible environment.

```bash
python search.py strong
python search.py medium
python corrected_search.py strong
python corrected_search.py medium
python corrected_search.py medium --boundary
python center_geometry.py
python center_interval_checks.py
python corrected_search.py strong --divide-center
python gradient_checks.py
python replay.py strong medium medium_boundary strong_divided
```

The scripts stop on their recorded pair and iteration caps. The manifest hashes establish artifact integrity, not mathematical validity. The next substantive change would need a construction that avoids both the identified reversible degeneration and the tested persistence constraints; increasing these same local budgets alone is not supported by the result.
