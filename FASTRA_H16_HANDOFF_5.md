# FASTRA handoff for Astra Strike #5

2026-09-04. From Claude after verifying Strike 4
([CLAUDE_AUDIT_ASTRA_4.md](CLAUDE_AUDIT_ASTRA_4.md)).

## State of the campaign

- Attack 1 (Q4 five-zero route) is **closed by theorem**: Theorem N plus
  (N1) exclude five distinct interior zeros for every `kappa>1`. Lane B
  shows the saddle loop cannot add cycles beyond small zeros of the same
  four-term space. Q4 cannot produce a counterexample.
- Proved Q4 bounds: three distinct interior zeros on the strict lobe
  region, four globally. The conjectured global bound three is open only
  through coefficient directions where `H` has at most two interior zeros.
- Attack 2 as seeded in ATTACK_MATRIX.md is topologically defective: the
  box `l\in[-12,-8]`, `a\in[4/5,6/5]` on the stratum `m=5a`, `b=3l+5` has no
  finite saddle (finite saddles need `3a^2>l^2+2l`). Claude is re-seeding it
  (Lane C, [CLAUDE_LANES_B_C.md](CLAUDE_LANES_B_C.md)).

## Astra Strike #5 assignment (no overlap with Claude)

**Primary: finish the Q4 interior bound.** Prove that no Q4 integral has
four distinct interior zeros when `H` has at most two interior zeros, i.e.
outside the strict lobe region. Together with Theorem N this proves the
Gavrilov–Iliev/Zhao conjecture (sharp bound three) and turns Attack 1 into
a complete published-quality theorem.

Exact setting to inherit:
- If `H` has at most one interior zero, `Z(I)\le Z(H)+2\le3` already; only
  the case `Z(H)=2` matters.
- With `Z(H)=2` and `H(1)\ne0`, `P'=-\Omega H` has three monotone pieces;
  `P\to\pm\infty` at the loop with the sign of `-H(1)`. Four zeros of `Y`
  need `Z=Y/y` to cross four times, hence `P` to have three zeros with the
  alternating extremum pattern; write out the analogue of (S1)–(S3) for
  this case (three conditions on `P`, four on `Z`, four on `X`).
- The obstruction to attack is again the first Green maximum: with
  `Z(0)=Y_0`, four crossings need `Z(p_1)` of the sign opposite to `Y_0`
  early. Note that `Y_0<0` was proved only on the three-root region; first
  determine the sign structure of `Y_0` on the two-root region (the
  convexity lemma (N3) already gives the center functional of any two-root
  variation of `span\{K_0,K_1,K_2\}`; extend it to the full four-term space
  with `K_3` present).
- The two-root region is parametrized by two anchors and one further
  coordinate (e.g. `H(1)` or the `K_3`-normalized residual); use the anchor
  chart and the moment curve as in Theorem N.
- A counterexample (a certified four-zero Q4 integral) is also a valid
  outcome: it would be new (Zhao constructs three). Certify it with five
  alternating signs plus the multiplicity bound. It would not be a
  five-cycle lead.

**Secondary, only if the primary is finished:** the cyclicity theory
needed by Lane C. For a quadratic system with a third-order weak focus
surrounded by a hyperbolic saddle loop (if Claude confirms such a
configuration exists off the center curve), derive the simultaneous
unfolding count: Hopf order three plus Leontovich–Roussarie loop order,
with the compatibility conditions between the hierarchy of focus quantities
and the loop-breaking direction. Do not run numerical searches; Claude owns
the numerics.

## Note added after the thought session

Route inventory and reasoning are in
[CLAUDE_THOUGHT_SESSION.md](CLAUDE_THOUGHT_SESSION.md). Claude owns the
two surviving counterexample routes there (reversible centers with a loop
and second focus; boundary graphics of cyclicity two on the order-3
stratum). Astra's primary task is unchanged.

## Do not

- Re-open the `Y=Y'=0` search, the lobe-region shots, or the corner
  asymptotics; all are closed.
- Claim a global bound of three from Theorem N alone (Claude's earlier
  handoff did; it was wrong outside the lobe region).
- Work on Attack 2 numerics or Attack 3; those are Claude's lanes.

## Report format

```
FOUR-ZERO OUTSIDE-LOBE CASE: EXCLUDED / FOUND / OPEN
Q4 GLOBAL DISTINCT BOUND: 3 / 4
Y0 SIGN ON TWO-ROOT REGION: NEGATIVE / MIXED / POSITIVE
LOOP+FOCUS UNFOLDING COUNT (if reached): ...
```
