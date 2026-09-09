# H16P coordination protocol — 2026-09-06 (Fable, auditor)

## Goal
A real planar quadratic vector field with >= 5 limit cycles. By Zegeling 2024 (only (n,0)/(n,1) distributions; two independent rotation
parameters, one per nest, Lemma 6.4) this is EQUIVALENT to: a quadratic nest with FOUR limit cycles in a field that also has a second
focus (the fifth cycle is then obtained by rotating only the other region through its Hopf angle). No (4,0) or (4,1) has ever been
exhibited (14 constructions since 1979, all (3,1)). No theorem bounds a nest by 3; Bautin's 3 is local. See LIT_A/B/C and H16P_SUMMARY in this folder.

## The instrument: the Andronov–Hopf curve of a uniformly rotated family
For a base field X=(P,Q), the rotated family is X_beta = (P cos b - Q sin b, P sin b + Q cos b) (still quadratic; all equilibria fixed).
Duff/Perko: the displacement D(s, beta) on any fixed transversal section from the focus is STRICTLY MONOTONE in beta at every s
(d_beta = -(w/|f|) int_0^T exp(-int div) |f|^2 dt). Hence for each s there is a unique beta*(s) closing the orbit through s.
Cycles of X_beta = level set {s : beta*(s) = beta}. Number of cycles in the nest at any beta = crossings of a horizontal line.
THREE cycles <=> beta* has 2 interior extrema (Cherkas–Artés–Llibre 2003 call it the Andronov–Hopf function AH(x) and publish one).
FOUR cycles in the nest <=> beta* has THREE interior extrema with overlapping height windows. This has never been searched.
Use the EXACT beta*(s) (root-solve in beta per s; monotone so bisection/Newton is safe), not the linearization -D/D_beta.

## How extrema of beta* can be born as the base field moves (Perko 1992 Thm 2,3; Perko 1995 Thm 4.1-4.3)
(a) at the focus end s->0: capped by Bautin (<=3 small cycles). Closed.
(b) at the outer end s->s_max (the nest's boundary graphic): needs a neutral graphic. Dead around weak foci of order 3 (Llibre–Schlomiuk 2004
    Thm 16) and for finite loops around any trace-zero focus (repo Proposition A). Open around STRONG foci and for order-2 focus + graphic
    through infinity (needs Artés–Llibre–Schlomiuk 2006, paywalled).
(c) in the interior: two extrema born together at a CUSP of beta* = a limit cycle of multiplicity 3 (D=D_s=D_ss=0). A multiplicity-4 cycle
    (D=D_s=D_ss=D_sss=0, swallowtail) with Perko's nondegeneracy conditions FORCES four simple cycles in the nest nearby (Perko 1995 Thm 4.3).
    Triple cycles are KNOWN to exist at small amplitude in the Bautin unfolding of a third-order weak focus:
    D(r) ~ V1 r + V3 r^3 + V5 r^5 + V7 r^7, and D = V7 r (r^2 - r0^2)^3 gives a triple cycle at r0 when
    V5 = -3 r0^2 V7, V3 = 3 r0^4 V7, V1 = -r0^6 V7. Nobody has continued this cusp manifold to normal amplitude.

## Evidence rules (binding; violations get the lane's results discarded)
1. A cycle is claimed ONLY from a sign change of the displacement D on a bracket [s1,s2] on a ray/section from the focus, with
   min(|D(s1)|,|D(s2)|) above a two-tolerance noise estimate (recompute at looser rtol; noise = 10*|difference| + 5e-12*s).
   |D| ~ 0 is NEVER evidence of a cycle. A near-fold is an organizer, not two cycles.
2. Every count must be reproduced by a second, independent integrator (different method or different section) before it is reported.
3. Anything with >= 4 sign changes in ONE nest is a TRIGGER: stop, recompute all bracket endpoints in binary128 or mpmath dps>=40,
   write TRIGGER_<lane>_<timestamp>.json with the exact rational coefficient vector, section, brackets, both engines' values, and push immediately.
   Do NOT announce a counterexample. The auditor decides.
4. Section completeness: cycles of quadratic systems are convex and enclose exactly one focus, so a ray from the focus meets each cycle
   once, transversally. Use rays from the focus (or the paper's own section for Cherkas seeds). Record the nest domain end s_max
   (where the return fails) and never count beyond it.
5. Every run is a ledger: append-only JSONL with coefficients (exact strings), parameters, results, engine name+hash, tolerances,
   wall time. Failures are recorded as UNRESOLVED, never as "no cycle".
6. Precision: double for sweeps; long double / binary128 (repo: audit/fable_engine, fastra_d1_2026_09_05/matching_quad.cpp,
   astra_afternoon_2026_09_05/full_return128.cpp on the branches named in H16P_SUMMARY.md §0) for anything that decides a trigger.
7. Validate the engine FIRST on the seeds in SEEDS.json: must reproduce Cherkas rows 1–8 cycle positions (to ~1e-2) and the KKL control.
   Row 4's published AH polynomial must be reproduced qualitatively (2 extrema on [0.6,0.9]). A lane whose engine fails validation must not sweep.
8. Push to your own branch every 30–45 minutes with a REPORT_<lane>.md at the branch root: what ran, ledger sizes, max extrema seen,
   best candidates, open problems, next step. The auditor reads only that file and the ledgers.
9. Do not spend effort on: re-reading the whole repo (the summaries here suffice), Q4/Abelian-integral routes, the reversible
   reseed, the KKL fold-surface continuation at radius > 1e10, or any near-integrable seed (D1) — all already exhausted or capped by Melnikov bounds.
10. Certification standard for a candidate (Lane 3 builds it): five disjoint positively/negatively invariant annuli bounded by polygons with
    rational vertices, each edge exactly transversal (sign of X·n along the edge decided in exact rational arithmetic), no equilibrium inside,
    => Poincaré–Bendixson gives five periodic orbits; a quadratic field without a center has only isolated periodic orbits.
