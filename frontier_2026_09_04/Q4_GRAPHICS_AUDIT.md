# Q4 endpoint, higher Melnikov and graphics audit

Audit date: 2026-09-04. Source canon: [commit 5bcfe11](https://github.com/git-df-scott/H16P/tree/5bcfe1172c124cb9a162a2e88c26aa5be26b40c1). This note preserves Theorem N and supplies an independent correction to the endpoint argument in [CLAUDE_LANES_B_C.md](https://github.com/git-df-scott/H16P/blob/5bcfe1172c124cb9a162a2e88c26aa5be26b40c1/CLAUDE_LANES_B_C.md).

## Decisions

- Q4 five-distinct-interior-zero mechanism: **CLOSED**, exactly as Theorem N states. Strict lobe: at most three distinct interior zeros; globally for finite kappa>1: at most four. No global-three or multiplicity extension is asserted.
- Generic Q4 “M1=0, use a larger M2 space”: **CLOSED as a distinct mechanism**. The four-function generating space already contains every first nonzero Melnikov function of an analytic quadratic arc at a fixed generic Q4 center.
- Generic Q4 endpoint/alien cycles and simultaneous interior+endpoint cyclicity: **UNKNOWN**. There is no justified closure in the canonical endpoint note. Neither 3+2 nor 4+1 has been constructed here.
- The genuinely quadratic Hamiltonian two-saddle mechanism is bounded by three for the closed annulus; standard generic Hamiltonian perturbations cannot supply five by adding separate local bounds.
- The accepted 2025 hemicycle paper is materially relevant, but its resonant a=-1 exception is a research gap, not evidence that an extra cycle exists. It should rank below Q4 endpoint compatibility and finite-distance continuation from a certified seed.

## 1. The original Q4 endpoint is not a finite one-saddle loop

The canonical original system is

    x' = y+(6+b)x²+2cxy−(2+b)y²,
    y' = −x+cx²+(8−2b)xy−cy²,
    b=2(1−rho²)/(1+rho²), c=4rho/(1+rho²), rho>0.

Gavrilov–Iliev explicitly identify a center and a node, a rational first integral, and a double ramified covering (not a birational projective coordinate change). Their Theorem 1 is for the **open** period annulus. Their section 2, pp.3–4 of the preprint, introduces the covering before the elliptic model. [GI 2009, JMAA357,69–76; arXiv0811.4602](https://arxiv.org/pdf/0811.4602).

The following compactification calculation is independent algebra, verified symbolically during this audit. Put v=y/x. If P2,Q2 are the degree-two homogeneous parts,

    R(v)=Q2(1,v)−vP2(1,v)
        =4/(1+rho²) (v−rho)(v²−2rho v−1).
    P2(1,v)=4[2+2rho²−(v−rho)²]/(1+rho²).

In the projective chart (u=1/x,v=y/x), with the usual desingularized time, the eigenvalues in radial/angular order at u=0 are:

| Direction | (−P2,R') | Type |
|---|---|---|
| v=rho | (−8,−4) | node |
| v=rho+sqrt(1+rho²) | (−4,8) | saddle |
| v=rho−sqrt(1+rho²) | (−4,8) | saddle |

Time orientation reverses between appropriate antipodal sphere charts; the saddle classification and ratios are unaffected. These are projective directions, not a count of six independent plane equilibria.

The identification of the apparent Hamiltonian saddle is also explicit. Write (U,V) for the final elliptic chart and kappa=1+rho². The inverse map on the chosen covering branch is

    y=[U²+kappa(V−U)²−1]/(8U²),
    x=[(V−U)/U+(2+b)y]/c.

The Hamiltonian saddle is (U,V)=(0,1/sqrt(kappa)). Its two nodal branch slopes are

    V=1/sqrt(kappa)+m U+O(U²), m=±rho/sqrt(kappa).

Substitution into the inverse map gives

    lim y/x = rho−sqrt(kappa) or rho+sqrt(kappa).

Thus the two branches meeting at one saddle in the elliptic chart map to **two distinct saddles at infinity in the original quadratic field**. The original closure must be studied with its infinity connection and the finite connecting orbit; the transformed one-saddle theorem does not transport automatically through this singular map. Even without identifying every orbit in the compactified boundary, this calculation decisively invalidates the asserted one-finite-saddle justification.

An arbitrary quadratic perturbation in the Hamiltonian chart is likewise not an arbitrary quadratic perturbation of the original Q4 family. Consequently neither the Hamiltonian single-loop bound nor its alien exclusion can be imported by changing names of coordinates.

## 2. Exact defects in the canonical Lane B closure

1. It reads Theorem N as a global-three bound. Canonical Theorem N section7 explicitly proves only global-four and describes the outside-lobe gap.
2. It treats the elliptic-chart saddle as an original finite saddle. The calculation above shows why that fails.
3. Its rank-two numerical check for the two endpoint coefficient functionals at kappa=2,4,9 does not establish a uniform normal form or a versal displacement family.
4. Killing H(1) defines an affine coefficient hyperplane. It does not force all three primitive roots to collide at the single corner (94/77,−17/77,1). A root at the boundary and a triple corner are different degeneracies.
5. Sampled two-functional slices with at most one interior sign change are not a theorem over all kappa and coefficient directions.
6. A surjective map to two leading scalar coefficients does not imply that every zero pattern of the full two-saddle displacement, including alien cycles, is reproduced by one Abelian integral.

The defensible replacement is **Q4 endpoint UNKNOWN**, not “LIVE with an explicit five-cycle construction” and not “CLOSED.”

The canonical chart expansion `I=c0+c1 w log(w)+c2 w+...` and its reported `c1 proportional H(1)` remain useful diagnostics if independently verified with the exact coordinate and perturbation transports. They do not count original endpoint cycles by themselves. In the strict lobe H(1)<0, hence the logarithmic coefficient is nonzero; this rules out the particular two-shadowed-zero degeneration c0=c1=0 inside that lobe. It does not rule out an alien two-cycle unfolding of the original two-saddle graphic. The four-interior case has the additional obstacle that no four-zero Q4 integral has yet been supplied.

## 3. Higher-order Melnikov: exact finite-dimensional statement

Buica–Gine–Grau prove a Bautin-function representation for a first nonzero Mk (Theorem4), then treat every quadratic center stratum (Lemma5 and Theorem6). Generic Darboux/Q4 has essential order1 and four essential parameters. Symmetric Q4 has essential order2; symmetric Hamiltonian and symmetric LV also require second order; Hamiltonian-LV has essential order3; the Hamiltonian triangle order4. These are essential orders for representing bifurcation functions, not upper bounds on the first nonzero order of an arbitrarily reparametrized arc. [CPAA14 (2015),1073–1095, arXiv1406.7612, pp.7–16](https://arxiv.org/pdf/1406.7612).

Their standard Bautin form is

    x'=lambda1*x−y−lambda3*x²+(2lambda2+lambda5)xy+lambda6*y²,
    y'=x+lambda1*y+lambda2*x²+(2lambda3+lambda4)xy−lambda2*y².

Its four generating focus quantities are

    v1=lambda1,
    v3=lambda5(lambda3−lambda6),
    v5=lambda2 lambda4(lambda3−lambda6)(lambda4+5lambda3−5lambda6),
    v7=lambda2 lambda4(lambda3−lambda6)²(lambda3 lambda6−2lambda6²−lambda2²).

For Mj=0, j<k,

    Mk(h)=v1,k h B1(h)+v3,k h³ B3(h)+v5,k h⁵ B5(h)+v7,k h⁷ B7(h).

Generic Q4 imposes v1=v3=v5=v7=0 with
`lambda2 lambda4(lambda3−lambda6) != 0`; the four coefficient directions have full rank. The canonical elliptic integral is the same generating space in a different chart (GI p.4, citing Iliev 1998 Thm2(iii); Zhao2011 section2 explicitly discusses the first nonzero Mk).

**Inference:** Merely choosing M1=0 cannot create a fifth independent generic-Q4 Abelian function. Any proposed higher-order escape must identify (i) a nongeneric center intersection not covered by fixed finite kappa; (ii) multiplicity splitting not ruled out by the distinct-zero theorem; or (iii) endpoint/nonuniform bifurcation. These cannot be relabelled as an independent fifth-zero M2 lane.

Theorem N bounds distinct real zeros. A multiple zero of Mk can split under higher corrections, so one must prove a bound with multiplicity, or a versal real-zero theorem, before converting it to a uniform bound on all nearby cycles.

An exact transverse family useful for original-field calculations (derived directly from the quantities above) is

    lambda=(tau,1+u,3,−10+v,w,1).

At the generic Q4 base lambda=(0,1,3,−10,0,1),

    (dv1,dv3,dv5,dv7)=(d tau,2d w,−20d v,80d u).

It gives four independent normal controls using actual quadratic vector fields, avoiding a perturbation-degree error in the elliptic chart.

## 4. Quadratic Hamiltonian graphics: the additive loophole is explicitly obstructed

Gavrilov–Iliev2015 prove cyclicity <=3 for every nondegenerate Hamiltonian two-saddle loop under arbitrary quadratic deformation (Thm1), and <=3 for its closed annulus (Thm3). Their Thm4 gives mutually exclusive interior/endpoint alternatives, including higher Mk. Thm5 bounds nearby systems by3 when the cubic Hamiltonian has four distinct critical points and three distinct critical values. These hypotheses matter. [Ann.IHP C32 (2015),307–324, arXiv1306.2340, Thms1,3–5, pp.3,16–19](https://arxiv.org/pdf/1306.2340).

For their normal form,

    Mk(h)= integral_delta(h) [alpha+beta x+gamma/x] y dx, k>=2.

If gamma!=0, endpoint cycles are absent and the closed-annulus bound is2. If gamma=0 and Mk(0)=0, the open annulus supplies none and the loop at most2. The remaining pure beta case can have upper bound3 at the loop, with none from the open annulus. Therefore adding “two interior plus three endpoint” to get five contradicts their theorem in its exact family.

The paper's AppendixA gives a truly quadratic alien-cycle control:

    H=y(x²+y²/12−1),
    x'=Hy,
    y'=−Hx−epsilon[(16+c x−pi sqrt(3)y)y+mu1+mu2 y], c>16.

The upper segment-plus-half-ellipse loop has two cycles although its first Melnikov function has at most one zero near the endpoint. Proposition8 establishes cyclicity2. This is an excellent **method-validation control**, not a fifth-cycle candidate: Thm5 blocks five sufficiently near it.

For a finite Hamiltonian homoclinic loop, the established cyclicity is2 (Horozov–Iliev1994; see also GI2015 references16,19,21). Han1997 proves2 for specified integrable non-Hamiltonian homoclinic cases with an exception; the primary PDF was found but retrieval returned429, so this audit does not resolve that exception from the abstract. That unresolved bibliographic point cannot rescue applying Han to the wrong original graph.

## 5. Accepted 2025 infinity result and its exact exceptional lane

Marin–Villadelprat study

    x'=(b−2)/4+(1−b)y+a x²+b y²,
    y'=−2xy,             (a,b) in (−2,0)x(0,2).

The centers are (0,1/2) and (0,(b−2)/(2b)); both annuli have infinity hemicycles. ThmB: each hemicycle has cyclicity2 when a!=−1, at least2 when a=−1. ThmC: simultaneous cyclicity is3 on K1\{a=−1}, 2 on K2, and at least3 on a=−1, where

    K1: a+b<=0 OR a−b+2<=0,
    K2: a+b>0 AND a−b+2>0.

ThmD proves alien bifurcations at (-1,1),(-1/2,1/2),(-1/2,3/2). The two individual bounds cannot be added. [JDE433 (2025),113281; arXiv2501.16924, ThmsB–D pp.6–9](https://arxiv.org/pdf/2501.16924), [publication record](https://portalrecerca.uab.cat/en/persons/jordi-villadelprat-yague/).

The exact full local quadratic quotient family (Lemma3.1, p.14) is

    x'=(b−2)/4+epsilon1*x+(1−b)y+a*x²+epsilon2*x*y+b*y²,
    y'=epsilon0−2*x*y.

The upper difference map satisfies (Prop3.2)

    Du(s;mu)=delta_u+Delta0_u*s^lambda+remainder,
    lambda=−(a+2)/a,
    delta_u=rho0(mu)*epsilon0, rho0>0,
    Delta0_u=−kappa01(mu)[2sqrt(b(a+2)/(a(b−2)))*epsilon1+epsilon2]
             +kappa02(mu)*delta_u, kappa01>0.

The secondary exponents lambda+1 and 2lambda coalesce at a=−1. The upper-bound/division theorem3.5 excludes that line; this is the precise technical gap, not a general assertion of arbitrarily high cyclicity.

**Independent compatibility obstruction:** At the upper equilibrium the linearized trace is

    tr = 2(a−1)*epsilon0+epsilon1+epsilon2/2+O(||epsilon||²).

On a=−1 and b=1/5 (the intersection with the maximal center-cyclicity line3a+5b+2=0), the leading endpoint degeneracy delta_u=Delta0_u=0 gives

    epsilon0=0, epsilon2=−(2/3)epsilon1,
    tr=(2/3)epsilon1.

Thus the naive “three small cycles at this center plus two resonant endpoint cycles” cannot use a nonzero first normal direction satisfying both degeneracies. An integrable leading direction or higher parameter blow-up would have to be handled explicitly; this is an obstruction, not a proof that all arcs are impossible.

At b=1 the trace obstruction vanishes, but this is the quadratic LV double center. Francoise–Gavrilov2022 prove total compact-annulus cyclicity2, allowing only (i,j) with i+j<=2; their conclusion explicitly leaves boundary cycles untreated. Their double Bautin ideals are

    B1=<lambda1,lambda3,lambda2 lambda5>,
    B2=<lambda1+lambda3+lambda1 lambda2,lambda5,lambda3 lambda4>.

[CCM24 (2022),2150064, arXiv2011.08316, Thms1,11 and conclusion](https://arxiv.org/pdf/2011.08316). Hence a putative five-cycle mechanism there needs compatibility of two compact cycles with three endpoint cycles, not a sum of separately attained upper bounds. No such compatible arc was found.

## 6. Graphic ledger: which numbers are actually available

| Graphic / perturbation class | Rigorous statement retrieved | Can it by itself justify fifth cycle? |
|---|---|---|
| Quadratic weak focus/nondegenerate center | local cyclicity <=3 | No; need independent compatible outer mechanism |
| Finite Hamiltonian saddle loop, quadratic perturbations | cyclicity2 | No; closed-annulus restrictions apply |
| Nondegenerate Hamiltonian two-saddle loop | <=3; explicit alien2 | No five in associated Thm5 neighborhoods |
| Hamiltonian triangle open annulus | cyclicity3; higher essential order4 | Does not by itself bound every boundary cycle; verify triangle unfolding separately |
| Two infinity hemicycles, 2025 family, a!=−1 | individual2; simultaneous2 or3 | Not additive4; interior compatibility still needed |
| Same, a=−1 | lower bounds2 individual,3 simultaneous; upper bound omitted | Unknown; resonance alone is insufficient evidence |
| Reversible infinity bicycle | Swirszcz1999 Thm1 gives2 or3 on a curve | Not generic Q4; applies within reversible family base stratum |
| Original generic Q4 infinity endpoint | explicit two-saddle geometry; no applicable sharp joint cyclicity theorem retrieved | UNKNOWN |
| H3_4,H3_5 generic elementary hemicycles | cyclicity2 under stated regular-transition/saddle-quantity hypotheses | Specific finite bound, not all graphs |
| H3_6 and listed saddle-node/nilpotent graphics | finiteness under exact hypotheses | Existential finite cyclicity is not <=4 |

For the elementary graphics, source is Dumortier–Guzman–Rousseau2002, QTDS3,123–154, and Dumortier–Ilyashenko–Rousseau, Cor4.2 for I2_10a. The latter supplies the exact family

    x'=delta*x−y+A*x², y'=x+gamma*y+x*y,
    1/2<A<1, (delta−gamma*A)²−4A=0

at the relevant finite double singularity. It proves finite cyclicity, not an explicit numerical maximum. [Author PDF, section4.2](https://dms.umontreal.ca/~rousseac/DIR.pdf).

Swirszcz's title is *Cyclicity of Infinite Contour around Certain Reversible Quadratic Center*, JDE154(2) (1999),239–266, DOI10.1006/jdeq.1998.3538. The 2025 primary paper pp.6–7 explicitly states its hypotheses and distinguishes the bicycle from its own hemicycle. It identifies the Q4 intersection (-4,2), not the entire generic Q4 component. Original1999 text not retrieved; no stronger claim is made.

## 7. Other modern developments and degree gates

- Buzzi–Gasull–Santana, JDE429 (2025),646–677, establishes lower cyclicity conditions for hyperbolic polygons. Its polynomial realization generally raises degree substantially, so it does **not** manufacture additional quadratic cycles. [Primary preprint2407.20721](https://arxiv.org/abs/2407.20721), [accepted publication](https://doi.org/10.1016/j.jde.2025.02.061).
- Arakaki–Santana, *On the cyclicity of persistent hyperbolic polycycles*, arXiv2504.07225, gives explicit leading return-map asymptotics for persistent polycycles. It supplies tools, not a quadratic five example. [Primary source](https://arxiv.org/abs/2504.07225); DOI10.1007/s10884-025-10469-9 (acceptance metadata should be checked against journal if included as accepted).
- Marin–Villadelprat, *Derivatives of the Separation Function of Generalized Saddle Connections*, QTDS24 (2025),227, provides general separation derivatives and quadratic unbounded-polycycle applications in Props6.2,6.4,6.5. These are suitable for the missing original Q4 transport calculation. [Journal primary text](https://link.springer.com/article/10.1007/s12346-025-01379-8).
- Haibo Lu, arXiv2607.13785 (July2026), claims local uniform finite cyclicity of H3_14 in the full quadratic family. This is a **PREPRINT/UNVERIFIED** contribution to one graphic, with an existential bound. It neither proves H(2) finite nor gives five cycles. [Primary abstract/version history](https://arxiv.org/abs/2607.13785). Full proof was not audited here.

## 8. Bounded Q4 endpoint experiment worth doing next

This is an **eligibility test for a construction**, not a claim of an already viable five-cycle example.

Use the actual Bautin-form Q4 base and four transverse controls above. If a base modulus is needed, use

    lambda3=2+r², lambda6=1, lambda2=r+u,
    lambda4=−5(1+r²)+v, lambda1=tau, lambda5=w,
    r in [1/2,2], ||(tau,u,v,w)||_infty<=2^−8.

All coefficients are quadratic and rational whenever parameters are rational. The normal rank is nonzero at every r>0. This bounded r interval is an experiment, not a global exclusion.

1. Derive the compactified vector field, select the two saddle branches of the **original** period-annulus boundary, and construct two Dulac transition maps with one common orientation and section coordinate. Use the 2025 separation-derivative formulas where their hypotheses hold. Confirm which infinity connection persists under all polynomial perturbations.
2. Calculate the linear maps from (tau,u,v,w) to the four elliptic generating coefficients and to the genuine separatrix displacement/saddle ratios. The gate is exact or interval full rank on the proposed slice; no approximate surjectivity inference.
3. Choose three interior section locations in a fixed compact subannulus, initially t=1/4,1/2,3/4. Solve the three homogeneous first-Melnikov equations for the projective perturbation direction. Certify that the three roots are simple and disjoint. This uniquely determines a direction if the 3x4 evaluation matrix has rank3.
4. Allow only the base r and the three root locations within t1 in[1/8,3/8],t2 in[3/8,5/8],t3 in[5/8,7/8], with gaps>=1/16. Continue the actual finite connection splitting through zero while retaining those three simple interior roots. Stop after this finite box, reporting failure as a box failure only.
5. At any connection-zero point, compute the two-saddle return map and determine whether an unfolding produces **two additional** simple zeros near the original endpoint while retaining the three interior zeros. Do not require c1=0 if the second cycle is alien; conversely do not infer alien2 merely from c0=0.
6. Certification gate: one exact rational quadratic field, five disjoint transverse return intervals, validated complete returns, opposite displacement signs at each interval's ends, and interval derivative separation from1 (or interval Newton on P(s)−s). Every return must remain on the asserted section itinerary in compactified charts. Failure of the loop unfolding at step5 terminates this lane before expensive certification.

This is stronger than the old Lane B test because it uses the original two-saddle return and exact admissible perturbation rank. The unresolved work is step2 and the simultaneous two-endpoint unfolding at step5. These are concrete mathematical obligations, not details already proved by Theorem N.

## 9. Suggested relative assessment

Scores are judgment, not mathematical evidence. Scale0–5: plausibility, parameter freedom, open status, numerical tractability, certificate tractability, proximity to certified four.

| Mechanism | Six scores | Audit judgment |
|---|---|---|
| Original Q4 endpoint +3 simple interior | 3,3,4,2,2,1 | Best new analytic compatibility target among lanes audited here; UNKNOWN |
| Resonant a=-1 hemicycle at generic b | 2,3,4,3,2,0 | Upper-bound gap plus explicit inner/outer obstruction; exploratory |
| Resonant LV b=1, compact2+endpoint3 | 2,3,4,3,2,0 | Needs genuine joint three-endpoint degeneracy; no candidate |
| Generic Q4 higher-order larger-space claim | 0,0,0,3,3,0 | Closed in claimed form |
| Quadratic Hamiltonian two-saddle | 0,2,0,4,4,0 | Useful alien control; five locally obstructed |
| Finite-distance continuation from certified4 | 3,4,4,4,4,5 | Best construction/certification proximity, but local protected box may exclude a fifth |

A focused Q4 endpoint compatibility strike can now be justified by newly corrected geometry, not by preexisting repository code. If root prioritizes producing an actual fifth-cycle candidate quickly, finite-distance continuation from a rigorously known four still has a better certificate starting point. Neither route currently has a fifth-cycle candidate.
