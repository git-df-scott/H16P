# Four-cycle seed ledger — September 4, 2026

All systems below are real, smooth, autonomous, planar and polynomial of degree two. Values written as fractions or decimal powers designate exact rationals. Numerical rows are freshly computed controls, not interval proofs. Stability U/S refers to the nontrivial periodic-orbit multiplier greater/less than one.

## Families, status and remaining freedom

| ID | Construction / exact family below | Configuration | Evidence | Remaining controls |
|---|---|---|---|---|
| SHI | Order-three focus plus a remote trapping cycle | (3,1) | Accepted analytical existence; original hierarchy | Three local unfolding controls, two normalized shape controls; local extra-loop route obstructed |
| GT | Exact Shi/Songling coefficients | (3,1) | Galias–Tucker accepted computer-assisted exact count four | A single fixed rational point certified, not a published parameter box; full surrounding quadratic family remains available |
| CW | Order-two focus + one preexisting intermediate + remote cycle | (3,1) | Accepted analytical parameter-hierarchy construction; visual coefficients numerical only | Two local controls; shape changes must preserve both old cycles |
| KKL | Leonov/Kuznetsov small/large-cycle construction | (3,1) | Published sufficient-domain theorem; explicit decimal seed independently numerical | Five normalized coefficients; beta is Hopf control |
| YH/YZ | Reversible two-center perturbation; local and finite-amplitude Melnikov cycles | (3,1) or reversed labeling | Yu–Han accepted analytical mechanism; Yu–Zeng concrete visual point numerical | Two center-shape controls and three perturbation coefficients; joint constraints matter |

Shi and Songling are the same example. “Four normal size” does not mean four in one nest. No accepted (4,0) or (2,2) seed is added to this ledger. Affine/time-equivalent presentations and whole parameter continua are grouped rather than falsely counted as new mechanisms.

## Exact coefficient records

**SHI / GT / visual Shi**

    x'=lambda x−y−10x²+(5+delta)xy+y²,
    y'=x+x²+(−25+8epsilon−9delta)xy.

Original Shi choice as reproduced in the primary numerical reconstruction: delta=−10^-13, epsilon=−10^-52, lambda=−10^-250. GT uses **lambda=−10^-200**, with the same delta and epsilon. Do not confuse these exact systems. Visual Shi: (delta,epsilon,lambda)=(−1/10,−1/1000,−1/50000000).

**Chen–Wang**

    x'=−delta2 x−y−3x²+(1−delta1)xy+y²,
    y'=x+(2/9)x²−3xy.

Original existence uses sufficiently small positive parameters; it does not specify a universal rational box. Visual point: delta1=1/100, delta2=1/50000.

**KKL**

    x'=y+x²+xy,
    y'=a x²+bxy+c y²+alpha x+beta y,
    (a,b,c,alpha,beta)=(−10,11/5,7/10,−363889/5000,3/2000).

The published domain template is 1<b<3, 1/3<c<1, 4a(c−1)>(b−1)², bc>1, alpha just above a(b+2)/(bc−1), and positive beta on a smaller scale. Its asymptotic qualifiers do not certify a rectangular numerical parameter box. On the beta0 target, l1 has sign K=−alpha(bc−1)−10(b+2); the incumbent K<0, whereas the recommended precursor requires K>0.

**Yu–Han / Yu–Zeng**

    x'=y(1+a1 x)+eta a10 x,
    y'=−x+x²+a4 y²+eta(b01 y+b11 xy).

Concrete Yu–Zeng equation (20):

    a1=−30/7, a4=−671/210, eta=1/100,
    a10=1/200, b01=−500001/100000000, b11=49182857/968100000.

The broader Yu–Han construction has two classes, up to exchanging nests: three local plus one remote; or two local plus one intermediate in the same nest and one remote. Representative base-shape choices in the latter are (a1,a4)=(−4,−18/5) and (−4/3,−6/5), with the perturbation hierarchies specified in the source. These are not additional verified fixed numerical points.

Primary references: [Shi original DOI](https://doi.org/10.1360/YA1980-23-2-153); [Chen–Wang original journal](https://doi.org/10.12386/A1979sxxb0069); [GT equation (1), Lemmas 2–3](https://doi.org/10.1016/j.amc.2021.126691); [KKL 2013 publication](https://doi.org/10.1007/s12591-012-0118-6); [Yu–Han 2012](https://arxiv.org/abs/1002.1055); [Yu–Zeng exact numerical reconstructions](https://arxiv.org/abs/2002.09987). Original Chinese full texts were not all retrieved; reconstructed coefficients are explicitly attributed to the later primary reconstruction.

## Fresh numerical return ledger

The script integrates a half-crossing and then the next crossing with the original orientation. It does not approach cycles by integrating spirals for trillions of steps. Multipliers below use the divergence integral on the numerically closed orbit. Tightened tolerances were checked for KKL and Yu–Zeng; the latter's weakest root moves by about 2e-8, so its listed coordinates are deliberately rounded. Cycle coordinates are section coordinates, not normal-form amplitudes.

### Visual Shi

Section: x=0, y in (0,1).

| Cycle | Coordinate | Period | Multiplier | Stability |
|---|---:|---:|---:|---|
| 1 | 0.00345363458 | 6.28228199 | 1.00000008526 | U |
| 2 | 0.00587720041 | 6.28035986 | 0.999999734865 | S |
| 3 | 0.0305712494 | 6.02185887 | 1.00415562974 | U |
| 4 | 0.0444271421 | 4.67804657 | 0.000116991858083 | S |

### Visual Chen–Wang

Section: x=0, y in (0,1).

| Cycle | Coordinate | Period | Multiplier | Stability |
|---|---:|---:|---:|---|
| 1 | 0.0575431026 | 6.26390134 | 1.0001056106 | U |
| 2 | 0.119517518 | 6.14587839 | 0.999356340185 | S |
| 3 | 0.185600079 | 5.53969064 | 1.02457589472 | U |
| 4 | 0.426900992 | 3.86060262 | 1.58718952967e-06 | S |

### KKL

Section: y=0, downward.

| Cycle | Coordinate | Period | Multiplier | Stability |
|---|---:|---:|---:|---|
| 1 | 0.683210217 | 0.728551432 | 0.999226902945 | S |
| 2 | 2.18369982 | 0.70498681 | 1.00242005513 | U |
| 3 | 15.962784 | 0.622436541 | 0.962020809942 | S |
| 4 | -3711.56081 | 0.434622845 | 11.4622677285 | U |

### Yu–Zeng

Section: y=0, initial crossing orientation.

| Cycle | Coordinate | Period | Multiplier | Stability |
|---|---:|---:|---:|---|
| 1 | -0.0288242755 | 6.27826881 | 1.00000000054 | U |
| 2 | -0.0985020073 | 6.24288679 | 0.999999996453 | S |
| 3 | -0.306152849 | 6.09644113 | 1.00000008436 | U |
| 4 | 245.330539 | 3.96554344 | 0.998996046218 | S |

### Galias–Tucker remote cycle

Section: x=0, y in (0,1).

| Cycle | Coordinate | Period | Multiplier | Stability |
|---|---:|---:|---:|---|
| 1 | 0.0426896039 | 4.76038718 | 9.11846740377e-05 | S |

### Original Shi replay

The original lambda=−10^-250 point was also independently replayed at1200 bits, Taylor orders128/144: opposite signs at7.0e-100 and7.2e-100 bracket its innermost cycle; the same middle/outer brackets and remote return were recomputed. All35 printed significant digits agree across orders. This remains numerical evidence without a remainder enclosure.

### Exact GT tiny cycles

Positive x=0 section, inner to outer. Published certified enclosures (slightly coarsened outward from Lemma 2) are:

    [7.071067811865475244, 7.071067811865475245] × 10^-75,
    [2.247805947, 2.247805948] × 10^-21,
    [6.6666601481, 6.6666601482] × 10^-8.

The remote certified bracket is [0.042689603882075,0.042689603882085]. Their stabilities are U/S/U and remote S. The exact coefficient “box” is the singleton specified above; no width in parameter space was reported in the inspected theorem.

Our MPFR Taylor replay uses scaled (x,y)=r(X,Y), 900 bits and independent Taylor orders112/128. The following displacement signs agree to 35 printed significant digits; no interval remainder is enclosed:

| r | D=P(r)−r |
|---|---:|
| 7.0e-75 | −4.3982297150e-276 |
| 7.2e-75 | +8.3239638950e-276 |
| 2.2e-21 | +2.8155162801e-115 |
| 2.3e-21 | −3.5914405781e-115 |
| 6.5e-8 | −7.1239592604e-48 |
| 6.8e-8 | +7.3048579059e-48 |

The GT Lemma2 printed sign expressions appear reversed relative to its declared f=y−P convention and Lemma3 multipliers. This does not change the bracketing proof. Stability here follows Lemma3 and the independently consistent sign pattern, not the mislabeled expressions.

## All finite equilibria and infinity directions at the reproduced seeds

Each reproduced field has exactly two real finite equilibria; its other pair is nonreal. Thus none has a finite saddle. Values below are numerical classifications backed by the exact elimination polynomials stored in each JSON, not newly interval-certified boxes.

### shi_visual

| Equilibrium (x,y) | Trace | Determinant | Type |
|---|---:|---:|---|
| (0, 0) | -2e-08 | 1 | attracting focus |
| (0, 1) | 4.89999998 | 23.108 | repelling focus |

Infinity: slope 0.06919531706; chart eigenvalues [-14.800478082879739, 9.65615495451153]; a saddle, with antipodal counterpart.

### chen_visual

| Equilibrium (x,y) | Trace | Determinant | Type |
|---|---:|---:|---|
| (0, 0) | -2e-05 | 1 | attracting focus |
| (0, 1) | 0.98998 | 2 | repelling focus |

Infinity: slope 0.399860204; chart eigenvalues [-1.2713877519350603, 2.4442502153844172]; a saddle, with antipodal counterpart.

### kkl

| Equilibrium (x,y) | Trace | Determinant | Type |
|---|---:|---:|---|
| (-6.259641416, 7.449768446) | -8.409549678 | 378.8194543 | attracting focus |
| (0, 0) | 0.0015 | 72.7778 | repelling focus |

Infinity: vertical direction; chart eigenvalues [0.3, -0.7]; a saddle, with antipodal counterpart.

### yu_zeng

| Equilibrium (x,y) | Trace | Determinant | Type |
|---|---:|---:|---|
| (0, 0) | -1e-10 | 0.9999999975 | attracting focus |
| (0.9999999938, 1.521739133e-05) | 0.0003455710112 | 3.285714238 | repelling focus |

Infinity: vertical direction; chart eigenvalues [-1.0904761904761904, 3.1952380952380954]; a saddle, with antipodal counterpart.

### gt_remote

| Equilibrium (x,y) | Trace | Determinant | Type |
|---|---:|---:|---|
| (0, 0) | -1e-200 | 1 | attracting focus |
| (0, 1) | 5 | 24 | repelling focus |

At the origin the exact trace is −10^-200; the printed ordinary floating trace/eigenvalues may round to zero. Classify it from the exact coefficient, not that rounding. At the unperturbed Shi base the remote eigenvalues are (5±i sqrt(71))/2.

Infinity: slope 0.06522985221; chart eigenvalues [-15.665063322998371, 9.669595805310292]; a saddle, with antipodal counterpart.

The Chen–Wang remote orbit reaches coordinate scale about 1.5e11 in this replay; the KKL remote reaches about1.3e4, visual Shi about170 and Yu–Zeng about245. These are sampled extents, not certified extrema. They explain the different conditioning. The three tiny GT cycles range over roughly67 orders in section scale.

## Nearby surfaces and separatrix geometry

For the Shi chart P=lambda x−y+l x²+mxy+y², Q=x(1+ax+by), the finite-equilibrium gate uses

    A2=l b²−a m b+a²,
    A1=lambda b²+ab−mb+2a,
    A0=b+1,
    Delta_finite=A1²−4A2A0.

Delta_finite=0 is the generic finite saddle-node condition; A2=0 signals a root escaping to infinity. Hopf at origin is lambda=0; Hopf at (0,1) is lambda+m=0 when its determinant is positive. On QW3 the extra-pair discriminant simplifies and the two-focus condition forces it negative. Infinity roots are those of a+(b−l)u−m u²−u³; repeated roots or vanishing radial eigenvalue identify singularity changes. Exact focus quantities are in the [Shi appendix](frontier_2026_09_04/SHI_TOPOLOGY_AUDIT.md).

For GT, the origin's trace distance to Hopf is exactly10^-200 in the lambda coordinate; remote trace is about5, finite-pair discriminant has a strict nonzero margin, and the sole infinity saddle is hyperbolic. This is a coordinate-specific surface value, **not a norm-minimized distance to every bifurcation**. The paper brackets an escape initial condition at y in [0.03689093,0.03689096]. Its return discontinuity separates itineraries and must not be counted as a cycle. A full nearest saddle-loop/heteroclinic surface was not computed.

For KKL the origin Hopf surface is beta=0, first focus coefficient vanishes at alpha=−10(b+2)/(bc−1), and all finite singularities are given by T_beta in the handoff. At b=11/5 the infinity strata are c=241/250 and c=1. The global line x=−1 is without contact (P=1), so its antipodal vertical infinity separatrices cannot be glued into the proposed extra loop. A 24-state asymptotic branch probe found the same sign of the splitting wherever its selected branches returned; this is consistent with the exact obstruction, not an interval proof. Other infinity portraits require separate analysis.

For Yu–Zeng the origin Hopf surface is a10+b01=0 (with nonzero eta), and the remote trace must be evaluated at its shifted equilibrium, not at (1,0). Its vertical infinity direction is a saddle. For each seed, cycle-fold surfaces require D=D_r=0 on a valid full-return itinerary. No exhaustive nearest-fold or separatrix atlas has been computed; treating the explicit surface equations as measured distances would be incorrect.

## Executed GT continuation experiment

Use the exact one-parameter path

    delta=−s, epsilon=−s^4, lambda=−10^8 s^16,
    10^-13 <= s <= 10^-2.

It passes through the exact GT point at its left endpoint. MPFR return signs at s=10^-8,10^-4,10^-2 bracket three local cycles. At s=10^-2 the remote return root is about0.04285849204, multiplier9.3547408e-5, so all four have numerical support at that endpoint. The smaller local scales are near7e-9 and7e-5, enormously easier than the original e-75/e-21 scales. Full-path validated continuation remains to be done. This is the explicit conditioning experiment requested for the GT seed; no fifth was found or inferred.

The cheapest local attempted addition at GT (another cycle from the maximal focus or a finite saddle loop) is obstructed. An ordinary fold elsewhere could add two cycles, but must be studied in the larger nest; splitting the remote singleton while the three persist gives forbidden (3,2). The recommended different construction is the first-order Hopf precursor, where one new local cycle would be enough.

## Candidate ledger

| Candidate | Present evidence | Next exact gate | Status |
|---|---|---|---|
| Exact GT point | Four published interval-certified cycles | Replay as precision benchmark | VERIFIED FOUR, not five |
| Relaxed GT s=1/100 | Four numerical sign/return controls | One common interval parameter/cycle certificate | NUMERICAL FOUR |
| KKL incumbent | Four robust numerical returns | Certify exact field; trace continuation control | NUMERICAL FOUR |
| KKL beta0, c=.7, alpha=−80 | One origin S root near64.5554; remote U return near−5391.14116 numerically found | Reach three origin roots while preserving remote | LIVE PRECURSOR SEARCH; not a five candidate |
| Original Q4 endpoint | Two-saddle original geometry identified | Admissible perturbation transport + joint endpoint/interior displacement | UNKNOWN |
| Resonant 2025 hemicycle | Accepted alien examples; no sharp bound on resonance | Joint center/endpoint compatibility | UNKNOWN |
| 2022/2026 six-cycle families | Explicit coefficients, incompatible local parameter sums | Any single-field return count; 3:3 route already excluded | CLAIM INVALID AS WRITTEN |
| 2026 seven-cycle claim | Abstract only; no coefficients accessible | Obtain full manuscript before numerical work | UNVERIFIED / RETRIEVAL GAP |

All numerical source, exact equations, return records, unsuccessful samples and replay instructions are in [frontier_2026_09_04](frontier_2026_09_04/README.md). No large optimization was launched.
