# Global H(2) status and claims audit — 2026-09-04

Scope: smooth real autonomous planar polynomial vector fields of total degree at most two. The notes below use primary manuscripts, publishers, institutional repositories and author-posted manuscripts. “Published” and “mathematically accepted” are separate predicates. The audit found a 2026 journal six-cycle claim with an elementary fatal error and an August 2026 seven-cycle preprint with inaccessible full text; neither supplies an accepted record.

## Status supported by the retrieved primary literature

* **Best accepted lower bound: 4.** No accepted five-cycle quadratic example was located. This is a literature status conclusion, not a proof that such a system does not exist.
* **Uniform finiteness H(2)<infinity remains open in the accepted literature located.** Finiteness for each individual quadratic field is a different theorem.
* Gasull–Santana, *A note on Hilbert 16th problem*, **Proc. AMS 153(2) (2025), 669–677**, Introduction p.1 of postprint, explicitly states unknown finiteness for H(2) and best lower bound 4. Their Theorem 2 says maximal finite counts can be realized with hyperbolic limit cycles; when H(n)=infinity one can obtain arbitrarily many hyperbolic cycles. This supports a rational-coefficient hyperbolic certificate target: after obtaining a hyperbolic finite collection, density of rational coefficients and persistence produce a rational example. [Accepted postprint](https://ddd.uab.cat/pub/artpub/2025/309367/GasSan24-Postprint.pdf), [institutional publication record](https://ddd.uab.cat/record/309367?ln=en), [arXiv v2](https://arxiv.org/html/2407.13465v2), [DOI](https://doi.org/10.1090/proc/17116). The journal year is 2025; the preprint revision is 2024-10-01. Do not describe it solely as a 2024 journal article.
* Craciun–Erban, *Planar chemical reaction systems with algebraic and non-algebraic limit cycles*, **J. Math. Biology 90:64, 2025-05-22**, Section 3/Lemma 7, again records H(2)>=4 and uses the finite/infinite dichotomy, without supplying finiteness. [Publisher](https://link.springer.com/article/10.1007/s00285-025-02221-0).
* Artés–Cairó–Llibre, *Phase Portraits of the Families V to X for the Quadratic Polynomial Differential Systems*, **QTDS 25:145, published 2026-08-07**, DOI **10.1007/s12346-026-01563-4**, Introduction Section 1, explicitly states that the maximum remains unknown, examples reach four, and four remains a conjecture. This is a particularly recent expert primary status statement, not merely a pre-2025 survey. [Publisher full text](https://link.springer.com/article/10.1007/s12346-026-01563-4).
* Arakaki, *New lower bounds for limit tori in 3D polynomial vector fields*, **V Symposium on Planar Vector Fields 2026**, poster Introduction, explicitly lists the best global planar lower bound H(2)>=4. This is current author conference material, not a proof of nonexistence of a newer preprint. [Poster](https://www.gsd.uab.cat/symp2026/data/uploads/slides/lqueiroz.pdf).
* Bamón, **Publ. Math. IHÉS 64 (1986), 111–142**, proves finite numbers for individual quadratic fields; do not exchange `for every X there exists N_X` with `there exists N for every X`. The Introduction also records the Petrovskii–Landis <=3 claim, withdrawal in 1959, and contradiction by Shi's four-cycle construction. [Primary archived paper](https://www.numdam.org/item/PMIHES_1986__64__111_0.pdf), [record](https://numdam.org/articles/10.1007/BF02699193/).

## A substantive 2026 six-cycle claim, and its direct failure

**L. Malyarets, O. Dorokhov, A. Voronin, I. Lebedeva, S. Lebedev, T. Denysova**, *Topology of Quadratic Systems of Differential Equations with Six Limit Cycles*, **Mathematica Montisnigri 65 (2026), 23–35**, received 2025-11-12. Status: **PUBLISHED CLAIM; proof invalid as written; no independent six-cycle certification or expert acceptance found**. Its source family is

\[
\dot x=\lambda x-y-\frac73 n x^2+5a xy+n y^2,\qquad
\dot y=x+a x^2-2nxy,\quad n\ne0.
\]

The paper gives this family on pp.29–31 from equations (23)–(25), and computes a third focus quantity in (27). Its p.30 count explicitly combines the origin at lambda=0 with the second equilibrium at lambda=-5a/n. The p.31–32 alternative calculation discards the same trace term after translating. [Full journal PDF](https://www.montis.pmf.ac.me/allissues/65/Mathematica-Montisnigri-65-3.pdf), [institutional publication record](https://repository.hneu.edu.ua/handle/123456789/39936?locale=en). Downloaded locally and all 13 pages text-inspected.

### Independent exact algebraic check

This is our calculation, not an attributed external refutation. Let E0=(0,0), E1=(0,1/n). Their Jacobians are

\[
J(E_0)=\begin{pmatrix}\lambda&-1\\1&0\end{pmatrix},\quad
J(E_1)=\begin{pmatrix}\lambda+5a/n&1\\-1&0\end{pmatrix}.
\]

Thus their traces vanish at **different** values of the same global coefficient lambda. One cannot add cycle counts computed in neighborhoods of these two different parameter values. Simultaneous trace-zero requires a=0. At that value the family is reversible under (x,y,t)->(-x,y,-t); both elementary monodromic singularities are centers, not order-three weak foci. The source's own displayed third focus expression also has a factor a. This defeats the stated simultaneous double-Bautin construction.

There is a second independent logical error: a nonzero third Lyapunov quantity with the first two zero describes an order-three weak focus and its potential cyclicity under an adequate unfolding, not three present periodic orbits at the degenerate parameter itself. No concrete common perturbation giving six disjoint return-map roots is supplied. Fixing the first error alone would not repair this omission.

For potential numerical reuse, all other finite equilibria lie on y=(1+a x)/(2n), where

\[
(33a^2-28n^2)x^2+(12n\lambda+30a)x-3=0.
\]

At a=1,n=3 and either lambda=0 or lambda=-5/3, this quadratic has negative discriminant (-1728), so only E0,E1 are real. At lambda=0 the second focus has eigenvalues (5 +/- i sqrt(11))/6; at lambda=-5/3 the origin has eigenvalues (-5 +/- i sqrt(11))/6. The family is explicit and interval-testable after a **single** rational lambda is chosen, but its claimed 3:3 architecture is closed by the distribution theorem below and its base points are not six-cycle seeds.

## Antecedents of that six-cycle claim

**Voronin–Lebedev 2022**, *Example of a quadratic system of two differential equations with six limit cycles*, author-posted preprint June 2022, four pages. Equation (1), p.1:

\[
\dot x=\lambda x-y-(7+\varepsilon_1)x^2+(5+\varepsilon_2)xy+3y^2,\qquad
\dot y=x+x^2-6xy.
\]

The two half-traces are mu0=lambda/2 and mu1=(lambda+(5+epsilon2)/3)/2. At epsilon1=epsilon2=0, the p.3 conclusion adds the local calculations at lambda=0 and lambda=-5/3. Hence the incompatibility predates the journal version. The asserted third focus quantities are +/-775/64, with no common fixed parameter triple, cycle brackets or interval proof. Status **PREPRINT; invalid simultaneous-count argument**. Exact coefficient family survives and can be tested. [Author landing page](https://www.researchgate.net/publication/361398861_Example_of_a_quadratic_system_of_two_differential_equations_with_six_limit_cycles), [actual four-page PDF inspected](https://www.researchgate.net/profile/Anatoliy-Voronin-2/publication/361398861_Example_of_a_quadratic_system_of_two_differential_equations_with_six_limit_cycles/links/62ae13f8a920e8693efeff76/Example-of-a-quadratic-system-of-two-differential-equations-with-six-limit-cycles.pdf).

**Voronin–Lebedev 2023**, *Quadratic System of Two Differential Equations with Six Limit Cycles: Two Approaches to Problem Analysis*, February 2023, DOI **10.13140/RG.2.2.32572.92808**. Author abstract claims the same 3:3 count by two methods. Indexed institutional first page reproduces the same epsilon family; RG full-text extraction is font-garbled, and direct institutional PDF retrieval timed out. Status **PREPRINT; antecedent of 2026 invalid proof**, not a separately verified six-cycle example. [Author record](https://www.researchgate.net/publication/368570270_Quadratic_System_of_Two_Differential_Equations_with_Six_Limit_Cycles_Two_Approaches_to_Problem_Analysis_AV_Voronin_SS_Lebedev), [institutional PDF location](https://repository.hneu.edu.ua/bitstream/123456789/29182/1/Voronin_A.V._Quadratic_system_ResearchQate.pdf).

**Malyarets–Voronin–Lebedev 2024**, *Mathematical model of the dynamics of the interaction of goods prices on adjacent markets: the problem of limit cycles*, Riga conference proceedings pp.47–51, DOI **10.30525/978-9934-26-398-9-13**, is explicitly cited as an antecedent by the 2026 article and its first page claims six cycles for a two-market model. Status **CONFERENCE CLAIM**, not independent confirmation. [Institutional primary PDF](https://repository.hneu.edu.ua/jspui/bitstream/123456789/34521/1/1.Malyarets_%D0%9B%D0%B8%D1%82%D0%B2%D0%B0_2024.pdf). Exact interior equations not fully retrieved in this audit; do not assert they are identical without checking.

## August 2026 seven-cycle claim: real lead, unavailable candidate

**Manuel Hernandez Rosales**, *A Multiscale Hilbert-Fourier Variety Approach to Hilbert's 16th Problem: Establishing the Upper Bound H(2)=7 for Quadratic Systems*, ResearchGate **PREPRINT**, August 2026. The author-posted abstract claims both an upper bound and attainment by a specific field, via a Fourier variety organized into 4+2+1 scales. The accessible page says **No file available**, offers an author request, and lists zero citations. No arXiv manuscript, journal, explicit coefficients, independent review or certificate was located by author/title searches. Status **CLAIMED/UNVERIFIED, FULL TEXT UNAVAILABLE**. Do not call it refuted from the abstract. It supplies no reproducible candidate yet. A finite Fourier truncation alone would not certify existence or exhaust an infinite Fourier variety; full proof must establish convergence, isolation modulo phase, minimal-period distinctions and complete scale coverage. [Author preprint record](https://www.researchgate.net/publication/412197587_A_Multiscale_Hilbert-Fourier_Variety_Approach_to_Hilbert%27s_16th_Problem_Establishing_the_Upper_Bound_H2_7_for_Quadratic_Systems).

## The distribution theorem that changes the attack map

**André Zegeling**, *Nests of limit cycles in quadratic systems*, **Advances in Nonlinear Analysis 13 (2024), 20240012**, accepted 2024-04-10, published 2024-06-12, DOI **10.1515/anona-2024-0012**. **Theorem 1.2, printed p.3, proof Section 7**, states that distributions are (n,0) or (n,1), filling essential omissions in Zhang Pingguang 2002. Every quadratic cycle surrounds exactly one focus and there are at most two nests. This closes (2,2), (3,2), (3,3), and three-nest five-cycle arrangements; the remaining five-cycle distributions are **(5,0) and (4,1)**. This is a published repaired proof, not merely the historically unchecked Zhang claim. Full PDF obtained and exact theorem confirmed. [Publisher full text](https://www.degruyterbrill.com/_language/de?uri=%2Fde%2Fdocument%2Fdoi%2F10.1515%2Fanona-2024-0012%2Fhtml), [archived full PDF](https://d-nb.info/1332906729/34), [2023 author talk](https://www.gsd.uab.cat/aqtde2023/data/uploads/slides/3-zegeling.pdf).

## Global upper-bound claims and actual criticisms

### Information geometry (2024) and the user-supplied rebuttal

Da Silva–Vieira–Leonel, **Entropy 26(9) (2024), 745**, Theorem 4 claims H(n)=2(n-1)(4(n-1)-2), giving H(2)=4. Buzzi–Novaes, **arXiv:2411.09594v1, 2024-11-14**, refutes the all-degree formula using established superquadratic lower growth and identifies its faulty curvature-based replacement for the definition of a limit cycle. Sections 2–3 provide counterexamples. Especially their quadratic Example 4 is x'=-y+x², y'=x+xy: a center with no limit cycles but satisfying the proposed curvature criterion. Their cubic Example 1 has exactly the unit-circle cycle despite failing that criterion. These are objections to the proof and all-degree theorem, not a five-cycle quadratic example and not a refutation of the conjecture H(2)=4. Status **PUBLISHED CLAIM / EXPLICITLY REFUTED BY PRIMARY MATHEMATICAL NOTE**. No retraction notice was found; mathematical refutation is not publisher withdrawal. [Original DOI](https://doi.org/10.3390/e26090745), [rebuttal full text](https://arxiv.org/html/2411.09594v1), [rebuttal metadata](https://arxiv.org/abs/2411.09594). Rebuttal metadata has no journal reference found.

### Pedregal's variational program

* **Llibre–Pedregal arXiv:1411.6814**, initial 2014 manuscript, current v3 is **WITHDRAWN**. Author note expressly says the cycle-counting method contained a mistake. Initial Theorem 4 claimed <=6 for quadratic systems; this was an upper bound, not a six-cycle construction. [Withdrawal](https://arxiv.org/abs/1411.6814), [initial full-text mirror](https://citeseerx.ist.psu.edu/document?doi=6cb103cbfb15155a8b75d20c248a171e99b4bca1&repid=rep1&type=pdf).
* **arXiv:1904.01292**, Llibre–Pedregal Part I, current v4 2020-10-08 is **WITHDRAWN**, comment says an improvement of the bound was found. **arXiv:1904.01300**, Pedregal Part II, v2 2020-10-19 also **WITHDRAWN**, comment says some N=3 bounds can be improved. These comments should not be replaced by the 2014 error explanation. [Part I](https://arxiv.org/abs/1904.01292), [Part II](https://arxiv.org/abs/1904.01300).
* **Pablo Pedregal arXiv:2103.07193v3**, **2024-08-30**, claims degree-four global bounds (Theorem 1.2) and **H(2)=4 (Corollary 1.5)**. Central Theorem 1.1 uses a bound 1+(n-1)^2(M+N) based on components of divergence=0 and contact points. v2 was withdrawn one day earlier; **v3 is active**, so calling the whole current claim withdrawn is incorrect. No journal reference or independent acceptance was located. **PREPRINT / UNVERIFIED**, not an obstruction to a construction. [Version history](https://arxiv.org/abs/2103.07193), [full v3](https://arxiv.org/html/2103.07193v3).

Independent preliminary check on v3: Theorem 1.4 as printed gives zero as the maximum number of components for a degree-one real algebraic curve and gives one for degree-two. A real affine line has one component and a hyperbola has two. The application in the proof of Theorem 1.2 similarly writes M<=0 at n=3. This is a concrete erroneous supporting statement as written, not an independent full refutation of central Theorem 1.1 or of the quadratic corollary. The deeper variational counting proof was not fully verified in this audit. Subsequent Gasull–Santana 2025 and the 2026 symposium still treat the global question as open.

### Gaiko and earlier published <=4 attempts

**Valery A. Gaiko, arXiv:math/0611142 (2006)**, *Geometry of Planar Quadratic Systems*, Theorem 3.2 claims at most four; Theorem 4.3 claims at most three around any focus by field rotation and Wintner–Perko termination. **CLAIMED GLOBAL UPPER BOUND; not treated as accepted by later field literature**. No dedicated primary rebuttal found. The proof's own alternatives allow termination at a separatrix cycle, while the ensuing contradiction invokes Bautin's bound for a focus; these are different objects. It also moves from four distinct cycles to a multiplicity-four swallow-tail surface without establishing that merger. These are visible unproved bridges in that argument; do not invoke Gaiko as a certified global obstruction. [Primary full text](https://arxiv.org/html/math/0611142), specifically Section 4, Theorem 4.3 and proof.

**Qin Yuanxun's published <=4 proof** was explicitly attacked by **Shi Songling, Bull. London Math. Soc. 20 (1988), 597–599**, *A Counterexample to a Proposed Solution of Hilbert's Sixteenth Problem for Quadratic Systems*. The author's abstract says two results used in Qin's proof are contradicted. This is not a fifth-cycle example. Full three-page coefficients were not retrieved, so the precise counterexample system remains a retrieval item. Shi's 1990 retrospective, reference [19], recovers the original bibliography: Y. Chin, *On surfaces defined by ordinary differential equations*, LNM 1151 (1985), 115–131. [Primary rebuttal record](https://londmathsoc.onlinelibrary.wiley.com/doi/abs/10.1112/blms/20.6.597).

## Other relevant 2025–26 developments and false positives

1. **Haibo Lu arXiv:2607.13785v3 (2026-08-26), 99 pages**, claims local uniform finite cyclicity of the H^3_14 semihyperbolic hemicycle in the full 12-dimensional quadratic coefficient neighborhood, via an analytic slice and finite stopped transition maps. Bound existential, nonsharp; no global H(2) bound claimed. **PREPRINT, not independently verified here.** [Current record](https://arxiv.org/abs/2607.13785).
2. **Lu arXiv:2607.27464v4 (2026-08-10)** is now titled *Neutral Entry–Exit Cycles with Quadratic Grazing...* and concerns **continuous piecewise-smooth slow-fast** systems. “Quadratic grazing” does not make it a smooth quadratic polynomial field. It includes CAPD and dual-backend ancillary certificates, useful methodology but no H(2) advance. [Current record](https://arxiv.org/abs/2607.27464).
3. **Eshkobilov–Kadyrov–Mamayusupov arXiv:2604.12883v1 (2026-04-14)** gives replication H(nm+m-1)>=m²H(n). It raises the degree, hence cannot add cycles at degree two. New lower estimates listed are at degrees 14,29,31,39. **PREPRINT.** [Record](https://arxiv.org/abs/2604.12883).
4. **Gasull–Santana arXiv:2510.11705**, *Limit cycles and invariant algebraic curves*, has journal DOI **10.3934/cpaa.2026012**. Relevant algebraic-curve restrictions, not a proof of uniform quadratic finiteness. [Record](https://arxiv.org/abs/2510.11705).
5. **Artés–Cairó–Llibre 2025**, *Phase Portraits of Family IV*, QTDS 24:66, classifies x'=y, y'=d+ax+by+l x²+mxy+n y²; Theorem 1.1 reports at least 85 portraits, not every quadratic family. [Primary publisher](https://link.springer.com/article/10.1007/s12346-025-01227-9). It links a **2026-08-07** follow-up for Families V–X. **Artés–Chen–Ferrer–Jia, Dynamical Systems 40(2) (2025), 191–222**, treats the Chinese Class I, with 261 parameter subsets, 49 portraits and 7 containing cycles. “Class I” depends on which classification and must not be conflated with the ten-class x' normal forms. [Institutional source](https://portalrecerca.uab.cat/en/publications/quadratic-vector-fields-in-class-i/).
6. **Zhao arXiv:1011.2253** says **at most five** cycles from Q4 period annulus, and exhibits **three**, not five. It is a local upper bound. [Record](https://arxiv.org/abs/1011.2253).
7. **Yang–Liang–Zhang 2010**, *Five Limit Cycles for a Class of Quadratic Systems in R3*, is **three-dimensional** with cycles on invariant cones, outside planar H(2). [Journal primary record](https://camath.fudan.edu.cn/camaen/ch/reader/view_abstract.aspx?file_no=31A412&flag=1).
8. **Tian–Yu 2014**, *Seven Limit Cycles Around a Focus Point in a 3-D Quadratic System*, is also **3D**, not H(2). [Author's published PDF](https://publish.uwo.ca/~pyu/pub/preprints/TY_IJBC2014.pdf).
9. **Da Cruz–Torregrosa et al.** modern six/twelve-plus piecewise quadratic results are a different object class. A smooth quadratic Bautin upper bound cannot be applied to their discontinuity, and their lower counts cannot be imported into H(2). [2024 piecewise paper](https://doi.org/10.1016/j.nonrwa.2024.104124), [2025 piecewise Kolmogorov primary manuscript](https://repositorio.usp.br/directbitstream/52e2f003-979e-460c-b384-1ba7b3ba374d/3270082.pdf).

## Search coverage and remaining uncertainty

Searches included title/author searches and citation leads for both user-started arXiv records; H(2)=4/5/6/7/global finite bounds; smooth quadratic five/six/seven cycles with discontinuous/cubic/R3 exclusions; 2025–26 quadratic classification and hemicycle work; Voronin/Lebedev and Malyarets antecedents; Gaiko and Qin criticism; and small exploratory Chinese/Russian phrases for five/six quadratic limit cycles. Primary arXiv records and current version histories were opened, not inferred from dated snippets. A secondary 2026 problem registry was used only as a lead list; it missed the six-cycle journal claim and conflated some recency details, and is not status evidence here.

This is **not an exhaustive bibliography of all quadratic literature**. No MathSciNet full reviews or comprehensive zbMATH citation audit, no exhaustive Chinese/Russian journal scan, and no private correspondence were available. Hernandez Rosales's unavailable full manuscript is a concrete open retrieval item; exact systems in the Qin rebuttal and the 2024 economics conference contribution are additional retrieval gaps. No author was contacted. No independent published mathematical rebuttal of the 2026 six-cycle article was located; its defects above are transparent algebra checked here. Nothing found warrants changing the accepted lower bound from four.


## Final historical source checks

Shi's own [1990 retrospective, pp.202–203](https://www.labmath.uqam.ca/~annales/volumes/14-2/PDF/193-206.pdf) records his corrected December 1978 five-cycle claim and preserves the attempted family and inequalities. The source scan says **V3**, not V5. It also distinguishes the Petrovskii–Landis 1967 correction letter from earlier abandonment. The exact family and correction are integrated in [the claims ledger](../HISTORICAL_FIVE_CYCLE_CLAIMS.md).

Yeung's *Dulac's Theorem Revisited* was [published in QTDS 24:57 on January 27, 2025](https://doi.org/10.1007/s12346-025-01220-2). It challenges a closure assertion in Ilyashenko's 1991 proof route; it does not construct a quadratic counterexample. Bamón's 1986 individual-quadratic finiteness result is independent. Scope and remaining verification limits are recorded in the same claims ledger.
