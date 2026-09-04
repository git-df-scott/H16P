# Shi/Chen–Wang/Bautin and topology audit — 2026-09-04

This is an audit contribution, not a new global upper-bound theorem. The algebra below was independently recomputed with exact rational symbolic arithmetic in `shi_algebra.py`; numerical eigenvalues are labelled as such. This appendix is incorporated into the frontier audit.

## Decisive conclusions

* Shi Songling, Song-ling Shi, S. Shi, and the “Songling system” refer to the same researcher/example. They are not independent four-cycle mechanisms.
* Three local cycles, not four, are the maximum available from an elementary quadratic weak focus. Shi produces three at one focus plus a distinct cycle at the other.
* Chen–Wang uses an order-two weak focus with a preexisting surrounding cycle, adds two local cycles, and retains a cycle at the other focus. This is a different construction mechanism, but again the final distribution is (3,1).
* Four around one focus is not among the accepted examples located. It must not be advertised as known. A four-cycle distribution (2,2) is ruled out by the published distribution theorem.
* The exact two-focus Shi weak-focus stratum algebraically excludes every additional real finite equilibrium. There is no finite saddle to seed a loop.
* More strongly, in the generic FOUR-REAL-equilibrium family, Zegeling’s Theorem 5.4 permits only (1,1) or (n,0). Thus making a finite saddle+node while retaining a hyperbolic (3,1) configuration cannot work even after leaving the weak-focus stratum.
* The audited theorems do not rule out an order-ONE weak focus with three preexisting nonzero cycles at that same focus and one at the other. This is an exact four-cycle precursor which, if found with simple cycles, would yield five by an ordinary Hopf bifurcation. It is not known to exist.

## Primary source ledger and proof status

1. **Shi Songling**, “A concrete example of the existence of four limit cycles for plane quadratic systems,” Scientia Sinica 23 (1980), 153–158; Chinese original 1979. DOI [10.1360/YA1980-23-2-153](https://doi.org/10.1360/YA1980-23-2-153). Accepted analytical example. Exact family and original small parameters reproduced in Yu–Zeng, below.
2. **Chen Lansun and Wang Mingshu**, “The relative position and number of limit cycles of the quadratic differential system,” Acta Math. Sinica 22 (1979), 751–758. [Original journal record](https://actamath.cjoe.ac.cn/Jwk_sxxb_cn/EN/10.12386/A1979sxxb0069). Accepted analytical construction. Original unspecified sufficiently small parameters must not be confused with subsequent visual numerical values.
3. **Pei Yu and Yanni Zeng**, “Visualization of Four Limit Cycles in Near-Integrable Quadratic Polynomial Systems,” [arXiv:2002.09987](https://arxiv.org/html/2002.09987), 2020. Primary numerical reconstruction; exact equations (10)–(15), original-vs-visual parameters, and independent explanation of the two historical mechanisms. Its plotted trajectories are numerical evidence, not interval certificates. The displayed radial equation (14) has evident notational slips (missing radial factor, `v4` in place of `v3`); do not copy it literally as a correct normal form.
4. **Chengzhi Li**, “Non-existence of limit cycle around a weak focus of order three for any quadratic system,” Chinese Ann. Math. B 7 (1986), 174–190. [Original journal abstract](https://camath.fudan.edu.cn/cambcn/ch/reader/view_abstract.aspx?file_no=7B205&flag=1). Analytical nonexistence theorem for a preexisting cycle around a third-order weak focus. This does not say every nearby strong focus has at most three global cycles.
5. **Jaume Llibre and Dana Schlomiuk**, “The Geometry of Quadratic Differential Systems with a Weak Focus of Third Order,” Canadian J. Math. 56 (2004), 310–343. [Primary PDF](https://www.cambridge.org/core/services/aop-cambridge-core/content/view/FCA46BEF322C4B8A02C2C47270177DF5/S0008414X00031758a.pdf/the-geometry-of-quadratic-differential-systems-with-a-weak-focus-of-third-order.pdf). DOI 10.4153/CJM-2004-015-2. Exact normal form and classification. **Important:** Theorem 12 explicitly leaves uniqueness/shape of the nonalgebraic connection curve numerical; Proposition 17(b) and Theorem 18(c,d) explicitly say numerical evidence. They are not proved global neighborhood bounds.
6. **André Zegeling**, “Nests of limit cycles in quadratic systems,” Advances in Nonlinear Analysis 13 (2024), 20240012. Accepted April 10, 2024. [Primary PDF](https://d-nb.info/1332906729/34), [DOI](https://doi.org/10.1515/anona-2024-0012). Theorem 1.2: only (n,0) or (n,1). Theorem 5.4: four real finite singularities allow (1,1) or (n,0). The paper supplies details missing from Zhang’s 2002 proof. These are publication-level rigorous results; no retraction or contradiction located.
7. **Zhang Pingguang and Cai Suilin**, “Quadratic systems with a weak focus,” Bull. Austral. Math. Soc. 44 (1991), 511–526. [Primary PDF](https://www.cambridge.org/core/services/aop-cambridge-core/content/view/F0A3CF508A37A97F71E9D37F1837F974/S0004972700030008a.pdf/quadratic_systems_with_a_weak_focus.pdf). DOI 10.1017/S0004972700030008. Restricted distribution/remote-cycle results; no bound of two cycles around every first-order weak focus stated.
8. **Zhang Pingguang and Zhao Shenqi**, “On number of limit cycles for the quadratic systems with a weak focus and a strong focus,” 2001, [DOI10.1007/s11766-001-0018-y](https://doi.org/10.1007/s11766-001-0018-y). Restricts the strong-focus nest; order-two/order-three cases have at most (1,1)/(0,1) before unfolding. Does not rule out the first-order precursor discussed below.
9. **Yu–Han**, “Four limit cycles in quadratic near-integrable systems,” JAAC 1(2) (2011), 291–298 [journal](https://www.jaac-online.com/article/doi/10.11948/2011020); longer IJBC22(10) (2012), 1250254, DOI10.1142/S0218127412502549, [arXiv1002.1055](https://arxiv.org/abs/1002.1055). Accepted reversible-center/Melnikov mechanism for (3,1), distinct from simply counting Bautin focus quantities. Some citations mistype the arXiv identifier as 1002.1005.

## Exact historical seeds

### Shi

\[
\dot x=\lambda x-y-10x^2+(5+\delta)xy+y^2,\quad
\dot y=x+x^2+(-25+8\epsilon-9\delta)xy.
\]

Original analytical values: \(\lambda=-10^{-250},\epsilon=-10^{-52},\delta=-10^{-13}\).

Yu–Zeng visual seed: \(\lambda=-2\cdot10^{-8},\epsilon=-1/1000,\delta=-1/10\).

Cycle distribution (3,1). At origin: stable focus; cycles, inner to outer, unstable/stable/unstable. At (0,1): unstable focus, stable cycle. The original existence argument constructs trapping regions and invokes Poincaré–Bendixson. The visual seed is not by itself rigorously certified by that numerical paper.

### Chen–Wang

\[
\dot x=-\delta_2x-y-3x^2+(1-\delta_1)xy+y^2,\quad
\dot y=x+\frac29x^2-3xy.
\]

Original: sufficiently small positive \(\delta_1,\delta_2\), unspecified. Visual seed: \(\delta_1=1/100,\delta_2=1/50000\). At the unperturbed origin the first nonzero focus coefficient in Yu–Zeng’s convention is \(-77/972\), order two. Two cycles are preexisting: an intermediate one around origin and a very large stable one around (0,1). Two local cycles are then born. Final stability is the same as Shi. The enormous remote scale in the Chen–Wang visual example makes it a poorer first numerical seed than KKL/Galias–Tucker.

### Broader mechanisms

Kuznetsov–Leonov and related “normal-size” examples deform these two-focus mechanisms and improve scale separation. “Three large cycles with a first-order weak focus” in that literature counts two around this focus and one remotely, not three around the weak focus. [Leonov2010 DOI10.1134/S1064562410060220](https://doi.org/10.1134/S1064562410060220).

The reversible Yu–Han family is
\[
\dot x=y(1+a_1x)+\varepsilon a_{10}x,\qquad
\dot y=-x+x^2+a_4y^2+\varepsilon(b_{01}y+b_{11}xy).
\]
Its (3,1) cycles come from center-annulus/Melnikov analysis, with combinations of small and finite-amplitude cycles. Neither “two small plus two large” nor “four normal-size cycles” means distribution (2,2) or four around one.

## Independently reproduced focus algebra

Work in the normalized local chart
\[
\dot x=-y+lx^2+mxy+y^2,\qquad \dot y=x+ax^2+bxy.
\]
Find homogeneous polynomials \(F_j\) such that
\[
F=(x^2+y^2)/2+\sum_{j\ge3}F_j,\qquad
\dot F=V_1r^4+V_2r^6+V_3r^8+\cdots.
\]
At each even degree our gauge sets the coefficient of \(x^j\) in \(F_j\) to zero; linear coefficient equations were solved exactly. Focus coefficients off lower-order zero loci depend on convention, so compare on the appropriate strata.

\[
V_1={m(l+1)-a(b+2l)\over8}.
\]
On \(V_1=0\), with \(a\ne0\), set \(b=(l+1)m/a-2l\), and define
\[
G=(l+1)^2(1+b)-a^2(b+2l+1).
\]
Then the exact computation gives
\[
V_2={m(5a-m)G\over48a}.
\]
A third-order weak focus therefore has
\[
m=5a,\qquad b=3l+5,
\]
with
\[
V_3={25a\over64}(2a^2+l+2)
\{3(l+1)^2(l+2)-a^2(5l+6)\}\ne0.
\]
The excluded factors identify center/lower-degeneracy branches; they are not extra free cycles. For Shi (a,l)=(1,-10),
\[
V_1=-\epsilon,\qquad
V_2\big|_{\epsilon=0}={5\delta(\delta+5)(36\delta+95)\over12},\qquad
V_3\big|_{\epsilon=\delta=0}={35625\over8}.
\]
The linear trace is \(\lambda\), so the three independent small unfolding controls are trace, \(V_1\), \(V_2\). This is codimension three; the two shape coordinates remaining in the normalized five-dimensional family are a,l. Bautin supplies cyclicity three, never four.

An exact five-parameter chart extending the Shi family is
\[
\dot x=\lambda x-y+lx^2+(5a+\delta)xy+y^2,
\]
\[
\dot y=x+ax^2+\left[3l+5+{(l+1)\delta+8\epsilon\over a}\right]xy,
\qquad a\ne0.
\]
It is useful for continuation, but does not override the geometric obstructions below.

## Finite-equilibrium obstruction, proved directly

On the exact third-order stratum, let b=3l+5. The equilibria are O=(0,0), N=(0,1), and any roots of the restriction to \(1+ax+by=0\). For b≠0 they satisfy
\[
A_2x^2+A_1x+A_0=0,
\]
\[
A_2=l(3l+5)^2-a^2(15l+24),\quad
A_1=-6a(2l+3),\quad A_0=3(l+2),
\]
with \(y=-(1+ax)/b\). The exact discriminant is
\[
A_1^2-4A_2A_0=12(3l+5)^2[3a^2-l(l+2)].
\]
At N the trace and determinant are
\[
T_N=5a,\qquad D_N=-3(l+2).
\]
For N to be a focus,
\[
25a^2+12(l+2)<0.
\]
Consequently \(l<-2-25a^2/12\), and
\[
3a^2-l(l+2)<-\frac76a^2-\frac{625}{144}a^4<0.
\]
Hence **having the second focus forces the extra pair to be nonreal**. This is an exact algebraic obstruction on the relevant Shi chart, with no numerical premise.

At the Shi base, O has eigenvalues ±i and N has \(5/2\pm i\sqrt{71}/2\). The extra-root discriminant bracket is -77. The extra pair cannot become real by an arbitrarily small perturbation. For illustration, a=1,l=-3 has a finite saddle-node at (1/3,1/3), but N is then a node, not the required second focus. At a=1,l=-5/2 there is a finite saddle and another node; the original two-focus geometry is gone.

Beyond the exact stratum, Zegeling’s four-real-equilibrium theorem forbids retaining (3,1) after creation of a real saddle+node pair. Persistence also excludes a hyperbolic (3,1) system exactly on a generic saddle-node boundary that can be perturbed to four real equilibria. This is not a prohibition on one-nest (5,0) families with saddles.

For arbitrary coefficients in this chart, useful continuation gates are
\[
A_2=lb^2-amb+a^2,\quad A_1=\lambda b^2+ab-mb+2a,\quad A_0=b+1.
\]
A generic finite saddle-node has \(A_1^2-4A_2A_0=0\); an equilibrium escaping to infinity has \(A_2=0\). At N the focus condition is \((\lambda+m)^2+4(1+b)<0\). At O it is \(\lambda^2<4\).

## Infinity, graphic compatibility, and proof limits

In the chart x=1/z,y=u/z, after polynomial desingularization, the infinity polynomial is
\[
C(u)=a+(b-l)u-mu^2-u^3.
\]
At a simple zero, eigenvalues of the desingularized local system are \(C'(u)\) tangentially and \(-[l+mu+u^2]\) transversely. On QW3,
\[
C(u)=a+(2l+5)u-5au^2-u^3,
\]
\[
\operatorname{disc}C=4[125a^4+25a^2l^2+170a^2l+262a^2+8l^3+60l^2+150l+125].
\]
At Shi, \(C(u)=1-15u-5u^2-u^3\), with discriminant -8752. Its only real projective direction is approximately u=0.065229852214. The eigenvalues are approximately -15.665063323 and 9.669595805: a hyperbolic saddle, with an antipodal copy on the Poincaré circle. The other two directions are nonreal. No saddle-node at infinity is infinitesimally close in coefficient space.

Llibre–Schlomiuk’s class QW3 has three relevant kinds of portraits:

| Base portraits | Base geometry | Published neighborhood/cyclicity status |
|---|---|---|
| W1–W12 | no cycle and no graphic | Theorem18(a): nearby quadratic systems ≤3 cycles |
| W13 | one infinity graphic around strong focus | Prop17(a)/Thm18(b): graphic cyclicity1, nearby ≤4 |
| W15,W18 | other infinity graphics around strong focus | Prop17(b)/Thm18(c): corresponding ≤4 claims explicitly numerical |
| W14,W16,W17 | one cycle around strong focus, no graphic | Thm18(d) conditional on hyperbolicity, explicitly numerical except specified subregion |

All of these graphics involve infinity and surround the strong focus. The weak-focus nest has no preexisting surrounding cycle (Li) and no surrounding graphic in this classification. A hypothetical two-cycle splitting of the remote graphic **while three local cycles survive** would have distribution (3,2), excluded by Zegeling. Thus the numerical gaps in the 2004 neighborhood results do not furnish a preserved Shi 3+2 construction.

This is a focused obstruction to the proposed construction, not a new proof of every numerical clause of Theorem18. In particular, do not label the entire quadratic neighborhood of every exceptional QW3 member rigorously ≤4 merely by citing that paper.

The exact finite-saddle variant of “Shi4+1” is closed; a broader five-cycle continuation must be reformulated away from a maximal-focus unfolding. Infinity can still matter to nonperturbative continuation, but must add to the larger nest to reach (4,1), not split the remote singleton into two.

## Configuration map

| Five-cycle distribution | Status from audited results |
|---|---|
| (5,0), often written (5) | open globally; necessarily one focus inside every cycle |
| (4,1) | open globally; natural continuation target from (3,1) |
| (3,2) | excluded by Zegeling2024 |
| (3,1,1),(2,2,1), etc. | excluded: at most two nests |
| a cycle surrounding both foci or any extra saddle | excluded for quadratic limit cycles: each surrounds exactly one focus |

Four-cycle examples located all have (3,1). “Four around one” is an open configuration, not a known building block. Four-real-equilibrium phase portraits can support a prospective (5,0), but cannot support a robust (4,1).

## First-order Hopf precursor: strongest surviving exact reformulation

The root audit proposes the KKL family
\[
\dot x=y+x^2+xy,\qquad
\dot y=-10x^2+bxy+cy^2+\alpha x+\beta y.
\]
At \(\beta=0,\alpha<0\), independently normalizing the linear rotation gives
\[
\ell_1={K\over8(-\alpha)^{3/2}},\qquad
K=-\alpha(bc-1)-10(b+2).
\]
A possible bounded target is
\[
b\in[3/2,3],\quad c\in[1/2,19/20],\quad \alpha\in[-200,-10],\quad\beta=0,
\]
retaining two real foci and a nonreal additional pair, with K bounded strictly away from zero. The task is to find **three NONZERO simple return-map roots around origin**, plus the remote simple cycle. A validated ordinary Hopf unfolding would then supply cycle five.

For the incumbent KKL cycle stability S/U/S around origin and remote U, the natural precursor sector is **K>0**: origin is repelling at beta0, allowing an innermost stable cycle. A small negative beta produces a new unstable cycle inside the three existing ones. Setting beta to zero at incumbent b,c,alpha with K<0 instead causes the innermost stable cycle to collapse and leaves only two nonzero roots; this is not a fifth-cycle construction.

Red-team findings:

* Bautin’s local cyclicity does not forbid this: the three nonzero roots are preexisting finite-amplitude cycles, while only one new cycle is created locally.
* Li’s order-three theorem and the order-two uniqueness theorem do not apply if K≠0.
* The resulting distribution is (4,1), permitted by the audited topology theorems. No theorem bounding every first-order weak focus to at most two surrounding cycles was located after targeted searches of Zhang–Cai1991, Zhang–Zhao2001, and later classification/citation chains.
* A fourth nonzero root should not be silently assumed from three nearby roots in a truncated focus polynomial. Numerical integration must locate actual full return-map roots and then validation must isolate each one away from zero.
* The remote cycle must be certified at the SAME coefficients; parameter freedom does not imply survival.
* No primary theorem imposing opposite stability on every pair of outermost nests was located. The proposed S/U/S plus remoteU sector already retains the incumbent opposite outer stabilities.
* This is a specific open precursor search with a sufficient certification gate, not an existing five-cycle candidate or a proof that the bounded box contains one.

## Modern classification leads retained for completeness

* Artés–Trullàs, IJBC36(6), 2026, “Quadratic Differential Systems with a Weak Focus of First-Order and a Finite Saddle-Node.” The [2023 thesis precursor](https://upcommons.upc.edu/bitstreams/c1f65358-cbbd-494b-903e-e64cb81e204f/download), §5.4, discusses numerically undetected crossings that could produce four around one focus; absence is conjectured. This supplies an independent one-nest frontier, not a ready five-cycle seed. The thesis is a prepublication source; the 2026 journal text needs a separate version comparison.
* [2026 EJDE codimension paper](https://ejde.math.txstate.edu/Volumes/2026/28/abstr.html) is relevant to organizing phase portraits modulo cycles. “Modulo limit cycles” explicitly leaves cycle counts unresolved and must never be read as H(2)=4.

### Final root-selected bounded slice and no-contact red team

For a smaller first experiment, fix \(b=11/5\) and use
\[
(c,\alpha)\in[1/2,3/2]\times[-200,-10],\qquad\beta=0,
\]
with exactly two finite foci, the nonreal extra pair, and
\[
K=-\alpha(11c/5-1)-42>0.
\]
There is no dimension-count objection to two shape parameters: once the Hopf equality is imposed, coexistence of three simple nonzero cycles and the remote simple cycle consists of open conditions. Nevertheless, fixing b can miss a viable component entirely. A negative numerical search closes neither the full KKL family nor the first-order Hopf route. Expand b only after useful evidence or a certified exclusion of the selected slice.

At beta0 the nonzero finite equilibria of KKL are the roots of
\[
(c-b-10)x^3+(\alpha-b-20)x^2+(2\alpha-10)x+\alpha=0,
\qquad y=-x^2/(1+x).
\]
The strict two-real-equilibrium gate is that this cubic has exactly one real root (negative), and its Jacobian has negative discriminant and positive determinant. Together with alpha<0 this supplies the two foci. A cubic discriminant <0 gives the first condition, except excluded degree-drop boundaries; the full Jacobian provides the second. Root x=-1 is impossible since substituting it gives -c≠0 in this selected c-box.

The line x=-1 is globally without contact because \(\dot x=1\) everywhere on it, independently of b,c,alpha,beta. Every periodic orbit remains entirely in one of its two open half-planes. The origin nest lies in x>-1 and the remote nest in x<-1 when its equilibrium has x<-1. This is a useful automatic return-map validation gate.

The homogeneous angular polynomial is
\[
XQ_2-YP_2=X[-10X^2+(b-1)XY+(c-1)Y^2].
\]
There is always the vertical projective direction X=0. Other infinite directions exist when
\[
(b-1)^2+40(c-1)\ge0.
\]
At b=11/5 their collision occurs at c=241/250=0.964. At the vertical direction, using v=x/y, w=1/y, the desingularized eigenvalues are 1-c tangentially and -c transversely: a saddle for c<1, a node for c>1, and degenerate at c=1. Thus this broadened c-box deliberately spans changes at infinity; retain this in the numerical ledger rather than pretending it preserves the incumbent complete portrait.

The classical quadratic antipodal-saddle connection restriction (Llibre–Schlomiuk2004 §7.1(v), citing its primary antecedents) says a finite separatrix connecting antipodal infinity saddles must be a straight invariant line. A connection between the two vertical infinity saddles in KKL would therefore require an invariant vertical line x=k. But
\[
P(k,y)=k^2+(1+k)y
\]
can never vanish identically in y: it requires both k=0 and k=-1. Hence **the obvious vertical antipodal infinity connection is impossible for the whole KKL family**. A numerical separatrix probe that never changes sign is consistent with this exact obstruction. This does not rule out all other graphics involving distinct infinity directions; it does rule out treating that particular probe as the cheapest fifth-cycle source.
