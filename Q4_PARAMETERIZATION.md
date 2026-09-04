# Q4 parameterization and corrected necessary conditions

Updated 2026-09-04. This supersedes the earlier sign-reversed beta strip
and inferred kappa cutoff. Historical numerical records are retained;
their old filter exclusions have no mathematical force.

## Exact family and coordinates

For \(\rho>0\), set
\[
b=2(1-\rho^2)/(1+\rho^2),\quad c=4\rho/(1+\rho^2),\quad
\kappa=1+\rho^2>1.
\]
Rational rho gives rational original Q4 coefficients. The other sign of c
is equivalent by reflection/time reversal. The raw domain remains
\((1,\infty)_\kappa\times\mathbb{RP}^3_{[\mu]}\).
Use the four-term basis of Q4_THEORY.md. Normalize
\(\|\mu\|_\infty=1\), first nonzero entry positive, for storage.
\[
h=-\frac23\sqrt{s/\kappa},\quad 1<s<\kappa,\quad
r=\frac{s-1}{\kappa-1},\quad t=1-r=\frac{\kappa-s}{\kappa-1}.
\]
The center is t=0 (s=kappa), the homoclinic boundary t=1 (s=1).

## Coefficient transport

Zhao uses \(g=P_2+(\beta_1s-\beta_0)w\), with
\(P_2=\alpha_2s^2+\alpha_1s+\alpha_0\), \(w=J_2/J_1\), and
\[
\alpha_0=\beta_0-\kappa\beta_1-\kappa\alpha_1-\kappa^2\alpha_2.
\]
The map in q4/q4_integrals.py correctly first applies
\[
\widetilde\mu_2=-\tfrac23\mu_2-
\tfrac{2(\kappa-1)}{3\kappa}\mu_3,\qquad
\widetilde\mu_3=-\tfrac{2}{3\kappa}\mu_3,
\quad\widetilde\mu_1=\mu_1,\quad\widetilde\mu_4=\mu_4.
\]
The source's chain, counting multiplicity, is
\(Z(I)\le Z(G)\le Z(\mathcal F)+2\le Z(g)+2\).

## Corrected necessary strip

If beta1 is nonzero, divide all alpha/beta coefficients by beta1.
The necessary condition for five zeros is
\[
\boxed{\frac{54-23\kappa}{31}<\beta_0<1.}\tag{P1-corrected}
\]
Indeed,
\[
g'''(\kappa)=3w''(\kappa)+(\kappa-\beta_0)w'''(\kappa)
=-\frac{25(23\kappa+31\beta_0-54)}{3888(\kappa-1)^3}.
\]
Zhao's Proposition 17 has the correct sign; Theorem 14 and Corollary 18
print its negative. The derivative identity and universal structure in
Q4_STRUCTURE.md resolve the inconsistency. The corrected strip has width
23(kappa-1)/31 for EVERY kappa>1. No kappa<85/23 cutoff follows.
This correction does not assert five-zero existence at any kappa.

Cases beta1=0, beta0>=1, and beta0<=(54-23kappa)/31 cannot give five
through this chain. The filter establishes no finite upper kappa bound.

## Corrected conservative curvature filter

Put beta=beta0<1 and P=P2(beta). Exact polynomial division gives
\[
\left(\frac{g(s)}{s-\beta}\right)''
=\frac{2P}{(s-\beta)^3}+w''(s).
\]
Since \(w'''>0\) and
\(w''(\kappa)=-25/[216(\kappa-1)^2]\), the condition
\[
P\le\frac{25(1-\beta)^3}{432(\kappa-1)^2}
\]
makes this second derivative strictly negative in the open interval.
Together with g(kappa)=0, strict concavity allows at most one interior
zero of g, hence at most three of I. A necessary survivor condition is
\[
\boxed{P_2(\beta_0)>\frac{25(1-\beta_0)^3}{432(\kappa-1)^2}.}
\tag{P2-safe}
\]
The earlier linear-power threshold is withdrawn: the source comments
omit the factor 2 and do not justify that threshold.
Floating-point rejection remains numerical near equality; exact or
interval evidence is needed for rigorous candidate-box exclusion.

## Universal coefficient chart

Let d=kappa-1 and, after beta1=1, define
\[
\eta=(\kappa-\beta_0)/d,\quad A=-(\alpha_1+2\kappa\alpha_2),
\quad B=d\alpha_2.
\]
Then, with the universal M of Q4_STRUCTURE.md,
\[
\frac{g(s)}{dt}=A+Bt-1+(t-\eta)M(t),\qquad1<\eta<54/31.
\]
The inverse chart is
\(\alpha_2=B/d,\ \alpha_1=-A-2\kappa B/d,\ \beta_0=\kappa-d\eta\),
followed by the alpha0 relation and inverse linear map to mu.
The subsequent lift to I still depends on kappa.

## Bounded work policy

No production sweep was launched. A finite numerical kappa window is a
chosen experiment, not an exhaustive analytic domain. The former
[1.01,3.69] window is historical and excludes none of its complement.
Three prescribed zeros generically determine projective mu; four require
a singular evaluation matrix. See Q4_ZERO_GEOMETRY.md.

Primary source: [Zhao, equations (14), (27)--(29), Lemmas 11, 15--16,
Proposition 17 and section 6](https://arxiv.org/html/1011.2253).
