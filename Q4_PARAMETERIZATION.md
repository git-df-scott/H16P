# Q4 parameterization and exact search domain

## Center parameter

Parameterize the circle \(b^2+c^2=4\) rationally by

\[
 b=2\frac{1-\rho^2}{1+\rho^2},\qquad
 c=\frac{4\rho}{1+\rho^2},\qquad \rho>0.
\]

Then \(\kappa=4/(2+b)=1+\rho^2\). The sign \(c<0\) is equivalent to
\(c>0\) by conjugation/reflection with time reversal, so it supplies no second
search chart. Rational \(\rho\) gives rational original Q4 coefficients.

## Raw exact domain

Before analytic pruning the zero-search object is

\[
 \mathcal D=(1,\infty)_\kappa\times\mathbb{RP}^3_{[\mu]},
\]

four real dimensions. Use the canonical representative
\(\|\mu\|_\infty=1\) with the first nonzero entry positive. For root location use

\[
 r=\frac{s-1}{\kappa-1}\in(0,1).
\]

The compactification \(u=\rho/(1+\rho)\in(0,1)\) is useful for storage, but
its two faces are singular degenerations, not uniformly regular boxes.

## Zhao reduction and rigorous pruning

Zhao derives

\[
 g(s)=P_2(s)+Q_1(s)w(s),\quad
 P_2=\alpha_2s^2+\alpha_1s+\alpha_0,\quad
 Q_1=\beta_1s-\beta_0,\quad w=J_2/J_1,
\]

where

\[
 \alpha_0=\beta_0-\kappa\beta_1-\kappa\alpha_1-\kappa^2\alpha_2.
\]

The explicit linear map from \(\mu\) is in q4/q4_integrals.py. One must first
apply Zhao's relabeling after his equation (20):

\[
 \widetilde\mu_2=-\frac23\mu_2-
 \frac{2(\kappa-1)}{3\kappa}\mu_3,\qquad
 \widetilde\mu_3=-\frac{2}{3\kappa}\mu_3,
\]

with \(\widetilde\mu_1=\mu_1,\widetilde\mu_4=\mu_4\). Omitting this step
silently applies the pruning theorem to the wrong coefficients.

Zhao's Rolle/Picard--Fuchs chain is

\[
 \#I\le\#G\le\#\mathcal F+2,\qquad
 \#\mathcal F\le\#\mathcal F'=\#g.
\]

After normalizing \(\beta_1=1\), three zeros of \(g\), hence the only route
through this bound to five zeros of \(I\), require

\[
 \frac{23\kappa-54}{31}<\beta_0<1. \tag{P1}
\]

The strip is nonempty only when

\[
 1<\kappa<\frac{85}{23},\qquad
 0<\rho<\sqrt{\frac{62}{23}}\approx1.641.
\]

The cases \(\beta_1=0\), \(\beta_0\ge1\), and
\(\beta_0\le(23\kappa-54)/31\) cannot support five through Zhao's chain.
His closing comments give a further necessary survivor condition:

\[
 P_2(\beta_0)\ge
 \frac{25(1-\beta_0)}{432(\kappa-1)^2}. \tag{P2}
\]

Boxes violating (P1) or strictly violating (P2) are rejected analytically.

## Final bounded search object

For endpoint margin \(\delta>0\), search

\[
\begin{aligned}
 \mathcal D_\delta=\{&(\kappa,[\mu]):
 1+\delta\le\kappa\le85/23-\delta,\ \|\mu\|_\infty=1,\\
 &\text{first nonzero }\mu_i>0,\ \beta_1\ne0,\ \text{(P1), (P2)}\},
\end{aligned}
\]

and seek five simple roots in \(r\in[\delta,1-\delta]\). This is explicit and
compact but still four-dimensional; the integral constraints are
transcendental. Each endpoint needs a separate asymptotic chart rather than
indefinitely shrinking \(\delta\).

The first Astra core chart fixes exact margins

\[
 \frac{101}{100}\le\kappa\le\frac{369}{100},\qquad
 2^{-16}\le r\le1-2^{-16},
\]

together with (P1), (P2), and the projective convention above. The omitted
slivers \(\kappa\downarrow1\), \(\kappa\uparrow85/23\), \(r\downarrow0\), and
\(r\uparrow1\) are four named boundary charts, not permission to extend a
floating grid.

The more effective candidate coordinates are
\((\kappa,r_1,r_2,r_3)\): solve \(I(s_i)=0\) for the unique generic
projective direction \([\mu]\), then ask whether that function acquires two
additional roots. This removes most root-poor directions without pretending
the search is exhaustive.
