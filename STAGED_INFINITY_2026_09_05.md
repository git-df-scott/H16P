# KKL infinity audit for the staged attack

Exact theory work, September 5, 2026. This note uses the existing beta-zero KKL field

\[
P=(1+x)y+x^2,\qquad Q=-10x^2+\tfrac{11}{5}xy+cy^2-mx,
\qquad m=\frac{K+42}{11c/5-1}.
\]

The inherited shape interval is \(1/2\le c\le3/2\), with the requested Stage 2 range \(0<K<6/5\). No ODE integrations are performed here. `staged_2026_09_05/infinity_check.py` checks the algebra independently with exact symbolic arithmetic.

**Main findings.** At the reported cutoff near \(c=0.9301\), the only projective infinity singularity is the vertical saddle. Its antipodal exponent product is one, but its required interior connection is impossible: the two separatrices lie on opposite sides of a line that can only be crossed in the wrong direction. Thus the proposed vertical infinity graphic cannot be the terminal mechanism there. At larger \(c\), a different pair of infinity saddles has a genuine *candidate eigenvalue-neutrality line* \(c=c_*\), independent of \(K\), where

\[
J(c_*)=0,\quad J(c)=305+634c-11c^2-1000c^3,
\quad c_*=0.968620633553494\ldots.
\]

A global separatrix connection and its first transition coefficient on this line remain uncomputed. These results neither exclude KKL folds nor certify a neutral graphic.

## 1. Exact charts and time orientation

In the signed vertical chart \(x=u/v,y=1/v\), with \(dt=v\,d\tau\),

\[
u_\tau=(1-c)u-\tfrac65u^2+10u^3+v(1+mu^2),\qquad
v_\tau=v(10u^2-\tfrac{11}{5}u-c)+muv^2.
\]

The north antipode has tangential/radial eigenvalues \((1-c,-c)\). At the south antipode the physically forward compactification reverses these signs, giving \((c-1,c)\). Using one signed time chart at both antipodes without reversing the southern time gives an incorrect itinerary.

In the signed finite-slope chart \(x=1/v,y=z/v\), with \(dt=v\,d\tau\),

\[
z_\tau=p(z)-v(m+z^2),\quad v_\tau=-v(1+z)-zv^2,
\quad p(z)=-10+\tfrac65z+(c-1)z^2.
\]

A simple root has radial/tangential eigenvalues \((-(1+z),p'(z))\) on the positive-x antipode, with both reversed at the other antipode. The discriminant is \((1000c-964)/25\). Lower-degree parameter \(m\), and hence \(K\), changes global transitions but none of these eigenvalues.

Here a saddle's Dulac exponent is \(\rho=|\lambda_{stable}|/\lambda_{unstable}\). An elementary graphic's exponent is the product of its saddle exponents.

## 2. The vertical graphic is impossible below 241/250

For \(c<241/250\), there are no real finite-slope infinity singularities. The vertical north saddle is stable transversely and unstable along the equator; the south saddle is unstable transversely and stable along the equator. Both equator semicircles run north to south. Closing either semicircle therefore requires an interior orbit from south to north.

The transverse invariant manifold has the exact local expansion

\[
u=-v-\frac{v^2}{1+c}+O(v^3),\qquad
x=-1-\frac{1}{(1+c)y}+O(y^{-2}).
\]

To verify it, substitute \(u=-v+av^2\) into \(u_\tau-u_vv_\tau\); its coefficient of \(v^2\) is \(a(1+c)+1\). Consequently the north stable separatrix lies in \(x<-1\), while the south unstable separatrix lies in \(x>-1\). But

\[
P(-1,y)=1
\]

everywhere. A forward orbit cannot cross from \(x>-1\) to \(x<-1\). The required interior connection is impossible for every \(K\) in the stated range (indeed every finite \(m\)). With no other infinity directions, there is no directed graphic through infinity assembled from these saddles and regular connecting orbits. The existing finite-equilibrium analysis gives two antisaddles, so no additional finite saddle can supply a missing connection in this parameter region.

The reciprocal exponents are nevertheless

\[
\rho_N=\frac{c}{1-c},\qquad \rho_S=\frac{1-c}{c},\qquad \rho_N\rho_S=1.
\]

This is a useful distinction: automatic local eigenvalue neutrality does not create a graphic. In particular, radius \(2^{20}\) near \(c=0.9301\) is an artificial numerical boundary, and is not evidence for the specific graphic proposed in Stage 2. The conclusion does not by itself identify the actual continuation or prove an a priori uniform radius bound.

## 3. The admissible candidate neutrality line above 241/250

For \(241/250<c<1\), let \(z_s\) be the smaller positive root of \(p\), with \(q_s=p'(z_s)>0\). This is a saddle; the larger root is a node. Equator arcs that avoid the node connect the positive-x finite saddle to the south vertical saddle, and the north vertical saddle to the negative-x finite saddle. A potential interior connection closes each pair within its corresponding half-plane. Unlike the vertical-to-vertical connection, these connections are not ruled out by \(x=-1\).

Their exponent products, for origin-side and remote-side itineraries respectively, are

\[
R(c)=\frac{(1-c)(1+z_s)}{cq_s},\qquad R(c)^{-1}.
\]

Neutrality requires \((1-c)(1+z_s)-cp'(z_s)=0\). Exact elimination gives

\[
\operatorname{Res}_z\left(p,(1-c)(1+z)-cp'\right)
=\frac{(c-1)J(c)}{25}.
\]

There is precisely one root of \(J\) on \([1/2,3/2]\), because \(J'<0\) there, and
\(J(241/250)>0>J(39/40)\). At that root,

\[
z_s=\frac{11c_*-5}{5(1+c_*-2c_*^2)}>0,
\quad p'(z_s)=\frac{(1-c_*)(1+z_s)}{c_*}>0.
\]

Thus no extraneous node root enters the elimination. The candidate eigenvalue-neutrality set in the requested parameter rectangle is exactly the vertical line \(c=c_*\) for these hyperbolic two-saddle itineraries. The existing `kkl/notes_local_unfolding.md` identifies \(J=K=0\) as a double-center organizer; that center fact does not establish a graphic at positive \(K\).

For \(1<c\le3/2\), the two finite-slope roots \(z_+>0>z_-\) are saddles and the vertical direction is a node. The node-free equator arcs pair the two finite saddles on the positive-x side, or their antipodes on the negative-x side. Since \(p'(z_-)=-p'(z_+)\), the corresponding exponent products are

\[
R_+(c)=-\frac{1+z_+}{1+z_-},\qquad R_+(c)^{-1}.
\]

Here \(z_-<-1\). The equation \(R_+=1\) requires \(z_++z_-=-2\), and Vieta's formula gives \(c=8/5\), outside the inherited box. Therefore these hyperbolic graphic itineraries are not neutral in \(1<c\le3/2\). This excludes their local neutrality, not all cycle folds. At \(c=241/250\) and \(c=1\), a nonhyperbolic infinity point requires separate analysis; hyperbolic ratio arguments do not apply.

## 4. What remains of the first graphic coefficient

Suppose an actual two-saddle graphic is established, with fixed positive transverse section coordinates. Write its Dulac maps and regular transition maps as

\[
d_i(s)=a_i s^{\rho_i}(1+o(1)),\qquad
g_i(s)=b_i s(1+o(1)),\quad a_i,b_i>0.
\]

For the itinerary \(g_2\circ d_2\circ g_1\circ d_1\),

\[
\Pi(s)=Cs^{\rho_1\rho_2}(1+o(1)),\qquad
C=b_2a_2(b_1a_1)^{\rho_2}.
\]

On the neutral line \(\rho_1\rho_2=1\), the first leading displacement coefficient is \(C-1\). If it vanishes, higher Dulac terms (which can include logarithms) must be examined. The constants include the global interior transition and cannot be recovered from the eigenvalue product. An unestablished connection also has a splitting parameter; the displayed return asymptotic assumes that splitting is zero. This audit supplies the expression, not a computed value for an unverified graphic.

The proposed Stage 2 kill test therefore remains unmet: the neutral eigenvalue line is present, the positive-K connection/splitting is uncomputed, and finite-amplitude folds have not been excluded. Lack of a sampled neutral graphic is not a theorem of absence over a parameter region.

## 5. Multiplier and physical-period diagnostic correction

A cycle approaching a hyperbolic graphic need not have multiplier tending to one. When the exponent product differs from one, the asymptotic derivative can tend to zero or infinity; at product one it can tend to the nonunit transition coefficient \(C\). Thus a multiplier bounded away from one does not exclude graphic absorption. Even genuine multiplier convergence to one is not a proof of a fold or of any particular boundary geometry.

Logarithmic residence times are standard in desingularized saddle time. The original physical time here obeys \(dt=v\,d\tau\) in the north charts and the corresponding positive radial factor in southern charts. This factor tends to zero at infinity. A logarithmically divergent compactified residence time can therefore have bounded physical duration; for example \(v\sim e^{-c\tau}\) gives a convergent physical-time integral. Period logging must distinguish these two times and cannot alone classify the terminal mechanism.

The useful next numerical task is continuation with the correct full-return section and compactified coordinates, explicitly logging the itinerary and connection splitting near \(c_*\). The exact exclusion near 0.9301 and the exact candidate line at \(c_*\) should guide that work.
