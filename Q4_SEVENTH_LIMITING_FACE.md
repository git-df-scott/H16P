# Q4 route 4: exact exclusion at the limiting lift

2026-09-05. **Route 4 remains open at finite lifts. No four-zero
counterexample is supplied.** This note proves a three-zero bound for the
limiting reconstruction at `a=1`, and a two-zero bound on an explicit
three-dimensional coefficient subspace at every finite lift. Neither
statement establishes the requested global four-zero exclusion.

The proofs use the Stieltjes representation and reconstruction already
established in [Q4_STRUCTURE.md](Q4_STRUCTURE.md) and
[Q4_RECONSTRUCTION_GEOMETRY.md](Q4_RECONSTRUCTION_GEOMETRY.md). The new
identities have exact symbolic replays. They have not received an independent
mathematical audit.

## 1. Notation and reconstruction

Put `F=2F1(1/6,5/6;1;t)`, `G=F'`, and

\[
M=1-6(1-t)G/F,\qquad
k=(K_0,K_1,K_2,K_3),\qquad
k'=tF(1,t,M,tM).
\]

The primitive vector vanishes at the center. For a coefficient vector `h`,
write `H_h=k h`. Its uniquely prescribed reconstruction satisfies

\[
L_aY_h=-\frac{H_h}{1152t^2(1-t)},\quad
L_a=(1-at)(1-t)D^2-\frac{1-a}{2}D+\frac{5a}{36},
\]

\[
Y_h(0)=\ell(h)=\frac9{3080}
\left(h_0+\frac{144}{221}h_1+\frac{11}{6}h_2+
\frac{204}{221}h_3\right),\quad
Y_h'(0)=-\frac32(1+a)\ell(h)+\frac{h_2}{192}.
\]

The original integral has the same interior zeros as `X_h`, where
`X_h(0)=0` and `X_h'=Y_h/(1-at)^(3/2)`. The value `a=1` is the
`kappa -> infinity` comparison limit, not a finite Q4 center parameter.

## 2. Three closed reconstructions at every lift

The following table gives exact coefficient vectors in the order
`(K0,K1,K2,K3)`.

| Reconstructed function | Coefficient vector |
|---|---|
| `F` | `(-176(9a+4)/3, 1088a, 16(54a+59)/3, -768a)` |
| `tF` | `(-96(12a-5), 272(21a-10)/3, 192, -16(18a+49)/3)` |
| `t(1-t)G` | `(-40(36a-25)/9, 680(3a-2)/9, 80/3, 40(a-5)/3)` |

Applying `L_a` and the Gauss equation verifies their forcing. The transported
center values and slopes also match exactly, so no homogeneous correction
is missing. This is checked in
[finite_basis.py](q4/seventh/finite_basis.py).

These three functions form an extended complete Chebyshev system on `(0,1)`.
Their first three Wronskians are

\[
F>0,\qquad F^2>0,\qquad
-\frac{F^3}{6}\bigl(2M'+tM''\bigr)<0.
\]

Here `M'>0` and `M''>0` follow from the positive Stieltjes density. Thus
every nonzero reconstruction in this subspace has at most two interior
zeros, counting multiplicity. Rolle's theorem, including the anchored
zero `X_h(0)=0`, gives the same bound for its original integral.

This excludes an explicit codimension-one coefficient subspace at every
`0<a<1`. It does not exclude the fourth coefficient direction.

If `h_F,h_tF,h_G` denote the three columns above, then

\[
\det[h_F,h_{tF},h_G,e_0]
=\frac{2263040}{81}
 (486a^3-441a^2-486a+236).
\]

Consequently these three functions together with `J=Y_K0` describe the
entire reconstructed space wherever the displayed polynomial is nonzero.
In particular it is nonzero at `a=1`, where its value is `-205`.

## 3. Exact determinant signs at a=1

Let `j=(Y_K0,...,Y_K3)` at `a=1`. All four entries have closed forms:

\[
\begin{aligned}
j_0&=\frac9{1185800}
 [(410t+385)F+(11448t^2-11653t)G],\\
j_1&=\frac9{52412360}
 [(24436t+11088)F+(420239t^2-426528t)G],\\
j_2&=\frac3{215600}
 [(60t+385)F+(6408t^2-6438t)G],\\
j_3&=\frac3{6166160}
 [(7213t+5544)F+(176586t^2-177228t)G].
\end{aligned}
\]

The replay verifies all four differential equations and all eight center
data, using rational arithmetic and the exact Gauss relation.

For the endpoint row use `e=(9061,6289,2431,1819)`, a positive multiple
of `k(1)`. Direct determinant factorization gives

\[
\det[k;e;j;j']
=\frac{27}{1600}F^3t^2(M-1)(t-1)
 \bigl[t(M^2+4M)-3M-2\bigr]<0.
\]

Indeed `0<M<1`, and the bracket equals
`(M+2)(M-1)-(1-t)(M^2+4M)<0`; all other displayed factors are positive.

For the confluent row,

\[
\det[k;k';j;j']=-\frac{81F^4t^3}{2620618000}Q(t,M),
\]

where

\[
\begin{aligned}
Q={}&M^4(5t^3-59t^2+9t)
 +M^3(40t^3-208t^2+240t)\\
&+M^2(80t^3+68t^2+86t-180)
 +M(-176t^2-28t+60)+68t-5.
\end{aligned}
\]

Here is a sign proof on an explicit region containing the actual moment
curve, with no numerical sampling. Write

\[
M(t)=\int_0^1\frac{\rho(u)}{1-tu}\,du,\qquad
d\nu(u)=\frac{\rho(u)}{1-u}\,du.
\]

The endpoint identities give `nu([0,1])=1` and `E_nu[u]=5/6`.
The function `(1-u)/(1-tu)` is strictly concave in `u` for `0<t<1`.
Jensen's inequality and strict increase of `M` therefore yield

\[
\frac16<M(t)<\frac1{6-5t}.
\]

Parameterize the entire intervening strip by

\[
t=\frac z{1+z},\qquad
M=\frac{1/6+v(1+z)/(6+z)}{1+v},\qquad z>0,\ v>0.
\]

Exact substitution gives

\[
Q=\frac{125z\,P(z,v)}
 {1296(1+v)^4(1+z)^3(6+z)^4},
\]

where **all 29 nonzero coefficients of `P` are positive integers**.
The full polynomial and a replay of this identity are in
[a1_determinants.json](q4/seventh/a1_determinants.json) and
[a1_determinants.py](q4/seventh/a1_determinants.py). In particular `Q>0`
throughout the strip, proving the second strict determinant sign.

## 4. Three-zero theorem for the limiting original family

The determinant proof gives more than a boundary diagnostic. Write
`s=-1/[1152t^2(1-t)^3]`. At `a=1`, the reconstruction equation implies

\[
W(Y_{K0},Y_{K1},Y_{K2},Y_{K3})
=s^2\det[k;k';j;j']<0.
\]

To check this identity, subtract the homogeneous terms from the third
Wronskian row to leave `s k`, then from the fourth to leave `s k'`.
Moving these two rows ahead of `j,j'` has positive permutation sign.

Change coefficient basis to `(h_F,h_tF,h_G,e0)`. Its determinant is nonzero
at `a=1`. Section 2 supplies three nonvanishing lower Wronskians; the last
display supplies the fourth. Hence the four-dimensional reconstructed
space at `a=1` is an extended complete Chebyshev space in this order.
Every nonzero `Y_h` has at most three interior zeros counted with
multiplicity. Since `X_h(0)=0`, four distinct interior zeros of `X_h`
would require four distinct interior zeros of `Y_h`, a contradiction.

**Thus the limiting original family at `a=1` has at most three distinct
interior zeros.** The two strict determinant signs also close every
two-anchor fibre on that limiting face by the Strike-6 boundary theorem.

## 5. What is still missing

Neither compact-interval continuity in `a` nor the three-dimensional
subspace proves the same statement for all `a<1`. Continuity gives a
neighborhood of each fixed compact interval; it supplies no uniform
control as the anchor tends to `1` simultaneously with `a`.

An exact parameter connection for the remaining function `J=Y_K0` is
derived and checked in [parameter_connection.py](q4/seventh/parameter_connection.py).
It has the form `J_a=c(a)J-(1-t)J_t/[a(1-a)]+U(a,t,F,G)`.
It is a usable identity, not a sign theorem. A comparison based on dropping
the negative double-integral term in the determinant is insufficient:
the resulting upper bound becomes positive near the joint corner while
the actual determinant remains negative. The frozen diagnostic records
this failure in [explore_upper.json](q4/seventh/explore_upper.json).

The finite-lift boundary inequalities in
[Q4_SIXTH_BOUNDARY_REDUCTION.md](Q4_SIXTH_BOUNDARY_REDUCTION.md) are still
unproved. No original integral with four verified interior zeros was
found. There is consequently no justified binary verdict of either
global closure or counterexample from this work, and no proof of `H(2)=4`.
