# Independent scope review of the KKL finite-fold candidate

September 5, 2026. Read-only review of the numerical artifacts plus exact rational algebra; **zero additional ODE evaluations**.

## Outcome

The new result is meaningful numerical evidence for an ordinary finite-amplitude fold of origin cycles at fixed \(K=1/512\), near

\[
c_f=0.9688884793906646,\qquad r_f=6.949087993605231.
\]

It is not an interval-certified fold, an M1 three-cycle point, or a five-cycle field. The perturbed field has numerical evidence for a stable/unstable pair. This is a materially stronger result than merely driving an existing cycle's multiplier toward one.

## Return and unfolding evidence

`fold_refined_candidate.json` records a full-return numerical root with \(D\approx-9.6\times10^{-14}\), \(D_r\approx6.7\times10^{-16}\), winding approximately \(-1\), period approximately 0.9152566646, and x-range approximately \([-0.552524,6.949088]\).

The independently integrated Cartesian replay in `final_verification.json` gives physical return displacement approximately \(-1.40\times10^{-10}\), logarithmic displacement approximately \(-2.02\times10^{-11}\), and return derivative approximately 0.999999999977. The difference between replays is a numerical error indicator; the smaller original residual is not a rigorous error bound.

For \(c=c_f+10^{-6}\), the independent Cartesian return signs are

| Section radius | Return displacement |
|---|---:|
| 4 | \(+6.81228\times10^{-6}\) |
| 6.949087993605231 | \(-1.13527\times10^{-6}\) |
| 12 | \(+6.21330\times10^{-5}\) |

Provided the same return is continuously defined across those intervals and the signs are correct, the intermediate value theorem gives a zero in each of the two intervals. The sign changes support an attracting inner cycle and repelling outer cycle; simplicity and exact existence are not yet certified. The middle sign becomes positive at \(c_f-10^{-6}\), consistent with an ordinary fold unfolding. A single positive middle sample does not prove absence of cycles on that side.

The positive return displacement at \(r=20000\) is useful outer-profile evidence. It neither supplies a third cycle nor excludes one beyond that radius or between unsampled radii.

## Why the Newton Jacobian is relevant

The solver uses \(z=\log r\), \(L(z,c)=\log(P(r,c)/r)\), and the residual vector \((L,L_z)/s\), with \(s=(1-q)^2>0\), \(q=1/(1+r)\). At an exact fold, \(L=L_z=0\), so dividing by \(s\) does not change the zero set or Jacobian rank. Also

\[
L_c=D_c/r,\qquad L_{zz}=rD_{rr}
\]

at that zero. Therefore a nonsingular two-dimensional Jacobian with nonzero \(L_c\) and \(L_{zz}\) is the correct ordinary-fold diagnostic, not an artifact of the log coordinate or removal of small-radius vanishing.

The last recorded finite-difference Jacobian is approximately

\[
\begin{pmatrix}
2.41\times10^{-9}&-0.2137458\\
2.89142\times10^{-5}&-0.2048610
\end{pmatrix},
\qquad \det\approx6.18\times10^{-6}.
\]

Converted to unnormalized derivatives, it suggests \(D_c\approx-1.13513\) and \(D_{rr}\approx3.18\times10^{-6}>0\). This agrees with a positive local minimum becoming negative as c increases. These are finite differences near, rather than rigorous derivative enclosures at, the final root. The recorded Jacobian supports ordinary nondegeneracy numerically; it does not certify it.

## Exact rational remote-focus gate

Interpret the decimal **pair field**, \(c_f+10^{-6}=0.9688894793906646\), as the exact rational

\[
c=\frac{4844447396953323}{5000000000000000},\quad
K=\frac1{512},\quad
m=\frac{1050048828125000000}{28288921366486553},\quad \alpha=-m,\quad\beta=0.
\]

Using the established exact threshold

\[
K_H(c)=-\frac{441J(c)}{125(16-10c)(1+2c)^2},
\]

direct rational evaluation gives \(J(c)<0\) and \(K<K_H(c)\), with \(K_H\approx0.0383581568420\). The value approximately 0.0382154 applies to the unperturbed fold field. Thus the remote trace is **positive**, not negative, at the pair field.

The focus property can also be checked without orbit integration. Set

\[
F(s)=10s+\frac{11s^2}{5(s-1)}-\frac{cs^3}{(s-1)^2},\quad
s_- = \frac{3324674}{10^6},\quad s_+=\frac{3324675}{10^6}.
\]

Exact rational comparisons give \(F(s_-)<m<F(s_+)\). The established strict monotonicity of \(F\) isolates the remote equilibrium \(x_*=-s_*\) inside this bracket. Its trace

\[
t(s)=s\left[\frac{(1+2c)s}{s-1}-\frac{21}{5}\right]
\]

is decreasing and obeys \(0<t(s_+)<t(s_-)<6/1000\). A direct rational lower bound is

\[
F_s(s)\ge10+\frac{11s_-(s_--2)}{5(s_+-1)^2}
-\frac{cs_+^2(s_+-3)}{(s_--1)^3},
\]

which yields \(\det=s_*(s_*-1)F_s(s_*)>88\). Hence \(t^2-4\det<0.006^2-352<0\). These exact comparisons prove an **unstable remote focus** for the rational pair field.

This fails the inherited stable-remote-focus/unstable-remote-cycle precursor gate. It does not exclude a stable remote cycle or any other finite-amplitude remote-cycle mechanism. No remote cycle at the pair field has been demonstrated here.

## Exact multiplier-density gate

With \(u=1+x\), \(d=16-10c\),

\[
W=m+(2m+10)x+(m+111/5)x^2+(61/5-c)x^3,
\]

the inherited exact multiplier density polynomial is

\[
N=\{du+(c+1)(21+dx)\}W-u(21+dx)W'.
\]

At the rational pair field, exact evaluation gives \(N(0)=5/512>0\) and \(N(1)<-180<0\). Both x-values lie inside the recorded candidate's x-range. Thus the necessary positive/negative density sampling for a multiplier-one cycle is allowed. This does not prove the weighted integral vanishes; the return computation supplies the numerical evidence for that balance.

## Relation to infinity continuation and the requested target

The independently replayed remote cycle at \(c=0.9683,K=1/64\) has section coordinate about \(-6.5372\times10^9\), relative return residual about \(3.55\times10^{-9}\), physical period about 0.56721, and numerical multiplier about 1.010887. This directly demonstrates that the old \(2^{20}\) section cap was not a disappearance certificate. Its large absolute return residual, about 23 at a coordinate of order \(10^9\), should be reported together with the relative residual. It supplies no completed infinity connection or graphic coefficient.

The finite-fold candidate and the large remote-cycle continuation occur at different parameter vectors. They cannot be combined into a coexisting cycle count. The new pair supplies two origin cycles numerically; M1 still requires three origin cycles at trace zero at one coefficient vector. Stage 4's coexistence and Stage 5's five-annulus certification remain unavailable from these data.
