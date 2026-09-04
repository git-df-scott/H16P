# Q4 reconstruction geometry: an exact Green lift and rigorous exclusions

2026-09-04. **PROVED:** five original Q4 zeros require both `κ>21636/19043` and a first universal primitive crossing after `5/11`. The complete certified weighted-lobe box fails the latter condition and is therefore excluded for every `κ>1`. Whether another lobe-region point can saturate the reconstruction allowance remains **UNKNOWN**.

Derivation and bounded replay: [q4_reconstruction.py](q4/q4_reconstruction.py).
The inherited original basis is
`(hI00,I10,I01,2I−1,0+3κhI−1,1)`, with positive area orientation.
The inherited zero bounds already exclude five zeros when `β1=0`.
Every possible five-zero direction can therefore be rescaled to `β1=1`,
which is the normalization used throughout. Rescaling changes no roots.
The universal lobe assumptions used below are necessary for five original
zeros by [Q4_ZERO_GEOMETRY.md](Q4_ZERO_GEOMETRY.md).

## Exact scalar reconstruction, including a newly resolved source error

Put

\[
k=\kappa,\ d=k-1,\quad a=d/k\in(0,1),\quad s=k-dt,
\quad C=J_1(k)=\pi/\sqrt d,
\]
\[
q(t)=A+Bt-1+(t-\eta)M(t),\quad
H(t)=\int_0^t uF(u)q(u)\,du,
\quad Y(t)=G(s(t))/C.
\]

**PROVED.** The correct scalar initial-value problem is

\[
\boxed{(1-at)(1-t)Y''-\frac{1-a}{2}Y'+\frac{5a}{36}Y
=-\frac{H(t)}{1152t^2(1-t)}}. \tag{R1}
\]

Its exact center data are

\[
\boxed{Y_0=\frac{3(1326A+864B-2431\eta-102)}{1361360}},
\qquad
\boxed{Y_1=-\frac32(1+a)Y_0-\frac\eta{192}}. \tag{R2}
\]

The quotient `H/t²` has the removable value `(A−1−η/6)/2` at zero.
Consequently the IVP starts directly at the center and needs no unstable
shooting from a nearby singular initial point.

The final original integral is reconstructed by one positive-weight
primitive and a negative nonvanishing multiplier:

\[
X(t)=\int_0^t\frac{Y(u)}{(1-au)^{3/2}}\,du,
\qquad
\boxed{I(s(t))=-\frac{aC}{2}\sqrt{1-at}\,X(t)}. \tag{R3}
\]

Thus the only remaining lift parameter is the bounded number `0<a<1`.
These formulas reconstruct the **original** four-integral normalization;
they do not define a surrogate auxiliary integral.

### Derivation and the sign error

The source operator and forcing in `h` are

\[
L_h=5kh-(9kh^2-8)D_h+h(9kh^2-4)D_h^2,
\qquad L_hG=\frac{2h\mathcal F}{(9h^2-4)^2(9kh^2-4)}.
\]

Since `s=9kh²/4`, direct differentiation gives

\[
L_h=-\frac{16s}{h}L_s,\qquad
L_s=s(1-s)D_s^2-\frac12D_s-\frac5{36}.
\]

Therefore

\[
\boxed{L_sG=-\frac{k\mathcal F}{1152(s-k)^2(s-1)}}. \tag{R4}
\]

The **minus sign** is required. The plus sign printed in [Zhao equation (24)](https://arxiv.org/html/1011.2253)
is inconsistent with the displayed `h` equations. It was harmless for the
inherited zero-count inequalities, but reverses the actual reconstruction
forcing. Substitution of `ℱ=−d²C H` and `D_s=−D_t/d` in (R4) gives (R1).
The independent numerical comparison below exposed and confirms this sign
correction; the proof is the exact variable conversion just displayed.

In the original coefficients the same center data are

\[
Y_0=\frac49\mu_1-\frac23\mu_2-\frac23\mu_3+\frac{4d}{3}\mu_4,
\]
\[
Y_1=\frac{36-31k}{81k}\mu_1+\frac1{54}\mu_2
+\frac{k-6}{54k}\mu_3-\frac{31d}{27}\mu_4.
\]

They follow by differentiating `G=hI_h−I` and the center PF series.
Transport through the audited invertible alpha/beta map, with
`α1=−A−2kB/d`, `α2=B/d`, `β0=k−dη`, `β1=1`, gives (R2).
An exact small symbolic solve verified these expressions. Finally
`G=−2sI_t/d−I`, together with `I(0)=0`, integrates to (R3).

## The homogeneous Green function is elementary and positive

Set `x=asinh(sqrt(s−1))`. The homogeneous `L_s` equation becomes

\[
y_{xx}-2\tanh x\,y_x+\frac59y=0.
\]

Two explicit solutions are

\[
E(x)=\frac{5\cosh(x/3)-\cosh(5x/3)}4,
\qquad
O(x)=\frac3{10}\sinh(5x/3)+\frac32\sinh(x/3).
\]

Substitution proves the formulas directly. Their `x` Wronskian is
`E O_x−E_x O=cosh²x`, and their `s` Wronskian is

\[
W_s(E,O)=\frac{\sqrt s}{2\sqrt{s-1}}>0.
\]

Write `E(t)=E(x(s(t)))` and `O(t)=O(x(s(t)))`. Then
`W_t=−d sqrt(s)/(2sqrt(s−1))<0`. The forward Green kernel for the
standard second-order equation is

\[
\boxed{K_a(t,u)=
\frac{O(t)E(u)-E(t)O(u)}{W_t(u)}>0},\qquad0\le u<t<1. \tag{R5}
\]

Indeed `O>0` for `t<1`, and `(E/O)_t=−W_t/O²>0`; the numerator and
denominator in (R5) are both negative. The kernel has
`K(u,u)=0`, `∂_tK(u,u)=1`.

If `Y_hom` is the homogeneous solution with the data (R2), then

\[
Y(t)=Y_{\rm hom}(t)
-\int_0^t K_a(t,u)
\frac{H(u)}{1152u^2(1-u)^2(1-au)}\,du. \tag{R6}
\]

This specifies the sign action exactly. The inhomogeneous response has a
negative kernel acting on `H`. The homogeneous term is fixed by the
original coefficient transport and cannot be freely chosen to saturate a
generic two-zero allowance.

### Exact limiting kernels

On compact subintervals of `0≤u<t<1`, the endpoints of the lift parameter
have finite limits. At `a=0`, the homogeneous solutions are `1` and
`sqrt(1−t)`, and

\[
K_0(t,u)=2\sqrt{1-u}\bigl(\sqrt{1-u}-\sqrt{1-t}\bigr)>0.
\]

At `a=1`, the homogeneous solutions are `(1−t)^(1/6)` and
`(1−t)^(5/6)`, and

\[
K_1(t,u)=\frac32\bigl[(1-u)^{5/6}(1-t)^{1/6}
-(1-u)^{1/6}(1-t)^{5/6}\bigr]>0.
\]

The corresponding forcings in the standard equations are respectively
`−H/[1152t²(1−t)²]` and `−H/[1152t²(1−t)³]`. These are compact-domain
parameter limits, not assertions of uniformity at the homoclinic endpoint.

## A rigorous additional exclusion for the remaining lift

Suppose the normalized coefficient point belongs to the strict lobe region:
`q` has three simple zeros with initial sign positive, and `H` has three
simple zeros and signs `+,−,+,−` on its four lobes.

First, **Y0 is strictly negative for every such point**. In fact three `q`
zeros imply `1<η<54/31`, `A+B<η`, and `B>−1`. The last inequality follows
at a zero of `q'` from
`B=−M+(η−t)M'>−M>−1`. Put
`cA=9/3080`, `cB=162/85085`, `cEta=3/560`, `c0=−9/40040`.
Since `cA<cEta` and `cB<cA`,

\[
Y_0=c_AA+c_BB-c_\eta\eta+c_0
<(c_A-c_\eta)+(c_A-c_B)+c_0<0.
\]

It follows from (R3) that the original integral starts positive at the center.

Normalize the positive homogeneous solution that vanishes at the loop:

\[
y(t)=\frac{O(x(k-dt))}{O(x(k))}>0,\quad y(0)=1,
\]
\[
r(a)=y'(0)
=-\frac{\sqrt a}{2}\frac{O_x(x(k))}{O(x(k))},
\quad O_x(x)=\frac{\cosh(5x/3)+\cosh(x/3)}2.
\]

An equivalent entirely algebraic parameterization is
`z=sinh²(x(k)/3)>0`, for which

\[
k=(1+z)(1+4z)^2,\qquad
r=-\frac56+\frac5{3(8z^2+10z+5)}.
\]

The identity follows from the triple-angle formulas. Rational `z` therefore
gives rational `k` and rational `r`, useful for exact lift certificates.

Define the explicit scalar

\[
\boxed{P_0(a,A,B,\eta)=Y_1-r(a)Y_0}. \tag{R7}
\]

**PROVED sign-chain exclusion.** If a point lies in the strict lobe region
and `P0≤0`, its original integral has at most three distinct simple interior
zeros. In particular, five simple zeros require the additional strict
inequality `P0>0`.

For a direct proof put `Z=Y/y`,
`p(t)=sqrt((1−t)/(1−at))`, and `P=p y² Z'`. The standard homogeneous
equation gives

\[
P'(t)=-\frac{p(t)y(t)H(t)}{1152t^2(1-t)^2(1-at)}.
\]

Here `P(0)=P0`, `Z(0)=Y0<0`. If `P0≤0`, then `P` starts nonpositive
and decreases strictly on the first `H` lobe. Its subsequent monotone
pieces have signs of derivative `+,−,+`; it can therefore change sign
at most three times, with initial sign negative. Since `Z` also starts
negative and `Z'` has the sign of `P`, the same monotone-lobe argument
gives at most three sign changes of `Z`, hence of `Y`, counting each
odd-order zero as one sign change regardless of its multiplicity.
Finally `X(0)=0`, `X'=Y/(1-at)^(3/2)` and its first lobe is negative,
so the anchored primitive `X`, and hence `I`, has at most three sign changes
and therefore at most three simple zeros. Touches at intermediate extrema do not create additional
sign crossings. This statement suffices to exclude the five-simple-zero
target; it makes no unsupported claim that every point of the lobe region
satisfies (R7).

Thus the exact next necessary condition is more restrictive than membership
in the universal lobe region. It measures whether the original fixed center
data even permit the first additional oscillation of the inverse operator.

### A uniform center-parameter interval is now excluded

**PROVED.** No nonzero original Q4 integral has five distinct interior zeros
when

\[
\boxed{1<\kappa\le\frac{21636}{19043}}.
\]

This is a new small-parameter exclusion, unrelated to the withdrawn upper
bound on `κ` from the old sign-reversed beta strip. It establishes at most
four **distinct** interior zeros in this range; no stronger global
multiplicity bound is claimed.

To prove it, every hypothetical five-zero point must lie in the strict
universal lobe region. There `A>1+η/6`, `B>−1`, and `η<54/31` imply

\[
1326A+864B-2431\eta-102>360-2210\eta,
\qquad
0<\frac{|Y_0|}{\eta}<\frac{601}{136136}.
\]

The homogeneous slope also satisfies `r(a)<−1/2`. Indeed the displayed
formula for `O` gives

\[
O(x)=\int_0^{\sinh x}\cosh\!\left(\frac23\operatorname{arsinh}v\right)dv.
\]

Its integrand is strictly increasing, so `O(x)/sinh(x)` is strictly
increasing. Thus `O_x/O>coth(x)`, and
`r=−tanh(x)O_x/(2O)<−1/2`. Combining these estimates,

\[
P_0<\eta\left[\left(1+\frac32a\right)\frac{601}{136136}
-\frac1{192}\right]\le0
\quad\text{if}\quad a\le\frac{2593}{21636}.
\]

This is exactly `κ≤21636/19043`. The sign-chain exclusion then contradicts
the supposed five simple zeros. If five distinct zeros existed, the inherited
multiplicity-five bound would make them simple, so they too are excluded.

## Exact shooting conditions for saturation of the inverse operator

The following reduction has been independently checked against (R1)–(R7), including its signs and endpoint limits. It isolates the next
construction problem without claiming its inequalities are satisfiable.

Let `0<tau1<tau2<tau3<1` be the three simple zeros of `H`. Define

\[
\Omega(t)=\frac{y(t)}{1152t^2(1-at)^{3/2}(1-t)^{3/2}}>0,
\quad V(t)=\int_0^t\Omega(u)H(u)\,du.
\]

Then `P=P0−V` and `P'=−ΩH`. Because `H(1)<0` and
`y(t)` is a positive constant times `sqrt(1−t)` at the loop,
`P(t)→+∞`. Therefore `P` has four distinct simple zeros precisely when

\[
\boxed{\max\{0,V(\tau_2)\}<P_0<\min\{V(\tau_1),V(\tau_3)\}}. \tag{S1}
\]

Indeed `P` decreases, increases, decreases, then increases on the four
successive `H` lobes. Formula (S1) is exactly the alternating-extremum
condition, including its initial sign. It is necessary for five `Y` zeros
and hence for five original `I` zeros. It can be extremely restrictive even
when `P0>0`; a positive initial slope alone is insufficient.

If (S1) holds, denote the four zeros of `P` by `p1<p2<p3<p4`.
Since `Z(0)=Y0<0` and `Z'=P/(p y²)`, the endpoint behavior just used also
implies `Z(t)→+∞`. Thus `Y=yZ` has five distinct simple zeros precisely if

\[
\boxed{Z(p_1)>0,\quad Z(p_2)<0,\quad Z(p_3)>0,\quad Z(p_4)<0}. \tag{S2}
\]

This explicitly requires the first positive lobe of `P` to overcome the
negative offset `Y0`; making `P0` arbitrarily small can defeat this condition.

Finally, if the resulting five zeros of `Y` are `v1<⋯<v5`, then its signs
are `−,+,−,+,−,+`. The final primitive `X` in (R3) has five distinct simple
zeros precisely when

\[
\boxed{X(v_2)>0,\ X(v_3)<0,\ X(v_4)>0,\ X(v_5)<0,\ X(1)>0}. \tag{S3}
\]

For every fixed `a<1`, the endpoint integral `X(1)` is finite because the
logarithmic divergence of `Y` is integrable and `1−at` stays positive.
Conditions (S1)–(S3) are sequential weighted-lobe criteria for actually
saturating the two-zero allowance. They remain exact conditions to solve,
not a construction. `green_coordinates` in the replay module evaluates
`Z,P` from the stable scalar IVP for a targeted diagnostic.

### A minimal reduced target for an ordinary original-integral fold

**PROVED conditional mechanism; no point satisfying it has been found.**
Assume (S1), (S2), and `X(1)>0`. Replace exactly one of the four strict
extremum inequalities in (S3), at `v_j` with `j∈{2,3,4,5}`, by

\[
X(v_j)=0,
\]

while retaining the other three strict signs. The monotone-lobe geometry
then gives exactly three simple original-integral roots and one ordinary
double root at `v_j`. Indeed `X'=Y/(1-at)^(3/2)` and `Y'(v_j)≠0`, so
`X''(v_j)=Y'(v_j)/(1-av_j)^(3/2)≠0`. Two neighboring simple crossings
have coalesced at this extremum; the other three remain on their strictly
monotone lobes. The nonvanishing multiplier in (R3) preserves these
multiplicities.

A transverse variation of this extremal height to the sign prescribed in
(S3) splits the double zero into two simple zeros. Such a transverse
coefficient direction is available: the original basis function
`hI00` never vanishes in the open annulus, so varying `μ1` changes the
value at the double root nontrivially. The other simple roots and strict
inequalities persist for a sufficiently small variation. Thus one scalar
extremal-height equality, after achieving (S1) and (S2), is an exact
reduced three-simple-plus-double target. This formulation does not assert
that the target is attainable.

## One exact lift exclusion and three bounded numerical lift tests

### Every live five-zero target must delay its first primitive crossing

**PROVED; independently reviewed.** For every `κ>1`, five
distinct original Q4 zeros require the first universal primitive zero to
satisfy

\[
\boxed{\tau_1>\frac5{11}}. \tag{E0}
\]

This is a necessary condition on the entire universal lobe region, not only
on the constructed candidate box. To prove it, suppose `tau1≤5/11` and that
five original zeros existed. The necessary four-crossing condition (S1)
gives `P0>0` and a first root `p1<tau1`. On `[0,p1]`, `P≤P0`, and the
positive-homogeneous bounds used below give

\[
\int_0^{p_1}\frac{dt}{p(t)y(t)^2}
<\int_0^{5/11}(1-t)^{-13/6}dt
=\frac67\left[\left(\frac{11}6\right)^{7/6}-1\right]<\frac89.
\]

The last strict inequality requires no floating-point exponent evaluation:
`(11/6)^(7/6)<55/27` is equivalent, after raising positive quantities to
the sixth power and cancelling, to
`11·3^11=1948617<2000000=2^7·5^6`.

Writing `m=|Y0|`, the universal estimates `P0<(9/4)m−η/192` and
`m/η<601/136136` therefore imply

\[
Z(p_1)<-m+\frac89P_0<m-\frac\eta{216}<0,
\]

since the exact rational gap is

\[
\frac1{216}-\frac{601}{136136}=\frac{395}{1837836}>0.
\]

This contradicts the first-peak condition (S2). Thus the first primitive
crossing must occur after `5/11`. No sampling or assertion about the later
primitive crossings is involved.

### The entire certified coefficient box is excluded for every kappa

**PROVED; independently audited.** Let `Bcert` be the closed box of radius
`10^−7` in the infinity norm about the frozen rational `(A,B,η)` below.
The [exact lobe certificate](q4/notes_certificate_second.md) places all of this
box in the strict lobe region and puts its first primitive zero below `3/8`.
Nevertheless **no point of this box produces five original Q4 zeros for
any `κ>1`**. This is a uniform hostile disproof of the constructed candidate
region, even though its three universal primitive zeros are certified.
It now follows immediately from (E0), since `3/8<5/11`. The independent
coefficient-margin proof below is retained as a separate exact check.

Here is a more general sufficient exclusion. Suppose a lobe-region point
has first `H` zero below `3/8` and

\[
\eta>307|Y_0|. \tag{E1}
\]

If five original zeros existed, the necessary shooting condition (S1) would
give `P0>0` and a first root `p1` of `P` before the first `H` zero. On
`[0,p1]`, `H>0` and hence `P≤P0`. The local logarithmic slope bounds
`−5/6<−(s−1)O_s/O<−1/2` integrate to

\[
(1-t)^{5/6}<y(t)<\sqrt{1-t},\qquad p(t)\ge\sqrt{1-t}.
\]

Thus, using `p1<3/8`,

\[
\int_0^{p_1}\frac{dt}{p(t)y(t)^2}
<\frac38\left(\frac85\right)^3=\frac{192}{125}.
\]

The fixed initial data also give the uniform estimate

\[
P_0=C_a|Y_0|-\frac\eta{192},\qquad
C_a=\frac32(1+a)+r(a)<\frac94. \tag{E2}
\]

For an exact algebraic proof of (E2), use
`z=sinh²(x(k)/3)>0`, `k=(1+z)(1+4z)²`, and `D=8z²+10z+5`.
Then

\[
\frac94-C_a=\frac{kD+18D-20k}{12kD},
\]
\[
kD+18D-20k=128z^5+352z^4+z(72z^2-118z+55)+75>0.
\]

The quadratic in the last expression is positive because its leading
coefficient is positive and its discriminant is `−1916`.
Consequently the first maximum of `Z` satisfies

\[
Z(p_1)<Y_0+\frac{192}{125}P_0
<\frac{307|Y_0|-\eta}{125}<0.
\]

This contradicts the necessary first-peak inequality in (S2). Therefore
five original zeros are impossible under (E1), for every `κ>1`.

The box verification is a particularly small exact calculation.
At its center,

\[
\eta+307Y_0=
\frac{14871489355525071}{272272000000000000}>0.
\]

The coefficient vector of this linear functional in `(A,B,η)` is
`(2763/3080, 49734/85085, −361/560)`. Subtracting `10^−7` times its
one-norm gives the valid bound throughout the closed box

\[
\boxed{\eta-307|Y_0|=\eta+307Y_0
\ge\frac{79526371464733}{1456000000000000}>0.}
\]

Here `Y0<0` throughout the lobe region was already proved. The primitive
sign certificate proves that its first zero is less than `3/8` throughout
the same box. This completes the uniform exclusion without any numerical
sampling of `κ` or the coefficient box.

The argument also bounds the number of original sign changes by three on
this box: when `P0≤0`, use (R7); when `P0>0`, `P` begins and ends positive
and has at most four sign changes. With at most two, `Z` has at most three;
with four, its first maximum is negative by the estimate above, again
allowing at most three. The center-anchored primitive `X` cannot increase
this number. No global multiplicity-three assertion is needed.

Use the independently certified rational universal point

\[
A=1243911778077/10^{12},\quad B=-86917392526/10^{12},
\quad\eta=1460428426173/10^{12}.
\]

**RIGOROUS COMPUTATION.** At

\[
k=\cosh^2(3\log(6/5))=\frac{3878922961}{2916000000},
\quad a=\frac{962922961}{3878922961},
\]

every hyperbolic value in (R7) is rational. Exact fraction arithmetic gives

\[
r=-\frac{27095705}{51954846},\qquad
P_0=-\frac{3056925605483331742344151782161}
{2151790655250172000064000000000000}<0.
\]

Together with the independently certified lobe membership, this proves that
this exact original Q4 lift has at most three simple interior zeros. It is a
failed five-zero construction with an analytic reason, not a numerical miss.

**NUMERICAL.** The script evaluated the same universal point at exactly three
chosen lift parameters, `a=1/4,1/2,3/4`, corresponding to `k=4/3,2,4`.
At the five specified points `t=1/8,3/8,5/8,7/8,0.99`, all original-integral
values were positive. These finite sign evaluations are diagnostics, not a
proof of absence of additional roots. The initial values were:

| a | Y0 | P0 |
|---|---:|---:|
| 1/4 | −0.004579180628428025 | −0.0014094080328338168 |
| 1/2 | −0.004579180628428025 | +0.00016316222665472592 |
| 3/4 | −0.004579180628428025 | +0.0016304797545134398 |

The condition (R7) alone leaves the latter two lifts unresolved. The stronger
first-peak theorem above now excludes five zeros for both of them and for
every other lift of this box; their positive sample values alone did not
establish that conclusion.

At `k=2`, the transported original coefficients were numerically

```
mu = (0.004338948182084654,
      0.001126145265353277,
      0.002593635392529454,
     -0.0030208112030745363).
```

An independent 40-digit original area evaluation at `s=1.5` gave
`0.0015985882612802535`; the scalar PF reconstruction gave
`0.0015985882612802537`. The absolute discrepancy was `2.17e−19`.
With the source's erroneous forcing sign it instead differed by `4.54e−9`.
That failure is preserved here as the reason the sign regression matters.

The entire diagnostic, including this independent area comparison and the
exact rational exclusion, is capped at ten CPU seconds with numerical
libraries restricted to one thread. The numerical part took about one second
in the observed run. No coefficient scan, optimization, or dense sweep was
performed.

### Frozen shooting reaches four P crossings and fails the next gate

**NUMERICAL ONLY.** A targeted one-dimensional shooting of the remaining
lift parameter selected the exact rational value

\[
a=\frac{189417314263391}{400000000000000},\qquad
k=\frac{400000000000000}{210582685736609}\approx1.899491397409134,
\]

using the same certified universal coefficient point. The frozen replay is
[q4_green_shoot.py](q4/q4_green_shoot.py), with results in
[second_green_shoot.json](q4/data/second_green_shoot.json).
It uses four predetermined scalar brackets and performs no new parameter
search. The reported CPU time was approximately `0.009` seconds.

Its initial weighted derivative was `P0≈2.1646953843×10^−8`. The diagnostic
found four crossings of `P`, meeting the numerical version of (S1), but all
four extrema of `Z` remained negative:

| Approximate root of P | Z at that root |
|---:|---:|
| 0.1444648725550713 | −0.004579179219335157 |
| 0.3821270921597670 | −0.004579180186171887 |
| 0.5838595137337097 | −0.004579178869461455 |
| 0.8100626745341928 | −0.004579200508688371 |

This is a failed original five-zero construction: the first positive
derivative lobe does not overcome the negative center offset `Y0`, so (S2)
fails. All five sampled original-integral values were positive as well.
The numerical crossings are not certified and do not show saturation of
the original reconstruction allowance. The separate analytic first-crossing
barrier already proves failure of five original zeros for this point and
its entire certified box, uniformly over all `k`.

## What remains open after the reconstruction analysis

The Green kernel has an exact sign, the original center data are explicit,
and the necessary extra inequality `P0>0` is now available. Nevertheless,
no case here realizes the generic `+2` allowance, and no proof rules that
out throughout the rest of the lobe region. The entire certified coefficient
box has now been rigorously excluded for all `κ`. Positivity of the Green kernel alone is insufficient
because the fixed homogeneous component can affect variation.

The next construction domain is precisely the certified universal lobe
region intersected with `P0>0`, with the one lift parameter
`2593/21636<a<1`, and with first primitive crossing strictly after `5/11`.
An additional weighted-lobe test on `P`, then `Z`, then the final primitive
`X`, or a rigorously controlled three-simple-plus-double configuration of
`I`, is required to force five original zeros.
