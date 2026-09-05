# Supplemental negative-K and actual-m review

The saved negative sheet remains consistent with a numerical fold pair.
Its inner cycle is unstable and its outer cycle stable; both finite
antisaddles are attracting at every audited field. The pole in the
(c,K) parametrization at c=5/11 is not a singularity of the actual field.
There is an analytic large-m reason to examine c=1/3, but no endpoint
or complete-component theorem has been established.

`theory_negative_checks.py` performs the exact identities and saved-data
checks below with **zero ODE calls**. Its JSON fixes file snapshots by
hash: 20 accepted fields plus one unresolved attempt in `events_negative`,
and six accepted fields plus one unresolved attempt in `events_m`.
It also checks eleven accepted fields plus one unresolved attempt in the
completed `events_logm` snapshot.
Every accepted field has two ordered opposite-sign root brackets, negative
G_z, and the recorded unstable/stable ordering. The attempted c=0.3
corrector failure supplies no endpoint evidence. The latest audited log-m
field has r approximately 1.19703159003e10, m approximately 2.76644661952e12,
and c approximately 0.33668600865535249. Later file updates are outside the
recorded hash snapshot. None of these saved profiles triggers three origin
cycles; the first requested (3,1) configuration has not been established.

Independent complete angular returns in `large_m_full_reproduction.json`
also reproduce both bracket endpoint signs at the same actual-m field:
c approximately 0.33668739232461736 and m approximately 2.76644661952e12.
These are complete returns, beyond the earlier two-half matching checks.

| Bracket radii, approximately | Left log displacement | Right log displacement | Endpoint derivative signs |
| --- | ---: | ---: | --- |
| 1.031660e10 to 1.053643e10 | -2.29867e-5 | +1.20739e-6 | Both positive |
| 1.355000e10 to 1.411362e10 | +5.01111e-7 | -6.41798e-5 | Both negative |

The first bracket is consistent with an unstable cycle and the second
with a stable cycle. The independent flux discrepancies are approximately
1e-24 or smaller. The audit checks exact equality of the saved parameter
strings across the four endpoints, their radius ordering, and recorded
signs. This is numerical reproduction of a two-cycle field, not interval
certification or a third-cycle result. Subsequent root polishing is not
needed to assert the recorded bracket signs and is outside this review's
scope. No A–D completion condition is established; the outcome remains
**PARTIAL**.

## A finite-topology gate covering all new fields

Use actual m>0 in

    xdot=(1+x)y+x²,
    ydot=-10x²+(11/5)xy+c y²-mx.

For 1/3<=c<=1 and m>=37, put u=1+x and

    W(u)=(61/5-c)u³+(m-72/5+3c)u²+(11/5-3c)u+c.

All audited fold and pair fields satisfy this box. Its quadratic
coefficient is at least 118/5. Thus for u>=0,

    W(u)>=(61/5-c)u³+(118/5)u²-(4/5)u+1/3>0.

The last quadratic has negative discriminant. For u=-v<0 write
W(-v)=-h v³+b v²-d v+c, with h>0, b>=118/5 and d=11/5-3c.
If d<=0, the coefficient signs give one positive root. If d>0,
d<=6/5, and at any positive critical point the cubic value is

    c+(b v²-2d v)/3 >= c-d²/(3b)>0.

Hence again the only positive root occurs on the final decreasing
branch. At this root W'(u)>0. As in the positive-sheet review,
determinant=x W'(u)/u>0. There are exactly two finite equilibria, the
origin and a single remote antisaddle x<-1; no finite equilibrium
collision or saddle is present. The transverse line x=-1 separates
any periodic orbit from the opposite nest.

The origin has determinant m>0 and trace zero. Its first Lyapunov sign
is the sign of K=m(11c-5)/5-42, which is negative at every audited field.
For the remote point x=-(v+1), the exact trace is

    (v+1)[(1+2c)/v-(16-10c)/5].

The unique remote v increases with m. Its Hopf trace threshold, expressed
without the singular (c,K) chart, is

    m_H(c)=21(1000c²+1021c+481)/[50(1+2c)²(8-5c)].

Exact rational comparisons at all stored folds and pair fields give
m>m_H(c), so the remote trace is negative. At c=5/11 the threshold is
1738/75, whereas the continued fold has m approximately 381.86. The
actual field and its derivatives are regular there. Every m on that
parameter line maps to K=-42; the inverse formula m=5(K+42)/(11c-5)
loses the parameter rather than detecting a dynamical boundary.

The sheet does cross an infinity-direction collision earlier, at
c=241/250: the discriminant 40c-964/25 of the finite projective-slope
quadratic vanishes there. The fold remains finite while crossing this
parameter line. For c below it, the remaining vertical infinity point
is hyperbolic for 0<c<1: in coordinates X=x/y,V=1/y its linearization
has eigenvalues 1-c and -c. In particular c=1/3 is not a finite-m loss
of hyperbolicity of that point. A putative endpoint there must involve
the singular m=infinity limit, not just the finite-field infinity chart.

## Large-m center scaling

The suggested x=mX,y=mY,time=mt limit is correct and produces the
nilpotent linearization [[0,0],[-1,0]] at the origin. A second scaling
exhibits an integrable leading system more directly. Put

    y=sqrt(m)Y, tau=sqrt(m)t, epsilon=m^(-1/2).

Then

    x'=(1+x)Y+epsilon x²,
    Y'=-x+cY²+(11/5)epsilon xY-10epsilon²x².

At epsilon=0 and 0<c<1/2 the origin has a center annulus in x>-1.
With u=1+x, its first integral, up to an additive constant, is

    H=u^(-2c)[Y²+a u+b], a=2/(1-2c), b=1/c.

The integrating factor is 2u^(-2c-1). Both end potentials tend to
infinity, so the annulus has ovals at arbitrarily large energy E.
The first Melnikov integral for the epsilon perturbation is

    M(E)=(4/5) integral[u^(-2c-2)(u-1)((16-10c)u+5+10c)
                         *sqrt(E u^(2c)-a u-b) du],

between the two positive turning roots. This formula and the first
integral are derived and checked algebraically in the accompanying script.

## Why one third is distinguished

The positive right-end contribution and negative left-end contribution
have leading forms C_R E^p_R and -C_L E^p_L, where

    p_R=(3-4c)/(2-4c), p_L=(2c+1)/(2c),

    C_R=(4/5)(16-10c) a^(1/2-p_R)
         *B((1-c)/(1-2c),3/2)/(1-2c),

    C_L=4(1+2c)b^(1/2-p_L)
         *B((c+1)/(2c),3/2)/(2c).

Here B is the beta integral. These endpoint asymptotics follow by setting
u=(E/a)^(1/(1-2c))v on the right and u=(b/E)^(1/(2c))v on the left.
The rescaled integrals converge for c in a neighborhood of 1/3. The
remaining compact region is of lower order. Their exponent difference is

    p_R-p_L=(3c-1)/[2c(1-2c)].

At c=1/3 both exponents equal 5/2, and the constants reduce exactly to

    C_R=152/675, C_L=8/27, C_L/C_R=25/19.

Thus the first Melnikov integral is negative at sufficiently large E
when c=1/3. For fixed c>1/3 close to it, the right term eventually
dominates positively; for fixed c<1/3 the left term dominates negatively.
This change makes one third a concrete organizing value for large
amplitude first-order balance. It does not establish uniqueness of a
Melnikov zero or imply that all actual folds have this limit.

If a joint c-down-to-1/3, E-to-infinity sequence satisfies cancellation
of these two leading Melnikov terms, then

    (c-1/3)log E -> (2/27)log(25/19).

For the positive horizontal radius r, E~a r^(1-2c), so equivalently

    (c-1/3)log r -> (2/9)log(25/19)
                  =0.060985965711502286487...

Applying this relation to actual cycles needs a **uniform joint
perturbation estimate**: the true displacement remainder must be small
relative to epsilon times the sum of the two displayed leading terms.
Compact-annulus Melnikov theory alone does not supply this estimate when
the oval escapes. Applying it to the double-cycle condition additionally
needs differentiated estimates and control of higher-order terms.
The observed r/m decreasing is a useful possible scaling condition,
but it has not been proved sufficient. It is therefore a conditional
asymptotic target for continuation, not a certified endpoint law.

The current c=0.35 point has m approximately 188826.674 and r approximately
3523.033; its (c-1/3)log r is about 0.136, still far from the conditional
limit above. At the latest audited log-m point that product is about
0.07780, closer to the conditional 0.06099 target but still not a proof
of convergence. Its size alone does not justify replacing the finite fold
by the limiting center or nilpotent graphic. No third root, exact K1
exclusion, or completed component boundary is established here.
