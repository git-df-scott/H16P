# Independent review of the new KKL center and remote-Hopf gates

2026-09-04. This is an exact analytic review of `notes_local_unfolding.md`
and `notes_other_strata.md`, following the published 202-call checkpoint.
No ODE solve, numerical return evaluation, coefficient sweep, or cycle
discovery was performed. Earlier checkpoint notes and data were not edited.

**Verdict: pass.** The original-coordinate reversal, remote Hopf sign,
local return factorization, and rational geometry illustration are
consistent. They do not exclude a separate finite-amplitude pair or prove
that any new proposed field has a cycle.

## 1. Original-coordinate reversing involution

Write `m=-alpha`, `J=305+634c-11c^2-1000c^3`, and impose
`K=0`, hence `m=210/(11c-5)`. The exact values

    J(241/250)=25281/2500>0,
    J(39/40)=-11333/800<0

and `J'<0` on the stated c interval isolate a unique root c*. It is
beyond the first infinity-direction transition, below c=1, and outside
the experimental K>=1/64 margin.

Independently expanding the original field in
`zeta=x-i y/sqrt(m)`, with positive time `tau=sqrt(m)t`, gives the A,B,C
coefficients in the local-unfolding note. On K=0,

    A=(8/5) conjugate(B),
    Im(conjugate(B)^3 C)=J/(35280 sqrt(m)).

For example, put `p=1/sqrt(m)` and `sigma=5(1+2c)/21`. Then
`B=(p-i sigma)/2`. Direct multiplication of `conjugate(B)^3 C` reduces
its imaginary numerator to

    250(1+2c)^2(1-c)-(11c-5)(c+11)=J.

At J=0, `lambda=-B/conjugate(B)` has modulus one. The three quadratic
reversal identities are

    A lambda=-conjugate(A),
    B=-lambda conjugate(B),
    C=-lambda^3 conjugate(C).

The first two follow on K=0; the third is exactly the imaginary-part
identity above. Translating `zeta -> lambda conjugate(zeta)` back gives

    M = [[m sigma^2-1, -2 sigma],
         [-2m sigma, 1-m sigma^2]]/(1+m sigma^2).

Thus `F(Mz)=-M F(z)` for the original quadratic field and original time,
at the exact organizer. Also `M^2=Id`, `det M=-1`, and its fixed line
is `x+sigma y=0`. The coordinate normalization used a positive time
scale; the minus sign is the actual reversing symmetry, not an accidental
time-orientation change in a Hopf coefficient.

A compact independent polynomial replay writes M=N/L, where

    L=1000c^2+1231c+145,
    N11=1000c^2+769c+355,
    N12=-10(1+2c)(11c-5),
    N21=-2100(1+2c), N22=-N11.

L is positive on the relevant interval. The replay verifies
`N^2=L^2 Id`, `det N=-L^2`, linear anticommutation, and every
coefficient of `F2(Nz)+L N F2(z)` modulo J. This checks the full
original-coordinate identity without substituting a floating approximation
to the algebraic root.

## 2. Why there are two centers

The independent nullcline variable `s=-x>1` gives the remote equation

    m=F(s,c)=10s+(11/5)s^2/(s-1)-c s^3/(s-1)^2.

For `w=s-1>0`,

    w^3 F_s=(61/5-c)w^3+(3c-11/5)w+2c>0.

The piecewise lower bounds in both source notes are valid on the stated
c interval. The endpoint limits give a unique simple remote root; the
origin-side positivity argument excludes any additional real equilibrium.
The remote determinant is `s(s-1)F_s>0`.

The reversing involution fixes the origin and must fix the unique other
real equilibrium. Intersecting `x+sigma y=0` with the first nullcline gives

    x*=-1/(1-sigma), y*=1/[sigma(1-sigma)].

Here 0<sigma<1, so x*<-1. Its trace is zero. Both equilibria have
nonzero imaginary linear eigenvalues and are fixed by the reversing
reflection. Nearby rotating half-arcs join their reversed reflections to
make periodic orbits, proving that both are centers. Calling the organizer
an order-three weak focus would be wrong: its return map is identically
the identity locally at the exact center.

## 3. Independent remote-Hopf sign and its normalization

In this section set `d=16-10c`, `e=1+2c`, both positive in the box.
The remote trace vanishes at `s_H=21/d`, `w_H=5e/d`. Substitution in
the preceding F gives the stated m_H, and direct simplification gives

    K_H=-441J/(125d e^2).

An independent useful form is `K_H=-s_H^3 J/(105 w_H^2)`.
Since `trace_s=-(d/5+e/w^2)<0` and F_s>0, increasing m, and hence K
at fixed c, decreases the trace. Therefore negative remote trace is
equivalent to K>K_H in the specified simple-remote regime. A negative
focus discriminant remains a separate gate.

For an original-coordinate coefficient check write

    xddot=A(x) xdot^2+B(x) xdot-C(x),
    A=(c+1)/(1+x), C=xW/(1+x).

At the remote equilibrium, with primes here denoting x derivatives,

    D=C'=s w F_s>0,
    C''=-2(2s-1)F_s-s w F_ss,
    B'=d/5+e/w^2, B''=2e/w^3.

Use the increasing local coordinate `dz/dx=rho=(-1-x)^(-c-1)>0`.
For the resulting Lienard equation the Hopf numerator is

    f_z g_zz-f_zz g_z
      =[B''C'-B'C''+2AB'C']/rho^2.

At Hopf its original-coordinate numerator reduces to

    B's[2(1-c)F_s+wF_ss]=882J/(625e^3),

using the independently checked polynomial identity

    5e^2(1-c)(61-5c)-c(c+1)d^2=J.

In `(U,V)=(z-z*,-zdot/sqrt(D))` and positive time `sqrt(D)t`,

    l1_z=441J/[5000e^3 rho^2 D^(3/2)].

If the leading amplitude coordinate is instead the original x-x*, the
positive factor rho^2 cancels and
`l1_x=441J/[5000e^3 D^(3/2)]`. Both formulas have the sign of J.
For J<0 the positive-K remote Hopf is supercritical: its small cycle is
stable, on the unstable-focus side K<K_H. It does not supply the target
remote U cycle on the stable-focus side. The J=0 case is the separate
center calculation, not an ordinary Hopf with a vanishing coefficient
silently ignored.

## 4. Local unfolding and stable-remote cone

The analytic-division argument in (U6) passes. `(K,J)` are valid local
coordinates, the common small return domain exists by uniform rotation,
and the exact center makes the K=0 displacement divisible by J. Thus

    D_return=r^3[K a(r,K,J)+J r^2 b(r,J)]

with positive analytic a,b in a sufficiently small common neighborhood.
For fixed nearby parameters `r^2 b/a` has positive derivative for small
r>0: its leading derivative is `2r b(0)/a(0)>0` and the remaining
terms are uniformly O(r^2). Consequently there is at most one small
isolated nonzero root, it requires KJ<0, and the sign of its displacement
derivative is the sign of J. No local S/U fold pair occurs on this
beta=0 slice. The center's nonisolated periodic ovals are not exceptions
to an assertion about limit cycles.

The remote stability cone is also consistent with these coordinates.
For K>0,J<0 it imposes a positive lower bound on K/|J|, whereas a
small root would require `r^2=(K/|J|)a/b` with a/b tending to
`75(11c*-5)/2>0`. Uniform lower bounds therefore exclude an origin
cycle shrinking to zero inside this cone. This is a local statement;
it supplies neither an explicit amplitude radius nor a finite-amplitude
cycle-count bound.

## 5. Other infinity strata and the rational illustration

The listed infinity types follow from the two eigenvalues
`(-(1+z),p'(z))` and the vertical pair `(-c,1-c)`. The changes at
c=241/250 and c=1 are real changes of singularity portrait. Their
existence neither invalidates all finite returns nor permits transporting
an old itinerary without checking it.

For the illustration `c=1001/1000,m=196/5`, the exact replay checks
the stated endpoint values of F at s=88/25 and 353/100, the K value,
and every coefficient and sign sample of the multiplier quartic N.
The lower endpoint lies above s_H. Hence the trace is negative; its
coarse lower bound -5 and the determinant bound >79 give the stated
strict focus discriminant bound. The slopes
`-600 +/- 100sqrt(37)` and their radial/angular eigenvalues are correct.

The N coefficient signs give two sign changes, while N(1)>0,N(4)<0,
N(64)>0 give two distinct positive roots. Descartes' rule then makes
them the only positive roots and simple. This is a sign chart for a
necessary multiplier density, not a return-displacement chart or a cycle
certificate. No remote U orbit, new stationary-return branch, or common
precursor is established at the rational illustration.

## 6. Replay and scope

Run `check_new_gate_review.py` in the existing q4 virtual environment.
It uses one thread, reduced priority, and a ten-CPU-second fuse. The first
implementation attempted nested rational-matrix expansion and hit that
fuse without completing. It was replaced by the polynomial numerator
check above. The final replay passed in **0.247154 CPU seconds**, checking
the original-coordinate involution modulo J, the remote trace/Hopf
identities, and the rational illustration's exact polynomial data.

This is an exact algebra replay, not an interval ODE calculation. The
hand proofs independently fix the Hopf orientation and positive scaling.
No further numerical work was run, no earlier checkpoint was rewritten,
and no global KKL exclusion or five-cycle claim follows.
