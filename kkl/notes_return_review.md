# KKL return sensitivities and fold review

Independent Astra analytic review, 2026-09-04. No orbit integrations or
numerical searches were run. This note supplies formulas and control
expectations for the authorized bounded pilot; none is an interval return
certificate.

## 1. Field, derivatives and fixed-time sensitivities

Use `b=11/5` and write the vector field as `F=(f,g)`, reserving R for the
scalar section return:

    f = y+x^2+xy,
    g = -10x^2+bxy+c y^2+alpha x+beta y.

The pilot keeps beta fixed at zero. All formulas below allow a fixed
nonzero beta as well. Along a nominal orbit z(t),

    J = [[2x+y,             1+x],
         [-20x+by+alpha, bx+2cy+beta]],
    div F = (2+b)x+(1+2c)y+beta.

For vectors v,w the symmetric state Hessian is

    Hess F[v,w] = (
      2 vx wx + vx wy + vy wx,
      -20 vx wx + b(vx wy+vy wx) + 2c vy wy).

The explicit parameter derivatives and their state Jacobians are

    F_c=(0,y^2),       J_c=[[0,0],[0,2y]],
    F_alpha=(0,x),     J_alpha=[[0,0],[1,0]].

All pure explicit second parameter derivatives vanish. If beta is later
unfolded, `F_beta=(0,y)` and `J_beta=[[0,0],[0,1]]`.

Start at `z(0)=(r,0)`. For j in `{r,c,alpha}`, let w_j denote the derivative
of the flow at **fixed physical time**, not at a moving return event:

    w_r' = J w_r,                    w_r(0)=(1,0),
    w_c' = J w_c + F_c,              w_c(0)=(0,0),
    w_alpha' = J w_alpha + F_alpha,  w_alpha(0)=(0,0).

Also integrate `I'=div F`, `I(0)=0`. The two orbit coordinates, these six
sensitivity coordinates and I give nine scalar state variables.

## 2. Moving-event correction

Let T be the selected full downward return to y=0, with endpoint
`z(T)=(R,0)` and `g_T != 0`. Differentiate the event equation:

    T_j = -w_j,y(T)/g_T,
    R_j = w_j,x(T) + f_T T_j
        = w_j,x(T) - (f_T/g_T) w_j,y(T).

Thus the projection matrix is

    Pi_T = Id - F_T (0,1)/g_T,

and `Pi_T w_j=(R_j,0)`. It applies equally to initial-coordinate and
parameter sensitivities. Using only the x component of w_j omits the
moving crossing time.

The independent radial identity is

    R_r = [g(r,0)/g(R,0)] exp(I(T)).

On one continuous downward itinerary both section speeds are negative,
so R_r is strictly positive. At an actual periodic point R=r the speed
ratio is one. Away from a fixed point the bare divergence exponential is
not the section derivative. Agreement of this identity with the projected
variational derivative is a meaningful same-orbit implementation check.

All expressions require a nonzero endpoint section speed and a differentiable
full-return itinerary. They do not extend across nonreturns or event-order
changes. A small numerical residual is not exact periodicity.

If integration uses an opposite crossing followed by the desired crossing,
one consistent convention is to continue the entire augmented **fixed-time**
state through the opposite crossing and apply Pi only at the final return.
This is equivalent to interrupting and resuming the same variational ODE.
Alternatively one can compose two projected half-return maps, carrying
their parameter and second-derivative chain rules. Do not project at the
first crossing and then treat those data as the original unprojected
fixed-time sensitivities. Discarding/reinitializing sensitivity columns at
the intermediate event is also incorrect.

## 3. Second sensitivities needed for folds

Let w_ij be a fixed-time second derivative. For i,j in `{r,c,alpha}`,

    w_ij' = J w_ij + Hess F[w_i,w_j]
            + J_i w_j + J_j w_i,
    w_ij(0)=0,

where `J_r=0`. For a fold evaluator only w_rr, w_rc and w_ralpha are
required, adding six scalar ODE coordinates. Parameter-parameter second
derivatives are optional. If they are used, the displayed equation includes
the factor two when i=j is a parameter.

At the final event put `F_r=0`, retain the explicit F_c,F_alpha above, and
define

    E_ij = w_ij + (J w_i+F_i) T_j + (J w_j+F_j) T_i
                   + (J F) T_i T_j.
    T_ij = -E_ij,y/g_T,
    R_ij = E_ij,x - (f_T/g_T) E_ij,y.

All quantities on the right are evaluated at T. This follows by twice
differentiating z(T(r,p);r,p): `z_tt=JF` and
`partial_t w_i=Jw_i+F_i`. The explicit F_i terms matter for mixed
initial-coordinate/parameter derivatives. Projecting w_ij alone is not
the correct second return derivative.

One can integrate directly with initial log coordinate u, but then
`w_u(0)=(r,0)` and **w_uu(0)=(r,0)**. The second initial derivative does
not vanish because r=exp(u). Integrating r derivatives and converting
after the final projection avoids this easy mistake.

## 4. Log-radius residuals and continuation

For an origin return require r>0 and R>0 on the selected itinerary. Set
`u=log r` and choose the dimensionless residual

    L(u,p) = log(R(exp(u),p)/exp(u)).

For near returns compute it as `log1p((R-r)/r)` if that improves numerical
accuracy; it does not remove the error already present in R. The exact
derivatives are

    L_u = r R_r/R - 1,
    L_p = R_p/R,
    L_uu = r R_r/R + r^2 R_rr/R - (r R_r/R)^2,
    L_up = r R_rp/R - r R_r R_p/R^2.

At a periodic point, `L_u=R_r-1`. At a fold additionally R_r=1, so

    L_uu = r R_rr,
    L_p = R_p/r,
    L_up = R_rp - R_p/r.

Do not use these last simplified expressions off the root/fold. The
alternative residual `G=(R-r)/r` has the same values of these derivatives
at a fold, but its off-fold derivatives differ. The unnormalized residual
`H=R-r` instead has `H_u=r(R_r-1)` and `H_uu=r^2 R_rr` at a fold.
These conventions cannot be mixed in a Newton corrector.

For a specified shape path p(s)=(c(s),alpha(s)), a regular root obeys

    du/ds = -(L_c c'(s)+L_alpha alpha'(s))/L_u.

This formula becomes singular near a fold. Use a pseudo-arclength
corrector there, with the path/slice stated explicitly. A single equation
L=0 in (u,c,alpha) defines a surface; it does not select a unique
continuation direction without an additional path or slice condition.

Once a fold is actually located, its equations are L=0,L_u=0. Their
Jacobian with columns (u,c,alpha) is

    [[L_u,  L_c,  L_alpha],
     [L_uu, L_uc, L_ualpha]].

At a nondegenerate fold `L_uu != 0` and `(L_c,L_alpha) != (0,0)`, this
matrix has rank two and defines a local curve. Choose a normalized tangent
in its nullspace, oriented consistently with the previous tangent, and
use a third pseudo-arclength equation for its corrector.

For the remote negative section one can instead take `u=log|r|`,
`r=-exp(u)` and `L=log|R|-log|r|`, requiring R<0 on the same itinerary.
The displayed derivative formulas remain valid because dr/du=r, but
increasing u now moves farther left in x. The origin ordering statements
below refer to positive radii only.

## 5. Which side of a fold adds which pair

At an origin fold (u0,p0), choose a shape direction v transverse to the fold
image. Write `p=p0+delta v`, `A=grad_p L dot v != 0`, and `B=L_uu != 0`.
The leading normal form is

    L = A delta + (B/2)(u-u0)^2 + higher-order terms.

The two-root side satisfies `A delta B < 0`; the leading offsets are
`u-u0 = +/- sqrt(-2 A delta/B)`. Because at each actual root
`R_r=1+L_u`, the lower-radius root is stable and the upper unstable when
B>0. For B<0 the order is U/S. These are local signs for a sufficiently
small generic crossing; a finite step still needs actual root checks.

Starting from one stable origin cycle, a new S/U pair inside it or a new
U/S pair outside it can supply the desired S/U/S triple. The original
stable cycle and remote unstable cycle must coexist at the SAME shape.
A fold that annihilates the only known stable cycle is not a precursor.
Nor does a fold equation locate an unseen disconnected pair automatically.

At beta=0, K>0 makes the local focus repelling. Therefore sufficiently
near r=0 the displacement R-r is positive. Three simple origin cycles in
the intended configuration have successive displacement signs `+,-,+,-`
and multipliers S/U/S. This is a consistency requirement, not a global
root-count theorem.

## 6. Meaningful controls within the charged pilot

Use existing controls and same-orbit identities before spending calls on
additional parameter paths. All finite-difference replays, failed returns
and derivative evaluations count toward the pilot budget.

1. **Linear-focus analytic check.** If the quadratic terms are suppressed,
   let `nu=sqrt(-alpha-beta^2/4)>0`. The exact downward full return is
   `R=exp(pi beta/nu) r`, with `T=2pi/nu`. Its derivative is
   `exp(pi beta/nu)`, whereas the bare divergence exponential is
   `exp(2pi beta/nu)`. The section-speed correction gives the missing
   factor. At beta=0, R_r=1, R_c=R_alpha=0 and
   `T_alpha=pi/(-alpha)^(3/2)`. This tests the event projection and the
   parameter-time derivative without relying on a numerical KKL cycle.
2. **Same-orbit derivative identity.** Compare projected R_r with
   `g(r,0)/g(R,0)*exp(I)` at a non-fixed section point as well as a root.
   A root-only test will not detect omission of the speed correction.
   Also require positive R_r whenever both crossing speeds are negative.
3. **KKL incumbent control.** At `(c,alpha,beta)=(7/10,-363889/5000,3/2000)`
   the saved approximate root/multiplier pairs are
   `(0.683210217,0.999226903)`, `(2.183699825,1.002420055)`,
   `(15.962783982,0.962020810)` and remote
   `(-3711.560806,11.46226773)`. They are numerical reference values,
   not exact expected equalities or interval bounds. The two weak
   multipliers require relative error well below their distance from one
   before inferring their stability numerically.
4. **K>0 actual starting control.** At `(7/10,-80,0)`, K=6/5. The saved
   origin root is about 64.555434341 with multiplier 0.809691136, and the
   separate remote control is about -5391.141160 with multiplier 12.168019.
   The pilot must reproduce both at one shape and carry a remote ledger.
   These are one-origin-plus-one-remote controls, not the desired precursor.
5. **Beta-zero negative control.** The incumbent shape at beta=0 has
   K=-674997/250000, with detected roots about 3.068845425 and
   15.064071451, multipliers about 1.005012827 and 0.967338367. A spurious
   persistent third root near the origin is a warning. The saved static
   probe is not an exhaustive two-root theorem.
6. **A small charged derivative spot-check.** At a safe point on the
   already verified itinerary, compare one variational parameter derivative
   with symmetric parameter differences at two step sizes. If second
   derivatives drive a putative fold, likewise compare a directional
   derivative of R_r with R_rr/R_rp. Near-escape or nearly tangent sections
   are unsuitable for this check. Agreement is numerical validation of the
   evaluator, not an interval proof; do not run a new coefficient grid.

An additional analytic local check follows from the cubic Hopf coefficient.
At beta=0 the full return has

    R-r = C r^3 + O(r^4),
    C = pi K/[4(-alpha)^(3/2)].

Thus `R_r=1+3C r^2+O(r^3)`, while the bare divergence exponential is
`1+4C r^2+O(r^3)` on non-fixed small returns. Their differing cubic-focus
coefficients provide another check on projection. Also
`R_c = [pi b/(4 sqrt(-alpha))] r^3+O(r^4)` and

    R_alpha = (pi/4)[-(bc-1)/omega^3 + 3K/(2omega^5)] r^3+O(r^4).

These are asymptotic expectations only. No explicit remainder radius is
proved here, and cancellation at very small r can defeat floating-point
verification. Do not count arbitrarily tiny numerical roots as cycles.

## 7. Scope of the proposed pilot

The 64-evaluation/16-step pilot can establish evaluator consistency,
continue a seeded branch over a small stated path, and record whether its
remote control survives. A fold candidate requires the corrected first
derivative, nonzero second derivative and an oriented shape crossing.
Only a common-parameter S/U/S plus remote-U ledger is precursor evidence.
No missing fold or failed return establishes an exclusion of this entire
box. Interval existence, derivative and persistence gates remain required
before a mathematical construction claim.

## 8. Independent checks and the selected geometry path

The PF agent independently checked sections 2–5 by differentiation, including
the explicit parameter terms in E_ij and the fold-only log simplifications;
no defect was found. No numerical or symbolic execution was used in that
cross-review.

The initial nine-state `return_map.py` implementation was also read:
its initial columns, variational sources, final projection, independent
section-speed formula and unprojected intermediate-event continuation match
sections 1–2. This is a code review of the first-derivative evaluator,
not a validation of a numerical result or a second-derivative fold solver.

For the independently checked constant-K geometry segment in
`notes_geometry.md`,

    alpha(c)=-216/(11c-5),
    alpha'(c)=2376/(11c-5)^2,

the regular-root predictor is

    dr/dc = -(R_c+alpha'(c) R_alpha)/(R_r-1),
    du/dc = (dr/dc)/r.

The predictor uses derivatives at the actual continued root. A small
denominator signals that this graph parameterization may fail and motivates
pseudo-arclength; it does not itself prove a fold. The exact geometry proof
preserves the equilibrium and infinity gates along this path, while the
origin and remote periodic orbits still need their numerical/validated
return ledger at each accepted common parameter.

## 9. Stable transverse determinants for large remote returns

The endpoint projection in section 2 can subtract two large numbers even
when the transverse derivative is moderate. An equivalent scalar evolution
avoids that subtraction. Define, at fixed time,

    h_j = det(F,w_j) = f w_j,y - g w_j,x.

For any 2x2 matrix J,

    det(JF,w_j)+det(F,Jw_j) = trace(J) det(F,w_j).

Since `F'=JF` along the nominal orbit and `w_j'=Jw_j+F_j`,

    h_j' = div(F) h_j + det(F,F_j).

In the KKL field the three equations and initial values are

    h_r'     = div(F) h_r,                h_r(0)=-g(r,0),
    h_c'     = div(F) h_c + f y^2,        h_c(0)=0,
    h_alpha' = div(F) h_alpha + f x,      h_alpha(0)=0.

There is no extra parameter term in `F'=JF`: parameters are constant while
the nominal orbit evolves in time. At the return event,
`w_j+F T_j=(R_j,0)`. Adding a multiple of F does not change its determinant
with F, so

    h_j(T) = det(F_T,(R_j,0)) = -g_T R_j,
    R_j = -h_j(T)/g_T.

This proves the determinant method gives the same first-return derivatives
as the moving-event projection, including both shape derivatives. It does
not assume periodicity and is valid on every regular selected return.
Continue the determinant variables through the intermediate crossing
without restarting them, exactly as for the fixed-time sensitivities.

The homogeneous radial equation gives the identity

    h_r(t) = -g(r,0) exp(I(t)).

For a downward start this is positive, and for a downward return it yields
positive R_r. Comparing this scalar evolution with the exponential identity
is a useful consistency check, but they express the same scalar variational
law and are not independent evidence for a correct orbit itinerary. The
original fixed-time projection remains a diagnostic for cancellation.
Tightened tolerances, meaningful control returns and, ultimately, validated
enclosures remain necessary; the reformulation does not justify relaxing a
fold or stability sign threshold.

The parameter equations can still experience cancellation in their signed
forcing integrals. This method specifically removes the potentially severe
final subtraction `w_x-(f/g)w_y`; it is not a general conditioning theorem.
No new formula for second return derivatives is asserted by this section.

## 10. Independent audit of the K=0 fifth-order return term

The following checks the PF lane's local calculation in `notes_lienard.md`
by hand; no numerical or symbolic execution was used. Here `m=-alpha`,
`u=1+x`, `d=16-10c`, and `QL=f_L/g_L` is the Liénard quotient, distinct
from the section return R above.

At K=0, `m=210/(11c-5)`. With W and N from the Liénard note,

    (QL)_x = -u^c N/(5W^2),
    N=N1 x+N2 x^2+O(x^3),
    N1=m(-10c^2+38c-10)+210c-4662/5,
    N2=m(-20c^2+43c-5)-100c^2+(4056/5)c-1430.

Direct expansion gives

    QL(x)=QL(0)-N1 x^2/(10m^2)
          -[N2+(c-4-20/m)N1]x^3/(15m^2)+O(x^4).

The harmonic-energy coordinate `q=sign(z)sqrt(2V(z))` satisfies

    q=sqrt(m) x[1+A x+O(x^2)],
    A=(10/m-2c-1)/3.

The positive time change `d tau/dt=g_L/q` gives
`q''+q QL(q)q'+q=0`. Consequently the q^3 coefficient of QL is

    -Delta/(15m^(7/2)),
    Delta=N2-(3-3c+30/m)N1=N2-d N1/7.

Substitution of m and expansion of the displayed N1,N2 yield

    Delta = 4J(c)/[5(11c-5)],
    J(c)=305+634c-11c^2-1000c^3.

At K=0 the lower terms of the damping `q QL(q)` are odd and hence define
a reversible center. Its first even term is h4 q^4 with
`h4=-Delta/(15m^(7/2))`. The first nonzero amplitude change is
`-(pi/8)h4 rho^5`, since
`integral_0^(2pi) cos^4(theta) sin^2(theta) dtheta=pi/8`.
Lower odd damping terms do not change this leading even-term coefficient;
their interactions with it occur at higher degree. The original section
has harmonic amplitude `rho=sqrt(m)r+O(r^2)`, so conjugating back gives

    R(r)-r = pi Delta r^5/(120m^(3/2)) + O(r^6),  K=0.

The exact values `J(33/40)=103619/400` and
`J(9/10)=13769/100` are positive. Thus the K=0 weak focus is repelling
of order two at each of these shapes.

The local no-collapse consequence is a one-sided neighborhood statement,
not an assertion for arbitrary finite r. Analytic dependence of the return
map gives near either shape

    R(r;c,K)-r = r^3 [K A0(c,K,r)+r^2 B0(c,r)],

where A0 and B0 are analytic and positive after shrinking the neighborhood.
At r=K=0 the first is the positive cubic coefficient derivative, and the
second is the positive fifth-order coefficient above. For K>=0 and
sufficiently small r>0 this displacement is strictly positive. Therefore
no nonzero cycle can shrink into the origin through this local K>=0
region. The statement does not exclude finite-amplitude folds, escape,
remote loss, or another K=0 shape with a different J sign. No explicit
numerical radius or parameter-neighborhood bound was certified here.

## 11. Final bounded review of the 202-evaluation checkpoint

Read-only code/ledger review found no blocking defect. The determinant
first derivatives and the fixed-time rr, rc and ralpha second-event
projection implement the formulas above. The exact Fraction equilibrium
gates are valid; the corrector uses the off-root log derivative, rather
than a root-only simplification. No ODE was rerun by this reviewer.

The ledger has 202 sequential charged evaluations. The first 64 reproduce
the saved pilot totals: 1.523849 seconds inside the evaluators and
20.881385836 seconds including subprocess wall time. Totals for all 202
are 6.683029 evaluator CPU seconds and 65.515520670 subprocess wall
seconds. The path files contain 24 attempted common-parameter records,
22 accepted. Two rejected path points record a derivative discrepancy
before the determinant upgrade and the remote section-range boundary.
These are not mathematical cycle-loss claims. All 202 evaluator outputs
themselves have NUMERICAL_ONLY status.

The first determinant-output record is 132. Historical records contain no
evaluator hash or timestamp fields; their absence is preserved. The
supervisor was subsequently extended to add such fields to future calls.
Frozen current-source hashes must not be represented as contemporaneous
hashes of the historical evaluator versions.

The last accepted exact field is
`(c,alpha,beta)=(9301/10000,-8403125/209244,0)`, with K=1/64.
Evaluations 183 and 186 support the reported stable origin and unstable
remote roots. The conservative saved numerical boundary sign pair uses
`c=0.9301046126889` and `c=0.9301046326889` at fixed `r=-2^20`,
with opposite displacements approximately -0.565318 and +0.565357.
This is numerical parameter bracketing, not an interval certificate or
proof that a periodic orbit ceases to exist at the arbitrary section cap.

Evaluation 202 has R_r approximately one, R_rr approximately -0.001122065,
and displacement approximately +0.242693252. It is a numerical local
maximum of the displacement, not a fold. Neither that observation nor
the sampled final profile gives an exact zero count. The incumbent and
K=0 auxiliary controls deliberately lie outside the precursor's K margin
and must remain labelled as controls, distinct from the in-box continuation.
