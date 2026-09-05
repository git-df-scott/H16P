# A complete rational curved section for the two KKL nests

Status: **elementary analytic theorem; exact identity replay**. This repairs
a coverage gap in a proposed future experiment. No return was evaluated
on the new section and no cycle-existence claim follows.

Consider the same family

\[
 P=y+x^2+xy,\qquad
 Q=-10x^2+\frac{11}{5}xy+cy^2+\alpha x+\beta y.
\]

Assume the two finite-focus gates hold, `c<61/5`, and the cubic

\[
 T_\beta(x)=(c-61/5)x^3+(\alpha-111/5-\beta)x^2
             +(2\alpha-10-\beta)x+\alpha
\]

has exactly one simple real root `x_*<-1`. These conditions hold at the
incumbent and at the accepted beta-zero branch points. The beta-zero
cubic/root geometry follows from the finite-equilibrium arguments already
recorded. Focus type and attracting remote trace remain separate gates;
they are not automatic throughout the beta-zero parameter box.

## Section theorem

Use the rationally parametrized nullcline

\[
 \sigma(r)=\left(r,-\frac{r^2}{1+r}\right).
\]

Every origin cycle intersects the branch `r>0` exactly once, transversely,
at its maximum x. Every remote cycle intersects the branch `r<x_*`
exactly once, transversely, at its minimum x. A rational value of r gives
an exact rational point of the section.

**Proof.** The barrier `x=-1` has P=1, so a periodic orbit cannot cross it.
On either half-plane the change of variables `(x,y)->(x,v=P(x,y))` is a
smooth diffeomorphism. Its Jacobian is `1+x`, never zero there. Direct
substitution gives

\[
 Q(\sigma(r))=\frac{rT_\beta(r)}{(1+r)^2},\qquad
 \dot P(\sigma(r))=(1+r)Q(\sigma(r)).                 \tag{S1}
\]

The negative leading coefficient and the single real root imply
`T_beta>0` on `r<x_*` and `T_beta<0` on `r>x_*`. Thus on the remote
branch Q<0 and Pdot>0: every intersection is a strict x minimum. On the
origin branch Q<0 and Pdot<0: every intersection is a strict x maximum.
Both branches are transverse, with no tangency away from the equilibrium.

A cycle enclosing the respective focus must have such an extremum. In
the `(x,v)` plane the section branch is a ray from that equilibrium on
the line v=0. All its crossings have the same orientation by (S1).
A Jordan curve enclosing the ray's endpoint has signed intersection
number of absolute value one with that ray. Therefore the transverse
crossing count is exactly one. This also prevents counting several
intersections of one cycle as different cycles. The orientation reversal
of the coordinate map on x<-1 does not change this absolute count. QED.

The theorem does not require a cycle to intersect the old horizontal
line y=0. It does not say every starting point on the new section has a
full return; escape, guards and changes of return domain still matter.

## Derivative and implementation consequences

For a full return with coordinate R(r), the usual planar determinant
identity gives

\[
 R'(r)=\frac{\det(F(\sigma(r)),\sigma'(r))}
              {\det(F(\sigma(R)),\sigma'(R))}
       \exp\!\int_0^{T(r)}\operatorname{div}F\,dt
       =\frac{Q(\sigma(r))}{Q(\sigma(R))}
        \exp\!\int_0^{T(r)}\operatorname{div}F\,dt .  \tag{S2}
\]

The last equality uses P=0 and `sigma'_x=1`, so each determinant is
exactly -Q. This is an off-root return derivative; at a fixed point it
reduces to the multiplier. The old engine's transverse-determinant
equations therefore remain applicable after changing the initial state,
initial radial sensitivity, and events.

Specifically, the radial initial sensitivity is
`sigma'(r)=(1,-r(r+2)/(1+r)^2)`, not `(1,0)`. Its initial second
derivative is `(0,-2/(1+r)^3)`. The section itself is independent of
`c,alpha,beta`, so their initial sensitivities vanish. Events are now
P=0. An origin maximum starts into P<0; the next minimum crosses into
P>0, then the complete return crosses into P<0. A remote minimum uses
the reverse sequence. A full-turn/branch check is still required.

At the final event the returned coordinate is x and P=0. Hence its
first event-time correction is `P*T_j=0`; the fixed-time x sensitivity
already equals the section-coordinate sensitivity there. The independent
determinant formula (S2) should still be checked. Second derivatives must
include the moving-event terms; this observation does not authorize
omitting them.

If an orbit meets both the old and the new sections, first-hit transport
locally conjugates their return maps and preserves its multiplier. That
transport must actually be computed and checked on the control cycles;
the theorem is not a substitute for an implementation benchmark.

## Scope of a revised experiment

The old handoff specifies downward y=0 sections and coordinate caps on
those sections. A new r here measures an x extremum; it is not the old
downward-crossing coordinate. The caps must be stated for the revised
chart or transported explicitly. Neither a section change nor the
existence theorem for the ray silently expands the allowed field box,
trajectory-size cap, or computation budget.

For a final certificate one can use rigorous returns to this smooth
rational curve, or transfer each isolated cycle to a small rational
straight section with a validated first-hit map. In either case the
field's coefficients remain exact rationals. No coefficients have been
changed and no experimental scope has been expanded by this note.

The identities, initial section derivatives and transversality formula are
checked by `check_nullcline_section.py`. The global intersection argument
and the beta-dependent extension were independently reviewed. See also
[COVERAGE_AUDIT.md](COVERAGE_AUDIT.md). The script does not validate return
parameter sensitivities or their moving-event implementation.
