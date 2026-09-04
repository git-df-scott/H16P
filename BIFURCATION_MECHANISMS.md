# Bifurcation mechanisms for a fifth cycle

## Strategic topology first

For quadratic systems with two foci, Zhang proved that one of the two nests has
at most one cycle. Thus a five-cycle target with two nests must be `4+1`, not
`3+2`. Since Bautin bounds the number born locally from a single quadratic
focus/center by three, the fourth member of the large nest must have a global
origin: a loop, period-annulus boundary, or other large-amplitude event.

This makes the best target a **simultaneous local/global construction**, not a
larger Hopf computation alone.

## Mechanism ledger

| Mechanism | Rigorous quadratic bound in its stated scope | Best relevant construction | Freedom left / fifth-cycle route | Verdict for rigorous numerics |
|---|---|---|---|---|
| Ordinary Hopf | One local cycle at generic Hopf | Standard | Cannot produce five alone | Good local seed only |
| Degenerate Hopf / Bautin | At most 3 small cycles from one quadratic focus or center; bound attained | Three cycles in Shi's weak-focus nest | A fourth in that nest must be nonlocal | Excellent analytic filter; poor standalone attack |
| Two foci / nested cycles | One nest contains at most one cycle | All known four-cycle examples are `3+1` | A five-cycle field must reach `4+1`; `3+2` is excluded | Strong topology pruning |
| Saddle homoclinic loop | Finite cyclicity proved in many generic/special quadratic cases, but no global bound covering every degenerate graphic | Loop-created cycles in near-Hamiltonian analyses | Add one global cycle outside a Bautin three-nest while retaining the remote cycle | High-value but difficult near saddles |
| Heteroclinic/polycycle | Numerous graph-specific cyclicity results; 121-graphics program remains incomplete as a global route | Dumortier--Roussarie--Rousseau program | Degenerate graphics can release cycles under several parameters | Candidate generation plausible; certification expensive |
| Hamiltonian perturbation | The quadratic infinitesimal problem is finite; exact counts depend on center class and perturbation order | Abelian/Melnikov analyses; `Q4` gives upper 5, lower 3 on the annulus | Find 4 or 5 simple zeros, or force first-order vanishing and use higher order | Best algebraically bounded opening |
| Reversible near-integrable two-center | Some subfamilies have exact annulus bounds (e.g. 2); Yu--Han attain four over both centers | Yu--Han `3+1` | One extra global zero in the three-cycle nest; cannot make `3+2` | Good narrow secondary attack |
| Saddle-node of limit cycles | Generic fold creates/destroys a pair | Continuation surfaces in five-parameter normal forms | A pair birth from a four-cycle field tends to give six unless another cycle is lost; codimension bookkeeping is essential | Useful continuation diagnostic, not an unconstrained objective |
| Cycle from infinity | No complete global quadratic bound | Poincare compactification/graphics analyses | A large cycle can enter from an infinite graphic while four finite cycles persist | Hardest enclosure; only pursue with a located graphic |
| Algebraic invariant oval | Strong uniqueness results in several quadratic infinity classes | Known algebraic quadratic cycles are sparse | Five algebraic cycles are not a credible degree-2 route | Low priority |

## Mechanisms that do not transfer

- Five cycles under **cubic perturbations** of a quadratic Hamiltonian system do
  not give a quadratic vector field.
- Five or more cycles in **piecewise quadratic** systems do not give a smooth
  quadratic vector field.
- High cyclicity in **three-dimensional quadratic** systems does not constrain
  or improve planar `H(2)`.
- Five zeros of a truncated radial or Melnikov series do not prove five cycles
  unless the remainder and parameter realization are controlled.

## Natural bifurcation surfaces

Within a normalized five-dimensional chart, the useful low-codimension sets
are:

- trace zero at an antisaddle (Hopf hypersurface);
- vanishing first, second, and third Lyapunov quantities (Bautin strata);
- stable/unstable separatrix matching (homoclinic or heteroclinic
  hypersurfaces, defined by a splitting function);
- double fixed points of a Poincare map (`D=0`, `D'=0`);
- degeneracy of a first Melnikov function and its endpoint values.

The intersections of these surfaces are thin but computable. The campaign
should solve or continue their defining functions, never optimize a generic
trajectory-count score over raw coefficients.
