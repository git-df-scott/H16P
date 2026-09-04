# Claude lanes B and C: first results

2026-09-04. Parallel to Astra Strike #4 (Theorem N). Nothing here overlaps
with Astra's assignment. Scripts: `audit/claude_laneB_loop_functionals.py`,
`audit/claude_laneC_shi_focus.py`, `audit/claude_laneC_focus_numeric.py`.

## Lane B: can the Q4 saddle loop supply the missing cycles?

**Question.** Zhao's bound and Theorem N concern zeros of the Melnikov
integral in the *open* annulus. Cycles can also be born at the homoclinic
loop (the `t=1` endpoint). If Theorem N holds (at most three interior
zeros), could three interior cycles plus two loop-born cycles give five?

**Literature.** Han (Sci. China A 40, 1997) proves that the cyclicity of a
homoclinic loop of a quadratic integrable non-Hamiltonian system is two,
with one exceptional case. Han, Yang, Tarta and Gao (JDDE 2008) give the
first four Dulac coefficients of the Melnikov function at a loop,
`I=c_0+c_1w\log w+c_2w+c_3w^2\log w+\dots`, and the Leontovich–Roussarie
theorem: a loop where `c_0=\dots=c_{k-1}=0\ne c_k` produces at most `k`
cycles. Alien limit cycles (cycles not detected by any Melnikov function)
are known to require a polycycle with at least two saddles; they do not
occur at a single saddle loop, a periodic orbit, or a non-degenerate
singularity (Caubergh–Dumortier–Roussarie; Gavrilov).

**Structure in the reconstruction variables.** Near `t=1`,
`Y=\Phi y+P y_2` with `y\to0`, `y_2(1)` finite, and
`P\sim-\Omega_0H(1)\log(1-t)`. Hence

- `c_1` (the `w\log w` coefficient) is proportional to `H(1)`, the exact
  affine functional `\frac{18}{85085\pi}(9061A+6289B-2431\eta-7242)`;
- `c_0=I(\text{loop})=-\frac{aC}{2}\sqrt{1-a}\,X(1)`, affine in `(A,B,\eta)`.

Measured: the map `(A,B,\eta)\mapsto(c_0,c_1)` has rank two at
`kappa=2,4,9` (the two loop functionals are independent, so both can be
killed by coefficient choice at fixed `kappa`).

**Consequences.**

1. On the strict lobe region `H(1)<0`, so `c_1\ne0` and the loop can
   produce at most **one** cycle at first order for any lobe point. The
   three-interior-plus-two-loop configuration would need `H(1)=0`, i.e. a
   primitive root *at* the loop, which is exactly the late-root corner
   `(94/77,-17/77,1)` that Astra's Strike 3 and the audit already analysed:
   there the interior roots have merged into the loop and (S1) fails.
2. The line `c_0=c_1=0` (two loop-born cycles) was computed at
   `kappa=2,4,9` for `\eta=1,1.2,1.5`. It lies far outside the lobe box
   (`A<0` or `B>1`, `q(0)<0`), and the corresponding primitive and integral
   have zero or one interior sign change. So two loop cycles come at the
   price of at most one interior zero: total three, not five.
3. Because alien cycles are excluded at a saddle loop, the cycles near the
   loop for parameters near a degenerate point are the small zeros of the
   first non-vanishing Melnikov function, whose loop expansion is
   controlled by the same functionals. Since `(c_0,c_1)` is surjective from
   the four-term space, any near-loop zero structure produced by
   `I_{\mu}+\epsilon M_2` is reproduced by some `I_{\mu''}` in the space,
   and the total count is again a zero count of one member of the space:
   bounded by Zhao (five with multiplicity) and, if Theorem N holds, by
   three distinct.

**Lane B verdict.** The loop is not an escape. If Theorem N is proved, the
closed Q4 annulus (loop included) cannot produce five cycles, and Q4 is
finished as a counterexample route. The one caveat, recorded for
completeness: Han's exceptional quadratic integrable case was not
identified from the abstract; if Q4's loop were that case, the loop
cyclicity bound would change but the surjectivity argument in item 3 would
not.

## Lane C: Attack 2 groundwork (Shi chart, third-order weak focus)

Chart: `x'=\lambda x-y+lx^2+mxy+y^2`, `y'=x+ax^2+bxy`.

**Exact focus quantities** (Lyapunov-function method, `sympy`, degree 8):

```text
eta_1 = -(ab+2al-lm-m)/8                       -> eta_1=0  <=>  m=a(b+2l)/(l+1)
eta_2 | eta_1=0 = a(b+2l)(b-3l-5)(a^2(b+2l+1)-(b+1)(l+1)^2) / (48(l+1)^2)
eta_3 | eta_1=eta_2=0, b=3l+5 = -25a(2a^2+l+2)(5a^2 l+6a^2-3l^3-12l^2-15l-6)/64
```

On the stratum `b=3l+5` the first condition gives `m=5a`, confirming the
ATTACK_MATRIX stratum. The factors `b=-2l` and `a^2(b+2l+1)=(b+1)(l+1)^2`
of `eta_2` are center strata (they also kill `eta_3`), not weak-focus
strata; an earlier version of the script picked `b=-2l` and reported
`eta_3\equiv0`, which was the wrong root, now fixed.

At Shi's seed `(l,a)=(-10,1)`: `eta_1=eta_2=0`, `eta_3=35625/8>0`, an
unstable third-order weak focus. Independent 30-digit return-map
integration gives `d(r)=4.42\cdot10^{-8},\,7.21\cdot10^{-6},\,1.87\cdot10^{-3}`
at `r=0.02,0.04,0.08`, i.e. `d(2r)/d(r)=163,\,259`, consistent with
`d\sim c\,r^7` (order three) with `c\approx2\pi\eta_3` plus higher-order
corrections. The `eta_3=0` locus on the stratum is
`a^2=3(l+2)/(5l+6)(l+1)^2`... explicitly `5a^2l+6a^2=3l^3+12l^2+15l+6`
or `2a^2+l+2=0`; the seed is away from both.

**Next Lane C steps** (not started): exact unfolding
`m=5a+\delta`, `b=3l+5-9\delta+8\epsilon`, `\lambda<0` with the dyadic
hierarchy of ATTACK_MATRIX.md; verify the sign pattern
`eta_1,eta_2,eta_3,\lambda` alternates to produce three small cycles;
locate the remote cycle about `(0,1)`; continue the outer separatrix
splitting in `(l,a)` toward a fourth cycle in the same nest. This is the
`4+1` target; it is independent of Q4 and of Theorem N.

## Lane C, second result: the loop around the weak focus exists only at a center

Scripts: `audit/claude_laneC_stratum_saddles.py`, `claude_laneC_saddle_region.py`,
`claude_laneC_splitting3.py`, `claude_laneC_splitting4.py`.

1. On the stratum `m=5a`, `b=3l+5` the finite equilibria other than `(0,0)`
   and `(0,1)` solve a quadratic in `x` with discriminant
   `12(3l+5)^2(3a^2-l^2-2l)`. Finite saddles exist only for
   `-1-\sqrt{1+3a^2}<l<-1+\sqrt{1+3a^2}`. The whole Attack-2 box
   (`l\in[-12,-8]`, `a\in[4/5,6/5]`) has none: only two foci. No finite
   separatrix loop can bound the origin nest there.
2. In the finite-saddle region a continuous signed splitting of the saddle
   whose stable branch comes from the origin nest was computed (crossing
   angles of the true stable and returning unstable branches on a circle of
   radius `0.05`). It changes sign exactly once in `l` for each
   `a\in\{1,1.5,2,3\}`, at

   | `a` | loop `l^*` | center curve root of `5a^2l+6a^2=3l^3+12l^2+15l+6` | difference | saddle trace |
   |---:|---:|---:|---:|---:|
   | 1.0 | `-1.183503419072` | `-1.183503419072` | `4e-15` | `1e-14` |
   | 1.5 | `-1.192053160605` | `-1.192053160605` | `2e-14` | `6e-14` |
   | 2.0 | `-1.195392237464` | `-1.195392237464` | `4e-14` | `7e-14` |
   | 3.0 | `-1.197905643736` | `-1.197905643736` | `6e-14` | `8e-14` |

   The homoclinic loop around the origin occurs precisely where
   `eta_3=0`, i.e. where `eta_1=eta_2=eta_3=0` and the origin is a **center**
   (Bautin), and the loop saddle then has zero trace. Off the center curve a
   genuine third-order weak focus is never surrounded by a finite saddle
   loop on this stratum (numerically, to machine precision at four values
   of `a`; a proof should follow from Li–Cherkas plus the fact that a loop
   is a limit of cycles in the unfolding).
3. Consequence for Attack 2. "Third-order weak focus plus a finite outer
   separatrix loop in the same nest" is not a configuration of quadratic
   systems; the only way to have both is at a quadratic center, and then
   every cycle in that nest is a perturbation-born cycle of a quadratic
   center's period annulus, governed by the Melnikov/Bautin-ideal cyclicity
   of that center component. The `4+1` target through Attack 2 therefore
   reduces to the integrable-perturbation program. Which center component
   the curve belongs to is determined in `audit/claude_center_identify.py`
   (see CLAUDE_THOUGHT_SESSION.md).
