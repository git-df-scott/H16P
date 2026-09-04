# Four-cycle frontier

## Evidence hierarchy

The frontier contains three different kinds of record:

1. an exact existence proof for a parameter hierarchy;
2. a decimal field supported by ordinary numerical integration;
3. a fixed exact field verified by interval arithmetic.

They must not be merged. The Galias--Tucker field is the cleanest baseline for
rigorous verification; the Kuznetsov field is the easiest baseline for a fast
floating-point pipeline.

## 1. Songling family

\[
\begin{aligned}
\dot x&=\lambda x-y-10x^2+(5+\delta)xy+y^2,\\
\dot y&=x+x^2+(-25+8\varepsilon-9\delta)xy.
\end{aligned}
\]

| Item | Record |
|---|---|
| Parameters | Shi: `0 < -lambda << -epsilon << -delta << 1`; a reported explicit choice is `delta=-10^-13`, `epsilon=-10^-52`, `lambda=-10^-250` |
| Configuration | Three nested cycles around `(0,0)` and one about `(0,1)` (`3+1`) |
| Mechanism | Unfold a third-order weak focus for three small cycles; use a global no-contact/trapping construction for the remote cycle |
| Proof | Shi (1980), Poincare--Bendixson and weak-focus analysis |
| Robustness | Strict trapping/focus inequalities give an open but fantastically thin parameter region; the convenient numeric scales span hundreds of decimal orders |
| Machine data | Coefficients are explicit; the original proof is not packaged as modern replay data |

### Certified Songling point

Galias--Tucker use the same equations with the exact decimal powers

```text
delta   = -10^-13
epsilon = -10^-52
lambda  = -10^-200
```

and prove by rigorous, adaptive-precision interval computation that the system
has exactly four limit cycles, with precise location bounds. This is the best
audited rigorous construction. The article is open access, but this audit did
not locate a maintained public source/data repository containing a one-command
replay; a future campaign should contact the authors or reimplement the proof
with current CAPD before using it as a bulk verifier.

The three tiny cycles make this a poor discovery benchmark. They make it an
excellent arbitrary-precision and enclosure stress test.

## 2. Chen--Wang family

\[
\begin{aligned}
\dot x&=-\delta_2x-y-3x^2+(1-\delta_1)xy+y^2,\\
\dot y&=x+\frac{2}{9}x^2-3xy,
\end{aligned}
\qquad 0<\delta_2\ll\delta_1\ll1.
\]

| Item | Record |
|---|---|
| Configuration | `3+1`, with one large and three small cycles |
| Mechanism | Two trapping regions plus unfolding of a second-order fine focus; an additional Poincare--Bendixson cycle completes the small nest |
| Proof | Analytic/bifurcation proof in Chen--Wang (1979) |
| Robustness | A parameter hierarchy, not a comfortably sized explicit box |
| Machine data | Original paper gives no modern machine-readable certificate; Yu--Zeng use `delta1=0.01`, `delta2=0.00002` for visualization |

## 3. Kuznetsov--Kuznetsova--Leonov five-coefficient family

\[
\begin{aligned}
\dot x&=y+x^2+xy,\\
\dot y&=a x^2+bxy+cy^2+\alpha x+\beta y.
\end{aligned}
\]

A published decimal instance is

```text
a = -10
b = 2.2
c = 0.7
alpha = -72.7778
beta = 0.0015
```

The field has equilibria near `(0,0)` and `(-6.2596,7.4498)`. The construction
is based on an explicit parameter domain and analytical-numerical localization;
the decimal point is visually much easier than the certified Songling point,
but the audited sources do not supply an interval certificate for that decimal
field.

The sufficient-domain template reported by Yu--Zeng is

\[
1<b<3,\quad \tfrac13<c<1,\quad
4a(c-1)>(b-1)^2,\quad bc>1,
\]

with `alpha` just above `a(b+2)/(bc-1)` and positive `beta` on a smaller scale.
The `just above`/`much smaller` qualifiers are exactly where a search must
replace asymptotic notation with explicit boxes.

## 4. Near-integrable two-center family

Yu--Han study

\[
\begin{aligned}
\dot x&=y(1+a_1x)+\eta a_{10}x,\\
\dot y&=-x+x^2+a_4y^2+\eta(b_{01}y+b_{11}xy).
\end{aligned}
\]

At `eta=0` there are two centers in suitable regions. Local focus expansions
and global Melnikov functions give four cycles with `3+1` or `1+3`
distribution for open asymptotic parameter regions. Yu--Zeng later give a
concrete numerical instance, including `a1=-30/7`, `a4=-671/210`,
`a10=1/200`, `b01=-500001/100000000`,
`b11=49182857/968100000`, and `eta=1/100`.

This family is machine-readable and structurally motivated, but the concrete
Yu--Zeng visualization is not itself an interval proof. It is useful for
Melnikov-guided searches rather than a global coefficient sweep.

## 5. Reproduced modern numerical control

[`controls/reproduce_four_cycle.py`](controls/reproduce_four_cycle.py) applies
DOP853 integration and Brent root finding to the Kuznetsov decimal field. On
the downward section `y=0` it finds:

| Cycle | Section `x` | Period | Nontrivial Floquet multiplier | Classification |
|---:|---:|---:|---:|---|
| 1 | 0.683210218760 | 0.728551432450 | 0.999226902941 | stable |
| 2 | 2.183699825305 | 0.704986809930 | 1.002420055138 | unstable |
| 3 | 15.962783983170 | 0.622436541326 | 0.962020809899 | stable |
| 4 | -3711.560806385 | 0.434622844906 | 11.4622677135 | unstable |

The first three surround the origin; the fourth surrounds the second focus.
The output plot is [`controls/four_cycle_control.png`](controls/four_cycle_control.png)
and the raw values are [`controls/four_cycle_control.json`](controls/four_cycle_control.json).

These are floating-point values. In particular, multipliers 1 and 2 are close
to one and the outer orbit crosses scales of thousands. The control is **not**
a proof and must never be cited as one.

## 6. What a rigorous replay would add

For each of the four numerical section roots it would need a rational interval
`I`, a validated return map `P(I)`, an interval-Newton inclusion for `P-id`, a
derivative enclosure excluding `1`, and proof of a unique oriented return.
Pairwise-disjoint flow tubes or isolating annuli would establish distinctness.
The exact protocol is in [`RIGOROUS_CERTIFICATION.md`](RIGOROUS_CERTIFICATION.md).
