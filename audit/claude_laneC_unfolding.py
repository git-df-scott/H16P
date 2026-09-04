#!/usr/bin/env python3
"""Lane C: first-order unfolding signs at Shi's seed. With m=5a+delta,
b=3l+5-9delta+8eps, lambda<0, check eta_1 ~ -eps, eta_2 ~ c*delta with c>0,
eta_3>0, so (lambda, eta_1, eta_2, eta_3) alternate for lambda<0<... i.e.
lambda<0, eps<0, delta<0 give (-,+,-,+): three small cycles in the hierarchy
|lambda|<<|eps|<<|delta|<<1."""
import sympy as S
exec(open(__file__.replace('unfolding','shi_focus')).read().split("for i, e in enumerate(etas, 1):")[0])
d, e = S.symbols('delta epsilon')
sub = {l: -10, a: 1, m: 5+d, b: -25-9*d+8*e}
E1 = S.expand(etas[0].subs(sub)); E2 = S.series(S.expand(etas[1].subs(sub)), d, 0, 2).removeO()
print("eta_1 (exact) =", E1)
print("eta_2 to first order in delta at eps=0:", S.expand(E2.subs(e, 0)))
print("eta_2 first-order coefficient in eps at delta=0:", S.expand(S.diff(etas[1].subs(sub), e)).subs({d: 0, e: 0}))
print("eta_3 at seed:", etas[2].subs(sub).subs({d: 0, e: 0}))
