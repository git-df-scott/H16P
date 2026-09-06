#!/usr/bin/env python3
"""Cross-check my eta_2 (V_6 reduced mod eta_1) against Astra/Fable's exact
eta_2 from FABLE_D2_ORDER_TWO_LOOP.md:

  eta_2 = a(b+2l)(b-3l-5)[a^2(b+2l+1)-(b+1)(l+1)^2] / (48 (l+1)^2)

theirs is written after substituting the first focal relation m = a(b+2l)/(l+1);
mine keeps m free.  They must agree on {eta_1 = 0}."""
import sympy as S
l, m, a, b = S.symbols('l m a b', real=True)

mine = a*(6*a**2*b + 12*a**2*l + 3*a*b*m - a*m + b**3 - b**2*l - 4*b**2
          - 6*b*l**2 - 11*b*l - b*m**2 - 5*b - 6*l**2 - 10*l + m**2)/48
theirs = a*(b+2*l)*(b-3*l-5)*(a**2*(b+2*l+1) - (b+1)*(l+1)**2)/(48*(l+1)**2)

msub = a*(b + 2*l)/(l + 1)                     # their eta_1 = 0
print("their first focal relation m =", msub)
print("my eta_1 = a*b+2*a*l-l*m-m at that m:",
      S.simplify((a*b + 2*a*l - l*m - m).subs(m, msub)))
mine_on = S.simplify(S.factor(mine.subs(m, msub)))
theirs_f = S.simplify(S.factor(theirs))
print("\nmy eta_2 on {eta_1=0}:")
print("  ", mine_on)
print("their eta_2:")
print("  ", theirs_f)
print("\ndifference:", S.simplify(mine_on - theirs_f))
print("ratio     :", S.simplify(mine_on/theirs_f))
