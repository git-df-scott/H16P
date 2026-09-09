"""Exact certificates for the KKL origin-cycle nonexistence gate.

Zero ODE calls, zero floating-point decisions.  Every assertion below is
decided in exact rational / symbolic arithmetic by SymPy.  Numerical values
are printed only as decimal renderings of exact rationals.

Objects (beta-zero KKL family, xdot=(1+x)y+x^2, ydot=-m x-10x^2+(11/5)xy+c y^2):

    u = 1+x,  e = 11c-5,  d = 16-10c,  m = 5(K+42)/e,
    W = m+(2m+10)x+(m+111/5)x^2+(61/5-c)x^3,
    N = {d u+(c+1)(21+d x)}W - u(21+d x)W',
    psi(x) = (21+d x) u^{c+1} / (5 W)        (= -f/g in the Lienard chart),
    lambda0 = psi(0) = 21/(5m),
    Theta(x) = 5 m W (psi - lambda0) = m(21+d x) u^{c+1} - 21 W.
"""
import json
from pathlib import Path
import sympy as sp

x, c, K, m, s, t = sp.symbols('x c K m s t')

e_ = 11*c - 5
d_ = 16 - 10*c
h_ = sp.Rational(61, 5) - c
m_ = 5*(K + 42)/e_
W_ = m + (2*m + 10)*x + (m + sp.Rational(111, 5))*x**2 + h_*x**3
N_ = sp.expand((d_*(1+x) + (c+1)*(21 + d_*x))*W_ - (1+x)*(21 + d_*x)*sp.diff(W_, x))
Theta_ = m*(21 + d_*x)*(1+x)**(c+1) - 21*W_

out = {'arithmetic': 'exact rational SymPy; zero ODE calls'}

# ---------------------------------------------------------------- identities
# psi' (in x) = u^c N / (5 W^2); dz/dx = u^{-c-1} so d psi/dz = u^{2c+1} N/(5W^2).
psi = (21 + d_*x)*(1+x)**(c+1)/(5*W_)
assert sp.simplify(sp.cancel(sp.expand(sp.diff(psi, x)/(1+x)**c)*5*W_**2) - N_) == 0
out['psi_derivative_identity'] = "d(psi)/dx = u^c N /(5 W^2); d(psi)/dz = u^{2c+1} N/(5 W^2)"

# Theta = 5 m W (psi - psi(0)) exactly.
assert sp.simplify(sp.expand(Theta_ - 5*m*W_*(psi - sp.Rational(21, 5)/m))) == 0
out['Theta_identity'] = "Theta = 5 m W (psi - lambda0), lambda0 = 21/(5m)"

# -f - lambda0 g = g (psi - lambda0), and g has the sign of x, W>0, m>0:
# the certificate integrand is positive off x=0 iff Theta has the sign of x.
out['certificate_integrand'] = "-f - lambda0*g = g*(psi-lambda0); sign(g)=sign(x)"

# Theta is nondecreasing in m at fixed (c,x>=0): dTheta/dm = (21+dx)u^{c+1} - 21 u^2.
assert sp.simplify(sp.expand(sp.diff(Theta_, m) - ((21 + d_*x)*(1+x)**(c+1) - 21*(1+x)**2))) == 0
out['Theta_m_monotone'] = ("dTheta/dm = (21+dx)u^{c+1} - 21 u^2 = u^2[(21+dx)u^{c-1} - 21] "
                           ">= 0 for x>=0, c>=1, d>=0, strict for x>0")

# ------------------------------------------------- closed-form rational gate
# (1+x)^{c-1} >= (1+x)/(1+(2-c)x) for 1<=c<=2, x>=0  (Bernoulli with exponent 2-c in [0,1]).
# Hence Theta >= Theta3/(1+(2-c)x) with Theta3 = m(21+dx)(1+x)^3 - 21 W (1+(2-c)x).
Theta3 = sp.expand(m*(21 + d_*x)*(1+x)**3 - 21*W_*(1 + (2-c)*x))
P3 = sp.Poly(Theta3, x)
assert P3.nth(0) == 0
a1 = sp.factor(P3.nth(1)); a2 = sp.factor(P3.nth(2))
a3 = sp.factor(P3.nth(3)); a4 = sp.factor(P3.nth(4))
# a1 = m e - 210 = 5K.
assert sp.simplify(a1 - (m*e_ - 210)) == 0
assert sp.simplify((a1 - 5*K).subs(m, m_)) == 0
thresholds = {
    'x^2': sp.nsimplify(sp.solve(a2, m)[0]),
    'x^3': sp.nsimplify(sp.solve(a3, m)[0]),
    'x^4': sp.nsimplify(sp.solve(a4, m)[0]),
}
out['Theta3_coefficients'] = {'x': str(a1), 'x^2': str(a2), 'x^3': str(a3), 'x^4': str(a4)}
out['Theta3_m_thresholds'] = {k: str(sp.simplify(v)) for k, v in thresholds.items()}
out['Theta3_criterion'] = ("K>0 and m >= max over the three thresholds implies Theta3>0 on x>0, "
                           "hence Theta>0 on x>0")
tab = []
for cv in ['1', '11/10', '6/5', '13/10', '7/5', '3/2', '31/20']:
    cq = sp.Rational(cv)
    mm = max(sp.simplify(v.subs(c, cq)) for v in thresholds.values())
    KK = sp.simplify(mm*(11*cq - 5)/5 - 42)
    tab.append({'c': cv, 'm_sufficient': str(mm), 'K_sufficient': str(KK),
                'K_sufficient_decimal': str(sp.N(KK, 8))})
out['Theta3_sufficient_points'] = tab

# ------------------------------------------- exact decision at rational c
def theta_positive_on_positive_x(cq, Kq):
    """Exact: is Theta(x) > 0 for every x > 0?  Sturm root count, no floats."""
    cq, Kq = sp.Rational(cq), sp.Rational(Kq)
    q = sp.denom(cq + 1)                      # u = t^q makes u^{c+1} polynomial
    d0, e0 = 16 - 10*cq, 11*cq - 5
    m0 = sp.Rational(5, 1)*(Kq + 42)/e0
    xx = t**q - 1
    WW = m0 + (2*m0 + 10)*xx + (m0 + sp.Rational(111, 5))*xx**2 + (sp.Rational(61, 5) - cq)*xx**3
    P = sp.Poly(sp.expand(m0*(21 + d0*xx)*t**(q*(cq + 1)) - 21*WW), t)
    Q = sp.Poly(sp.expand(P.as_expr().subs(t, 1 + s)), s)   # x>0 <=> t>1 <=> s>0
    k = min(mo[0] for mo in Q.monoms())
    if k:
        Q = Q.quo(sp.Poly(s**k, s))
    if Q.eval(0) <= 0 or sp.LC(Q) <= 0:
        return False
    return sp.polys.polytools.count_roots(Q, 0, sp.oo) == 0

# c=1 is exactly solvable: u^{c+1}=u^2 makes Theta a cubic in x.
Th1 = sp.expand(Theta_.subs(c, 1))
P1 = sp.Poly(Th1, x)
assert P1.nth(0) == 0
assert sp.simplify(P1.nth(1) - (6*m - 210)) == 0
assert sp.simplify(P1.nth(2) - (12*m - sp.Rational(2331, 5))) == 0
assert sp.simplify(P1.nth(3) - (6*m - sp.Rational(1176, 5))) == 0
# cubic coefficient vanishes at m=196/5, i.e. K=126/25; the other two stay positive.
assert sp.simplify((6*m - sp.Rational(1176, 5)).subs(m, sp.Rational(196, 5))) == 0
assert sp.simplify(m_.subs([(c, 1), (K, sp.Rational(126, 25))]) - sp.Rational(196, 5)) == 0
assert theta_positive_on_positive_x(1, sp.Rational(126, 25))
assert not theta_positive_on_positive_x(1, sp.Rational(126, 25) - sp.Rational(1, 1000))
out['c_equal_1_exact_threshold'] = 'K* (c=1) = 126/25 exactly; Theta>0 on x>0 iff K >= 126/25'

def bracket(cq, hi, steps=14):
    lo, hi = sp.Rational(0), sp.Rational(hi)
    assert theta_positive_on_positive_x(cq, hi)
    for _ in range(steps):
        mid = sp.nsimplify((lo + hi)/2)
        if theta_positive_on_positive_x(cq, mid):
            hi = mid
        else:
            lo = mid
    return lo, hi

rows = []
for cv, hi in [('1', 10), ('11/10', 10), ('6/5', 15), ('13/10', 20), ('7/5', 20), ('3/2', 50)]:
    lo, up = bracket(sp.Rational(cv), hi)
    rows.append({'c': cv, 'K_star_lower_open': str(lo), 'K_star_upper_certified': str(up),
                 'decimal': f'({sp.N(lo, 8)}, {sp.N(up, 8)}]'})
out['exact_K_star_brackets'] = rows
# c = 8/5 has d = 0 and leading behaviour m u^{c+1} vs (61/5-c)x^3: never certified.
assert not theta_positive_on_positive_x(sp.Rational(8, 5), 10**6)
out['c_equal_8_5'] = 'd=0; Theta<0 for large x at every K>0 (checked exactly at K=10^6)'

Path(__file__).with_suffix('.json').write_text(json.dumps(out, indent=2) + '\n')
print(json.dumps(out, indent=2))
