"""Exact amplitude gate in the residual window {0 < K < K*(c)}.

If Theta > 0 on (0,X] then every origin-surrounding closed orbit whose
x-range is contained in (-1,X] has log M > 0, hence is hyperbolic repelling.
Contrapositive: every stable, semistable or multiplier-one origin cycle has
x_max > X.  Zero ODE calls; all decisions by exact Sturm root counting.
"""
import json
from pathlib import Path
import sympy as sp

s = sp.symbols('s')

def theta_shifted(cq, Kq):
    """Q(s) with s = u^{1/q}-1 >= 0; sign of Q on s>0 is the sign of Theta on x>0."""
    cq, Kq = sp.Rational(cq), sp.Rational(Kq)
    q = sp.denom(cq + 1)
    d0, e0 = 16 - 10*cq, 11*cq - 5
    m0 = 5*(Kq + 42)/e0
    xx = (1 + s)**q - 1
    WW = m0 + (2*m0 + 10)*xx + (m0 + sp.Rational(111, 5))*xx**2 + (sp.Rational(61, 5) - cq)*xx**3
    Q = sp.Poly(sp.expand(m0*(21 + d0*xx)*(1 + s)**(q*(cq + 1)) - 21*WW), s)
    k = min(mo[0] for mo in Q.monoms())
    return (Q.quo(sp.Poly(s**k, s)) if k else Q), q

def positive_upto(Q, sigma):
    """Exact: Q > 0 on (0, sigma]?"""
    if Q.eval(0) <= 0 or Q.eval(sp.Rational(sigma)) <= 0:
        return False
    return sp.polys.polytools.count_roots(Q, 0, sp.Rational(sigma)) == 0

def gate(cq, Kq, sigma_max=sp.Rational(64), steps=18):
    """Largest dyadic sigma with Theta>0 on (0,(1+sigma)^q-1]; returns that x bound."""
    Q, q = theta_shifted(cq, Kq)
    if positive_upto(Q, sigma_max):
        return None, q            # certified (no cycle) or gate beyond the probe range
    lo, hi = sp.Rational(0), sigma_max
    for _ in range(steps):
        mid = sp.nsimplify((lo + hi)/2)
        if positive_upto(Q, mid):
            lo = mid
        else:
            hi = mid
    return (1 + lo)**q - 1, q

rows = []
for cv, Kv in [('1', '1'), ('1', '3'), ('1', '5'),
               ('11/10', '1'), ('11/10', '3'), ('11/10', '6'),
               ('6/5', '1'), ('6/5', '4'), ('6/5', '8'),
               ('13/10', '2'), ('13/10', '10'),
               ('7/5', '2'), ('7/5', '13'),
               ('3/2', '2'), ('3/2', '7'), ('3/2', '21'),
               ('8/5', '7')]:
    X, q = gate(sp.Rational(cv), sp.Rational(Kv))
    rows.append({'c': cv, 'K': Kv, 'q': int(q),
                 'x_max_strict_lower_bound': (str(X) if X is not None else 'certified/none'),
                 'decimal': (str(sp.N(X, 10)) if X is not None else '-')})
    print(rows[-1])

out = {'arithmetic': 'exact rational SymPy Sturm counts; zero ODE calls',
       'statement': ('Theta>0 on (0,X] implies every closed orbit with x-range in (-1,X] '
                     'is hyperbolic repelling; hence any stable, semistable or '
                     'multiplier-one origin cycle has x_max > X.'),
       'rows': rows}
Path(__file__).with_suffix('.json').write_text(json.dumps(out, indent=2) + '\n')
