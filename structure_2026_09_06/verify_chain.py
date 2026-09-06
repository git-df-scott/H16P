#!/usr/bin/env python3
"""Is my V_8 (reduced modulo eta_1, eta_2) proportional to the repository's
eta_3 on the order-3 stratum?  A constant ratio means the identity
N = 640*eta_3 is independent of sign convention."""
import subprocess, sympy as S
l, m, a, b = S.symbols('l m a b', real=True)
loc = {'l': l, 'm': m, 'a': a, 'b': b}

out = subprocess.run(['python3', 'lyapunov.py'], capture_output=True, text=True).stdout
V8 = S.sympify([ln for ln in out.split('\n') if ln.startswith('V_8 =')][0][6:], locals=loc)
V6 = S.sympify([ln for ln in out.split('\n') if ln.startswith('V_6 =')][0][6:], locals=loc)
V4 = S.sympify([ln for ln in out.split('\n') if ln.startswith('V_4 =')][0][6:], locals=loc)
print("free symbols of V_8:", V8.free_symbols)

sub3 = {m: 5*a, b: 3*l + 5}
print("\nOn the order-3 stratum m=5a, b=3l+5:")
print("  V_4 ->", S.simplify(V4.subs(sub3)))
print("  V_6 ->", S.simplify(V6.subs(sub3)))
V8s = S.factor(S.simplify(V8.subs(sub3)))
print("  V_8 ->", V8s)
eta3 = S.factor(-25*a*(2*a**2+l+2)*(5*a**2*l+6*a**2-3*l**3-12*l**2-15*l-6)/64)
print("\nrepository eta_3 =", eta3)
ratio = S.simplify(V8s/eta3)
print("\nratio V_8/eta_3 =", ratio, "   constant?", ratio.free_symbols == set())
