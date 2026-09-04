#!/usr/bin/env python3
"""Tiny exact check of the NEW Theorem-N kernel transformation/residual.

No shooting, tangency, numerical parameter scan, or old proof replay.
"""
import os
for name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
             "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[name] = "1"
os.nice(10)
import resource
resource.setrlimit(resource.RLIMIT_CPU, (10, 10))
import time
import sympy as S


def check():
    started = time.process_time()
    # u=1-t; use u>0 to make fractional-power identities unambiguous.
    u, a = S.symbols("u a", positive=True)
    h = 1-a+a*u
    v = S.Function("v")(u)
    y = h**S.Rational(3,2)*v
    original = h*u*S.diff(y,u,2)+(1-a)*S.diff(y,u)/2+5*a*y/36
    transformed = h*u*S.diff(v,u,2)+(1-a+6*a*u)*S.diff(v,u)/2+8*a*v/9
    assert S.simplify(original-h**S.Rational(3,2)*transformed) == 0
    limiting = S.Rational(3,2)*(u**(-S.Rational(4,3))-u**(-S.Rational(2,3)))
    assert limiting.subs(u,1) == 0
    assert -S.diff(limiting,u).subs(u,1) == 1
    residual = (h*u*S.diff(limiting,u,2)
                +(1-a+6*a*u)*S.diff(limiting,u)/2+8*a*limiting/9)
    expected = (1-a)*(22-7*u**S.Rational(2,3))/(6*u**S.Rational(7,3))
    assert S.simplify(residual-expected) == 0
    assert S.simplify(residual.subs(a,1)) == 0
    print("EXACT PASS: homogeneous conjugation, center normalization, limiting solution, positive residual.")
    print("Residual: (1-a)*(22-7*u^(2/3))/(6*u^(7/3)); positive for 0<a<1,0<u<1.")
    print("CPU seconds:", round(time.process_time()-started,6))


if __name__ == "__main__":
    check()
