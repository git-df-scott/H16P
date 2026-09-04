#!/usr/bin/env python3
"""Exact Strike-5 Green identities; no parameter evaluation or search."""
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
    PB,PC,ZB,ZC,B,C,omega,Rprime = S.symbols(
        "PB PC ZB ZC B C omega Rprime")
    determinant = PB*ZC-ZB*PC
    derivatives = {PB:-omega*B, PC:-omega*C,
                   ZB:Rprime*PB, ZC:Rprime*PC}
    differentiate = lambda expression: sum(
        S.diff(expression,x)*dx for x,dx in derivatives.items())
    assert S.expand(differentiate(determinant)
                    -omega*(C*ZB-B*ZC)) == 0
    theta = PB/(PB-PC)
    height = (1-theta)*ZB+theta*ZC
    assert S.cancel(height-determinant/(PB-PC)) == 0
    assert S.cancel(differentiate(theta)
                    -omega*(B*PC-PB*C)/(PB-PC)**2) == 0
    print("EXACT PASS: determinant derivative, critical mixture height, mixture-coordinate derivative.")

    X,Z,P,u,j,R,udot,Rdot,Ut = S.symbols("X Z P u j R udot Rdot Ut")
    height_derivatives = {X:udot*Z, Z:Rdot*P, u:udot, j:R*udot, R:Rdot}
    height_diff = lambda expression: sum(
        S.diff(expression,x)*dx for x,dx in height_derivatives.items())
    U = u-j/R
    original_height = X-j*Z/R
    assert S.cancel(height_diff(U)-j*Rdot/R**2) == 0
    assert S.cancel(height_diff(original_height)-j*Rdot/R**2*(Z-R*P)) == 0
    kernel = R*(Ut-u)+j
    assert S.cancel(kernel-R*(Ut-U)) == 0
    assert S.cancel(height_diff(kernel)-Rdot*(Ut-u)) == 0
    print("EXACT PASS: original-height derivative, positive-kernel identity, kernel derivative.")
    print("Residual endpoint determinant sign remains OPEN; no parameter signs evaluated.")
    print("CPU seconds:",round(time.process_time()-started,6))


if __name__ == "__main__":
    check()
