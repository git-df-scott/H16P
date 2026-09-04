#!/usr/bin/env python3
"""Exact original-loop identities only; no numerical evaluation or search."""
import os
import resource

resource.setrlimit(resource.RLIMIT_CPU, (10, 10))
for variable in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
                 "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[variable] = "1"
try:
    os.nice(10)
except OSError:
    pass

import sympy as S


def run_checks():
    k = S.symbols("k", positive=True)
    x, m = S.symbols("x m")
    d, c = k-1, 1/S.sqrt(k)
    h = -2*c/3
    denominator = k*m**3-3*d*m+2*d
    y = c+m*x
    hamiltonian = 2*d*x**3/3-d*x*x*y+k*y**3/3-y
    factorization = x*x*(c*(k*m*m-d)+x*denominator/3)
    assert S.simplify(hamiltonian-h-factorization) == 0
    assert S.Matrix([[1, 0], [m, x]]).det() == x
    assert S.expand(denominator-2*d*(1-m)-m*(k*m*m-d)) == 0
    assert S.diff(denominator, m) == 3*k*m*m-3*d

    mu1, mu2, mu3, mu4 = S.symbols("mu1 mu2 mu3 mu4")
    numerator = (-mu1/3+mu3/2-k*mu4*m)*denominator+(mu2+mu3*m)*(d-k*m*m)
    expanded = (
        -k*k*mu4*m**4+k*(-mu1/3-mu3/2)*m**3
        +k*(3*d*mu4-mu2)*m**2
        +d*(mu1-mu3/2-2*k*mu4)*m+d*(mu2-2*mu1/3+mu3)
    )
    assert S.expand(numerator-expanded) == 0

    length = 3*c*(d-k*m*m)/denominator
    i00_density = length**2/2
    i10_density = length**3/3
    i01_density = c*length**2/2+m*length**3/3
    im10_density = length
    im11_density = c*length+m*length**2/2
    original_density = (
        mu1*h*i00_density+mu2*i10_density+mu3*i01_density
        +mu4*(2*im10_density+3*k*h*im11_density)
    )
    proposed_density = 9*c**3*(d-k*m*m)**2*numerator/denominator**3
    assert S.factor(original_density-proposed_density) == 0

    print("Exact saddle-ray factorization and positive-area Jacobian: PASS")
    print("Exact denominator derivative and tangent-endpoint identity: PASS")
    print("Exact affine quartic numerator: PASS")
    print("Exact reconstruction from all four ORIGINAL area basis terms: PASS")
    print("No sign on the residual fibre or endpoint-cycle claim is asserted.")


if __name__ == "__main__":
    run_checks()
