#!/usr/bin/env python3
"""Exact rational checks for the bounded K1 theory review; no orbit evaluations.

Run: python staged_2026_09_05/theory_exact_checks.py
The Bernstein enclosure is a polynomial identity on a closed unit cube.
"""
import itertools
import json
import math
from pathlib import Path
import sympy as s


def bernstein_coefficients(expr, variables):
    p = s.Poly(s.expand(expr), *variables)
    degrees, monomials = p.degree_list(), p.as_dict()
    coefficients = []
    for index in itertools.product(*(range(n + 1) for n in degrees)):
        value = sum(
            a * s.prod(s.Rational(math.comb(i, j), math.comb(n, j))
                       for i, j, n in zip(index, powers, degrees))
            for powers, a in monomials.items()
            if all(j <= i for i, j in zip(index, powers))
        )
        coefficients.append(value)
    return {"degrees": list(degrees), "count": len(coefficients),
            "min": str(min(coefficients)), "max": str(max(coefficients))}, coefficients


def main():
    x, c, K, C, T, X = s.symbols("x c K C T X")
    m = 5 * (K + 42) / (11 * c - 5)
    u, d, h = 1 + x, 16 - 10 * c, s.Rational(61, 5) - c
    W = m + (2*m + 10)*x + (m + s.Rational(111, 5))*x**2 + h*x**3
    N = s.cancel((d*u + (c + 1)*(21 + d*x))*W - u*(21 + d*x)*s.diff(W, x))
    assert s.cancel(N.subs(x, 0) - 5*K) == 0
    assert s.cancel(N.subs(x, -1) - 5*c*(c + 1)*(2*c + 1)) == 0
    assert s.cancel(s.Poly(N, x).coeff_monomial(x**4) - d*(c-1)*h) == 0
    P = s.Poly(s.cancel(5*(11*c-5)*N), x)
    substitutions = {c: s.Rational(9, 10) + C/10, K: s.Rational(6, 5)*T}
    output = {"status": "EXACT_CHECKS_PASS", "orbit_evaluations": 0,
              "domain": {"c": "[9/10,1]", "K": "[0,6/5]"},
              "positive_denominator": "5*(11*c-5)", "coefficients": {}}
    for power in (1, 2, 3):
        bounds, values = bernstein_coefficients(
            P.coeff_monomial(x**power).subs(substitutions), (C, T))
        assert max(values) < 0
        output["coefficients"][str(power)] = bounds
    bounds, values = bernstein_coefficients(
        P.as_expr().subs({**substitutions, x: -X}), (C, T, X))
    assert min(values) >= 0
    output["left_interval_N"] = bounds
    # Generic smooth-planar control: weak focus and three simple cycles,
    # without any double cycle. This field has degree 9, NOT degree 2.
    r = s.symbols("r")
    radial = r**3*(1-r**2)*(4-r**2)*(9-r**2)
    derivatives = [s.diff(radial, r).subs(r, a) for a in (1, 2, 3)]
    assert derivatives == [-48, 480, -6480]
    output["nonquadratic_logic_control"] = {
        "degree": 9, "radii": [1, 2, 3],
        "radial_derivatives": list(map(str, derivatives)),
        "focus_cubic_coefficient": str(s.expand(radial).coeff(r, 3)),
        "claim": "Trace zero alone does not require a double cycle."
    }
    destination = Path(__file__).with_name("theory_exact_checks.json")
    destination.write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
