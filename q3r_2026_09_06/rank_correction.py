#!/usr/bin/env python3
"""Astra's correction: the four displayed generators are LINEARLY DEPENDENT.

R(y) = h y^{-a} - A - B y - C y^2,  A=(b-2)/(4a), B=(1-b)/(a+1), C=b/(a+2).

Key identity:   y R'(y) + a R(y) = -[ (b-2)/4 + (1-b) y + b y^2 ]
(the constants A,B,C are exactly what makes this hold), hence

  d/dy ( y^{a-1} R^{3/2} ) = -(1/2) y^{a-2} sqrt(R) [ (a+2) R + 3((b-2)/4+(1-b)y+by^2) ]

and integrating over an oval (R vanishes at both turning points) gives

  (a+2) U + 3[ (b-2)/4 T_{a-2} + (1-b) T_{a-1} + b T_a ] = 0.
"""
import sympy as S
y, h, a, b = S.symbols('y h a b', positive=True)
A = (b-2)/(4*a); B = (1-b)/(a+1); C = b/(a+2)
R = h*y**(-a) - A - B*y - C*y**2

lhs = S.simplify(S.expand(y*S.diff(R, y) + a*R))
rhs = -((b-2)/4 + (1-b)*y + b*y**2)
print("y R' + a R           =", lhs)
print("-[(b-2)/4+(1-b)y+by^2] =", S.expand(rhs))
print("difference           =", S.simplify(lhs - rhs))

d = S.simplify(S.diff(y**(a-1)*R**S.Rational(3,2), y))
claim = -S.Rational(1,2)*y**(a-2)*S.sqrt(R)*((a+2)*R + 3*((b-2)/4 + (1-b)*y + b*y**2))
print("\nd/dy(y^(a-1) R^(3/2)) - claimed form =", S.simplify(S.expand(d - claim)))
