"""Exact council gates only: no orbit integration or parameter search."""
import os
import resource

for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[name] = "1"
resource.setrlimit(resource.RLIMIT_CPU, (10, 10))

import sympy as s

x, y, k, p, q = s.symbols("x y k p q", nonzero=True)
P, Q = -y * (1 + k*x), x + p*x*x + q*y*y
J = s.Matrix([P, Q]).jacobian([x, y])
off = s.simplify(J.subs(x, -1/k))
assert s.simplify(off-s.Matrix([[-k*y, 0], [1-2*p/k, 2*q*y]])) == s.zeros(2)
assert s.simplify(off.det().subs(y*y, (k-p)/(k*k*q))+2*(k-p)/k) == 0
assert s.simplify(J.subs({x: -1/p, y: 0}).det()-(k/p-1)) == 0
# H_X=mu Q, H_Y=-mu P; cross-partial compatibility proves local exactness.
u = 1+k*x
mu = 2*u**(2*q/k-1)
assert s.simplify(s.diff(mu*Q, y) + s.diff(mu*P, x)) == 0
weight = (1+k*x)**(2*q/k-2)
divergence = s.diff(weight*y*P,x)+s.diff(weight*y*Q,y)
assert s.simplify(divergence-weight*(x+p*x*x+(q+k)*y*y)) == 0
exponent, ratio = s.symbols("exponent ratio", nonzero=True)
us = 1-1/ratio
barrier_difference = ratio*us**2/(exponent+2)+(1-2*ratio)*us/(exponent+1)-(1-ratio)/exponent
assert s.factor(barrier_difference-us*(exponent+2*ratio)/(exponent*(exponent+1)*(exponent+2))) == 0
print("PASS reversible equilibrium types, integrating factor, moment relation and loop-energy identity")

c, alpha, beta, z = s.symbols("c alpha beta z")
P, Q = y+x*x+x*y, -10*x*x+s.Rational(11,5)*x*y+c*y*y+alpha*x+beta*y
K = -alpha*(s.Rational(11,5)*c-1)-42
assert K.subs({c:s.Rational(7,10), alpha:-80}) == s.Rational(6,5)
assert K.subs({c:s.Rational(7,10), alpha:-s.Rational(363889,5000)}) == -s.Rational(674997,250000)
T = (c-s.Rational(61,5))*x**3+(alpha-s.Rational(111,5)-beta)*x*x+(2*alpha-10-beta)*x+alpha
assert s.factor(Q.subs(y,-x*x/(1+x))*(1+x)**2-x*T) == 0
assert P.subs(x,-1) == 1
angular = -10+s.Rational(6,5)*z+(c-1)*z*z
assert s.discriminant(angular,z).subs(c,s.Rational(241,250)) == 0
print("PASS KKL start K=6/5, incumbent K=-674997/250000, cubic and infinity gate")

r, tau, u, v, w = s.symbols("r tau u v w")
l1,l2,l3,l4,l5,l6 = s.symbols("l1:7")
focus = s.Matrix([l1,l5*(l3-l6),l2*l4*(l3-l6)*(l4+5*l3-5*l6),l2*l4*(l3-l6)**2*(l3*l6-2*l6*l6-l2*l2)])
family = {l1:tau,l2:r+u,l3:2+r*r,l4:-5*(1+r*r)+v,l5:w,l6:1}
base = {tau:0,u:0,v:0,w:0}
F = focus.subs(family)
normal = s.simplify(F.jacobian([tau,w,v,u]).subs(base))
assert normal == s.diag(1,1+r*r,-5*r*(1+r*r)**2,10*r*r*(1+r*r)**3)
assert normal.subs(r,0).rank() == 2
P = -y-(2+r*r)*x*x+2*r*x*y+y*y
Q = x+r*x*x+(-1-3*r*r)*x*y-r*y*y
P2 = P+y
Q2 = Q-x
angular = s.factor(Q2.subs({x:1,y:z})-z*P2.subs({x:1,y:z}))
assert s.factor(angular+(z+r)*(z*z+2*r*z-1)) == 0
for direction in [-r+s.sqrt(1+r*r),-r-s.sqrt(1+r*r)]:
    assert s.simplify(-P2.subs({x:1,y:direction})-(1+r*r)) == 0
    assert s.simplify(s.diff(angular,z).subs(z,direction)+2*(1+r*r)) == 0
print("PASS generic Q4 four normal directions; rank two at symmetric limit; two infinity saddles")
print("NO numerical returns, candidate, cyclicity bound, or full attack computed")
