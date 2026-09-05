"""Ten-CPU-second exact replay; no ODEs or cycle-existence assertion."""
import os
import resource
for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
            'VECLIB_MAXIMUM_THREADS'):
    os.environ[key] = '1'
resource.setrlimit(resource.RLIMIT_CPU, (10, 10))
import sympy as s

x, y, c, alpha, beta, r = s.symbols('x y c alpha beta r')
P = y+x*x+x*y
Q = -10*x*x+s.Rational(11,5)*x*y+c*y*y+alpha*x+beta*y
T = (c-s.Rational(61,5))*r**3+(alpha-s.Rational(111,5)-beta)*r*r+(2*alpha-10-beta)*r+alpha
sigma_y = -r*r/(1+r)
at_section = {x:r, y:sigma_y}
Qs = s.cancel(Q.subs(at_section))
assert s.cancel(P.subs(at_section)) == 0
assert s.cancel(Qs-r*T/(1+r)**2) == 0
Pdot = s.diff(P,x)*P+s.diff(P,y)*Q
assert s.cancel(Pdot.subs(at_section)-(1+r)*Qs) == 0
assert s.cancel(s.diff(sigma_y,r)+r*(r+2)/(1+r)**2) == 0
assert s.cancel(s.diff(sigma_y,r,2)+2/(1+r)**3) == 0
F = s.Matrix([0,Qs]); tangent = s.Matrix([1,s.diff(sigma_y,r)])
assert s.cancel(s.Matrix.hstack(F,tangent).det()+Qs) == 0
assert s.det(s.Matrix([x,P]).jacobian([x,y])) == 1+x
print('PASS: rational section, cubic identity, transversality, sensitivities, determinant ratio')
print('No ODE evaluations or periodic-orbit certificates were produced.')
