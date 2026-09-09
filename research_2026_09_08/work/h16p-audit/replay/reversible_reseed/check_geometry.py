"""Exact algebra for the reversible two-center re-seed. No ODE sampling."""
from pathlib import Path
import json
import sympy as S

ROOT = Path(__file__).resolve().parent
x, y, a, b = S.symbols('x y a b', real=True)
X, Y, k, p, q = S.symbols('X Y k p q', real=True)
P = (b-2)/4 + (1-b)*y + a*x*x + b*y*y
Q = -2*x*y
V = S.Matrix([P,Q])
variables = S.Matrix([x,y])
checks = {}
def zero(name, expression):
    result = S.factor(expression)
    assert result == 0, (name, result)
    checks[name] = True

# First scale (X,Y) by k. Then x=kY/2, y=(1+kX)/2.
oldP, oldQ = -Y*(1+k*X), X+p*X*X+q*Y*Y
subs = {X:(2*y-1)/k, Y:2*x/k}
zero('affine_P', k*oldQ.subs(subs)/2-P.subs({a:2*q/k,b:2*p/k}))
zero('affine_Q', k*oldP.subs(subs)/2-Q)
J = V.jacobian(variables)
checks['upper_J'] = str(J.subs({x:0,y:S.Rational(1,2)}))
checks['lower_J'] = str(S.simplify(J.subs({x:0,y:(b-2)/(2*b)})))
checks['upper_det'] = str(S.factor(J.det().subs({x:0,y:S.Rational(1,2)})))
checks['lower_det'] = str(S.factor(J.det().subs({x:0,y:(b-2)/(2*b)})))
checks['off_axis_x_squared'] = str((2-b)/(4*a))

# H = y^a (x^2 - R(y)), y>0; special a=-1,-2 use logarithms.
R = -b*y*y/(a+2) - (1-b)*y/(a+1) + (2-b)/(4*a)
zero('orbit_equation', y*S.diff(R,y)+a*R+((b-2)/4+(1-b)*y+b*y*y))
zero('H_derivative', (a*(x*x-R)-y*S.diff(R,y))*Q + 2*x*y*P)
zero('integrating_factor', y*S.diff(P,x)+y*S.diff(Q,y)+(a-1)*Q)
e0,e1,e2 = S.symbols('epsilon0 epsilon1 epsilon2')
zero('normal_weighted_divergence',
     y*S.diff(e1*x+e2*x*y,x)+y*S.diff(e0,y)+(a-1)*e0
     -(e1*y+e2*y*y+(a-1)*e0))

# General quadratic perturbations modulo affine changes and time rescaling.
monomials=[1,x,y,x*x,x*y,y*y]
def coeffs(W):
    return S.Matrix([S.Poly(S.expand(f),x,y).coeff_monomial(m) for f in W for m in monomials])
columns=[]
for g in [S.Matrix([1,0]),S.Matrix([0,1]),S.Matrix([x,0]),S.Matrix([y,0]),S.Matrix([0,x]),S.Matrix([0,y])]:
    columns.append(coeffs(J*g-g.jacobian(variables)*V))
columns.append(coeffs(V))
columns += [coeffs(V.diff(a)),coeffs(V.diff(b)),coeffs(S.Matrix([0,1])),coeffs(S.Matrix([x,0])),coeffs(S.Matrix([x*y,0]))]
chartdet = S.factor(S.Matrix.hstack(*columns).det())
checks['full_chart_determinant'] = str(chartdet)
gamma_column=coeffs(S.Matrix([0,x*x]))
swapped=columns[:]
swapped[-1]=gamma_column
checks['boundary_replacement_determinant'] = str(S.factor(S.Matrix.hstack(*swapped).subs(a,-2).det()))
assert checks['boundary_replacement_determinant']=='-48*b'

# Flux identity eliminates x^2/y at a != -2, but loses that term at -2.
# div(y^(a-2) x V) / y^(a-2) = (a+2)x^2 + f(y).
zero('area_flux_identity',
    x*(S.diff(P,x)+S.diff(Q,y)+(a-2)*Q/y)+P
    -((a+2)*x*x+(b-2)/4+(1-b)*y+b*y*y))
checks['boundary_invisible_old_direction'] = {
    'epsilon0':'(2-b)/(12*b)', 'epsilon1':'(1-b)/b', 'epsilon2':'1'}

# At a=-2, b=1 the reciprocal chart is a logarithmic Hamiltonian.
U,W=S.symbols('U W',positive=True)
Hlog=x*x/y**2+S.log(y)+1/(8*y*y)
zero('boundary_log_first_integral',S.diff(Hlog,x)*P.subs({a:-2,b:1})+S.diff(Hlog,y)*Q)
zero('reciprocal_U_equation',((P.subs({a:-2,b:1})*y-x*Q)/y**2).subs({x:U/W,y:1/W})-(1/W-W/4))
zero('reciprocal_W_equation',(-Q/y**2).subs({x:U/W,y:1/W})-2*U)

u,v = S.symbols('u v')
up = -u*(a+b*v*v+(1-b)*u*v+(b-2)*u*u/4)
vp = -v*(a+2+b*v*v)-(1-b)*u*v*v-(b-2)*u*u*v/4
checks['infinity_horizontal_eigenvalues'] = [str(-a),str(-a-2)]
zero('bicycle_radial_eigenvalue',(-a-b*v*v-2).subs(v*v,-(a+2)/b))
zero('bicycle_angular_eigenvalue',(-a-2-3*b*v*v-2*(a+2)).subs(v*v,-(a+2)/b))

# Exact published Yu-Zeng shape, with reversal of time in their chart.
a1,a4 = -S.Rational(30,7),-S.Rational(671,210)
ayz,byz=-2*a4/a1,-2/a1
checks['yu_zeng_shape']={'a':str(ayz),'b':str(byz),'lambda':str(-(ayz+2)/ayz)}
assert -2 < ayz < 0 and 0 < byz < 2 and ayz != -1

# A concrete genuinely different two-center seed, with bicycle boundaries.
seed={a:S.Rational(-7,3),b:S.Integer(1)}
checks['bicycle_seed']={str(z):str(w) for z,w in seed.items()}
checks['bicycle_seed_R'] = str(S.factor(R.subs(seed)))
checks['bicycle_seed_h_center'] = str(S.factor((-R*y**a).subs(seed).subs(y,S.Rational(1,2))))
checks['versions']={'sympy':S.__version__}
out=ROOT/'data'/'geometry.json'
out.parent.mkdir(parents=True,exist_ok=True)
out.write_text(json.dumps(checks,indent=2)+'\n')
print(json.dumps(checks,indent=2))
