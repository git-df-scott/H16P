"""Exact certificates only: no ODE calls and no floating-point sampling."""
import json
from pathlib import Path
import sympy as s

c, K, x, X, C, t = s.symbols('c K x X C t')
u = 1+x
e = 11*c-5
d = 16-10*c
m = 5*(K+42)/e
W = m+(2*m+10)*x+(m+s.Rational(111,5))*x*x+(s.Rational(61,5)-c)*x**3
N = s.cancel((d*u+(c+1)*(21+d*x))*W-u*(21+d*x)*s.diff(W,x))
assert s.factor(N.subs(x,0)-5*K) == 0
assert s.factor(s.diff(N,K)-5*u*u*(c*d*u+5*(2*c+1)*(c-1))/e) == 0

def bernstein(poly, variables):
    p = s.Poly(poly, *variables)
    degrees = [p.degree(v) for v in variables]
    import itertools
    answer = {}
    for index in itertools.product(*(range(n+1) for n in degrees)):
        value = 0
        for powers in itertools.product(*(range(i+1) for i in index)):
            monomial = s.prod(v**j for v,j in zip(variables,powers))
            weight = s.prod(s.binomial(i,j)/s.binomial(n,j) for i,j,n in zip(index,powers,degrees))
            value += p.coeff_monomial(monomial)*weight
        answer[index] = s.factor(value)
    return degrees, answer

left = s.cancel(5*e*N.subs(K,0)).subs({x:-X,c:1+s.Rational(3,5)*C})
degrees, left_certificate = bernstein(left, [X,C])
assert degrees == [4,4]
assert all(v >= 0 for v in left_certificate.values())
assert all(v > 0 for (i,j),v in left_certificate.items() if i > 0)

coefficients = [s.factor(s.Poly(N,x).nth(i)) for i in range(5)]
thresholds = [s.solve(q,K)[0] for q in coefficients[1:4]]
A = 21300*c**4-100925*c**3+160280*c*c-102937*c+22570
B = 2500*c**4-19280*c**3+46119*c*c-53272*c+24095
assert s.factor(thresholds[1]-thresholds[2]-e*A/(50*c*(5*c-8)*(20*c*c-43*c+5))) == 0
assert s.factor(thresholds[1]-thresholds[0]+e*B/(25*(5*c*c-19*c+5)*(20*c*c-43*c+5))) == 0
_, cert_A = bernstein(A.subs(c,1+t/10), [t])
_, cert_minus_B = bernstein((-B).subs(c,s.Rational(11,10)+t/2), [t])
assert min(cert_A.values()) > 0
assert min(cert_minus_B.values()) > 0
for i in [1,2]:
    numerator = s.cancel(s.diff(coefficients[i],K)*5*e)
    _, cert = bernstein(numerator.subs(c,1+s.Rational(3,5)*t), [t])
    assert min(cert.values()) > 0
assert s.factor(s.diff(coefficients[3],K)-10*c*(8-5*c)/e) == 0
assert coefficients[3].subs(c,s.Rational(8,5)) == -s.Rational(2226,25)

# The local first focal numerator in the Lienard chart is exactly 2K.
fz = -s.Rational(21,5)
fzz = -(11+c)/5
gz = m
gzz = 20+m*(1-c)
assert s.factor(fz*gzz-fzz*gz-2*K) == 0

# Lowest radial resonance for an analytic density of divergence weight q.
xx,yy,eta,omega,beta_nf,q,n = s.symbols('xx yy eta omega beta_nf q n')
rho2 = xx*xx+yy*yy
Xnf = -omega*yy+eta*xx*rho2-beta_nf*yy*rho2
Ynf = omega*xx+eta*yy*rho2+beta_nf*xx*rho2
assert s.expand(s.diff(Xnf,xx)+s.diff(Ynf,yy)-4*eta*rho2) == 0
assert s.expand(Xnf*2*xx+Ynf*2*yy-2*eta*rho2*rho2) == 0
assert s.factor(16*c/(2*c+1)-4-4*(2*c-1)/(2*c+1)) == 0

# Verify the quartic generalized Dulac recurrence with arbitrary functions.
z,T,kappa,gamma = s.symbols('z T kappa gamma')
f,h,R = [s.Function(name)(z) for name in ['f','h','R']]
A0,A1,A2 = [s.Function(name)(z) for name in ['A0','A1','A2']]
C3 = s.symbols('C3')
Psi = T**4+C3*T**3+A2*T*T+A1*T+A0
Phi = s.expand((T-h)*s.diff(Psi,z)+(-kappa*f*T-R)*s.diff(Psi,T)+4*kappa*f*Psi)
relations = {s.diff(A2,z):4*R-kappa*f*C3,
             s.diff(A1,z):h*s.diff(A2,z)-2*kappa*f*A2+3*R*C3,
             s.diff(A0,z):h*s.diff(A1,z)-3*kappa*f*A1+2*R*A2}
for degree in range(1,5):
    q = s.Poly(Phi,T).nth(degree)
    for _ in range(3):
        q = q.subs(relations)
    assert s.expand(q) == 0
assert s.factor(s.Poly(Phi,T).nth(0)-(-R*A1-h*s.diff(A0,z)+4*kappa*f*A0)) == 0

output = {
    'arithmetic':'exact rational SymPy; zero ODE calls',
    'N_left_bernstein_degrees':degrees,
    'N_left_bernstein_rows':[[str(left_certificate[i,j]) for j in range(5)] for i in range(5)],
    'N_K_derivative':str(s.factor(s.diff(N,K))),
    'positive_root_coefficient_thresholds':[str(q) for q in thresholds],
    'A_bernstein_c_1_to_11_over_10':[str(v) for v in cert_A.values()],
    'minus_B_bernstein_c_11_over_10_to_8_over_5':[str(v) for v in cert_minus_B.values()],
    'quartic_recurrence_verified':True,
    'first_focal_numerator':'f_z*g_zz-f_zz*g_z = 2K',
    'analytic_density_necessary_order':'n = 4q = 16c/(2c+1)',
    'monic_quartic_one_sign_scalar_certificate_at_fold':'impossible for c>1/2, K!=0, m>0',
    'cycle_count_bound_proved':False,
    'fold_component_excluded':False,
}
Path(__file__).with_suffix('.json').write_text(json.dumps(output,indent=2)+'\n')
print(json.dumps(output,indent=2))
