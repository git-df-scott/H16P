import sympy as s
x,y,l,m,a,b=s.symbols('x y l m a b')

def lyapunov(P2,Q2,maxdeg):
 prev=(x*x+y*y)/2
 vals=[]
 for d in range(3,maxdeg+1):
  cs=s.symbols(f'c0:{d+1}')
  H=sum(cs[k]*x**(d-k)*y**k for k in range(d+1))
  base=s.Poly(s.expand(-y*s.diff(H,x)+x*s.diff(H,y)+P2*s.diff(prev,x)+Q2*s.diff(prev,y)),x,y)
  V=s.Symbol('V') if d%2==0 else s.S.Zero
  eq=s.Poly(base.as_expr()-V*(x*x+y*y)**(d//2),x,y)
  equations=[eq.coeff_monomial(x**(d-k)*y**k) for k in range(d+1)]
  variables=list(cs)
  if d%2==0: equations.append(cs[0]);variables.append(V)
  sol=s.solve(equations,variables)
  prev=s.Poly(s.expand(H.subs(sol)),x,y).as_expr()
  if V: vals.append(s.factor(sol[V])); print(d, vals[-1],flush=True)
 return vals
print('FIRST GENERAL',flush=True)
lyapunov(l*x*x+m*x*y+y*y,a*x*x+b*x*y,4)
print('SECOND L1=0',flush=True)
lyapunov(l*x*x+m*x*y+y*y,a*x*x+((l+1)*m/a-2*l)*x*y,6)
print('THIRD WEAK3',flush=True)
vs=lyapunov(l*x*x+5*a*x*y+y*y,a*x*x+(3*l+5)*x*y,8)
print('SHI', [v.subs({a:1,l:-10}) for v in vs],flush=True)
P=l*x*x+5*a*x*y+y*y-y;Q=x+a*x*x+(3*l+5)*x*y
b0=3*l+5
F=s.factor(P.subs(y,(-1-a*x)/b0)*b0*b0)
print('EXTRA FINITE POLY',F)
print('DISCR',s.factor(s.discriminant(F,x)))
u=s.symbols('u');C=s.expand((a+b0*u)-u*(l+5*a*u+u*u))
print('INFINITY',C,'DISCR',s.factor(s.discriminant(C,u)))
print('RESULTANT',s.factor(s.resultant(C,l+5*a*u+u*u,u)))
print('SHI_INFINITY',s.nroots(C.subs({a:1,l:-10})))
for r in s.nroots(C.subs({a:1,l:-10})):
 if abs(s.im(r))<1e-12:
  print('INF_EIGS',r,s.N(s.diff(C,u).subs({a:1,l:-10,u:r})),s.N(-(l+5*a*u+u*u).subs({a:1,l:-10,u:r})))
