"""Independent exact elimination and reversibility certificates. No Fable import."""
import json
from pathlib import Path
import sympy as s
x,y,a,b,l,m=s.symbols('x y a b l m',real=True)
V=s.Matrix([-y+l*x*x+m*x*y+y*y,x*(1+a*x+b*y)])
T=s.diff(V[0],x)+s.diff(V[1],y)
M=a*(b+2*l)/(l+1); W=V.subs(m,M)
H=a*a-b*(l+1); C=a*a*(b+2*l+1)-(b+1)*(l+1)**2
# Solve line and trace simultaneously, rather than substituting into Fable's calculation.
sol=s.solve([1+a*x+b*y,(l+1)*x+a*y],[x,y]);assert s.factor(sol[x]+a/H)==0;assert s.factor(sol[y]-(l+1)/H)==0
assert s.factor(W[0].subs(sol)+C/H**2)==0
assert s.factor(T.subs(m,M).subs(y,-(1+a*x)/b)+(b+2*l)*(H*x+a)/(b*(l+1)))==0
assert T.subs({x:0,y:1})==m
# Reproduce the UNCANCELLED numerator, and distinguish it from P itself.
uncancelled=s.expand((b*b*(l+1)*W[0]).subs(y,-(1+a*x)/b))
assert s.factor(uncancelled.subs(x,-a/H)*H**2+b*b*(l+1)*C)==0

def reflection_residual(field,d):
 R=s.eye(2)-2*d*d.T/d.dot(d);z=R*s.Matrix([x,y])
 assert s.simplify(R*R)==s.eye(2)
 return [s.factor(v) for v in field.subs({x:z[0],y:z[1]},simultaneous=True)+R*field]
def divisible(expr,poly):
 n=s.fraction(s.cancel(expr))[0]
 return s.rem(s.Poly(n,x,y,a,b,l,m),s.Poly(poly,x,y,a,b,l,m)).as_expr()==0
res=reflection_residual(W,s.Matrix([l+1,a]));assert all(divisible(v,C) for v in res)
assert all(v==0 for v in reflection_residual(V.subs({a:0,m:0}),s.Matrix([1,0])))
assert s.expand(T.subs({b:-2*l,m:0}))==0
# b=0, a*l != 0: the only zero-trace candidate on x=-1/a.
b0={b:0,m:2*a*l/(l+1),x:-1/a,y:(l+1)/a**2}
assert s.factor(V[0].subs(b0)+(a*a*(2*l+1)-(l+1)**2)/a**4)==0
# l=-1: eta1=0 is a(b-2)=0; NEVER divide by l+1.
N=(b+1)*(b-2)**2-(b-1)*m*m
D=V.subs({l:-1,a:0})
assert s.factor(D[0].subs({x:m/(b*(b-2)),y:-1/b})-N/(b*b*(b-2)**2))==0
resdeg=reflection_residual(D,s.Matrix([b-2,m]));assert all(divisible(v,N) for v in resdeg)
assert s.expand(T.subs({l:-1,b:2,m:0}))==0
# b=2, m != 0: trace=m*y forces y=0; P=-x^2 forces x=0.
assert V[0].subs({l:-1,b:2,y:0})==-x*x
# Counterexample to BOTH universal and existential center -> neutral converses.
example={a:1,b:1,l:0,m:1}; E=s.solve(list(V.subs(example)),[x,y]);assert set(E)=={(0,0),(0,1),(-1,0)}
traces=[str(T.subs(example).subs({x:u,y:v})) for u,v in E]
assert all(v==0 for v in reflection_residual(V.subs(example),s.Matrix([1,1])))
out={'verdict_literal_proposition':'GAP: converse false; missing H=0 case',
 'corrected_implication':'VERIFIED for real Shi chart, including b=0 and l=-1',
 'line_trace':str(s.factor(T.subs(m,M).subs(y,-(1+a*x)/b))),
 'candidate':{str(k):str(s.factor(v)) for k,v in sol.items()},
 'P_at_candidate':str(-C/H**2),'uncancelled_numerator':str(-b*b*(l+1)*C),
 'generic_reflection_residual':[str(v) for v in res],
 'degenerate_l_minus_one_polynomial':str(s.expand(N)),
 'degenerate_reflection_residual':[str(v) for v in resdeg],
 'counterexample':{'vector':['0','0','-1','0','1','1','0','1','0','1','1','0'],'equilibria':[list(map(str,e)) for e in E],'traces':traces}}
p=Path(__file__).with_name('proposition_a_exact.json');p.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
