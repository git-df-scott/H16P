"""Exact audit of the leading determinant in Marín–Villadelprat eq.(23).
Set theta=2, remainder=0. A missing bounded factor does not invalidate the
small-positive-variable sign argument. This is not an H16P counterexample.
"""
import sympy as S,json
from pathlib import Path
x,y,z,c4,c5,k,kh=S.symbols('x y z c4 c5 k kh')
M=S.Matrix([[1,x,c5*k*x*x],[1,y,c5*k*y*y],[1,c4*kh*z*z,z]])
actual=S.factor(M.det()/(y-x))
correct=z*(1-c4*c5*kh*z*k*(x+y))+c5*k*x*y
printed=z*(1-c5*kh*z*k*(x+y))+c5*k*x*y
assert S.expand(actual-correct)==0
error=S.factor(actual-printed)
assert error!=0
example={x:S.Rational(1,100),y:S.Rational(2,100),z:S.Rational(3,100),c4:2,c5:3,k:1,kh:1}
out=dict(status='exact determinant identity passed; displayed coefficient mismatch verified',source='https://arxiv.org/html/2501.16924v1#S4.E23',actual=str(actual),displayed=str(printed),actual_minus_displayed=str(error),example={str(a):str(b) for a,b in example.items()},actual_example=str(actual.subs(example)),displayed_example=str(printed.subs(example)),interpretation='The missing nu4 factor remains bounded locally; for positive nu5 and positive coefficients, the original small-variable positivity argument survives. No theorem refutation and no H16P counterexample.')
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
