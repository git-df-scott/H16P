"""Exact affine identification of the known finite-Q4 reversible center locus."""
import json
import sympy as S
from pathlib import Path
r,e,x,y,X,Y=S.symbols('r e x y X Y',real=True)
A=2+r*r;B=-1-3*r*r-e
P=-y-A*x*x+y*y;Q=x*(1+B*y)
sub={x:-2*X/B,y:(2*Y-1)/B}
a=2*A/B;b=-2/B
PX=S.factor((-B*P/2).subs(sub));QY=S.factor((B*Q/2).subs(sub))
assert S.factor(PX-((b-2)/4+(1-b)*Y+a*X*X+b*Y*Y))==0
assert S.factor(QY+2*X*Y)==0
J=S.Matrix([P,Q]).jacobian([x,y]);eq=[]
for pt in [(0,0),(0,1)]:
    jac=J.subs(dict(zip([x,y],pt)));eq.append(dict(point=list(pt),trace=str(S.trace(jac)),determinant=str(S.factor(jac.det()))))
out=dict(status='exact polynomial identities passed',locus='tau=0,w=0,u=-r',condition='B=-1-3r^2-epsilon != 0 for affine map',map={'X':'-B*x/2','Y':'(1+B*y)/2','time':'unchanged'},normal_form={'a':str(S.factor(a)),'b':str(S.factor(b))},equilibria=eq,scope='Reversibility plus ellipticity gives local centers; two elliptic points require epsilon > -3r^2. No global cycle bound, no general classification of all center loci.')
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
