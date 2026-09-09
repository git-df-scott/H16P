"""Exact first nonzero center coefficient of the PF-reconstructed Wronskian."""
from pathlib import Path
import json
import sympy as S
P=Path(__file__).resolve().parent
a,t=S.symbols('a t');R=S.Rational
F=sum(S.rf(R(1,6),j)*S.rf(R(5,6),j)/S.factorial(j)**2*t**j for j in range(6))
M=S.series(1-6*(1-t)*S.diff(F,t)/F,t,0,4).removeO()
qs=[S.Integer(1),t,-M,t*M-1]
y0s=[R(3,1361360)*n for n in [1326,864,-2431,-102]]
cols=[]
for index,(q,y0) in enumerate(zip(qs,y0s)):
    H=S.integrate(S.series(t*F*q,t,0,5).removeO(),(t,0,t))
    ys=[y0,-R(3,2)*(1+a)*y0-(R(1,192) if index==2 else 0)]
    for n in [0,1]:
        z=S.symbols('z');Y=sum(c*t**j for j,c in enumerate(ys))+z*t**(n+2)
        residual=(1-a*t)*(1-t)*S.diff(Y,t,2)-(1-a)*S.diff(Y,t)/2+5*a*Y/36+S.series(H/(1152*t*t*(1-t)),t,0,n+1).removeO()
        value=S.solve(S.expand(residual).coeff(t,n),z)[0];ys.append(S.factor(value))
    Y=sum(c*t**j for j,c in enumerate(ys))
    X=S.integrate(S.series(Y*(1-a*t)**(-R(3,2)),t,0,4).removeO(),(t,0,t))
    cols.append([S.expand(X).coeff(t,j) for j in range(1,5)])
C=S.Matrix.hstack(*(S.Matrix(col) for col in cols));coefficient=S.factor(12*C.det(method='domain-ge'))
assert coefficient!=0
out={'scope':'Exact algebra conditional on inherited PF reconstruction normalization; local center leading term only','basis_order':['A','B','eta','c'],'X_coefficients_t1_to_t4':[[str(S.factor(v)) for v in row] for row in C.tolist()],'Wronskian_X_leading_power':4,'Wronskian_X_leading_coefficient':str(coefficient),'relation':'Original I has common nonzero scalar factor away from center; W_I is that scalar to the fourth power times W_X.'}
(P/'q4_center_wronskian.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
