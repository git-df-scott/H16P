"""Exact rational boundary integrands and nonvalidated high-precision values.

Integrates the first energy-variation form along the exact finite connection.
It does not establish a full parameter-dependent Dulac displacement theorem.
"""
from pathlib import Path
import json
import sympy as S
import mpmath as mp

HERE=Path(__file__).resolve().parent
s=S.symbols('s');g=2*s**3-3*s-2;k=2*s*s-1
X=(3+8*s+4*s**4)/(12*(1-2*s*s))
Y=(1-8*s-24*s*s-16*s**3-4*s**4)/(12*(1-2*s*s))
P=-Y-3*X*X+2*X*Y+Y*Y;Q=X+X*X-4*X*Y-Y*Y
entries=[]
for name,dp,dq in [('tau',X,Y),('u',2*X*Y,X*X-Y*Y),('v',0,X*Y),('w',X*Y,0)]:
    f=S.factor(S.cancel(36/g*(3*k/(S.sqrt(2)*g))**5*(Q*dp-P*dq)))
    numerator,denominator=S.fraction(f)
    assert S.rem(S.Poly(numerator,s,extension=S.sqrt(2)),S.Poly(k,s,extension=S.sqrt(2))).is_zero
    fn=S.lambdify(s,f,'mpmath')
    vals=[]
    for precision in [60,90]:
        with mp.workdps(precision):
            end=1/mp.sqrt(2)
            vals.append(-mp.quad(fn,[-end,-end/2,0,end/2,end]))
    with mp.workdps(90):
        change=abs(vals[1]-vals[0])
        assert change<mp.mpf('1e-55')
        entries.append({'control':name,'exact_integrand':str(f),'value':mp.nstr(vals[1],85),'precision_change':mp.nstr(change,6)})
with mp.workdps(90):
    direction=json.loads((HERE/'q4_original_probe.json').read_text())['tighter_direction']
    combined=sum(mp.mpf(e['value'])*mp.mpf(str(d)) for e,d in zip(entries,direction))
out={'status':'exact integrands; nonvalidated quadrature','orientation':'integrate s=+1/sqrt2 to -1/sqrt2',
     'entries':entries,'three_anchor_direction_boundary_integral':mp.nstr(combined,40),
     'scope':'The numerical three-anchor direction has nonzero boundary energy integral. No full splitting remainder, two-cycle boundary construction, or interval enclosure is supplied.',
     'ode_evaluations':0}
(HERE/'q4_boundary_integral.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
