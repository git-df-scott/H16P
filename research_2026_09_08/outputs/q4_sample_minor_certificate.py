"""Exact sign check for minors of frozen decimal tables, NOT exact integrals."""
from pathlib import Path
from itertools import combinations
from fractions import Fraction
from math import lcm
import json,hashlib
P=Path(__file__).resolve().parent

def integer_row(row):
    q=[Fraction(str(x)) for x in row];scale=lcm(*(a.denominator for a in q))
    return [a.numerator*(scale//a.denominator) for a in q]

def det(a):
    a=[row[:] for row in a];previous=1;sgn=1;n=len(a)
    for k in range(n-1):
        if not a[k][k]:
            pivot=next((j for j in range(k+1,n) if a[j][k]),None)
            if pivot is None:return 0
            a[k],a[pivot]=a[pivot],a[k];sgn=-sgn
        for i in range(k+1,n):
            for j in range(k+1,n):
                numerator=a[i][j]*a[k][k]-a[i][k]*a[k][j]
                assert numerator%previous==0
                a[i][j]=numerator//previous
        previous=a[k][k]
        for i in range(k+1,n):a[i][k]=0
    return sgn*a[-1][-1]

assert det([[1,2],[3,4]])==-2
assert det([[0,1],[2,3]])==-2
assert det([[1,2],[2,4]])==0
records=json.loads((P/'q4_modulus_boundary_constraints.json').read_text())['records']
refined={r['modulus']:r for r in json.loads((P/'q4_modulus_constraint_refine.json').read_text())['records']}
out={'scope':'EXACT for rational interpretations of frozen decimal table entries ONLY. These numbers are approximations to true integrals; no continuum or rigorous quadrature claim.', 'records':[]}
for record in records:
    r=record['modulus'];run=refined.get(r,record)['runs'][-1];rows=run['rows'];ints=[integer_row(row['M']) for row in rows]
    signs=[]
    for inds in combinations(range(len(rows)),4):
        value=det([ints[i] for i in inds]);signs.append((value>0)-(value<0))
    counts={str(s):signs.count(s) for s in [-1,0,1]}
    entry={'modulus':r,'order':run['order'],'minor_count':len(signs),'sign_counts':counts,'uniform_nonzero_orientation':len(set(signs))==1 and signs[0]!=0,'table_sha256':hashlib.sha256(json.dumps([row['M'] for row in rows]).encode()).hexdigest()}
    out['records'].append(entry);print(json.dumps(entry),flush=True)
(P/'q4_sample_minor_certificate.json').write_text(json.dumps(out,indent=2)+'\n')
