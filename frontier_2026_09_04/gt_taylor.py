"""MPFR noninterval Taylor shooting at the exact Galias–Tucker rational seed.
Scaled coordinates (x,y)=r(X,Y) keep small cycles resolvable. No certified
remainder: repeated precision/order agreement is numerical evidence only.
"""
import gmpy2 as g, json, sys, time
from pathlib import Path
def compute(radius,order=112,bits=900,svalue=None,lambda_value=None):
 g.get_context().precision=bits
 M=g.mpfr; r=M(radius); lam=-M(10)**-200; delta=-M(10)**-13; eps=-M(10)**-52
 if svalue is not None:
  s=M(svalue);delta=-s;eps=-s**4;lam=-M(10)**8*s**16
 if lambda_value is not None:lam=M(lambda_value)
 b=-25+8*eps-9*delta; a=5+delta
 def coeff(X,Y):
  xx=[X]; yy=[Y]
  for n in range(order):
   x2=sum(xx[j]*xx[n-j] for j in range(n+1)); xy=sum(xx[j]*yy[n-j] for j in range(n+1)); y2=sum(yy[j]*yy[n-j] for j in range(n+1))
   xx.append((lam*xx[n]-yy[n]+r*(-10*x2+a*xy+y2))/(n+1))
   yy.append((xx[n]+r*(x2+b*xy))/(n+1))
  return xx,yy
 def poly(c,h):
  v=c[-1]
  for a in c[-2::-1]:v=v*h+a
  return v
 X=M(0);Y=M(1);h=M(1)/10
 for _ in range(62):
  xx,yy=coeff(X,Y); X=poly(xx,h);Y=poly(yy,h)
 xx,yy=coeff(X,Y); dh=[(j+1)*xx[j+1] for j in range(order)]; tau=M('.0832')
 for _ in range(12): tau-=poly(xx,tau)/poly(dh,tau)
 Yr=poly(yy,tau); D=r*(Yr-1)
 return {'radius':radius,'order':order,'bits':bits,'displacement':format(D,'.35g'),'relative_displacement':format(Yr-1,'.35g'),'period':format(M('6.2')+tau,'.35g'),'section_residual_scaled':format(poly(xx,tau),'.5g')}
if __name__=='__main__':
 order=int(sys.argv[1]) if len(sys.argv)>1 else 112
 data=[]
 for r in ['7.0e-75','7.2e-75','2.2e-21','2.3e-21','6.5e-8','6.8e-8']:
  t=time.time();row=compute(r,order);row['seconds']=time.time()-t; data.append(row);print(json.dumps(row),flush=True)
 (Path(__file__).resolve().parent/'data'/('gt_taylor_'+str(order)+'.json')).write_text(json.dumps({'status':'NONRIGOROUS MPFR Taylor shooting; no remainder enclosure','samples':data},indent=2)+'\n')
