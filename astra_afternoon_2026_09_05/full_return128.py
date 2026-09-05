"""One call: full_returns(rational_vector, ray_angle, log_radius_grid, **options).
All precision-sensitive numbers are strings/Fractions, never binary64 internally.
CLI: python full_return128.py request.json > results.json
"""
from pathlib import Path
from fractions import Fraction as F
import json,subprocess,sys,hashlib,tempfile,os
HERE=Path(__file__).resolve().parent

def full_returns(coefficients,ray_angle,log_radius_grid,*,tolerance='2e-25',y_scale='auto',center=('0','0'),max_evaluations=500000):
 log_radius_grid=list(log_radius_grid)
 if len(coefficients)!=12:raise ValueError('Expected P,Q coefficients in order 1,x,y,x^2,xy,y^2')
 if any(isinstance(x,float) for x in [*coefficients,ray_angle,*log_radius_grid,tolerance,y_scale,*center]):raise TypeError('Use decimal/rational strings, integers, or Fractions; float inputs are forbidden')
 c=list(map(F,coefficients));X,Y=map(F,center);out=[]
 for i in [0,6]:
  a,b,d,e,f,g=c[i:i+6]
  translated=[a+b*X+d*Y+e*X*X+f*X*Y+g*Y*Y,b+2*e*X+f*Y,d+f*X+2*g*Y,e,f,g]
  if translated[0]!=0:raise ValueError('center must be an exact rational equilibrium')
  out.extend(map(str,translated))
 if y_scale=='auto':
  from decimal import Decimal,localcontext
  ratio=-F(out[7])/F(out[2]) if F(out[2]) else F(0)
  with localcontext() as ctx:
   ctx.prec=70
   y_scale=str((Decimal(ratio.numerator)/Decimal(ratio.denominator)).sqrt()) if ratio>0 else '1'
 src=HERE/'full_return128.cpp';key=hashlib.sha256(src.read_bytes()).hexdigest()[:16]
 exe=Path(tempfile.gettempdir())/('fastra-full-return128-'+key)
 if not exe.exists():
  tmp=exe.with_suffix('.'+str(os.getpid())+'.tmp')
  subprocess.run(['g++','-O2','-std=c++17','-fext-numeric-literals',str(src),'-o',str(tmp),'-lquadmath'],check=True)
  os.replace(tmp,exe)
 grid=list(map(str,log_radius_grid));payload='\n'.join(out+[str(ray_angle),str(tolerance),str(y_scale),str(max_evaluations),str(len(grid))]+grid)+'\n'
 run=subprocess.run([str(exe)],input=payload,text=True,capture_output=True,check=True)
 return {'method':'binary128 full angular return','certified':False,'coefficient_order':'P:1,x,y,x^2,xy,y^2; Q:same','coefficients':list(map(str,c)),'center':list(map(str,center)),'ray_angle':str(ray_angle),'tolerance':str(tolerance),'y_scale':str(y_scale),'points':[json.loads(line) for line in run.stdout.splitlines()]}
if __name__=='__main__':
 req=json.load(open(sys.argv[1])) if len(sys.argv)>1 else json.load(sys.stdin)
 print(json.dumps(full_returns(**req),indent=2))
