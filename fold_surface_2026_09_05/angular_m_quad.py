"""Build/run the archived long-double source; one process per charged call."""
from pathlib import Path
import hashlib,json,subprocess,sys,resource
import mpmath as mp
HERE=Path(__file__).resolve().parent;source=HERE/'angular_m_quad.cpp'
exe=Path('/tmp')/('h16_angular_m_quad_'+hashlib.sha256(source.read_bytes()).hexdigest()[:16])
if not exe.exists():subprocess.run(['g++','-O2','-std=c++17','-fext-numeric-literals',str(source),'-o',str(exe),'-lquadmath'],check=True,capture_output=True)
mp.mp.dps=50
def num(x):
    if isinstance(x,str) and '/' in x:
        a,b=x.split('/');return mp.mpf(a)/mp.mpf(b)
    return mp.mpf(str(x))
q=json.load(sys.stdin)
resource.setrlimit(resource.RLIMIT_CPU,(10,10))
data=' '.join(mp.nstr(num(q[k]),45) for k in ('r','c','m'))+' '+str(q.get('tol','2e-17'))+'\n'
p=subprocess.run([str(exe)],input=data,text=True,capture_output=True,timeout=15)
if p.returncode:print(json.dumps(dict(status='UNRESOLVED',error='long-double process exit '+str(p.returncode))))
else:print(p.stdout.strip())
