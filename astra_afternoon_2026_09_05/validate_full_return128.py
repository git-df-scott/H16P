"""Analytic precision controls; same interface used by the sign map."""
import json
from pathlib import Path
import mpmath as mp
from full_return128 import full_returns
mp.mp.dps=70
checks=[]
def linear(name,M,eps,scale,grid,angle='0'):
 vec=['0',eps,'1','0','0','0','0','-'+M,eps,'0','0','0']
 r=full_returns(vec,angle,grid,y_scale=scale,tolerance='2e-28')
 target=2*mp.pi*mp.mpf(eps)/mp.sqrt(mp.mpf(M))
 errors=[]
 for p in r['points']:
  assert p['status']=='OK_NUMERICAL',p
  error=abs(mp.mpf(p['log_displacement'])-target);errors.append(mp.nstr(error,8));assert error<mp.mpf('1e-24'),error
 checks.append({'control':name,'target':mp.nstr(target,65),'absolute_errors':errors,'result':r})
linear('near_integrable_1e_minus17','1','1e-17','1',['-20','0','3'])
linear('beyond_exp36','1','0.001','1',['37','44','60'],'1.234567890123456789')
linear('large_coefficient_1e14','100000000000000','0.001','10000000',['-20','0','44'],'0.7')
linear('exact_center_scaled_nonzero_ray','100','0','10',['-20','0','44'],'1.7')
# Exact nonlinear reversible center from the Proposition A counterexample.
vec=['0','0','-1','0','1','1','0','1','0','1','1','0']
r=full_returns(vec,'0',['-3','-2'],tolerance='2e-28')
assert all(p['status']=='OK_NUMERICAL' and abs(mp.mpf(p['log_displacement']))<mp.mpf('1e-24') for p in r['points'])
checks.append({'control':'nonlinear_exact_reversible_center','result':r})
# This field has a stable node, so the angular chart must fail: never return a fake zero.
r=full_returns(['0','-1','0','0','0','0','0','0','-2','0','0','0'],'0',['0'])
assert r['points'][0]['status']=='ANGULAR_CHART_UNRESOLVED' and r['points'][0]['log_displacement'] is None
checks.append({'control':'chart_failure_has_null_displacement','result':r})
try:full_returns([0.0]*12,'0',['0']);raise AssertionError('float silently accepted')
except TypeError:pass
p=Path(__file__).with_name('full_return128_validation.json');p.write_text(json.dumps(checks,indent=2)+'\n');print('PASS:',len(checks),'analytic/failure controls and binary64-input rejection')
