#!/usr/bin/env python3
"""Forty frozen boundary determinant diagnostics, not a global certificate.

Integrates the determinant itself in logarithmic time to avoid reconstructing
it solely by subtracting large products. An independent a=1 quadrature checks
one of the frozen values. The a=1 values are comparison limits, not finite lifts.
"""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[key]='1'
import sys,json,hashlib,resource,time
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import mpmath as mp
import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import hyp2f1
from q4_threshold_path import primitive_basis_closed


def basis_derivative(t):
    F=mp.hyp2f1(mp.mpf(1)/6,mp.mpf(5)/6,1,t)
    Fp=mp.mpf(5)/36*mp.hyp2f1(mp.mpf(7)/6,mp.mpf(11)/6,2,t)
    D=t*(F-6*(1-t)*Fp)
    return (t*F,t*t*F,D,t*D)


def anchored_pair(r,face):
    rows=[primitive_basis_closed(r),basis_derivative(r) if face=='confluent' else primitive_basis_closed(1)]
    mat=mp.matrix([[row[0],row[1]] for row in rows])
    pairs=[]
    for j in (3,2):
        ab=mp.lu_solve(mat,mp.matrix([-row[j] for row in rows]))
        coeff=[ab[0],ab[1],mp.mpf(0),mp.mpf(0)];coeff[j]=1
        pairs.append(tuple(coeff))
    return tuple(pairs)


def center(coeff):
    alpha,beta,gamma,delta=coeff
    return mp.mpf(9)/3080*(alpha+mp.mpf(144)/221*beta+mp.mpf(11)/6*gamma+mp.mpf(204)/221*delta)


def determinant(a,r,pair):
    co=np.array(pair,dtype=float)
    yc=np.array([float(center(c)) for c in pair])
    slope=-1.5*(1+a)*yc+co[:,2]/192
    initial_det=-yc[0]/192 # eta_E=0, eta_V=-1
    initial=np.r_[np.zeros(2),yc,slope,initial_det]
    def rhs(x,state):
        d=np.exp(-x);t=-np.expm1(-x);b=1-a*t
        H,Y,W=state[:2],state[2:4],state[4:6]
        F=hyp2f1(1/6,5/6,1,t)
        Fp=5/36*hyp2f1(7/6,11/6,2,t)
        M=1-6*d*Fp/F
        Hder=d*t*F*(co[:,0]+co[:,1]*t+(co[:,2]+co[:,3]*t)*M)
        hratio=(co[:,0]+co[:,2]/6)/2 if x==0 else H/(t*t)
        Wder=-W+(1-a)*W/(2*b)-5*a*d*Y/(36*b)-hratio/(1152*b)
        detder=(hratio[1]*Y[0]-hratio[0]*Y[1])/(1152*b**1.5*np.sqrt(d))
        return np.r_[Hder,W,Wder,detder]
    end=-np.log1p(-float(r))
    sol=solve_ivp(rhs,(0,end),initial,method='DOP853',rtol=3e-12,atol=2e-15,max_step=.08)
    if not sol.success:raise RuntimeError(sol.message)
    H,Y,W=sol.y[:2,-1],sol.y[2:4,-1],sol.y[4:6,-1]
    d=1-float(r);p=np.sqrt(d/(1-a*float(r)))
    subtraction=p/d*(W[0]*Y[1]-Y[0]*W[1])
    integrated=sol.y[6,-1]
    return integrated,subtraction,float(initial_det),H.tolist()


def independent_a1_quadrature(r,pair):
    """Direct high-precision integral of P and Phi, no numerical ODE."""
    values=[]
    with mp.workdps(45):
      for co in pair:
        y0=center(co);eta=-co[2]
        def hh(t):
            if t==0:return (co[0]+co[2]/6)/2
            row=primitive_basis_closed(t)
            return sum(c*k for c,k in zip(co,row))/(t*t)
        split=[mp.mpf(0),r/2,r]
        integ=mp.quad(lambda t:hh(t)*(1-t)**(-mp.mpf(13)/6)/1152,split)
        phii=mp.quad(lambda t:hh(t)*(mp.mpf(3)/2)*((1-t)**(-mp.mpf(2)/3)-1)*(1-t)**(-mp.mpf(13)/6)/1152,split)
        P=-mp.mpf(13)/6*y0-eta/192-integ
        Phi=y0+phii
        values.append((P,Phi))
      result=values[0][0]*values[1][1]-values[0][1]*values[1][0]
      return mp.nstr(result,35)


def main():
    resource.setrlimit(resource.RLIMIT_CPU,(60,60));os.nice(10)
    mp.mp.dps=55;started=time.process_time();rows=[];control=None
    for rt in ('.825','.95','.999','.99999'):
      r=mp.mpf(rt)
      for face in ('confluent','endpoint'):
        pair=anchored_pair(r,face)
        for a in (.66,.875,.99,.99999,1.):
          D,sub,D0,residual=determinant(a,r,pair)
          row={'r':rt,'a':a,'face':face,'D':D,'D_initial':D0,
               'D_via_product_subtraction':sub,'difference':D-sub,
               'H_anchor_residuals':residual}
          rows.append(row)
          if rt=='.95' and face=='endpoint' and a==1.:
            independent=independent_a1_quadrature(r,pair)
            difference=abs(float(independent)-D)
            assert difference<2e-12, 'independent quadrature disagrees'
            control={'r':rt,'a':a,'face':face,'mp_dps':45,'quadrature_D':independent,
                     'ODE_D':D,'absolute_difference':difference}
    result={'status':'NUMERICAL_ONLY','count':len(rows),
            'nonnegative_D_count':sum(int(row['D']>=0) for row in rows),
            'domain':'two limiting anchor faces; four first anchors and five lifts, including a=1',
            'warning':'Finite samples do not determine the sign between samples, at other lifts, or at unsampled joint limits.',
            'independent_control':control,'rows':rows,
            'cpu_seconds':time.process_time()-started,
            'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    Path(__file__).with_name('boundary_diagnostic.json').write_text(json.dumps(result,indent=2)+'\n')
    print('Frozen face diagnostics:',len(rows),'nonnegative',result['nonnegative_D_count'])
    print('Largest D:',max(rows,key=lambda row:row['D']))
    print('Independent quadrature:',control)
    print('CPU seconds:',result['cpu_seconds'])
if __name__=='__main__':main()
