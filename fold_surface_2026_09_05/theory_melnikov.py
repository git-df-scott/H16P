#!/usr/bin/env python3
"""Numerical quadrature of exact center Melnikov moments; NOT a certificate.

No ODE integration. The endpoint bisection and quadrature use 45-digit
mpmath point arithmetic. Two finite-difference steps are compared.
"""
import json
from pathlib import Path
import mpmath as p
p.mp.dps=45
c=p.findroot(lambda z:305+634*z-11*z*z-1000*z**3,p.mpf('.9686'))
m=210/(11*c-5); sigma=5*(2*c+1)/21; D=1+m*sigma**2
k=sigma*p.sqrt(m)
A=(10*c*sigma*m+5*sigma**3*m*m-sigma**2*m-5*sigma*m+100*sigma+11)/(5*D**2)
B=sigma*(5*c-50*sigma**2-6*sigma-5)/(5*D**2)
C=(5*c*sigma**3*m*m+11*sigma**2*m+5*sigma*m-50*sigma+5)/(5*D**2)
b=B*m/A; d=C/A
a0=(1+b)/d; a1=2*(1+2*b)/(2*d-1); a2=b/(d-1)
Einf=-(a0-a1+a2)
nu=2*d+1; am=1/(m*D); ac=p.sqrt(m)/(A*D*D)
Aco=2*D-nu*k*k/(d-1)
Bco=-2*D+nu-nu*k*k*b/(d-1)
BK=am*k*nu*5/(11*c-5)
CC1=ac*Aco-am*k*nu*11*m/(11*c-5)
CC2=ac*Bco
moment_calls=0

def energy(r):
    y=-r*r/(1+r)
    X=A/m*(y-m*sigma*r)
    Y=A/p.sqrt(m)*(r+sigma*y)
    ss=1-X
    return (Y*Y+a0-a1*ss+a2*ss*ss)*ss**(-2*d)+Einf

def moments(E):
    global moment_calls
    moment_calls+=1
    h=Einf-E
    def q(ss):return -a0+a1*ss-a2*ss*ss-h*ss**(2*d)
    left=p.mpf('.999'); right=p.mpf('1.001')
    while q(left)>0:left/=2
    while q(right)>0:right*=2
    def bisect(lo,hi):
        flo=q(lo)
        for _ in range(165):
            mid=(lo+hi)/2
            if q(mid)*flo>0:lo=mid
            else:hi=mid
        return (lo+hi)/2
    left=bisect(left,1); right=bisect(1,right)
    def integral(power):
        def integrand(t):
            ss=left+(right-left)*p.sin(t)**2
            ds=2*(right-left)*p.sin(t)*p.cos(t)
            # Point-roundoff clipping at an approximate endpoint; this is
            # another reason these results are NOT interval enclosures.
            return (1-ss)**power*ss**(-2*d-2)*p.sqrt(max(p.mpf(0),q(ss)))*ds
        return 4*p.quad(integrand,[0,p.pi/4,p.pi/2])
    return integral(1),integral(2)

def ratio(E):
    J1,J2=moments(E)
    return J2/J1

def fmt(z):return p.nstr(z,32)

rows=[]
for text in ('.0001','.00005'):
    step=p.mpf(text)
    def derivative(E):return (ratio(E+step)-ratio(E-step))/(2*step)
    E=p.findroot(derivative,(p.mpf('1.20'),p.mpf('1.21')),tol=p.mpf('1e-30'))
    N=ratio(E)
    slope=-BK/(CC1+CC2*N)
    r=p.findroot(lambda z:energy(z)-E,7)
    second=(ratio(E+step)-2*N+ratio(E-step))/step**2
    rows.append(dict(finite_difference_step=text,energy=fmt(E),
                     section_r=fmt(r),moment_ratio=fmt(N),
                     limiting_dc_dK=fmt(slope),ratio_second_derivative=fmt(second)))
output=dict(status='NUMERICAL_QUADRATURE_ONLY',orbit_evaluations=0,
            precision_digits=p.mp.dps,moment_calls=moment_calls,
            center_c=fmt(c),center_m=fmt(m),normalized_b=fmt(b),
            normalized_d=fmt(d),annulus_outer_energy=fmt(Einf),rows=rows,
            limitations='Local stationary ratio only; no uniqueness, interval enclosure, or complete-component theorem.')
Path(__file__).with_name('theory_melnikov.json').write_text(json.dumps(output,indent=2)+'\n')
print(json.dumps(output,indent=2))
