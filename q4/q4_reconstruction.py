#!/usr/bin/env python3
"""Q4 Green/PF reconstruction, with a bounded three-parameter diagnostic.

Importing provides the evaluator; invoking this script tests only the frozen
universal coefficient point and three specified lift parameters. No search.
"""
import os
for name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
             "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[name] = "1"

import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import hyp2f1

from q4_integrals import alpha_beta_from_mu, basis_mp


def mu_from_universal(kappa, A, B, eta):
    """Numerical transport into the ORIGINAL four-integral coefficient basis."""
    k, d = float(kappa), float(kappa)-1
    mapping = np.column_stack([
        alpha_beta_from_mu(k, row)[1:] for row in np.eye(4)
    ])
    return np.linalg.solve(mapping, [-A-2*k*B/d, B/d, k-d*eta, 1.])


def center_data(a, A, B, eta):
    """Exact rational expressions, evaluated in the caller's arithmetic."""
    y0 = 3*(1326*A+864*B-2431*eta-102)/1361360
    y1 = -3*(1+a)*y0/2-eta/192
    return y0, y1


def loop_homogeneous_slope(a):
    """Initial logarithmic slope of the positive homogeneous solution
    that vanishes at t=1, normalized to one at t=0.
    """
    k = 1/(1-a)
    x = np.arcsinh(np.sqrt(k-1))
    odd = 3/10*np.sinh(5*x/3)+3/2*np.sinh(x/3)
    odd_x = (np.cosh(5*x/3)+np.cosh(x/3))/2
    return -np.sqrt(a)*odd_x/(2*odd)


def initial_weighted_derivative(a,A,B,eta):
    """If q lies in the certified lobe region and this value is <=0,
    the sign-chain theorem in the reconstruction notes excludes five I zeros.
    Floating evaluation alone is not a rigorous parameter certificate.
    """
    y0,y1=center_data(a,A,B,eta)
    return y1-loop_homogeneous_slope(a)*y0


def reconstruct(a, A, B, eta, t_end=0.99):
    """Return dense numerical PF solution (H,Y,Y',X); NOT interval rigorous.

    Y=G/J1(k), X=int_0^t Y(u)/(1-a*u)^(3/2)du,
    I=-a*J1(k)*sqrt(1-a*t)*X/2 and k=1/(1-a).
    """
    if not 0 < a < 1 or not 0 < t_end < 1:
        raise ValueError("require 0<a<1 and 0<t_end<1")
    y0, y1 = center_data(a,A,B,eta)
    q0 = A-1-eta/6

    def rhs(t, state):
        H, Y, V, X = state
        F = hyp2f1(1/6,5/6,1,t)
        Fp = 5/36*hyp2f1(7/6,11/6,2,t)
        M = 1-6*(1-t)*Fp/F
        q = A+B*t-1+(t-eta)*M
        h_over_t2 = q0/2 if t == 0 else H/(t*t)
        forcing = -h_over_t2/(1152*(1-t))
        vprime = (forcing+(1-a)*V/2-5*a*Y/36)/((1-a*t)*(1-t))
        return [t*F*q,V,vprime,Y/(1-a*t)**1.5]

    solution = solve_ivp(rhs,(0,t_end),[0,y0,y1,0],method="DOP853",
                         rtol=3e-12,atol=1e-14,dense_output=True,max_step=.025)
    if not solution.success:
        raise RuntimeError(solution.message)
    return solution


def original_values(a, solution, t):
    k = 1/(1-a)
    C = np.pi/np.sqrt(k-1)
    return -a*C*np.sqrt(1-a*np.asarray(t))*solution.sol(t)[3]/2


def green_coordinates(a, solution, t):
    """Return (Z=Y/y, P=p*y^2*Z') for the exact positive Green factor.
    Numerical evaluation of the scalar solution; t must remain below one.
    """
    k=1/(1-a)
    d=k-1
    t=np.asarray(t)
    s=k-d*t
    x=np.arcsinh(np.sqrt(s-1))
    x0=np.arcsinh(np.sqrt(d))
    odd=lambda xx:3/10*np.sinh(5*xx/3)+3/2*np.sinh(xx/3)
    odd_x=(np.cosh(5*x/3)+np.cosh(x/3))/2
    normalization=odd(x0)
    y=odd(x)/normalization
    yt=-d*odd_x/(2*np.sqrt(s)*np.sqrt(s-1)*normalization)
    state=solution.sol(t)
    p=np.sqrt((1-t)/(1-a*t))
    return state[1]/y,p*(y*state[2]-yt*state[1])


def bounded_diagnostic():
    import resource
    resource.setrlimit(resource.RLIMIT_CPU,(10,10))
    import sympy as S
    z=S.symbols('z',positive=True)
    derivative=lambda f:z*S.diff(f,z)/3  # z=exp(x/3)
    ch=lambda n:(z**n+z**(-n))/2
    sh=lambda n:(z**n-z**(-n))/2
    even=(5*ch(1)-ch(5))/4
    odd=S.Rational(3,10)*sh(5)+S.Rational(3,2)*sh(1)
    for homogeneous in (even,odd):
        assert S.factor(ch(3)*derivative(derivative(homogeneous))
                        -2*sh(3)*derivative(homogeneous)
                        +S.Rational(5,9)*ch(3)*homogeneous) == 0
    assert S.factor(even*derivative(odd)-derivative(even)*odd-ch(3)**2) == 0
    k,s=S.symbols('k s',positive=True)
    h_squared=4*s/(9*k)
    transformed_forcing=-2*h_squared*k*k/(16*s*64*(s-k)**2*(s-1))
    assert S.factor(transformed_forcing+k/(1152*(s-k)**2*(s-1))) == 0
    print("Exact elementary Green basis and forcing-sign conversion passed.")
    A = 1243911778077/10**12
    B = -86917392526/10**12
    eta = 1460428426173/10**12
    witnesses = np.array([1/8,3/8,5/8,7/8,.99])
    print("Numerical diagnostic only; frozen universal point:",A,B,eta)
    for a in (.25,.5,.75):
        solution = reconstruct(a,A,B,eta)
        values = original_values(a,solution,witnesses)
        print("a=",a,"kappa=",1/(1-a),"center Y,Y'=",center_data(a,A,B,eta))
        print("positive-homogeneous slope r=",loop_homogeneous_slope(a),
              "initial weighted derivative=",initial_weighted_derivative(a,A,B,eta))
        print("t=",witnesses.tolist())
        print("I=",values.tolist())
        if a == .5:
            mu = mu_from_universal(2,A,B,eta)
            independent = float(sum(m*v for m,v in zip(mu,basis_mp(2,1.5,dps=40))))
            reconstructed = float(original_values(a,solution,.5))
            difference = abs(independent-reconstructed)
            print("original mu=",mu.tolist())
            print("independent area at k=2,s=1.5:",independent)
            print("PF reconstruction:",reconstructed,"absolute difference:",difference)
            assert difference < 2e-11, "independent original-normalization check failed"
    # An exact negative initial weighted derivative: no transcendental bounds
    # are required at k=cosh(3*log(6/5))^2, because all hyperbolics are rational.
    from fractions import Fraction as Q
    Aq,Bq,eq=Q(1243911778077,10**12),-Q(86917392526,10**12),Q(1460428426173,10**12)
    v=Q(6,5)
    sinh=lambda n:(v**n-v**(-n))/2
    cosh=lambda n:(v**n+v**(-n))/2
    kq=cosh(3)**2
    aq=1-1/kq
    odd=Q(3,10)*sinh(5)+Q(3,2)*sinh(1)
    odd_x=(cosh(5)+cosh(1))/2
    rq=-sinh(3)/cosh(3)*odd_x/(2*odd)
    y0q,y1q=center_data(aq,Aq,Bq,eq)
    p0q=y1q-rq*y0q
    assert p0q == -Q(3056925605483331742344151782161,2151790655250172000064000000000000)
    print("Exact exclusion lift kappa=",kq,"r=",rq,"P0=",p0q)
    ratio_bound=Q(601,136136)
    a_max=Q(2,3)*(1/(192*ratio_bound)-1)
    assert a_max == Q(2593,21636)
    assert 1/(1-a_max) == Q(21636,19043)
    assert Q(1,192)-Q(23,20)*ratio_bound == Q(2147,16336320)
    print("Exact uniform exclusion endpoint kappa=",1/(1-a_max))
    point_margin=eq+307*y0q
    margin_coefficients=[Q(2763,3080),Q(49734,85085),-Q(361,560)]
    box_margin=point_margin-Q(1,10**7)*sum(map(abs,margin_coefficients))
    assert point_margin == Q(14871489355525071,272272000000000000)
    assert box_margin == Q(79526371464733,1456000000000000)
    assert box_margin > 0
    k_poly=(1+z)*(1+4*z)**2
    D_poly=8*z*z+10*z+5
    numerator=S.expand(k_poly*D_poly+18*D_poly-20*k_poly)
    assert S.expand(numerator-(128*z**5+352*z**4+z*(72*z*z-118*z+55)+75)) == 0
    assert 118**2-4*72*55 == -1916
    assert Q(3,8)*Q(8,5)**3 == Q(192,125)
    assert -1+Q(192,125)*Q(9,4) == Q(307,125)
    assert (Q(5,3)**2-1)/2 == Q(8,9)
    assert Q(1,216)-ratio_bound == Q(395,1837836)
    assert 11*3**11 == 1948617 < 2000000 == 2**7*5**6
    print("Exact all-kappa candidate-box exclusion margin:",box_margin)
    print("Exact universal first-primitive-zero exclusion threshold: 5/11")
    print("Bounded diagnostic finished; no five-sign candidate asserted.")


if __name__ == "__main__":
    bounded_diagnostic()
