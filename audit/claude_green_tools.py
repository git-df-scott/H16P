"""Claude's independent Green-reconstruction toolkit (quadrature based, no ODE).
P(t)=P0-int_0^t Omega H,  Z(t)=Y0+int_0^t P/(p y^2),  Phi(t)=Y0+int_0^t Rcal Omega H,
Y = Phi*y + P*y2  (variation of parameters, y2=y*Rcal)."""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "q4"))
import mpmath as mp
from q4_threshold_path import primitive_basis_closed, coefficients_from_r, threshold_anchors, from_primitive_anchors_closed

one6, five6 = mp.mpf(1)/6, mp.mpf(5)/6
def Fh(t): return mp.hyp2f1(one6, five6, 1, t)
def Fph(t): return mp.mpf(5)/36*mp.hyp2f1(one6+1, five6+1, 2, t)
def Mh(t): return 1-6*(1-t)*Fph(t)/Fh(t)

def Hval(co, t):
    A, B, eta = co
    r = primitive_basis_closed(t)
    return (A-1)*r[0]+B*r[1]-eta*r[2]+r[3]
def Hstar(t):
    return 6*t*(1-t)**2*Fh(t)*(6*Mh(t)-1)/77

class Lift:
    def __init__(self, k):
        self.k = mp.mpf(k); self.a = 1-1/self.k; self.d = self.k-1
        self.xk = mp.asinh(mp.sqrt(self.d))
        self.Ok = self.O(self.xk)
    @staticmethod
    def O(x): return mp.mpf(3)/10*mp.sinh(5*x/3)+mp.mpf(3)/2*mp.sinh(x/3)
    @staticmethod
    def E(x): return (5*mp.cosh(x/3)-mp.cosh(5*x/3))/4
    @staticmethod
    def Ox(x): return (mp.cosh(5*x/3)+mp.cosh(x/3))/2
    def x_of_t(self, t): return mp.asinh(mp.sqrt(self.d*(1-t)))
    def y(self, t): return self.O(self.x_of_t(t))/self.Ok
    def r(self): return -mp.sqrt(self.a)*self.Ox(self.xk)/(2*self.Ok)
    def Omega(self, t):
        return self.y(t)/(1152*t*t*(1-self.a*t)**mp.mpf(1.5)*(1-t)**mp.mpf(1.5))
    def Rcal(self, t):
        x = self.x_of_t(t)
        return 2*self.Ok**2/mp.sqrt(self.k*self.d)*(self.E(x)/self.O(x)-self.E(self.xk)/self.Ok)
    def center(self, co):
        A, B, eta = co
        Y0 = 3*(1326*A+864*B-2431*eta-102)/1361360
        Y1 = -mp.mpf(3)/2*(1+self.a)*Y0-eta/192
        return Y0, Y1-self.r()*Y0

def breakpoints(t_end, eps_scale=None):
    pts = [mp.mpf('1e-12'), mp.mpf('0.25'), mp.mpf('0.5'), mp.mpf('0.75')]
    g = mp.mpf('0.875')
    while 1-g > 4*(1-t_end) and g < t_end:
        pts.append(g); g = 1-(1-g)/2
    pts.append(t_end)
    return pts

def P_Phi_at(lift, co, t_end, Hfun=None):
    Hf = (lambda u: Hval(co, u)) if Hfun is None else Hfun
    Y0, P0 = lift.center(co)
    pts = breakpoints(t_end)
    V = mp.quad(lambda u: lift.Omega(u)*Hf(u), pts)
    W = mp.quad(lambda u: lift.Rcal(u)*lift.Omega(u)*Hf(u), pts)
    return P0-V, Y0+W, Y0, P0
