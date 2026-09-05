/* Compactified return map for quadratic fields in log-polar coordinates about the origin (the
   caller shifts the field so the focus is at the origin).  x = e^u cos th, y = e^u sin th.
   With P = P0 + P1(th) e^u + P2(th) e^{2u} (P0 constant terms vanish after the shift, P1 linear,
   P2 quadratic parts) and the same for Q:
     du/dt  = e^{-u} (x P + y Q) e^{-u} ... -> du/dt = (cos P + sin Q)/e^u,  dth/dt = (cos Q - sin P)/e^u,
   and we rescale time by ds = e^{-u} ... so that the right-hand side stays bounded as u -> +inf:
     du/ds = (c P + s Q) e^{-2u} * e^{u}?  Concretely with P = A1 e^u + A2 e^{2u}:
     du/ds = (c(A1 + A2 e^u) + s(B1 + B2 e^u)) * e^{-u} * e^{u}/(1+e^u)   (bounded)
     dth/ds = (c(B1 + B2 e^u) - s(A1 + A2 e^u)) * e^{-u} * e^{u}/(1+e^u)
   i.e. multiply the true rates by e^{u}/(1+e^{u}) (positive: same orbits).
   Return: start at (u0, th0), integrate until th has advanced by 2 pi (either sense) and report u1.
   status: 0 ok, 1 u > umax (escape), 2 time cap, 3 step cap, 4 stalled. */
#include <math.h>
#include <stdlib.h>
#include <omp.h>
static inline void rates(const double *c, double u, double th, double *du, double *dth){
  double co = cos(th), si = sin(th), eu = exp(u);
  double A1 = c[1]*co + c[2]*si, A2 = c[3]*co*co + c[4]*co*si + c[5]*si*si;
  double B1 = c[7]*co + c[8]*si, B2 = c[9]*co*co + c[10]*co*si + c[11]*si*si;
  double fac = 1.0/(1.0+eu);                   /* = e^{-u} * e^{u}/(1+e^{u}) */
  double Pn = A1 + A2*eu, Qn = B1 + B2*eu;     /* P e^{-u}, Q e^{-u} */
  *du  = (co*Pn + si*Qn)*fac;
  *dth = (co*Qn - si*Pn)*fac;
}
static double dp_step(const double *c, double u, double th, double h, double *un, double *thn, double rtol, double atol){
  static const double a21=1.0/5, a31=3.0/40, a32=9.0/40, a41=44.0/45, a42=-56.0/15, a43=32.0/9,
  a51=19372.0/6561, a52=-25360.0/2187, a53=64448.0/6561, a54=-212.0/729,
  a61=9017.0/3168, a62=-355.0/33, a63=46732.0/5247, a64=49.0/176, a65=-5103.0/18656,
  a71=35.0/384, a73=500.0/1113, a74=125.0/192, a75=-2187.0/6784, a76=11.0/84,
  e1=71.0/57600, e3=-71.0/16695, e4=71.0/1920, e5=-17253.0/339200, e6=22.0/525, e7=-1.0/40;
  double k1u,k1t,k2u,k2t,k3u,k3t,k4u,k4t,k5u,k5t,k6u,k6t,k7u,k7t;
  rates(c,u,th,&k1u,&k1t);
  rates(c,u+h*a21*k1u, th+h*a21*k1t,&k2u,&k2t);
  rates(c,u+h*(a31*k1u+a32*k2u), th+h*(a31*k1t+a32*k2t),&k3u,&k3t);
  rates(c,u+h*(a41*k1u+a42*k2u+a43*k3u), th+h*(a41*k1t+a42*k2t+a43*k3t),&k4u,&k4t);
  rates(c,u+h*(a51*k1u+a52*k2u+a53*k3u+a54*k4u), th+h*(a51*k1t+a52*k2t+a53*k3t+a54*k4t),&k5u,&k5t);
  rates(c,u+h*(a61*k1u+a62*k2u+a63*k3u+a64*k4u+a65*k5u), th+h*(a61*k1t+a62*k2t+a63*k3t+a64*k4t+a65*k5t),&k6u,&k6t);
  *un = u+h*(a71*k1u+a73*k3u+a74*k4u+a75*k5u+a76*k6u);
  *thn = th+h*(a71*k1t+a73*k3t+a74*k4t+a75*k5t+a76*k6t);
  rates(c,*un,*thn,&k7u,&k7t);
  double eu_ = h*(e1*k1u+e3*k3u+e4*k4u+e5*k5u+e6*k6u+e7*k7u);
  double et_ = h*(e1*k1t+e3*k3t+e4*k4t+e5*k5t+e6*k6t+e7*k7t);
  double su = atol + rtol*fmax(1.0,fabs(u)), st = atol + rtol*1.0;
  return sqrt(0.5*((eu_/su)*(eu_/su)+(et_/st)*(et_/st)));
}
static int full_return_log(const double *c, double u0, double th0, double rtol, double umax, double Smax, long maxsteps, double *u1, double *S1){
  double u = u0, th = th0, s = 0.0, h = 1e-3;
  double du, dth; rates(c,u,th,&du,&dth);
  if (fabs(dth) < 1e-300) return 4;
  double sense = dth > 0 ? 1.0 : -1.0, target = th0 + sense*2*M_PI;
  double atol = 1e-15; long steps = 0;
  while (steps < maxsteps){
    double un, thn; double err = dp_step(c,u,th,h,&un,&thn,rtol,atol);
    if (err > 1.0){ h *= fmax(0.2, 0.9*pow(err,-0.2)); steps++; continue; }
    if ((sense>0 && thn >= target) || (sense<0 && thn <= target)){
      double lo=0, hi=h, um=un, tm=thn;
      for (int it=0; it<60; ++it){ double mid=0.5*(lo+hi); dp_step(c,u,th,mid,&um,&tm,rtol,atol);
        if ((sense>0 && tm>=target) || (sense<0 && tm<=target)) hi=mid; else lo=mid; if (hi-lo<1e-16*h) break; }
      *u1 = um; *S1 = s+0.5*(lo+hi); return 0;
    }
    u = un; th = thn; s += h; steps++;
    if (u > umax) return 1;
    if (s > Smax) return 2;
    rates(c,u,th,&du,&dth);
    if (fabs(du)+fabs(dth) < 1e-14 && s > 1.0) return 4;
    h *= fmin(5.0, 0.9*pow(fmax(err,1e-10),-0.2));
  }
  return 3;
}
void returns_log(int nsets, const double *coef, int nr, const double *u0s, double th0,
                 double *u1s, double *S1s, int *status, double rtol, double umax, double Smax, long maxsteps){
  #pragma omp parallel for schedule(dynamic,1)
  for (int k=0;k<nsets;++k){
    for (int i=0;i<nr;++i){ double u1=NAN,S1=NAN; int st=full_return_log(coef+12*k, u0s[k*nr+i], th0, rtol, umax, Smax, maxsteps, &u1,&S1);
      u1s[k*nr+i]=u1; S1s[k*nr+i]=S1; status[k*nr+i]=st; }
  }
}
