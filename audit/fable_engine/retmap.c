/* Fable return-map engine for planar quadratic fields.
   P = c0 + c1 x + c2 y + c3 x^2 + c4 x y + c5 y^2
   Q = c6 + c7 x + c8 y + c9 x^2 + c10 x y + c11 y^2
   For each parameter set and each radius r, start at focus + r*dir and
   integrate (Dormand-Prince 5(4), adaptive) until the orbit has turned by
   2*pi about the focus and crosses the ray again; report the crossing
   radius R (displacement D = R - r computed by the caller).
   status: 0 ok, 1 escaped (|p-focus|>Rmax), 2 time cap, 3 step cap,
           4 stalled (approaching an equilibrium), 5 bad start.          */
#include <math.h>
#include <stdlib.h>
#include <omp.h>

static inline void field(const double *c, double x, double y, double *fx, double *fy){
  *fx = c[0] + c[1]*x + c[2]*y + c[3]*x*x + c[4]*x*y + c[5]*y*y;
  *fy = c[6] + c[7]*x + c[8]*y + c[9]*x*x + c[10]*x*y + c[11]*y*y;
}

/* one DP45 step; returns error estimate (scaled) and writes new state */
static double dp_step(const double *c, double x, double y, double h,
                      double *xn, double *yn, double rtol, double atol){
  static const double
  a21=1.0/5,
  a31=3.0/40, a32=9.0/40,
  a41=44.0/45, a42=-56.0/15, a43=32.0/9,
  a51=19372.0/6561, a52=-25360.0/2187, a53=64448.0/6561, a54=-212.0/729,
  a61=9017.0/3168, a62=-355.0/33, a63=46732.0/5247, a64=49.0/176, a65=-5103.0/18656,
  a71=35.0/384, a73=500.0/1113, a74=125.0/192, a75=-2187.0/6784, a76=11.0/84,
  e1=71.0/57600, e3=-71.0/16695, e4=71.0/1920, e5=-17253.0/339200, e6=22.0/525, e7=-1.0/40;
  double k1x,k1y,k2x,k2y,k3x,k3y,k4x,k4y,k5x,k5y,k6x,k6y,k7x,k7y;
  field(c,x,y,&k1x,&k1y);
  field(c,x+h*a21*k1x, y+h*a21*k1y,&k2x,&k2y);
  field(c,x+h*(a31*k1x+a32*k2x), y+h*(a31*k1y+a32*k2y),&k3x,&k3y);
  field(c,x+h*(a41*k1x+a42*k2x+a43*k3x), y+h*(a41*k1y+a42*k2y+a43*k3y),&k4x,&k4y);
  field(c,x+h*(a51*k1x+a52*k2x+a53*k3x+a54*k4x), y+h*(a51*k1y+a52*k2y+a53*k3y+a54*k4y),&k5x,&k5y);
  field(c,x+h*(a61*k1x+a62*k2x+a63*k3x+a64*k4x+a65*k5x), y+h*(a61*k1y+a62*k2y+a63*k3y+a64*k4y+a65*k5y),&k6x,&k6y);
  *xn = x+h*(a71*k1x+a73*k3x+a74*k4x+a75*k5x+a76*k6x);
  *yn = y+h*(a71*k1y+a73*k3y+a74*k4y+a75*k5y+a76*k6y);
  field(c,*xn,*yn,&k7x,&k7y);
  double ex = h*(e1*k1x+e3*k3x+e4*k4x+e5*k5x+e6*k6x+e7*k7x);
  double ey = h*(e1*k1y+e3*k3y+e4*k4y+e5*k5y+e6*k6y+e7*k7y);
  double sx = atol + rtol*fmax(fabs(x),fabs(*xn));
  double sy = atol + rtol*fmax(fabs(y),fabs(*yn));
  return sqrt(0.5*((ex/sx)*(ex/sx)+(ey/sy)*(ey/sy)));
}

/* full return of the point focus + r*dir. Returns status; *Rout crossing radius,
   *Tout time, *turn = +1/-1 rotation sense. */
static int full_return(const double *c, double fx0, double fy0, double dx, double dy,
                       double r, double rtol, double Rmax, double Tmax, long maxsteps,
                       double *Rout, double *Tout){
  double x = fx0 + r*dx, y = fy0 + r*dy, t = 0;
  double px = -dy, py = dx;           /* perpendicular to the ray */
  double vx, vy; field(c,x,y,&vx,&vy);
  double sp = vx*px + vy*py;          /* initial transverse velocity */
  if (fabs(sp) < 1e-300) return 5;
  double sense = sp > 0 ? 1.0 : -1.0;
  double cum = 0.0;                   /* cumulative angle (signed) */
  double atol = 1e-16*r + 1e-300;
  double speed = sqrt(vx*vx+vy*vy);
  double h = 1e-3*r/(speed+1e-300); if (h<=0||!isfinite(h)) h=1e-6;
  long steps = 0;
  double gprev = 0.0; /* signed perpendicular distance at previous accepted state = 0 at start */
  double relx = x-fx0, rely = y-fy0;
  while (steps < maxsteps){
    double xn, yn;
    double err = dp_step(c,x,y,h,&xn,&yn,rtol,atol);
    if (err > 1.0){ h *= fmax(0.2, 0.9*pow(err,-0.2)); steps++; continue; }
    /* accepted */
    double nrelx = xn-fx0, nrely = yn-fy0;
    double dang = atan2(relx*nrely - rely*nrelx, relx*nrelx + rely*nrely);
    double ncum = cum + dang;
    double gnew = nrelx*px + nrely*py;
    double dot = nrelx*dx + nrely*dy;
    /* crossing of the ray in the right sense after (nearly) a full turn */
    if (sense*ncum >= 2*M_PI - 1e-9 && dot > 0 && ((gprev<=0 && gnew>=0 && sense>0) || (gprev>=0 && gnew<=0 && sense<0) || (steps>0 && gprev*gnew<=0))){
      /* bisect on step length to hit g = 0 */
      double lo = 0.0, hi = h, xm=xn, ym=yn;
      for (int it=0; it<60; ++it){
        double mid = 0.5*(lo+hi);
        dp_step(c,x,y,mid,&xm,&ym,rtol,atol);
        double gm = (xm-fx0)*px + (ym-fy0)*py;
        if ((gm>0) == (gnew>0)) hi = mid; else lo = mid;
        if (hi-lo < 1e-15*h) break;
      }
      *Rout = (xm-fx0)*dx + (ym-fy0)*dy;
      *Tout = t + 0.5*(lo+hi);
      return 0;
    }
    x = xn; y = yn; t += h; cum = ncum; gprev = gnew; relx = nrelx; rely = nrely;
    double dist = sqrt(relx*relx+rely*rely);
    if (dist > Rmax) return 1;
    if (t > Tmax) return 2;
    /* stall: tiny speed relative to radius (approaching a point) */
    field(c,x,y,&vx,&vy); speed = sqrt(vx*vx+vy*vy);
    if (speed < 1e-13*(1.0+dist) && t > 1.0) return 4;
    h *= fmin(5.0, 0.9*pow(fmax(err,1e-10),-0.2));
    steps++;
  }
  return 3;
}

void returns(int nsets, const double *coef, const double *foc, const double *dir,
             int nr, const double *radii, double *Rout, double *Tout, int *status,
             double rtol, double Rmax, double Tmax, long maxsteps){
  #pragma omp parallel for schedule(dynamic,1)
  for (int s=0; s<nsets; ++s){
    const double *c = coef + 12*s;
    for (int i=0;i<nr;++i){
      double R=NAN, T=NAN;
      int st = full_return(c, foc[2*s], foc[2*s+1], dir[2*s], dir[2*s+1],
                           radii[s*nr+i], rtol, Rmax, Tmax, maxsteps, &R, &T);
      Rout[s*nr+i]=R; Tout[s*nr+i]=T; status[s*nr+i]=st;
    }
  }
}
