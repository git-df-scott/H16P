/* lane1_ahcurve/ahcurve.c
 * Andronov-Hopf curve engine for planar quadratic vector fields.
 *
 * Field is supplied ALREADY EXPANDED ABOUT THE FOCUS (local coordinates u,v):
 *     P = c0 + c1 u + c2 v + c3 u^2 + c4 u v + c5 v^2
 *     Q = c6 + c7 u + c8 v + c9 u^2 + c10 u v + c11 v^2
 * with c0 = c6 = 0 enforced by the caller (this is REVIEW_engine.md fix A2:
 * never integrate in a global coordinate whose magnitude dwarfs the orbit).
 *
 * Section: the ray  s -> s*(dx,dy),  s > 0,  from the focus at the origin.
 * By PROTOCOL rule 4 a quadratic limit cycle is convex and encloses exactly one
 * focus, so it meets this ray exactly once, transversally.
 *
 * Return map: integrate Dormand-Prince 5(4) accumulating the SIGNED ANGLE about
 * the origin and stop when |angle| = 2*pi, bisecting on the last step.  At that
 * instant the orbit is on the ray by construction, so no separate crossing
 * predicate is needed (this removes REVIEW_engine.md bug A3 entirely, rather
 * than patching its dead disjunct).
 *
 * Displacement D(s) = R(s) - s.
 *
 * Two one-parameter families are supported, both with the strict-monotonicity
 * property of Duff/Perko:
 *   rot: X_b = (P cos b - Q sin b, P sin b + Q cos b)         (uniform rotation)
 *   lin: X_t = X + t*E   for a caller-supplied coefficient direction E
 *        (used for the Cherkas rotating parameter a11, where E has no constant
 *         term so the focus is fixed).
 * For each s the engine root-solves the family parameter so that D = 0; that
 * value is beta*(s), the Andronov-Hopf curve.
 *
 * status codes
 *   return map: 0 ok, 1 escaped, 2 time cap, 3 step cap, 4 stalled, 5 bad start
 *   ah curve  : 0 ok, 10 no sign change found inside the span (UNRESOLVED),
 *               11 base evaluation failed, 12 root solve did not converge
 */
#include <math.h>
#include <stdlib.h>
#include <string.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static inline void field(const double *c, double x, double y, double *fx, double *fy){
  *fx = c[0] + c[1]*x + c[2]*y + c[3]*x*x + c[4]*x*y + c[5]*y*y;
  *fy = c[6] + c[7]*x + c[8]*y + c[9]*x*x + c[10]*x*y + c[11]*y*y;
}

/* build the rotated coefficient vector */
static void build_rot(const double *c, double b, double *out){
  double cb = cos(b), sb = sin(b);
  for (int k=0;k<6;++k){
    out[k]   = c[k]*cb - c[6+k]*sb;
    out[6+k] = c[k]*sb + c[6+k]*cb;
  }
}
/* build c + t*e */
static void build_lin(const double *c, const double *e, double t, double *out){
  for (int k=0;k<12;++k) out[k] = c[k] + t*e[k];
}

/* one DP45 step */
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

/* Full return of the point s*(dx,dy) on the ray.  Stops at signed angle 2*pi. */
static int full_return(const double *c, double dx, double dy, double s,
                       double rtol, double Rmax, double Tmax, long maxsteps,
                       double *Rout, double *Tout){
  double x = s*dx, y = s*dy, t = 0.0;
  double vx, vy; field(c,x,y,&vx,&vy);
  /* transverse component of the velocity at the start: the section must be
     transversal there, otherwise the "one turn" notion is meaningless */
  double px = -dy, py = dx;
  double sp = vx*px + vy*py;
  if (!(fabs(sp) > 0.0)) return 5;
  double sense = sp > 0 ? 1.0 : -1.0;
  double atol = 1e-16*s + 1e-300;                    /* REVIEW_engine.md fix A1 */
  double speed = sqrt(vx*vx+vy*vy);
  double h = 1e-3*s/(speed+1e-300);
  if (!(h > 0.0) || !isfinite(h)) h = 1e-9;
  double cum = 0.0, rx = x, ry = y;
  long steps = 0;
  const double TARGET = 2.0*M_PI;
  while (steps < maxsteps){
    double xn, yn;
    double err = dp_step(c,x,y,h,&xn,&yn,rtol,atol);
    if (!(err <= 1.0)){                              /* NaN-safe rejection */
      h *= isfinite(err) ? fmax(0.2, 0.9*pow(err,-0.2)) : 0.1;
      steps++;
      if (!(fabs(h) > 1e-290)) return 3;
      continue;
    }
    double nrx = xn, nry = yn;
    double dang = atan2(rx*nry - ry*nrx, rx*nrx + ry*nry);
    double ncum = cum + dang;
    if (sense*ncum >= TARGET){
      /* bisect on the step length so that sense*(cum+dang(hm)) == 2*pi */
      double lo = 0.0, hi = h, xm = xn, ym = yn;
      for (int it=0; it<80; ++it){
        double mid = 0.5*(lo+hi);
        dp_step(c,x,y,mid,&xm,&ym,rtol,atol);
        double dm = atan2(rx*ym - ry*xm, rx*xm + ry*ym);
        if (sense*(cum+dm) >= TARGET) hi = mid; else lo = mid;
        if (hi-lo <= 1e-16*h) break;
      }
      double R = sqrt(xm*xm + ym*ym);
      if (!(R > 0.0) || !isfinite(R)) return 3;
      *Rout = R;
      *Tout = t + 0.5*(lo+hi);
      return 0;
    }
    x = xn; y = yn; t += h; cum = ncum; rx = nrx; ry = nry;
    double dist = sqrt(x*x + y*y);
    if (!isfinite(dist)) return 1;
    if (dist > Rmax) return 1;
    if (t > Tmax) return 2;
    field(c,x,y,&vx,&vy); speed = sqrt(vx*vx+vy*vy);
    if (speed < 1e-14*(1.0+dist)) return 4;          /* approaching an equilibrium */
    h *= fmin(5.0, 0.9*pow(fmax(err,1e-10),-0.2));
    steps++;
  }
  return 3;
}

/* ---------------- exported: raw displacement on a fixed member ------------- */

void returns_rot(int nsets, const double *coef, const double *dir, const double *bval,
                 int ns, const double *svals,
                 double *Rout, double *Tout, int *status,
                 double rtol, double Rmax, double Tmax, long maxsteps){
#pragma omp parallel for schedule(dynamic,1)
  for (int i=0;i<nsets*ns;++i){
    int sIdx = i/ns, j = i%ns;
    double cr[12];
    build_rot(coef+12*sIdx, bval[sIdx], cr);
    double R = NAN, T = NAN;
    int st = full_return(cr, dir[2*sIdx], dir[2*sIdx+1], svals[i],
                         rtol, Rmax, Tmax, maxsteps, &R, &T);
    Rout[i]=R; Tout[i]=T; status[i]=st;
  }
}

void returns_lin(int nsets, const double *coef, const double *evec, const double *dir,
                 const double *tval, int ns, const double *svals,
                 double *Rout, double *Tout, int *status,
                 double rtol, double Rmax, double Tmax, long maxsteps){
#pragma omp parallel for schedule(dynamic,1)
  for (int i=0;i<nsets*ns;++i){
    int sIdx = i/ns, j = i%ns;
    double cl[12];
    build_lin(coef+12*sIdx, evec+12*sIdx, tval[sIdx], cl);
    double R = NAN, T = NAN;
    int st = full_return(cl, dir[2*sIdx], dir[2*sIdx+1], svals[i],
                         rtol, Rmax, Tmax, maxsteps, &R, &T);
    Rout[i]=R; Tout[i]=T; status[i]=st;
  }
}

/* ---------------- exported: the Andronov-Hopf curve ----------------------- */

typedef struct {
  const double *c;      /* base coefficients, local */
  const double *e;      /* linear direction (lin family only) */
  int mode;             /* 0 = rotation, 1 = linear */
  double dx, dy, s;
  double rtol, Rmax, Tmax;
  long maxsteps;
} Ctx;

/* D(param) = R - s; returns 0 on success */
static int disp(const Ctx *k, double p, double *D){
  double cc[12], R, T;
  if (k->mode == 0) build_rot(k->c, p, cc); else build_lin(k->c, k->e, p, cc);
  int st = full_return(cc, k->dx, k->dy, k->s, k->rtol, k->Rmax, k->Tmax, k->maxsteps, &R, &T);
  if (st != 0) return st;
  *D = R - k->s;
  return 0;
}

/* monotone root solve: bracket outward from p0, then bisection + secant */
static int solve_param(const Ctx *k, double p0, double span, double ptol,
                       double *pout, int *nev){
  double D0;
  int st = disp(k, p0, &D0);
  *nev = 1;
  if (st != 0) return 11;
  if (D0 == 0.0){ *pout = p0; return 0; }

  /* probe both directions with a small step, then commit to the descending one */
  double h0 = fmax(1e-4, span*1e-3);
  double Dp, Dm; int stp, stm;
  stp = disp(k, p0+h0, &Dp); (*nev)++;
  stm = disp(k, p0-h0, &Dm); (*nev)++;

  double lo=0, hi=0, Dlo=0, Dhi=0; int have=0;
  if (stp == 0 && D0*Dp < 0){ lo=p0; hi=p0+h0; Dlo=D0; Dhi=Dp; have=1; }
  else if (stm == 0 && D0*Dm < 0){ lo=p0-h0; hi=p0; Dlo=Dm; Dhi=D0; have=1; }

  if (!have){
    double sgn;
    if (stp == 0 && stm == 0) sgn = (fabs(Dp) < fabs(Dm)) ? 1.0 : -1.0;
    else if (stp == 0)        sgn = 1.0;
    else if (stm == 0)        sgn = -1.0;
    else                      return 10;
    double pa = p0, Da = D0;
    double pb = p0 + sgn*h0, Db = (sgn>0)?Dp:Dm;
    double h = h0;
    while (fabs(pb-p0) <= span){
      if (Da*Db < 0){
        if (pa < pb){ lo=pa; hi=pb; Dlo=Da; Dhi=Db; } else { lo=pb; hi=pa; Dlo=Db; Dhi=Da; }
        have=1; break;
      }
      pa = pb; Da = Db;
      h *= 1.7;
      pb = pa + sgn*h;
      double Dn; int stn = disp(k, pb, &Dn); (*nev)++;
      if (stn != 0) return 10;                 /* family member has no return: UNRESOLVED */
      Db = Dn;
    }
    if (!have) return 10;
  }

  /* bisection with secant acceleration; D is strictly monotone in the parameter */
  for (int it=0; it<200; ++it){
    if (hi-lo <= ptol){ *pout = 0.5*(lo+hi); return 0; }
    double mid = 0.5*(lo+hi);
    double sec = lo + (hi-lo)*(-Dlo)/(Dhi-Dlo);
    double p = sec;
    if (!isfinite(p) || p <= lo + 0.02*(hi-lo) || p >= hi - 0.02*(hi-lo)) p = mid;
    double Dn; int stn = disp(k, p, &Dn); (*nev)++;
    if (stn != 0){ p = mid; stn = disp(k, p, &Dn); (*nev)++; if (stn != 0) return 12; }
    if (Dn == 0.0){ *pout = p; return 0; }
    if (Dlo*Dn < 0){ hi = p; Dhi = Dn; } else { lo = p; Dlo = Dn; }
  }
  *pout = 0.5*(lo+hi);
  return 12;
}

void ahcurve(int nsets, const double *coef, const double *evec, const double *dir,
             const double *p0, int mode,
             int ns, const double *svals,
             double *pout, int *status, int *nev,
             double ptol, double span,
             double rtol, double Rmax, double Tmax, long maxsteps){
#pragma omp parallel for schedule(dynamic,1)
  for (int i=0;i<nsets*ns;++i){
    int sIdx = i/ns;
    Ctx k;
    k.c = coef + 12*sIdx;
    k.e = evec ? (evec + 12*sIdx) : NULL;
    k.mode = mode;
    k.dx = dir[2*sIdx]; k.dy = dir[2*sIdx+1];
    k.s = svals[i];
    k.rtol = rtol; k.Rmax = Rmax; k.Tmax = Tmax; k.maxsteps = maxsteps;
    double p = NAN; int ne = 0;
    int st = solve_param(&k, p0[sIdx], span, ptol, &p, &ne);
    pout[i] = (st==0) ? p : NAN;
    status[i] = st;
    nev[i] = ne;
  }
}
