/* lane1/retmap1.c  --  return-map / Andronov-Hopf-curve engine for Lane 1.
 *
 * Quadratic planar field X = (P,Q) given by 12 coefficients in the order
 *     P: 1, x, y, x^2, xy, y^2      Q: 1, x, y, x^2, xy, y^2
 * The driver re-expands the field about the focus F and passes the LOCAL
 * 10-vector (p1,p2,p3,p4,p5,q1,q2,q3,q4,q5) with
 *     Pl(u,v) = p1 u + p2 v + p3 u^2 + p4 u v + p5 v^2
 *     Ql(u,v) = q1 u + q2 v + q3 u^2 + q4 u v + q5 v^2
 * so the constant term is exactly zero at the focus (REVIEW_engine bug A2:
 * never integrate in a global chart whose origin is far from the orbit).
 *
 * Uniformly rotated family:
 *     f = ( Pl cos b - Ql sin b ,  Pl sin b + Ql cos b )
 *
 * Section: the ray from the focus with direction e = (cos phi, sin phi),
 * parametrised by s > 0.  h(u) = u.n with n = (-sin phi, cos phi).
 * A return is the first crossing h = 0 with u.e > 0 after at least half a
 * turn; the winding angle theta is carried as a third state and gated.
 *
 * Integrator: Dormand-Prince 5(4) with FSAL, PI-free standard step control,
 * atol = ATOL_REL * s (bug A1: never an absolute floor).  Crossing located by
 * Newton on the step length using full DP5 steps, so the crossing point is
 * accurate to the local truncation error, not to an interpolant.
 *
 * status codes: 0 ok, 1 escaped (|u|>Rmax), 2 t>Tmax, 3 step limit,
 *               4 winding gate failed, 5 non-finite, 6 shrunk to zero,
 *               7 bracket failed (beta solve), 8 root solve failed
 */

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

/* ---- precision layer -------------------------------------------------- */
#ifdef LANE1_QUAD
#include <quadmath.h>
typedef __float128 real;
#define RSQRT sqrtq
#define RFABS fabsq
#define RCOS  cosq
#define RSIN  sinq
#define RPOW  powq
#define RFMAX fmaxq
#define RISFIN(x) (!isnanq(x) && !isinfq(x))
#define REPS  FLT128_EPSILON
#define RTINY 1e-4000Q
#define ATOL_REL 1e-30Q
#define XTOL 1e-30Q
#define SUF(n) n##_q
#else
typedef double real;
#define RSQRT sqrt
#define RFABS fabs
#define RCOS  cos
#define RSIN  sin
#define RPOW  pow
#define RFMAX fmax
#define RISFIN(x) isfinite(x)
#define REPS  2.220446049250313e-16
#define RTINY 1e-300
#define ATOL_REL 1e-16
#define XTOL 1e-15
#define SUF(n) n
#endif


#define TWO_PI ((real)6.283185307179586476925286766559)

typedef struct {
    real p1,p2,p3,p4,p5;
    real q1,q2,q3,q4,q5;
    real cb, sb;
} Fld;

static inline void rhs(const Fld *F, const real *y, real *dy)
{
    const real u = y[0], v = y[1];
    const real P = F->p1*u + F->p2*v + F->p3*u*u + F->p4*u*v + F->p5*v*v;
    const real Q = F->q1*u + F->q2*v + F->q3*u*u + F->q4*u*v + F->q5*v*v;
    const real fx = P*F->cb - Q*F->sb;
    const real fy = P*F->sb + Q*F->cb;
    dy[0] = fx;
    dy[1] = fy;
    dy[2] = (u*fy - v*fx) / (u*u + v*v);
}

/* ---- Dormand-Prince 5(4) ---- */
static const real A21= 1.0/5.0;
static const real A31= 3.0/40.0,        A32= 9.0/40.0;
static const real A41= 44.0/45.0,       A42= -56.0/15.0,     A43= 32.0/9.0;
static const real A51= 19372.0/6561.0,  A52= -25360.0/2187.0,
                    A53= 64448.0/6561.0,  A54= -212.0/729.0;
static const real A61= 9017.0/3168.0,   A62= -355.0/33.0,
                    A63= 46732.0/5247.0,  A64= 49.0/176.0,     A65= -5103.0/18656.0;
static const real B1 = 35.0/384.0,      B3 = 500.0/1113.0,   B4 = 125.0/192.0,
                    B5 = -2187.0/6784.0,  B6 = 11.0/84.0;
static const real E1 = 35.0/384.0    - 5179.0/57600.0,
                    E3 = 500.0/1113.0  - 7571.0/16695.0,
                    E4 = 125.0/192.0   - 393.0/640.0,
                    E5 = -2187.0/6784.0+ 92097.0/339200.0,
                    E6 = 11.0/84.0     - 187.0/2100.0,
                    E7 = -1.0/40.0;

/* one DP5 step; k1 must hold f(y).  Writes yn (5th order) and err vector. */
static void dp5_step(const Fld *F, const real *y, const real *k1, real h,
                     real *yn, real *k7, real *errv)
{
    real t[3], k2[3], k3[3], k4[3], k5[3], k6[3];
    int i;
    for (i=0;i<3;i++) t[i] = y[i] + h*(A21*k1[i]);
    rhs(F,t,k2);
    for (i=0;i<3;i++) t[i] = y[i] + h*(A31*k1[i]+A32*k2[i]);
    rhs(F,t,k3);
    for (i=0;i<3;i++) t[i] = y[i] + h*(A41*k1[i]+A42*k2[i]+A43*k3[i]);
    rhs(F,t,k4);
    for (i=0;i<3;i++) t[i] = y[i] + h*(A51*k1[i]+A52*k2[i]+A53*k3[i]+A54*k4[i]);
    rhs(F,t,k5);
    for (i=0;i<3;i++) t[i] = y[i] + h*(A61*k1[i]+A62*k2[i]+A63*k3[i]+A64*k4[i]+A65*k5[i]);
    rhs(F,t,k6);
    for (i=0;i<3;i++) yn[i] = y[i] + h*(B1*k1[i]+B3*k3[i]+B4*k4[i]+B5*k5[i]+B6*k6[i]);
    rhs(F,yn,k7);
    for (i=0;i<3;i++)
        errv[i] = h*(E1*k1[i]+E3*k3[i]+E4*k4[i]+E5*k5[i]+E6*k6[i]+E7*k7[i]);
}

typedef struct {
    real phi, rtol, Tmax, Rmax;
    long   nstep;
} Opt;

/* Full return from s along the ray.  Returns status; *R = return coordinate. */
static int ret_once(const Fld *F, const Opt *O, real s, real *R, real *Tp)
{
    const real ce = RCOS(O->phi), se = RSIN(O->phi);
    /* n = (-se, ce) ; h = -se*u + ce*v */
    real y[3], k1[3], k7[3], yn[3], errv[3];
    real t = 0.0, h;
    long steps = 0;
    const real atol = ATOL_REL * s + RTINY;

    y[0] = s*ce; y[1] = s*se; y[2] = 0.0;
    rhs(F, y, k1);
    {
        real sp = RSQRT(k1[0]*k1[0] + k1[1]*k1[1]);
        if (!(sp > 0.0)) return 6;
        h = (real)1e-4 * s / sp;
        if (!(h > 0.0)) return 6;
    }

    real hprev = -se*y[0] + ce*y[1];   /* == 0 at start */
    /* first sub-step: take the sign of h just after departure */
    int have_sign = 0;
    real sgn = 0.0;

    while (1) {
        if (++steps > O->nstep) return 3;
        if (t > O->Tmax) return 2;
        if (h > O->Tmax - t) h = O->Tmax - t + 1e-30;

        dp5_step(F, y, k1, h, yn, k7, errv);

        real err = 0.0;
        for (int i=0;i<3;i++) {
            real sc = atol + O->rtol * RFMAX(RFABS(y[i]), RFABS(yn[i]));
            if (i==2) sc = (real)1e-12 + O->rtol*RFMAX((real)1.0, RFABS(yn[2]));
            real e = errv[i]/sc;
            err += e*e;
        }
        err = RSQRT(err/3.0);

        if (!(err <= 1.0)) {                       /* NaN-safe reject */
            h *= RISFIN(err) ? RFMAX((real)0.2, (real)0.9*RPOW(err, (real)-0.2)) : 0.1;
            if (!(RFABS(h) > REPS*((real)1.0+t))) return 3;
            continue;
        }

        /* accepted */
        real hn = -se*yn[0] + ce*yn[1];
        real rr = RSQRT(yn[0]*yn[0] + yn[1]*yn[1]);
        if (!RISFIN(rr) || !RISFIN(hn)) return 5;
        if (rr > O->Rmax) return 1;

        if (!have_sign) {
            if (RFABS(hn) > (real)1e-13*s) { sgn = hn > 0 ? 1.0 : -1.0; have_sign = 1; }
        } else if (hprev*hn < 0.0 || hn == 0.0) {
            /* candidate crossing inside [t, t+h]; refine by Newton on tau */
            real tau = h * hprev/(hprev - hn);
            real yc[3], kc[3], ec[3], k7c[3];
            int ok = 0;
            for (int it=0; it<40; it++) {
                if (!(tau > 0.0)) tau = 0.5*h;
                if (tau > h) tau = h;
                dp5_step(F, y, k1, tau, yc, k7c, ec);
                real hc = -se*yc[0] + ce*yc[1];
                rhs(F, yc, kc);
                real dh = -se*kc[0] + ce*kc[1];
                if (!RISFIN(hc) || !RISFIN(dh) || dh == 0.0) break;
                real d = hc/dh;
                tau -= d;
                if (RFABS(hc) <= XTOL*((real)1.0 + RSQRT(yc[0]*yc[0]+yc[1]*yc[1]))) { ok = 1; break; }
                if (RFABS(d) < XTOL*(real)0.01*((real)1.0+RFABS(tau))) { ok = 1; break; }
            }
            if (ok) {
                real proj = yc[0]*ce + yc[1]*se;
                real th = yc[2];
                if (proj > 0.0 && RFABS(th) > 3.0) {
                    if (RFABS(RFABS(th) - TWO_PI) > 0.7) return 4;
                    *R = proj;
                    if (Tp) *Tp = t + tau;
                    return 0;
                }
            }
        }

        /* advance */
        t += h;
        memcpy(y, yn, 3*sizeof(real));
        memcpy(k1, k7, 3*sizeof(real));
        hprev = hn;
        if (RFABS(y[2]) > 4.0*TWO_PI) return 4;

        real fac = (real)0.9*RPOW(err > (real)1e-12 ? err : (real)1e-12, (real)-0.2);
        if (fac > 5.0) fac = 5.0;
        if (fac < 0.2) fac = 0.2;
        h *= fac;
    }
}

static void fill_fld(Fld *F, const double *loc10, double b)
{
    F->p1=loc10[0]; F->p2=loc10[1]; F->p3=loc10[2]; F->p4=loc10[3]; F->p5=loc10[4];
    F->q1=loc10[5]; F->q2=loc10[6]; F->q3=loc10[7]; F->q4=loc10[8]; F->q5=loc10[9];
    F->cb = RCOS(b); F->sb = RSIN(b);
}


/* ---------- exported: displacement over a grid of s at fixed b ---------- */
void SUF(d_curve)(const double *loc10, double phi, const double *sarr, int n,
             double b, double rtol, double Tmax, double Rmax, long nstep,
             double *Dout, int *stout, double *Tout)
{
    Fld F; fill_fld(&F, loc10, b);
    Opt O; O.phi=(real)phi; O.rtol=(real)rtol; O.Tmax=(real)Tmax;
    O.Rmax=(real)Rmax; O.nstep=nstep;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1)
#endif
    for (int i=0;i<n;i++) {
        real R = 0.0, T = 0.0;
        int st = ret_once(&F, &O, (real)sarr[i], &R, &T);
        stout[i] = st;
        Dout[i]  = (st==0) ? (double)(R - (real)sarr[i]) : NAN;
        if (Tout) Tout[i] = (st==0) ? (double)T : NAN;
    }
}

/* full return coordinate R itself, in the working precision, as a double */
void SUF(returns)(const double *loc10, double phi, const double *sarr, int n,
             double b, double rtol, double Tmax, double Rmax, long nstep,
             double *Rout, int *stout)
{
    Fld F; fill_fld(&F, loc10, b);
    Opt O; O.phi=(real)phi; O.rtol=(real)rtol; O.Tmax=(real)Tmax;
    O.Rmax=(real)Rmax; O.nstep=nstep;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1)
#endif
    for (int i=0;i<n;i++) {
        real R = 0.0, T = 0.0;
        int st = ret_once(&F, &O, (real)sarr[i], &R, &T);
        stout[i] = st;
        Rout[i]  = (st==0) ? (double)R : NAN;
    }
}
/* ---------- single displacement (serial, for the beta solve) ---------- */
static int Dval(const double *loc10, const Opt *O, real s, real b, real *D)
{
    Fld F; fill_fld(&F, loc10, b);
    real R;
    int st = ret_once(&F, O, s, &R, NULL);
    if (st) return st;
    *D = R - s;
    return 0;
}

/* Brent root of D(s,.) on a sign-changing bracket [ba,bb]. */
static int brent_b(const double *loc10, const Opt *O, real s,
                   real ba, real Da, real bb, real Db,
                   real btol, real *broot, real *nfev)
{
    real a=ba, fa=Da, bcur=bb, fb=Db, c=a, fc=fa, d=bcur-a, e=d;
    for (int it=0; it<200; it++) {
        if (fb*fc > 0.0) { c=a; fc=fa; d=bcur-a; e=d; }
        if (RFABS(fc) < RFABS(fb)) { a=bcur; bcur=c; c=a; fa=fb; fb=fc; fc=fa; }
        real tol1 = (real)2.0*REPS*RFABS(bcur) + (real)0.5*btol;
        real xm = 0.5*(c-bcur);
        if (RFABS(xm) <= tol1 || fb == 0.0) { *broot=bcur; return 0; }
        if (RFABS(e) >= tol1 && RFABS(fa) > RFABS(fb)) {
            real p,q,r,ss = fb/fa;
            if (a == c) { p = 2.0*xm*ss; q = 1.0-ss; }
            else {
                q = fa/fc; r = fb/fc;
                p = ss*(2.0*xm*q*(q-r) - (bcur-a)*(r-1.0));
                q = (q-1.0)*(r-1.0)*(ss-1.0);
            }
            if (p > 0.0) q = -q;
            p = RFABS(p);
            real min1 = 3.0*xm*q - RFABS(tol1*q), min2 = RFABS(e*q);
            if (2.0*p < (min1 < min2 ? min1 : min2)) { e=d; d=p/q; }
            else { d=xm; e=d; }
        } else { d=xm; e=d; }
        a=bcur; fa=fb;
        if (RFABS(d) > tol1) bcur += d;
        else bcur += (xm > 0.0 ? tol1 : -tol1);
        real fnew;
        int st = Dval(loc10, O, s, bcur, &fnew);
        *nfev += 1.0;
        if (st) return 8;
        fb = fnew;
    }
    return 8;
}

/* ---------- exported: beta*(s) over a grid ---------- */
/* dirhint: +1 if D increases with b, -1 if it decreases, 0 = probe.
 * bmax: largest |b| searched.  btol: tolerance in b. */
void SUF(betastar)(const double *loc10, double phi, const double *sarr, int n,
              double rtol, double Tmax, double Rmax, long nstep,
              double bmax, double btol, int dirhint,
              double *bout, int *stout, double *d0out, double *nfout)
{
    Opt O; O.phi=(real)phi; O.rtol=(real)rtol; O.Tmax=(real)Tmax; O.Rmax=(real)Rmax; O.nstep=nstep;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1)
#endif
    for (int i=0;i<n;i++) {
        real s = (real)sarr[i], D0, nfev = 1.0;
        int st = Dval(loc10, &O, s, 0.0, &D0);
        d0out[i] = NAN; bout[i] = NAN; stout[i] = st; nfout[i] = (double)nfev;
        if (st) continue;
        d0out[i] = (double)D0;
        if (D0 == 0.0) { bout[i] = 0.0; stout[i] = 0; continue; }

        int dir = dirhint;
        if (dir == 0) {
            real Dp, Dm; real eps = (real)1e-4;
            int s1 = Dval(loc10,&O,s, eps,&Dp);
            int s2 = Dval(loc10,&O,s,-eps,&Dm);
            nfev += 2.0;
            if (s1 || s2) { stout[i] = 7; nfout[i]=(double)nfev; continue; }
            dir = (Dp > Dm) ? 1 : -1;
        }
        /* march in the direction that drives D to zero */
        real sgnstep = (D0 > (real)0.0 ? (real)-1.0 : (real)1.0) * (real)dir;
        real ba = 0.0, Da = D0, bb, Db;
        real step = (real)1e-3;
        int bracketed = 0, fail = 0;
        while (step <= (real)bmax) {
            bb = sgnstep*step;
            int stb = Dval(loc10,&O,s,bb,&Db);
            nfev += 1.0;
            if (stb) { fail = 1; break; }
            if (Da*Db < 0.0) { bracketed = 1; break; }
            ba = bb; Da = Db;
            step *= 2.0;
        }
        nfout[i] = (double)nfev;
        if (!bracketed) { stout[i] = fail ? 7 : 7; continue; }
        real broot;
        int stR = brent_b(loc10,&O,s,ba,Da,bb,Db,(real)btol,&broot,&nfev);
        nfout[i] = (double)nfev;
        if (stR) { stout[i] = 8; continue; }
        bout[i] = (double)broot; stout[i] = 0;
    }
}
