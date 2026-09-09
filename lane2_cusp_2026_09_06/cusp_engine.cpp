// cusp_engine.cpp — Lane 2 (cusp/swallowtail) return-map engine for the Cherkas normal form
//
//   xdot = P = 1 + x y
//   ydot = Q = a00 + a10 x + a20 x^2 + a01 y + a11 x y + a y^2,   a00 = a01 + a11 - a10 - a20 - a
//
// Focus A = (1,-1).  Section: the ray { y = -1, x > 1 }.
// Displacement D(x0) = (x-coordinate of first return to the ray) - x0.
//
// Method
// ------
// Taylor series time-stepping (the vector field is polynomial, so the series
// recurrences are exact), with the STATE CARRIED AS A DEGREE-3 JET in the initial
// section coordinate eps:   x(0) = x0 + eps,  y(0) = -1.
// The jet of the return time is solved from the jet equation y(tau) = -1, and
// R(eps) = x(tau(eps)) is then a degree-3 jet, giving
//     D    = R.c0 - x0
//     D_x  = R.c1 - 1
//     D_xx = 2 R.c2
//     D_xxx= 6 R.c3
// exactly (to integration accuracy) -- NOT by finite differences.  This is what
// makes an honest D_xxx possible at all.
//
// Templated on the float type: `long double` (64-bit mantissa) for sweeps,
// `__float128` (113-bit mantissa) for anything that decides a trigger.
//
// Protocol: a failed/guard-stopped return is UNRESOLVED, never "no cycle".

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <string>
#include <vector>
#include <quadmath.h>

// ---------------------------------------------------------------- traits ----
template <typename T> struct FT;

template <> struct FT<long double> {
    static long double fromstr(const char* s) { return strtold(s, nullptr); }
    static void tostr(long double v, char* b, int n) { snprintf(b, n, "%.21Le", v); }
    static long double fabs_(long double v) { return fabsl(v); }
    static long double pow_(long double a, long double b) { return powl(a, b); }
    static long double log_(long double v) { return logl(v); }
    static long double exp_(long double v) { return expl(v); }
    static long double sqrt_(long double v) { return sqrtl(v); }
    static bool finite_(long double v) { return std::isfinite((double)(v == v ? v : 1.0/0.0)) && v == v && fabsl(v) < 1e4000L; }
    static const char* name() { return "long double (64-bit mantissa)"; }
    static int    order()  { return 20; }
    static long double tol() { return 1e-21L; }
};

template <> struct FT<__float128> {
    static __float128 fromstr(const char* s) { return strtoflt128(s, nullptr); }
    static void tostr(__float128 v, char* b, int n) { quadmath_snprintf(b, n, "%.36Qe", v); }
    static __float128 fabs_(__float128 v) { return fabsq(v); }
    static __float128 pow_(__float128 a, __float128 b) { return powq(a, b); }
    static __float128 log_(__float128 v) { return logq(v); }
    static __float128 exp_(__float128 v) { return expq(v); }
    static __float128 sqrt_(__float128 v) { return sqrtq(v); }
    static bool finite_(__float128 v) { return finiteq(v) != 0; }
    static const char* name() { return "binary128 (113-bit mantissa)"; }
    static int    order()  { return 26; }
    static __float128 tol() { return strtoflt128("1e-32", nullptr); }
};

// ------------------------------------------------------------------- jet ----
// Truncated polynomial in eps of degree JETDEG (3 by default; 4 gives D_xxxx,
// which makes the distance to a swallow-tail scale-free: nu = D_xxx/(D_xxxx r0)).
#ifndef JETDEG
#define JETDEG 3
#endif
#define NJ (JETDEG + 1)

template <typename T>
struct Jet {
    T c[NJ];
    Jet() { for (int i = 0; i < NJ; i++) c[i] = T(0); }
    explicit Jet(T v) { for (int i = 0; i < NJ; i++) c[i] = T(0); c[0] = v; }
};

template <typename T> inline Jet<T> operator+(const Jet<T>& a, const Jet<T>& b) {
    Jet<T> r; for (int i = 0; i < NJ; i++) r.c[i] = a.c[i] + b.c[i]; return r;
}
template <typename T> inline Jet<T> operator-(const Jet<T>& a, const Jet<T>& b) {
    Jet<T> r; for (int i = 0; i < NJ; i++) r.c[i] = a.c[i] - b.c[i]; return r;
}
template <typename T> inline Jet<T> operator*(T s, const Jet<T>& a) {
    Jet<T> r; for (int i = 0; i < NJ; i++) r.c[i] = s * a.c[i]; return r;
}
template <typename T> inline Jet<T> operator*(const Jet<T>& a, const Jet<T>& b) {
    Jet<T> r;
    for (int i = 0; i < NJ; i++) {
        T s = T(0);
        for (int j = 0; j <= i; j++) s += a.c[j] * b.c[i - j];
        r.c[i] = s;
    }
    return r;
}
// series inverse (constant term must be nonzero)
template <typename T> inline Jet<T> inv(const Jet<T>& a) {
    Jet<T> r;
    T i0 = T(1) / a.c[0];
    r.c[0] = i0;
    for (int i = 1; i < NJ; i++) {
        T s = T(0);
        for (int j = 1; j <= i; j++) s += a.c[j] * r.c[i - j];
        r.c[i] = -s * i0;
    }
    return r;
}
template <typename T> inline T jnorm(const Jet<T>& a) {
    T m = FT<T>::fabs_(a.c[0]);
    for (int i = 1; i < NJ; i++) { T v = FT<T>::fabs_(a.c[i]); if (v > m) m = v; }
    return m;
}

// ------------------------------------------------------------------ field ---
template <typename T>
struct Params { T a, a20, a11, a01, a10, a00; };

template <typename T>
Params<T> make_params(T a, T a20, T a11, T a01, T a10) {
    Params<T> p; p.a = a; p.a20 = a20; p.a11 = a11; p.a01 = a01; p.a10 = a10;
    p.a00 = a01 + a11 - a10 - a20 - a;
    return p;
}

// ------------------------------------------------------------- integrator ---
enum Status { OK = 0, GUARD_BOUND, GUARD_TIME, GUARD_STEPS, GUARD_STEPSIZE, GUARD_NAN, NO_RETURN };

template <typename T>
struct ReturnResult {
    Status st;
    Jet<T> R;      // jet of the returned x-coordinate
    T T_return;    // return time (eps = 0)
    long nsteps;
    T ymin_gap;    // smallest |Q| seen on the section crossing (transversality diagnostic)
};

template <typename T>
class Engine {
public:
    Params<T> p;
    int N;              // Taylor order in time
    T tol;              // local-truncation target
    T hmax;
    T bound;            // coordinate guard
    T tmax;             // physical-time guard
    long maxsteps;

    Engine(const Params<T>& pp) : p(pp) {
        N = FT<T>::order();
        tol = FT<T>::tol();
        if (const char* e = getenv("CUSP_ORDER")) N = atoi(e);
        if (const char* e = getenv("CUSP_TOL"))   tol = FT<T>::fromstr(e);
        hmax = T(0.25);
        if (const char* e = getenv("CUSP_HMAX"))  hmax = FT<T>::fromstr(e);
        bound = T(1e7);
        tmax = T(400);
        maxsteps = 400000;
    }

    // Taylor coefficients of the jet-valued solution on the current step.
    // X[k], Y[k] for k = 0..N ; solution ~ sum X[k] t^k.
    void taylor(const Jet<T>& x0, const Jet<T>& y0, std::vector<Jet<T>>& X, std::vector<Jet<T>>& Y) const {
        X.assign(N + 1, Jet<T>());
        Y.assign(N + 1, Jet<T>());
        X[0] = x0; Y[0] = y0;
        for (int k = 0; k <= N - 1; k++) {
            Jet<T> XY, X2, Y2;
            for (int j = 0; j <= k; j++) {
                XY = XY + X[j] * Y[k - j];
                X2 = X2 + X[j] * X[k - j];
                Y2 = Y2 + Y[j] * Y[k - j];
            }
            Jet<T> dx = XY;
            if (k == 0) dx.c[0] += T(1);
            Jet<T> dy = p.a10 * X[k] + p.a20 * X2 + p.a01 * Y[k] + p.a11 * XY + p.a * Y2;
            if (k == 0) dy.c[0] += p.a00;
            T ik = T(1) / T(k + 1);
            X[k + 1] = ik * dx;
            Y[k + 1] = ik * dy;
        }
    }

    T stepsize(const std::vector<Jet<T>>& X, const std::vector<Jet<T>>& Y) const {
        T mN  = jnorm(X[N]);     { T v = jnorm(Y[N]);     if (v > mN)  mN  = v; }
        T mN1 = jnorm(X[N - 1]); { T v = jnorm(Y[N - 1]); if (v > mN1) mN1 = v; }
        T h = hmax;
        if (mN > T(0))  { T r = FT<T>::pow_(tol / mN,  T(1) / T(N));     if (r < h) h = r; }
        if (mN1 > T(0)) { T r = FT<T>::pow_(tol / mN1, T(1) / T(N - 1)); if (r < h) h = r; }
        h *= T(0.85);
        if (h < T(1e-14)) h = T(1e-14);
        return h;
    }

    static Jet<T> horner(const std::vector<Jet<T>>& C, T t) {
        Jet<T> s = C.back();
        for (int k = (int)C.size() - 2; k >= 0; k--) s = (t * s) + C[k];
        return s;
    }
    // evaluate only the eps^0 component (cheap)
    static T horner0(const std::vector<Jet<T>>& C, T t) {
        T s = C.back().c[0];
        for (int k = (int)C.size() - 2; k >= 0; k--) s = t * s + C[k].c[0];
        return s;
    }

    // side = +1 : section is the ray { y = -1, x > 1 }
    // side = -1 : section is the ray { y = -1, x < 1 }  (the other side of the focus)
    int side = 1;

    // First return to the section, starting from (x0 + eps, -1).
    ReturnResult<T> ret(T x0) const {
        ReturnResult<T> res;
        res.st = NO_RETURN; res.nsteps = 0; res.T_return = T(0); res.ymin_gap = T(0);

        Jet<T> X0, Y0;
        X0.c[0] = x0; X0.c[1] = T(1);
        Y0.c[0] = T(-1);

        std::vector<Jet<T>> X, Y;
        T t = T(0);
        T prev = T(0);          // value of (y+1) at the start of the current segment
        bool armed = false;     // becomes true once we have left the section
        long steps = 0;

        while (true) {
            if (steps++ > maxsteps) { res.st = GUARD_STEPS; res.nsteps = steps; return res; }
            taylor(X0, Y0, X, Y);
            if (!FT<T>::finite_(X[N].c[0]) || !FT<T>::finite_(Y[N].c[0])) { res.st = GUARD_NAN; return res; }
            T h = stepsize(X, Y);
            if (h <= T(1e-14)) { res.st = GUARD_STEPSIZE; return res; }

            // scan the step for a crossing of y = -1 (on the eps=0 component)
            const int NSUB = 6;
            T tlo = T(0);
            if (!armed) {
                // leave the section first: skip a tiny initial slice
                tlo = h * T(1e-9);
                prev = horner0(Y, tlo) + T(1);
                if (FT<T>::fabs_(prev) > T(0)) armed = true;
            }
            for (int i = 1; i <= NSUB; i++) {
                T ta = tlo + (h - tlo) * (T(i - 1) / T(NSUB));
                T tb = tlo + (h - tlo) * (T(i)     / T(NSUB));
                T va = (i == 1) ? prev : horner0(Y, ta) + T(1);
                T vb = horner0(Y, tb) + T(1);
                if ((va > T(0) && vb < T(0)) || (va < T(0) && vb > T(0))) {
                    // bisect + Newton-polish on the eps=0 component
                    T lo = ta, hi = tb, flo = va;
                    for (int it = 0; it < 200; it++) {
                        T mid = T(0.5) * (lo + hi);
                        T fm = horner0(Y, mid) + T(1);
                        if ((flo > T(0)) == (fm > T(0))) { lo = mid; flo = fm; } else hi = mid;
                        if (hi - lo < FT<T>::fabs_(hi) * T(1e-34) + T(1e-300)) break;
                    }
                    T tr = T(0.5) * (lo + hi);
                    // Newton polish
                    for (int it = 0; it < 60; it++) {
                        T f = horner0(Y, tr) + T(1);
                        // derivative of the eps=0 series
                        T d = T(0);
                        for (int k = (int)Y.size() - 1; k >= 1; k--) d = tr * d + T(k) * Y[k].c[0];
                        if (d == T(0)) break;
                        T dt = f / d;
                        tr -= dt;
                        if (FT<T>::fabs_(dt) <= FT<T>::fabs_(tr) * T(1e-34)) break;
                    }
                    if (tr < tlo) tr = tlo;
                    if (tr > h) tr = h;
                    T xr = horner0(X, tr);
                    if (((side > 0 && xr > T(1)) || (side < 0 && xr < T(1))) && (t + tr) > T(1e-9)) {
                        // ---- solve the JET equation y(tau) = -1 ----
                        Jet<T> tau; tau.c[0] = tr;
                        for (int it = 0; it < 12; it++) {
                            // f(tau) and f'(tau) in jet arithmetic
                            Jet<T> f = Y[Y.size() - 1];
                            for (int k = (int)Y.size() - 2; k >= 0; k--) f = (tau * f) + Y[k];
                            f.c[0] += T(1);
                            Jet<T> d;  // derivative series
                            {
                                int n = (int)Y.size() - 1;
                                d = T(n) * Y[n];
                                for (int k = n - 1; k >= 1; k--) d = (tau * d) + (T(k) * Y[k]);
                            }
                            if (FT<T>::fabs_(d.c[0]) == T(0)) break;
                            Jet<T> dt = f * inv(d);
                            tau = tau - dt;
                            if (jnorm(dt) <= FT<T>::fabs_(tau.c[0]) * T(1e-34)) break;
                        }
                        Jet<T> Rj = X[X.size() - 1];
                        for (int k = (int)X.size() - 2; k >= 0; k--) Rj = (tau * Rj) + X[k];
                        res.st = OK;
                        res.R = Rj;
                        res.T_return = t + tau.c[0];
                        res.nsteps = steps;
                        // transversality diagnostic: |dy/dt| at the crossing
                        {
                            T d = T(0);
                            for (int k = (int)Y.size() - 1; k >= 1; k--) d = tau.c[0] * d + T(k) * Y[k].c[0];
                            res.ymin_gap = FT<T>::fabs_(d);
                        }
                        return res;
                    }
                }
                prev = vb;
            }

            X0 = horner(X, h);
            Y0 = horner(Y, h);
            t += h;
            armed = true;
            if (FT<T>::fabs_(X0.c[0]) > bound || FT<T>::fabs_(Y0.c[0]) > bound) { res.st = GUARD_BOUND; return res; }
            if (t > tmax) { res.st = GUARD_TIME; return res; }
            if (!FT<T>::finite_(X0.c[0]) || !FT<T>::finite_(Y0.c[0])) { res.st = GUARD_NAN; return res; }
        }
    }
};

static const char* stname(Status s) {
    switch (s) {
        case OK: return "OK";
        case GUARD_BOUND: return "GUARD_BOUND";
        case GUARD_TIME: return "GUARD_TIME";
        case GUARD_STEPS: return "GUARD_STEPS";
        case GUARD_STEPSIZE: return "GUARD_STEPSIZE";
        case GUARD_NAN: return "GUARD_NAN";
        default: return "NO_RETURN";
    }
}

// ------------------------------------------------------------------ main ----
// stdin protocol, one command per line:
//   D <a> <a20> <a11> <a01> <a10> <x0>
//     -> "OK <D> <Dx> <Dxx> <Dxxx> <T> <nsteps> <transv>"  or  "<STATUS>"
//   Q            quit
// Numbers are decimal strings, parsed at full working precision.

template <typename T>
int run() {
    char line[4096];
    setvbuf(stdout, nullptr, _IOLBF, 0);
    char buf[256];
    while (fgets(line, sizeof(line), stdin)) {
        if (line[0] == 'Q') break;
        if (line[0] == 'V') { printf("ENGINE cusp_engine %s order=%d jetdeg=%d\n", FT<T>::name(), FT<T>::order(), JETDEG); continue; }
        if (line[0] != 'D' && line[0] != 'L') { printf("ERR bad command\n"); continue; }
        int sd = (line[0] == 'L') ? -1 : 1;
        char sa[128], s20[128], s11[128], s01[128], s10[128], sx[128];
        if (sscanf(line + 1, "%127s %127s %127s %127s %127s %127s", sa, s20, s11, s01, s10, sx) != 6) {
            printf("ERR bad args\n"); continue;
        }
        Params<T> p = make_params<T>(FT<T>::fromstr(sa), FT<T>::fromstr(s20), FT<T>::fromstr(s11),
                                     FT<T>::fromstr(s01), FT<T>::fromstr(s10));
        T x0 = FT<T>::fromstr(sx);
        Engine<T> eng(p);
        eng.side = sd;
        ReturnResult<T> r = eng.ret(x0);
        if (r.st != OK) { printf("%s\n", stname(r.st)); continue; }
        T D    = r.R.c[0] - x0;
        T Dx   = r.R.c[1] - T(1);
        T Dxx  = T(2) * r.R.c[2];
        T Dxxx = T(6) * r.R.c[3];
#if JETDEG >= 4
        T Dxxxx = T(24) * r.R.c[4];
        const int NV = 7;
        T vals[7] = { D, Dx, Dxx, Dxxx, r.T_return, r.ymin_gap, Dxxxx };
#else
        const int NV = 6;
        T vals[6] = { D, Dx, Dxx, Dxxx, r.T_return, r.ymin_gap };
#endif
        std::string out = "OK";
        for (int i = 0; i < NV; i++) { FT<T>::tostr(vals[i], buf, sizeof(buf)); out += " "; out += buf; }
        out += " " + std::to_string(r.nsteps);
        printf("%s\n", out.c_str());
    }
    return 0;
}

int main(int argc, char** argv) {
    bool quad = true;
    for (int i = 1; i < argc; i++) if (!strcmp(argv[i], "--ld")) quad = false;
    if (quad) return run<__float128>();
    return run<long double>();
}
