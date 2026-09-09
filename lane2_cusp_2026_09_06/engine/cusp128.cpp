// cusp128.cpp -- Lane 2 (cusp / swallowtail) return-map engine for the
// Cherkas-Artes-Llibre normal form, in binary128 (__float128).
//
//   xdot = 1 + x y
//   ydot = a00 + a10 x + a20 x^2 + a01 y + a11 x y + a y^2,
//   a00  = a01 + a11 - a10 - a20 - a          (so A=(1,-1) is a singular point)
//
// Translate u = x-1, v = y+1:
//   u' = -u + v + u v
//   v' = A u + B v + a20 u^2 + a11 u v + a v^2,
//        A = a10 + 2 a20 - a11,  B = a01 + a11 - 2a
//   trace  T  = B - 1   ( = V1 = a11 + a01 - 2a - 1 )
//   det    L  = -B - A = 2a - a01 - a10 - 2 a20        (independent of a11)
//   omega  w  = sqrt(L - T^2/4)                        (focus iff L>0, T^2<4L)
//
// Linear-normalising chart  u = xi,  v = k1 xi - w eta,  k1 = 1 + T/2 :
//   xi'  = (T/2) xi - w eta + N_u
//   eta' =  w xi  + (T/2)eta + (k1 N_u - N_v)/w
// Polar xi = rho cos th, eta = rho sin th :
//   rho' = (T/2) rho + alpha(th) rho^2
//   th'  = w + beta(th) rho
//   ==>  drho/dth = ( (T/2) rho + alpha rho^2 ) / ( w + beta rho )
// with, writing Cu = cos th, Cv = k1 cos th - w sin th,
//   P2 = Cu Cv,  Q2 = ( k1 P2 - a20 Cu^2 - a11 P2 - a Cv^2 ) / w
//   alpha = Cu P2 + Cv_s... (see code: alpha = c*P2 + s*Q2, beta = c*Q2 - s*P2)
//
// The section {y = -1, x > 1} is {v = 0, u > 0} = the fixed ray
//   th = th0 = atan2(k1, w)  (cos th0 > 0),   u = rho * w / n,  n = hypot(w,k1)
// so x = 1 + rho * w / n.  Return map: integrate th from th0 to th0 + 2*pi.
//
// Derivatives in the initial radius are carried as an order-3 jet
//   rho(th) = r0 + r1 e + r2 e^2 + r3 e^3,   e = d rho0,
// which is exact (no finite differences anywhere).  Then
//   D = r0 - rho0,  D' = r1 - 1,  D'' = 2 r2,  D''' = 6 r3.
//
// Integrator: Gragg-Bulirsch-Stoer (modified midpoint + Neville extrapolation
// in h^2), adaptive step, tolerance settable.  No tabulated RK coefficients.
//
// Build: g++ -O2 -o cusp128 cusp128.cpp -lquadmath
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <quadmath.h>

typedef __float128 qd;

static const int NJ = 4;                 // jet order 3  -> 4 coefficients

struct Jet { qd c[NJ]; };

static inline void jset(Jet&A, qd v){ A.c[0]=v; A.c[1]=0; A.c[2]=0; A.c[3]=0; }
static inline Jet jmul(const Jet&A,const Jet&B){
  Jet C;
  C.c[0]=A.c[0]*B.c[0];
  C.c[1]=A.c[0]*B.c[1]+A.c[1]*B.c[0];
  C.c[2]=A.c[0]*B.c[2]+A.c[1]*B.c[1]+A.c[2]*B.c[0];
  C.c[3]=A.c[0]*B.c[3]+A.c[1]*B.c[2]+A.c[2]*B.c[1]+A.c[3]*B.c[0];
  return C;
}
static inline Jet jdiv(const Jet&A,const Jet&B){
  Jet C; qd b0=B.c[0];
  C.c[0]= A.c[0]/b0;
  C.c[1]=(A.c[1]-C.c[0]*B.c[1])/b0;
  C.c[2]=(A.c[2]-C.c[0]*B.c[2]-C.c[1]*B.c[1])/b0;
  C.c[3]=(A.c[3]-C.c[0]*B.c[3]-C.c[1]*B.c[2]-C.c[2]*B.c[1])/b0;
  return C;
}

struct Sys {
  qd a, a20, a11, a01, a10;      // moduli
  qd A, B, T, L, w, k1, th0, nrm; // derived
  bool ok;
};

static Sys make_sys(qd a, qd a20, qd a11, qd a01, qd a10){
  Sys S; S.a=a; S.a20=a20; S.a11=a11; S.a01=a01; S.a10=a10;
  S.A = a10 + 2*a20 - a11;
  S.B = a01 + a11 - 2*a;
  S.T = S.B - 1;
  S.L = -S.B - S.A;                     // = 2a - a01 - a10 - 2 a20
  qd disc = S.L - S.T*S.T/4;
  S.ok = (S.L > 0) && (disc > 0);
  if(!S.ok){ S.w=0; S.k1=0; S.th0=0; S.nrm=0; return S; }
  S.w = sqrtq(disc);
  S.k1 = 1 + S.T/2;
  S.nrm = sqrtq(S.w*S.w + S.k1*S.k1);
  S.th0 = atan2q(S.k1, S.w);            // cos th0 = w/nrm > 0
  return S;
}

// coefficients alpha(th), beta(th) of the polar field
static inline void ab_coeffs(const Sys&S, qd th, qd&alpha, qd&beta){
  qd c = cosq(th), s = sinq(th);
  qd Cu = c;
  qd Cv = S.k1*c - S.w*s;
  qd P2 = Cu*Cv;
  qd Q2 = (S.k1*P2 - S.a20*Cu*Cu - S.a11*P2 - S.a*Cv*Cv)/S.w;
  alpha = c*P2 + s*Q2;
  beta  = c*Q2 - s*P2;
}

// f(th, rho) = ((T/2) rho + alpha rho^2) / (w + beta rho), as a jet
static inline Jet rhs(const Sys&S, qd th, const Jet&R, qd&den0_out){
  qd alpha, beta; ab_coeffs(S, th, alpha, beta);
  Jet R2 = jmul(R,R);
  Jet num, den;
  qd h = S.T/2;
  for(int i=0;i<NJ;i++) num.c[i] = h*R.c[i] + alpha*R2.c[i];
  for(int i=0;i<NJ;i++) den.c[i] = beta*R.c[i];
  den.c[0] += S.w;
  den0_out = den.c[0];
  return jdiv(num,den);
}

struct Result {
  bool ok;
  const char* why;
  qd D, D1, D2, D3;
  qd th_reached;
  qd min_den;      // min over path of |w + beta rho| (leading jet order)
  qd min_rho;
  qd x_sec;        // x on the section for rho0
  qd steps;
};

// ---- modified midpoint over [th, th+H] with n substeps ----
static bool midpoint(const Sys&S, qd th, const Jet&y, qd H, int n, Jet&out,
                     qd&min_den, qd&min_rho){
  qd h = H/n;
  Jet z0=y, z1, zm, f;
  qd d0;
  f = rhs(S, th, z0, d0);
  if(!(fabsq(d0)>0) || z0.c[0]<=0) return false;
  if(fabsq(d0)<min_den) min_den=fabsq(d0);
  if(z0.c[0]<min_rho) min_rho=z0.c[0];
  for(int i=0;i<NJ;i++) z1.c[i] = z0.c[i] + h*f.c[i];
  for(int m=1;m<n;m++){
    f = rhs(S, th+m*h, z1, d0);
    if(!(fabsq(d0)>0) || z1.c[0]<=0) return false;
    if(fabsq(d0)<min_den) min_den=fabsq(d0);
    if(z1.c[0]<min_rho) min_rho=z1.c[0];
    for(int i=0;i<NJ;i++){ qd t = z0.c[i] + 2*h*f.c[i]; z0.c[i]=z1.c[i]; z1.c[i]=t; }
  }
  f = rhs(S, th+H, z1, d0);
  if(!(fabsq(d0)>0) || z1.c[0]<=0) return false;
  if(fabsq(d0)<min_den) min_den=fabsq(d0);
  if(z1.c[0]<min_rho) min_rho=z1.c[0];
  for(int i=0;i<NJ;i++) out.c[i] = (z0.c[i] + z1.c[i] + h*f.c[i])/2;
  return true;
}

static const int KMAX = 14;
static const int NSEQ[KMAX] = {2,4,6,8,10,12,14,16,18,20,22,24,26,28};

// one adaptive GBS step; returns achieved step (may be < H) or 0 on failure
static bool gbs_step(const Sys&S, qd th, const Jet&y, qd Htry, qd tol,
                     Jet&yout, qd&Hdone, qd&Hnext, qd&min_den, qd&min_rho){
  qd H = Htry;
  for(int attempt=0; attempt<60; attempt++){
    Jet tab[KMAX]; qd x[KMAX];
    bool bad=false; int kconv=-1; Jet best; qd errconv=0;
    qd md=min_den, mr=min_rho, md_try, mr_try;
    for(int k=0;k<KMAX;k++){
      md_try=md; mr_try=mr;
      Jet raw;
      if(!midpoint(S, th, y, H, NSEQ[k], raw, md_try, mr_try)){ bad=true; break; }
      x[k] = (H/NSEQ[k])*(H/NSEQ[k]);
      tab[k]=raw;
      // Neville in x -> x=0
      for(int j=k-1;j>=0;j--){
        qd f = x[k]/(x[j]-x[k]);
        for(int i=0;i<NJ;i++)
          tab[j].c[i] = tab[j+1].c[i] + f*(tab[j+1].c[i]-tab[j].c[i]);
      }
      if(k>=2){
        qd err=0;
        for(int i=0;i<NJ;i++){
          qd sc = fabsq(tab[0].c[i]) + fabsq(y.c[i]) + (qd)1e-30Q;
          qd e = fabsq(tab[0].c[i]-tab[1].c[i])/sc;
          if(e>err) err=e;
        }
        if(err<tol){ kconv=k; best=tab[0]; errconv=err;
                     min_den=md_try; min_rho=mr_try; break; }
      }
    }
    if(bad || kconv<0){ H *= 0.5Q; if(fabsq(H) < (qd)1e-13Q) return false; continue; }
    yout=best; Hdone=H;
    // work-optimal-ish step control: extrapolated order at level k is ~2k+2
    qd fac;
    if(errconv <= 0) fac = 4;
    else {
      qd ex = (qd)1/(2*kconv+1);
      fac = 0.9Q*powq(tol/errconv, ex);
      if(fac > 4) fac = 4;
      if(fac < 0.2Q) fac = 0.2Q;
    }
    Hnext = H*fac;
    return true;
  }
  return false;
}

static Result integrate(const Sys&S, qd rho0, qd tol, int side){
  Result R; memset(&R,0,sizeof(R));
  R.ok=false; R.why="none"; R.min_den=1e300Q; R.min_rho=1e300Q; R.steps=0;
  if(!S.ok){ R.why="not_a_focus"; return R; }
  if(!(rho0>0)){ R.why="rho0_nonpositive"; return R; }
  R.x_sec = (side==0) ? (1 + rho0*S.w/S.nrm) : (1 - rho0*S.w/S.nrm);
  Jet y; y.c[0]=rho0; y.c[1]=1; y.c[2]=0; y.c[3]=0;
  qd th_start = S.th0 + (side ? M_PIq : (qd)0);
  qd th = th_start;
  qd th_end = th_start + 2*M_PIq;
  qd H = (2*M_PIq)/24;
  int nsteps=0;
  while(th < th_end){
    if(th + H > th_end) H = th_end - th;
    Jet yn; qd Hd, Hn;
    if(!gbs_step(S, th, y, H, tol, yn, Hd, Hn, R.min_den, R.min_rho)){
      R.why="return_lost"; R.th_reached = th - th_start; R.steps=nsteps; return R;
    }
    y = yn; th += Hd; H = Hn; nsteps++;
    if(nsteps>200000){ R.why="step_limit"; R.th_reached=th-th_start; return R; }
  }
  R.ok=true; R.why="ok"; R.th_reached = 2*M_PIq; R.steps=nsteps;
  R.D  = y.c[0] - rho0;
  R.D1 = y.c[1] - 1;
  R.D2 = 2*y.c[2];
  R.D3 = 6*y.c[3];
  return R;
}

static void pq(char*b,size_t n,qd v){ quadmath_snprintf(b,n,"%.34Qe",v); }

int main(int argc, char**argv){
  qd tol = 1e-28Q;
  for(int i=1;i<argc;i++){
    if(!strcmp(argv[i],"--tol") && i+1<argc) tol = strtoflt128(argv[++i],0);
  }
  char line[4096], b[8][80];
  while(fgets(line,sizeof(line),stdin)){
    char cmd[32];
    if(sscanf(line,"%31s",cmd)!=1) continue;
    if(!strcmp(cmd,"QUIT")) break;
    if(!strcmp(cmd,"D")){
      char s1[80],s2[80],s3[80],s4[80],s5[80],s6[80]; int side=0;
      int nf = sscanf(line,"%*s %79s %79s %79s %79s %79s %79s %d",
                      s1,s2,s3,s4,s5,s6,&side);
      if(nf<6){ printf("ERR parse\n"); fflush(stdout); continue; }
      if(nf<7) side=0;
      Sys S = make_sys(strtoflt128(s1,0),strtoflt128(s2,0),strtoflt128(s3,0),
                       strtoflt128(s4,0),strtoflt128(s5,0));
      qd rho0 = strtoflt128(s6,0);
      Result R = integrate(S,rho0,tol,side);
      if(!R.ok){
        char t1[80],t2[80]; pq(t1,80,R.th_reached); pq(t2,80,R.min_den);
        printf("FAIL %s %s %s\n", R.why, t1, t2);
      } else {
        pq(b[0],80,R.D); pq(b[1],80,R.D1); pq(b[2],80,R.D2); pq(b[3],80,R.D3);
        pq(b[4],80,R.min_den); pq(b[5],80,R.min_rho); pq(b[6],80,R.x_sec);
        pq(b[7],80,S.w);
        printf("OK %s %s %s %s %s %s %s %s %d\n",
               b[0],b[1],b[2],b[3],b[4],b[5],b[6],b[7],(int)R.steps);
      }
      fflush(stdout); continue;
    }
    if(!strcmp(cmd,"INFO")){
      char s1[80],s2[80],s3[80],s4[80],s5[80];
      if(sscanf(line,"%*s %79s %79s %79s %79s %79s",s1,s2,s3,s4,s5)!=5){
        printf("ERR parse\n"); fflush(stdout); continue; }
      Sys S = make_sys(strtoflt128(s1,0),strtoflt128(s2,0),strtoflt128(s3,0),
                       strtoflt128(s4,0),strtoflt128(s5,0));
      pq(b[0],80,S.T); pq(b[1],80,S.L); pq(b[2],80,S.w); pq(b[3],80,S.th0);
      pq(b[4],80,S.nrm);
      printf("%s T=%s L=%s w=%s th0=%s nrm=%s\n", S.ok?"OK":"BAD",
             b[0],b[1],b[2],b[3],b[4]);
      fflush(stdout); continue;
    }
    printf("ERR cmd\n"); fflush(stdout);
  }
  return 0;
}
