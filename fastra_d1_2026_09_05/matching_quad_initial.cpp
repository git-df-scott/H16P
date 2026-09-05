// Binary128 two-sided angular matching, with y=sqrt(m)*Y to condition large m.
// Uses archived modified-midpoint extrapolation, here only w and its sensitivity.
// F=A-B has the same sign as B^{-1}(A)-z where both angular half maps exist.
// Numerical (not interval) evidence. An angular-chart failure is unresolved.
#include <quadmath.h>
#include <array>
#include <algorithm>
#include <stdexcept>
#include <string>
using R=__float128;using A=std::array<R,2>;
struct Engine {
 R c,m,beta,sm,dir,tol;long evaluations;
 A rhs(R s,const A&v){
  ++evaluations;if(evaluations>2000000)throw 3;
  R C=cosq(dir*s),S=-sinq(dir*s),E=expq(v[0]);
  R h1=beta/sm*S*S,g1=-1+beta/sm*C*S;
  R p2=C*C/sm+C*S,q2=-10/m*C*C+2.2Q/sm*C*S+c*S*S;
  R h2=C*p2+S*q2,g2=C*q2-S*p2,H=h1+E*h2,G=g1+E*g2;
  if(!(G<0))throw 4;
  return A{-dir*H/G,-dir*E*(h2*g1-h1*g2)/(G*G)};
 }
 A mid(R s,A y,R H,int n){
  R h=H/n;A old=y,cur=y,f=rhs(s,y);
  for(int j=0;j<2;j++)cur[j]+=h*f[j];
  for(int i=1;i<n;i++){f=rhs(s+i*h,cur);A nxt;for(int j=0;j<2;j++)nxt[j]=old[j]+2*h*f[j];old=cur;cur=nxt;}
  f=rhs(s+H,cur);A out;for(int j=0;j<2;j++)out[j]=(old[j]+cur[j]+h*f[j])/2;return out;
 }
 A integrate(R z,R direction){
  dir=direction;R s=0,end=acosq(-1.Q),H=.01Q;A y{z,0};int steps=0;
  while(s<end){
   if(++steps>100000)throw 3;H=std::min(H,end-s);if(H<1e-26Q)throw 4;
   std::array<A,9>tab;R err=1e100Q;A out;bool ok=false;
   try {
    for(int k=0;k<8;k++){
     int n=2*(k+1);tab[k]=mid(s,y,H,n);
     for(int j=k-1;j>=0;j--){R ratio=(R)n/(2*(j+1)),factor=ratio*ratio-1;for(int v=0;v<2;v++)tab[j][v]=tab[j+1][v]+(tab[j+1][v]-tab[j][v])/factor;}
     if(k>=3){err=0;for(int v=0;v<2;v++)err=std::max(err,fabsq(tab[0][v]-tab[1][v])/(tol*(1+fabsq(tab[0][v]))));if(err<=1){out=tab[0];ok=true;break;}}
    }
   }catch(int code){if(code==3)throw;ok=false;}
   if(ok){y=out;s+=H;H*=std::min(1.5Q,std::max(.5Q,.9Q*powq(std::max(err,1e-10Q),-1.Q/12)));H=std::min(H,.08Q);}
   else H*=.5Q;
  }
  return y;
 }
};
R parse(const char*s){std::string a(s);auto i=a.find('/');return i==std::string::npos?strtoflt128(s,0):strtoflt128(a.substr(0,i).c_str(),0)/strtoflt128(a.substr(i+1).c_str(),0);}
extern "C" void matching(const char*cs,const char*ms,const char*bs,const char*ts,int n,const double*z,double*F,double*G,int*status){
 Engine e;e.c=parse(cs);e.m=parse(ms);e.beta=parse(bs);e.sm=sqrtq(e.m);e.tol=parse(ts);
 for(int i=0;i<n;i++){
  e.evaluations=0;F[i]=G[i]=0;status[i]=0;
  try{auto a=e.integrate((R)z[i],1),b=e.integrate((R)z[i],-1);F[i]=(double)(a[0]-b[0]);G[i]=(double)(a[1]-b[1]);}
  catch(int code){status[i]=code;}
 }
}
