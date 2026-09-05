// Long-double independent angular shooting with modified-midpoint extrapolation.
// No interval arithmetic. A failed angular chart is unresolved.
#include <array>
#include <quadmath.h>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <algorithm>
#include <ctime>
using R=__float128;using A=std::array<R,9>;
R c,K,al,ac,ak,beta;int nfev=0;
A rhs(R s,const A&z){
 ++nfev;R w=z[0],M=z[1],C=cosq(s),S=-sinq(s),E=expq(w);
 R h1=(1+al)*C*S+beta*S*S,g1=al*C*C-S*S+beta*C*S,p2=C*C+C*S,q2=-10*C*C+2.2Q*C*S+c*S*S;
 R h2=C*p2+S*q2,g2=C*q2-S*p2,H=h1+E*h2,G=g1+E*g2;
 if(!(G<0))throw std::runtime_error("angular chart lost monotonicity");
 R h1c=ac*C*S,g1c=ac*C*C,h1K=ak*C*S,g1K=ak*C*C,h2c=S*S*S,g2c=C*S*S;
 R Hc=h1c+E*h2c,Gc=g1c+E*g2c,HK=h1K,GK=g1K;
 R N=h2*g1-h1*g2,Nc=h2c*g1+h2*g1c-h1c*g2-h1*g2c,NK=h2*g1K-h1K*g2;
 R fw=-E*N/(G*G),fww=fw*(1-2*E*g2/G);
 R fc=-(Hc*G-H*Gc)/(G*G),fk=-(HK*G-H*GK)/(G*G);
 R fwc=-E*Nc/(G*G)+2*E*N*Gc/(G*G*G),fwK=-E*NK/(G*G)+2*E*N*GK/(G*G*G);
 return A{-H/G,fw,fw*z[2]+fc,fw*z[3]+fk,fww*expq(M),fww*z[2]+fwc,fww*z[3]+fwK,
          -1/G,-(beta+E*(4.2Q*C+(1+2*c)*S))/G};
}
A mid(R s,const A&y,R H,int n){
 R h=H/n;A old=y,cur=y,f=rhs(s,y);
 for(int j=0;j<9;j++)cur[j]+=h*f[j];
 for(int i=1;i<n;i++){f=rhs(s+i*h,cur);A nxt;
  for(int j=0;j<9;j++)nxt[j]=old[j]+2*h*f[j];old=cur;cur=nxt;}
 f=rhs(s+H,cur);A out;for(int j=0;j<9;j++)out[j]=(old[j]+cur[j]+h*f[j])/2;return out;
}
int main(){
 R r,tol;std::string sr,sc,sk,sb,st;std::cin>>sr>>sc>>sk>>sb>>st;beta=strtoflt128(sb.c_str(),0);r=strtoflt128(sr.c_str(),0);c=strtoflt128(sc.c_str(),0);K=strtoflt128(sk.c_str(),0);tol=strtoflt128(st.c_str(),0);R start=std::clock()/(R)CLOCKS_PER_SEC;
 try{
 al=-K;ac=0;ak=-1;
 if(r<=0||al>=0)throw std::runtime_error("focus/section gate");
 R s=0,end=2*acosq(-1.L),H=.01Q;A y{logq(r),0,0,0,0,0,0,0,0};int steps=0,reject=0;
 while(s<end){
  if(++steps>2000000)throw std::runtime_error("step guard");
  H=std::min(H,end-s);if(H<1e-20Q)throw std::runtime_error("step resolution guard");
  std::array<A,9> tab;R err=1e100Q;A out;bool ok=false;
  try{
   for(int k=0;k<8;k++){
    int n=2*(k+1);tab[k]=mid(s,y,H,n);
    for(int j=k-1;j>=0;j--){R ratio=(R)n/(2*(j+1)),factor=ratio*ratio-1;
     for(int v=0;v<9;v++)tab[j][v]=tab[j+1][v]+(tab[j+1][v]-tab[j][v])/factor;}
    if(k>=3){err=0;for(int v=0;v<9;v++)err=std::max(err,fabsq(tab[0][v]-tab[1][v])/(tol*(1+fabsq(tab[0][v]))));
      if(err<=1){out=tab[0];ok=true;break;}}
   }
  }catch(const std::exception&){ok=false;}
  if(ok){y=out;s+=H;H*=std::min(1.5Q,std::max(.5Q,.9Q*powq(std::max(err,1e-10Q),-1.L/12)));H=std::min(H,.08Q);}
  else {H*=.5Q;if(++reject>200000)throw std::runtime_error("rejection guard");}
 }
 R RR=expq(y[0]),aa=expq(y[1]),Q0=-10*r*r+al*r,Qf=-10*RR*RR+al*RR;
 R fluxA=r/RR*Q0/Qf*expq(y[8]);
 std::cout<<std::setprecision(23)<<"{\"status\":\"NUMERICAL_ONLY\",\"section\":\"positive horizontal ray; clockwise\"";
 auto emit=[&](const char*k,R v){char buf[100];quadmath_snprintf(buf,sizeof(buf),"%.36Qg",v);std::cout<<",\""<<k<<"\":\""<<buf<<"\"";};
 emit("beta",beta);emit("r",r);emit("c",c);emit("m",K);emit("K",K*(11*c-5)/5-42);emit("alpha",al);emit("return_coordinate",RR);
 emit("L",y[0]-logq(r));emit("L_z",expm1q(y[1]));emit("L_zz",aa*y[4]);
 emit("L_c",y[2]);emit("L_m",y[3]);emit("L_zc",aa*y[5]);emit("L_zm",aa*y[6]);
 emit("multiplier",expq(y[8]));emit("period",y[7]);emit("first_derivative_discrepancy",aa-fluxA);
 emit("cpu_seconds",std::clock()/(R)CLOCKS_PER_SEC-start);std::cout<<",\"nfev\":"<<nfev<<",\"steps\":"<<steps<<"}\n";
 }catch(const std::exception&e){std::cout<<"{\"status\":\"UNRESOLVED\",\"error\":\""<<e.what()<<"\",\"nfev\":"<<nfev<<"}\n";}
}
