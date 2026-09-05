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
R c,K,al,ac,ak;int nfev=0;R direction=1;
A rhs(R s,const A&z){
 ++nfev;R w=z[0],M=z[1],C=cosq(direction*s),S=-sinq(direction*s),E=expq(w);
 R h1=(1+al)*C*S,g1=al*C*C-S*S,p2=C*C+C*S,q2=-10*C*C+2.2Q*C*S+c*S*S;
 R h2=C*p2+S*q2,g2=C*q2-S*p2,H=h1+E*h2,G=g1+E*g2;
 if(!(G<0))throw std::runtime_error("angular chart lost monotonicity");
 R h1c=ac*C*S,g1c=ac*C*C,h1K=ak*C*S,g1K=ak*C*C,h2c=S*S*S,g2c=C*S*S;
 R Hc=h1c+E*h2c,Gc=g1c+E*g2c,HK=h1K,GK=g1K;
 R N=h2*g1-h1*g2,Nc=h2c*g1+h2*g1c-h1c*g2-h1*g2c,NK=h2*g1K-h1K*g2;
 R fw=-E*N/(G*G),fww=fw*(1-2*E*g2/G);
 R fc=-(Hc*G-H*Gc)/(G*G),fk=-(HK*G-H*GK)/(G*G);
 R fwc=-E*Nc/(G*G)+2*E*N*Gc/(G*G*G),fwK=-E*NK/(G*G)+2*E*N*GK/(G*G*G);
 A ret{-H/G,fw,fw*z[2]+fc,fw*z[3]+fk,fww*expq(M),fww*z[2]+fwc,fww*z[3]+fwK,
          -1/G,-E*(4.2Q*C+(1+2*c)*S)/G};
 for(auto &v:ret)v*=direction;return ret;
}
A mid(R s,const A&y,R H,int n){
 R h=H/n;A old=y,cur=y,f=rhs(s,y);
 for(int j=0;j<9;j++)cur[j]+=h*f[j];
 for(int i=1;i<n;i++){f=rhs(s+i*h,cur);A nxt;
  for(int j=0;j<9;j++)nxt[j]=old[j]+2*h*f[j];old=cur;cur=nxt;}
 f=rhs(s+H,cur);A out;for(int j=0;j<9;j++)out[j]=(old[j]+cur[j]+h*f[j])/2;return out;
}
A integrate(R r,R tol,R dir){direction=dir;
 R s=0,end=acosq(-1.L),H=.01Q;A y{logq(r),0,0,0,0,0,0,0,0};int steps=0,reject=0;
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
 return y;
}

int main(){R r,tol;std::string sr,sc,sk,st;std::cin>>sr>>sc>>sk>>st;r=strtoflt128(sr.c_str(),0);c=strtoflt128(sc.c_str(),0);K=strtoflt128(sk.c_str(),0);tol=strtoflt128(st.c_str(),0);R start=std::clock()/(R)CLOCKS_PER_SEC;
try{
 al=-K;ac=0;ak=-1;
 if(r<=0||al>=0)throw std::runtime_error("focus/section gate");
 A a=integrate(r,tol,1),b=integrate(r,tol,-1);R ea=expq(a[1]),eb=expq(b[1]);
 std::cout<<std::setprecision(23)<<"{\"status\":\"NUMERICAL_TWO_HALF_PASSAGES\",\"section\":\"positive horizontal ray matched on negative ray\"";
 auto emit=[&](const char*k,R v){char buf[100];quadmath_snprintf(buf,sizeof(buf),"%.36Qg",v);std::cout<<",\""<<k<<"\":\""<<buf<<"\"";};
 emit("r",r);emit("c",c);emit("m",K);emit("K",K*(11*c-5)/5-42);emit("alpha",al);
 emit("F",a[0]-b[0]);emit("G",a[1]-b[1]);
 emit("F_z",ea-eb);emit("F_c",a[2]-b[2]);emit("F_m",a[3]-b[3]);
 emit("G_z",a[4]-b[4]);emit("G_c",a[5]-b[5]);emit("G_m",a[6]-b[6]);
 emit("multiplier_at_match",expq(a[1]-b[1]));emit("period_at_match",a[7]-b[7]);
 emit("negative_forward",-expq(a[0]));emit("negative_backward",-expq(b[0]));
 emit("forward_log_sensitivity",a[1]);emit("backward_log_sensitivity",b[1]);
 emit("divergence_multiplier_at_match",expq(a[8]-b[8]));
 emit("cpu_seconds",std::clock()/(R)CLOCKS_PER_SEC-start);std::cout<<",\"nfev\":"<<nfev<<"}\n";
}catch(const std::exception&e){std::cout<<"{\"status\":\"UNRESOLVED\",\"error\":\""<<e.what()<<"\",\"nfev\":"<<nfev<<"}\n";}}
