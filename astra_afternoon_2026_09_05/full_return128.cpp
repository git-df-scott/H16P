// Binary128 full angular return, modified-midpoint extrapolation.
// Numerical evidence only. Non-monotone angular charts remain UNRESOLVED.
#include <quadmath.h>
#include <array>
#include <algorithm>
#include <iostream>
#include <string>
#include <stdexcept>
using R=__float128;
static_assert(__FLT128_MANT_DIG__ == 113, "IEEE binary128 required");
R parse(std::string s){auto k=s.find('/');if(k!=std::string::npos){R den=parse(s.substr(k+1));if(den==0)throw std::runtime_error("zero denominator");return parse(s.substr(0,k))/den;}char*end;R q=strtoflt128(s.c_str(),&end);if(*end||!finiteq(q))throw std::runtime_error("invalid number");return q;}
std::string fmt(R x){char b[128];quadmath_snprintf(b,sizeof(b),"%.36Qg",x);return b;}
struct Engine{
 std::array<R,12> c;R tol,scale,theta,sense;long eval=0,maxeval;R minG=1e4900Q;int steps=0,reject=0;
 R rhs(R s,R w){
  if(++eval>maxeval)throw 3;
  if(!finiteq(w)||fabsq(w)>2000)throw 5;
  R C=cosq(theta+sense*s),S=sinq(theta+sense*s),E=expq(w);
  // y=scale*Y. Coefficients were translated exactly in Python first.
  R p1=c[1]*C+c[2]*scale*S,q1=c[7]/scale*C+c[8]*S;
  R p2=c[3]*C*C+c[4]*scale*C*S+c[5]*scale*scale*S*S;
  R q2=c[9]/scale*C*C+c[10]*C*S+c[11]*scale*S*S;
  R G=C*q1-S*p1+E*(C*q2-S*p2),H=C*p1+S*q1+E*(C*p2+S*q2);
  if(!finiteq(G)||!finiteq(H))throw 5;
  if(!(sense*G>0))throw 4;
  minG=std::min(minG,sense*G);R f=sense*H/G;if(!finiteq(f))throw 5;return f;
 }
 R mid(R s,R y,R H,int n){R h=H/n,old=y,cur=y+h*rhs(s,y);for(int i=1;i<n;i++){R nxt=old+2*h*rhs(s+i*h,cur);old=cur;cur=nxt;}return (old+cur+h*rhs(s+H,cur))/2;}
 R integrate(R u,R angle){
  R C=cosq(angle),S=sinq(angle)/scale;theta=atan2q(S,C);R z=u+logq(sqrtq(C*C+S*S));
  R p=c[1]*C+c[2]*scale*S+expq(u)*(c[3]*C*C+c[4]*scale*C*S+c[5]*scale*scale*S*S);
  R q=c[7]/scale*C+c[8]*S+expq(u)*(c[9]/scale*C*C+c[10]*C*S+c[11]*scale*S*S);
  R G=C*q-S*p;if(!finiteq(G)||G==0)throw 4;sense=G>0?1:-1;
  R t=0,end=2*acosq(-1.Q),h=.01Q,y=z;
  while(t<end){
   ++steps;h=std::min(h,end-t);if(h<1e-28Q)throw 4;
   std::array<R,8>tab;R err=1e100Q,out=0;bool ok=false;
   try{for(int k=0;k<8;k++){
    int n=2*(k+1);tab[k]=mid(t,y,h,n);if(!finiteq(tab[k]))throw 5;
    for(int j=k-1;j>=0;j--){R ratio=(R)n/(2*(j+1));tab[j]=tab[j+1]+(tab[j+1]-tab[j])/(ratio*ratio-1);if(!finiteq(tab[j]))throw 5;}
    if(k>=3){err=fabsq(tab[0]-tab[1])/(tol*(1+fabsq(tab[0])));if(!finiteq(err))throw 5;if(err<=1){out=tab[0];ok=true;break;}}
   }}catch(int code){if(code==3)throw;ok=false;}
   if(ok){y=out;t+=h;h*=std::min(1.5Q,std::max(.5Q,.9Q*powq(std::max(err,1e-10Q),-1.Q/12)));h=std::min(h,.08Q);}
   else{++reject;h*=.5Q;}
  }
  return y-z;
 }
};
int main(){try{
 Engine base;std::string s;for(auto&v:base.c){std::cin>>s;v=parse(s);}std::cin>>s;R angle=parse(s);std::cin>>s;base.tol=parse(s);std::cin>>s;base.scale=parse(s);std::cin>>base.maxeval;int n;std::cin>>n;
 if(base.c[0]!=0||base.c[6]!=0||base.scale<=0||base.tol<1e-32Q||base.tol>1e-3Q)throw std::runtime_error("invalid origin/scale/tolerance");
 for(int i=0;i<n;i++){std::cin>>s;R u=parse(s);Engine e=base;
  try{R D=e.integrate(u,angle);std::cout<<"{\"status\":\"OK_NUMERICAL\",\"log_radius\":\""<<fmt(u)<<"\",\"log_displacement\":\""<<fmt(D)<<"\",\"radial_displacement\":\""<<fmt(expq(u)*expm1q(D))<<"\",\"direction\":"<<(int)e.sense;}
  catch(int code){std::cout<<"{\"status\":\""<<(code==3?"EVALUATION_LIMIT":code==4?"ANGULAR_CHART_UNRESOLVED":"NONFINITE_OR_RANGE")<<"\",\"log_radius\":\""<<fmt(u)<<"\",\"log_displacement\":null";}
  std::cout<<",\"evaluations\":"<<e.eval<<",\"steps\":"<<e.steps<<",\"rejections\":"<<e.reject<<"}\n";
 }
 }catch(const std::exception&e){std::cerr<<e.what()<<"\n";return 2;}}
