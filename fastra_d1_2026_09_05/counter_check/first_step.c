/* Instrument the unchanged Fable first step; no production counter changes. */
#include "../../audit/fable_engine/retmap_log.c"
#include <stdio.h>
int main(void){
 double c[12],u,theta;
 for(int i=0;i<12;i++) if(scanf("%lf",c+i)!=1)return 2;
 if(scanf("%lf %lf",&u,&theta)!=2)return 2;
 double du,dt,un,tn;rates(c,u,theta,&du,&dt);
 double err=dp_step(c,u,theta,1e-3,&un,&tn,1e-12,1e-15);
 printf("{\"initial_du\":%.17g,\"initial_dtheta\":%.17g,\"h\":0.001,\"error_finite\":%s,\"next_u_finite\":%s,\"next_theta_finite\":%s,\"original_reject_condition\":%s}\n",du,dt,isfinite(err)?"true":"false",isfinite(un)?"true":"false",isfinite(tn)?"true":"false",err>1?"true":"false");
 return 0;
}
