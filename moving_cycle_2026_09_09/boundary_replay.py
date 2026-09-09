"""Resolve the two observed outer-domain failures; no exhaustive root claim."""
from engine import *

def line_diagnostic(s,q,rtol=2e-13):
    a,b,e0,e1,e2=coefficients(q)
    def rhs(t,Y):return field(Y,q,1,regularized=True)[0]
    def section(t,Y):return Y[0]
    section.terminal=True;section.direction=-1
    def switch(t,Y):return Y[1]+14
    switch.terminal=True;switch.direction=-1
    kick=solve_ivp(rhs,[0,1e-7],[0.,s],method='DOP853',rtol=rtol,atol=rtol*.01)
    sol=solve_ivp(rhs,[1e-7,300],kick.y[:,-1],method='DOP853',rtol=rtol,atol=rtol*.01,
                  events=[section,switch],max_step=.2)
    if len(sol.t_events[0]):return {'status':'returned_before_chart_switch','z':float(sol.y_events[0][0][1])}
    if not len(sol.t_events[1]):return {'status':'unresolved','message':sol.message}
    v,z=sol.y_events[1][0];initial=np.array([np.sinh(v)*np.exp(z),np.exp(z)])
    def cart(t,Y):
        x,y=Y
        return [(b-2)/4+(1-b)*y+a*x*x+b*y*y+e1*x+e2*x*y,e0-2*x*y]
    def line(t,Y):return Y[1]
    line.terminal=True;line.direction=-1
    def xs(t,Y):return Y[0]
    xs.terminal=True;xs.direction=-1
    cr=solve_ivp(cart,[0,10],initial,method='DOP853',rtol=rtol,atol=rtol*.001,events=[line,xs],max_step=.01)
    return {'status':'crossed_y_zero' if len(cr.t_events[0]) else 'other','switch_xy':initial.tolist(),
            'events':[x.tolist() for x in cr.y_events],'e0':float(e0),
            'message':cr.message,'rtol':rtol,'nfev':kick.nfev+sol.nfev+cr.nfev}

def main(kind):
    start=time.perf_counter();r=json.loads((ROOT/f'{kind}_result.json').read_text());last=r['accepted'][-1]
    q=np.array(last['q']);s3=last['roots'][2];out={'kind':kind,'evidence':'NUM; finite boundary probes only','q':q.tolist(),
         'rescaled_checks':[],'bisection':[],'domain_side_samples':[]}
    # Check the extra time factor on successful returns before using it near failure.
    for gap in [1,10]:
        p=pair(s3+gap,q,rtol=8e-14,regularized=True)
        old=next(x for x in last['outer'] if x['status']=='passed' and abs(x['s']-(s3+gap))<1e-10)
        out['rescaled_checks'].append({'new':p,'old_D':old['D'],'D_change':p.get('D',0)-old['D'],
              'derivative_change':float(p['derivatives'][0]-old['derivatives'][0])})
    out['failed_height_replay']=pair(s3+20,q,rtol=3e-13,regularized=True)
    out['line_crossing']=line_diagnostic(s3+20,q)
    lo,hi=s3+10,s3+20
    for _ in range(12):
        mid=(lo+hi)/2;p=pair(mid,q,rtol=3e-13,regularized=True)
        out['bisection'].append(p)
        if p['status']=='passed':lo=mid
        else:hi=mid
        save(f'{kind}_boundary_checkpoint.json',out)
    out['transition_interval']=[lo,hi]
    for delta in [.1,.03,.01,.003]:
        p=pair(lo-delta,q,rtol=8e-14,regularized=True)
        out['domain_side_samples'].append(p)
    out['last_resolved_line_diagnostic']=line_diagnostic(lo,q)
    out['first_unresolved_line_diagnostic']=line_diagnostic(hi,q)
    out['counts']=COUNTS.copy();out['wall_seconds']=time.perf_counter()-start
    out['total_branch_pairs_including_boundary']=r['calls']+out['counts']['pairs']
    save(f'{kind}_boundary.json',out)
    print(kind,'transition',out['transition_interval'],'line',out['line_crossing'],
          'inside',[(p['s'],p.get('D'),p.get('derivatives',[None])[0]) for p in out['domain_side_samples']],
          'counts',COUNTS,flush=True)

if __name__=='__main__':
    import sys
    main(sys.argv[1])
