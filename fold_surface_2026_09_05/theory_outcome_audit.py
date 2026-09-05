#!/usr/bin/env python3
"""Read-only numerical-event audit. No flow integration or external calls."""
from pathlib import Path
from fractions import Fraction as Q
from collections import Counter
import hashlib,json,math
HERE=Path(__file__).resolve().parent
def number(x): return Q(str(x))
def KH(c):
    J=305+634*c-11*c*c-1000*c**3
    return -441*J/(125*(16-10*c)*(1+2*c)**2)
def remote_enclosure(c,K):
    m=5*(K+42)/(11*c-5)
    h=Q(61,5)-c; b=m-Q(72,5)+3*c; aa=3*c-Q(11,5)
    def polynomial(v):return -h*v**3+b*v*v+aa*v+c
    lo,hi=Q(0),Q(1)
    while polynomial(hi)>0:hi*=2
    for _ in range(64):
        mid=(lo+hi)/2
        if polynomial(mid)>0:lo=mid
        else:hi=mid
    # trace=(v+1)*[(1+2c)/v-(16-10c)/5], decreasing in v here.
    def trace(v):return (v+1)*((1+2*c)/v-(16-10*c)/5)
    return {'x_interval':[str(-hi-1),str(-lo-1)],
            'trace_interval':[str(trace(hi)),str(trace(lo))],
            'trace_approximation':float((trace(lo)+trace(hi))/2),
            'trace_positive_exact':trace(hi)>0}
result={'status':'SAVED_DATA_AUDIT_ONLY','orbit_evaluations':0,'files':{}}
for name in ['increasing','decreasing','arclength','angular_ld','half','quad']:
    path=HERE/f'events_{name}.json'; raw=path.read_bytes(); events=json.loads(raw)
    accepted=[]; failures=Counter(); roots=Counter(); raw_roots=Counter(); problems=[]
    remote_gaps=[]; curvatures=[]; residuals=[]; flux=[]; stability=[]
    for i,event in enumerate(events):
        a=event.get('fold',{}); ok='ACCEPTED' in a.get('status','') or event.get('status')=='ACCEPTED'
        if not ok:
            failures[a.get('status',event.get('status','unknown'))]+=1;continue
        accepted.append(a); ret=a.get('return_data',a)
        c,K=number(a['c']),number(a['K'])
        assert Q(24,25)<=c<Q(8,5) and K>0
        remote_gaps.append(float(KH(c)-K))
        curvatures.append(float(ret.get('G_z',ret.get('L_zz'))))
        residuals.append([abs(float(ret[k])) for k in (('F','G') if 'F' in ret else ('L','L_z'))])
        if 'first_derivative_discrepancy' in ret:flux.append(abs(float(ret['first_derivative_discrepancy'])))
        if 'divergence_multiplier_at_match' in ret:
            flux.append(abs(float(number(ret['multiplier_at_match'])-number(ret['divergence_multiplier_at_match']))))
        pair=event.get('pair',event.get('pair_profile',{}))
        if not pair:continue
        cp,kp=number(pair['c']),number(pair.get('K',K));remote_gaps.append(float(KH(cp)-kp))
        brackets=pair.get('root_brackets',pair.get('root_sign_brackets',[])); roots[len(brackets)]+=1
        for j,b in enumerate(brackets):
            if 'signs' in b:
                good=number(b['signs'][0])*number(b['signs'][1])<0 and b['log_bracket'][0]<b['log_bracket'][1]
                root=b.get('root',next((x for x in reversed(b.get('history',[])) if x.get('status')=='NUMERICAL_ONLY'),{}))
                st='stable' if root.get('log_displacement_derivative',0)<0 else 'unstable'
            else:
                l,r=b['left'],b['right']
                key='F' if 'F' in l else 'L'
                good=number(l[key])*number(r[key])<0 and number(l['r'])<number(r['r'])
                st=b.get('stability','unassigned')
                if st in ['stable','unstable'] and 'G' in b.get('approximation',{}):
                    if (number(b['approximation']['G'])<0)!=(st=='stable'):problems.append([i,j,'stability label mismatch'])
            if not good:problems.append([i,j,'invalid saved sign bracket'])
            stability.append(st)
        samples=pair.get('profile',pair.get('samples',[])); count=0
        for l,r in zip(samples,samples[1:]):
            l,r=l['result'],r['result']; key='F' if 'F' in l else ('L' if 'L' in l else 'log_displacement')
            if key in l and key in r and number(l[key])*number(r[key])<0:count+=1
        raw_roots[count]+=1
    result['files'][path.name]={
        'sha256':hashlib.sha256(raw).hexdigest(),'events':len(events),'accepted':len(accepted),
        'failures':dict(failures),'saved_root_bracket_count_histogram':dict(roots),
        'adjacent_profile_sign_change_histogram':dict(raw_roots),'bracket_or_label_problems':problems,
        'minimum_positive_fold_curvature':min(curvatures) if curvatures else None,
        'max_abs_residual_components':[max(x[j] for x in residuals) for j in (0,1)] if residuals else None,
        'maximum_recorded_flux_discrepancy':max(flux) if flux else None,
        'minimum_exact_remote_Hopf_gap_KH_minus_K':min(remote_gaps) if remote_gaps else None,
        'all_checked_remote_traces_positive':all(x>0 for x in remote_gaps),
        'root_stability_labels':dict(Counter(stability)),
        'last_accepted':{k:accepted[-1][k] for k in ('r','c','K')} if accepted else None,
        'last_remote_exact_enclosure':remote_enclosure(number(accepted[-1]['c']),number(accepted[-1]['K'])) if accepted else None}
(HERE/'theory_outcome_audit.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
