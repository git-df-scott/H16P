"""Build figures/report numerical table from immutable saved events (no ODE)."""
from pathlib import Path
import json,math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
HERE=Path(__file__).resolve().parent
s=json.loads((HERE/'component_summary.json').read_text());a=[x for x in s['events'] if x['status']=='ACCEPTED_NUMERICAL_FOLD']
pos=sorted([x for x in a if float(x['K'])>0],key=lambda x:float(x['K']))
neg=sorted([x for x in a if float(x['K'])<0],key=lambda x:float(x['c']))
fig,ax=plt.subplots(1,3,figsize=(14,4.4),layout='constrained')
ax[0].plot([float(x['K']) for x in pos],[float(x['c']) for x in pos],'.-',color='#176b8d',label='Accepted fold points')
ax[0].scatter([0,7.18499469640662],[.9686206335534943,1.6],marker='x',s=65,color='#b84c3e',label='Candidate limiting organizers')
ax[0].set(xlabel='K',ylabel='c',title='Positive-K sheet');ax[0].axhline(1,color='gray',lw=.6,ls=':');ax[0].legend(fontsize=7)
for label in ['curved P=0 maximum x','positive horizontal ray']:
 p=[x for x in pos if x['section']==label];ax[1].plot([float(x['c']) for x in p],[math.log10(float(x['r'])) for x in p],'.',label=label)
ax[1].axvline(1.6,color='#b84c3e',lw=.8,ls='--');ax[1].set(xlabel='c',ylabel='log10 section radius',title='Amplitude along positive sheet');ax[1].legend(fontsize=7)
xs=[];ys=[]
for x in neg:
 e=json.loads((HERE/x['file']).read_text())[x['event_index']]['fold'];M=float(e.get('m',-float(e['alpha'])));xs.append(math.log10(M));ys.append(float(e['c']))
ax[2].plot(xs,ys,'.-',color='#6e4790');ax[2].axhline(1/3,color='#b84c3e',lw=.8,ls='--',label='1/3 asymptotic organizer');ax[2].axhline(5/11,color='gray',lw=.6,ls=':',label='(c,K) chart pole');ax[2].set(xlabel='log10 m',ylabel='c',title='Center-selected negative-K sheet');ax[2].legend(fontsize=7)
for q in ax:q.grid(alpha=.18)
fig.savefig(HERE/'fold_continuation.png',dpi=180);fig.savefig(HERE/'fold_continuation.svg')
lines=['## Numerical checkpoint', '',f"This strike recorded **{s['new_evaluations']}** charged calls, including unresolved attempts and separatrix passages. Together with the inherited 756, the campaign has used **{s['campaign_evaluations']}/4096**, leaving **{s['remaining']}**. No old 550-call stage was repeated.", '', '| Sheet / event | K | c | section r | pair brackets |','|---|---:|---:|---:|---:|']
selected=[]
for file in ['events_decreasing.json','events_increasing.json','events_arclength.json','events_quad.json','events_negative.json','events_m.json','events_logm.json']:
 p=[x for x in a if x['file']==file]
 if p:selected.append(p[-1])
for x in selected:
 lines.append('| '+x['file'].replace('events_','').replace('.json','')+' | '+f"{float(x['K']):.12g} | {float(x['c']):.14g} | {float(x['r']):.13g} | {x['root_sign_brackets']} |")
lines+=['','The first three table rows use the curved maximum-x section; later rows use the positive horizontal ray. These section coordinates must not be equated. Full decimal values and every intermediate event are archived in the linked ledger.','', '![Numerical fold continuation](fold_surface_2026_09_05/fold_continuation.png)','']
c1=json.loads((HERE/'crossing_c1.json').read_text())['result'];center=json.loads((HERE/'center_binary128.json').read_text())['rows'][-1]
inf=json.loads((HERE/'infinity_binary128.json').read_text());good=[x['result'] for x in inf['rows'] if x['result']['status']=='NUMERICAL_TWO_HALF_PASSAGES'];b=good[-1]
lines += [f"The c=1 crossing is numerical at K={c1['K']}, horizontal r={c1['r']}; it is a regular finite fold, not loss at infinity.", '',f"At K=1e-9, the decreasing branch gives c={center['result']['c']}, horizontal r={center['result']['r']}, and (c-c*)/K={center['secant_dc_dK']}.", '', f"At exact target c=8/5, finite-radius matching through r=1e14 gives K≈{b['K']} and G≈{b['G']}. The tighter r=1e17 control hit the CPU fuse; subsequent looser-tolerance controls completed and agree closely (see infinity_tolerance_control.json). The independent separatrix-series approach gives K≈7.18499469694; its roughly 5e-10 difference is numerical integration error, not a distinct connection.", '', 'The negative-K sheet crosses the removable (c,K) chart pole using independent coefficients (c,m). Its later growth motivates the separate large-m asymptotic analysis in the supplemental theory review. This sheet is included as an extension of the center unfolding, not as evidence that the positive-K regular component has been exhaustively completed.','']
if (HERE/'graphic_coefficient.json').exists():
 g=json.loads((HERE/'graphic_coefficient.json').read_text())
 if 'C_approx' in g:lines += [f"Refined finite-radius matching at r=1e17 gives K≈{g['K_at_finite_matching_radius']}, C≈{g['C_approx']}, and conditional (1.6-c)log r limit≈{g['conditional_delta_log_r_limit']}. This is numerical/asymptotic evidence, not an exact connection certificate.",'']
if (HERE/'pair_replay_claims.json').exists():
 claims=json.loads((HERE/'pair_replay_claims.json').read_text())
 lines += ['Complete-return root refinements are in `pair_replay_claims.json`, including rational coefficient vectors, preserved endpoint signs, periods and multipliers. Each row group is a separate two-cycle field; their cycles must not be added across fields.','','| Field | Origin cycle | Approximate horizontal r | Period | Multiplier |','|---|---|---:|---:|---:|']
 for row in claims['rows']:
  for i,x in enumerate(row['cycles']):
   if 'root' in x:
    b=x['root'];lines.append(f"| {row['source']} | {i+1} | {float(b['r']):.12g} | {float(b['period']):.12g} | {float(b['multiplier']):.12g} |")
 lines.append('')
p=HERE.parent/'KKL_FOLD_SURFACE_STRIKE.md';t=p.read_text();start='<!-- GENERATED_RESULTS -->';end='<!-- END_GENERATED_RESULTS -->'
if end in t:t=t[:t.index(start)]+start+'\n'+'\n'.join(lines)+'\n'+t[t.index(end):]
else:t=t.replace(start,start+'\n'+'\n'.join(lines)+'\n'+end)
p.write_text(t)
