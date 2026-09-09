"""Target previously unresolved outer return boundary using saddle-section radius."""
import json
import numpy as np
from q4_finite_continuation import P,difference,COUNTS

def main():
    fields={(a['r'],a['epsilon']):a for a in json.loads((P/'q4_finite_extended.json').read_text())['records']}
    source=json.loads((P/'q4_saddle_separatrices.json').read_text())['records']
    out={'scope':'Nonvalidated half-return sampling up to numerically shot stable saddle boundary; no exclusion between samples','records':[]}
    for a in source:
        f=fields[(a['r'],a['epsilon'])]
        stable=[q for q in a['rows'] if q['manifold']=='stable' and q['status']=='negative_section']
        unstable=[q for q in a['rows'] if q['manifold']=='unstable' and q['status']=='negative_section']
        rs=stable[-1]['radius'];ru=unstable[-1]['radius'];lower=max(f['anchors'])*1.05
        radii=rs-(rs-lower)*np.geomspace(1,1e-8,17)
        rows=[difference(a['r'],a['epsilon'],f['normalized_controls'],R) for R in radii]
        changes=[[x['radius'],y['radius']] for x,y in zip(rows,rows[1:]) if 'difference' in x and 'difference' in y and x['difference']*y['difference']<0]
        item=dict(r=a['r'],epsilon=a['epsilon'],stable_radius=rs,unstable_radius=ru,section_splitting=ru-rs,stable_offset_change=abs(stable[0]['radius']-rs),unstable_offset_change=abs(unstable[0]['radius']-ru),rows=rows,sign_change_brackets=changes)
        out['records'].append(item);out['counts']=COUNTS.copy()
        print(json.dumps({k:v for k,v in item.items() if k!='rows'}|{'resolved':sum('difference' in q for q in rows),'last_difference':rows[-1].get('difference')}),flush=True)
        (P/'q4_saddle_boundary_profile.json').write_text(json.dumps(out,indent=2)+'\n')
if __name__=='__main__':main()
