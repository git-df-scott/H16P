from engine import *
seed=json.loads((ROOT/'seed_replay.json').read_text())
s3=seed['roots'][2]['s'];out={'evidence':'NUM','q':SEED.tolist(),'rows':[]}
for gap in [.02,.1,.25,.5,1,1.5,2,3,4,6,8,10,12,16,20]:
    p=pair(s3+gap,SEED)
    out['rows'].append(p);save('outer_seed_profile.json',out)
    print(gap,p['status'],p.get('D'),p.get('derivatives',[None])[0],flush=True)
out['counts']=COUNTS.copy();save('outer_seed_profile.json',out)
