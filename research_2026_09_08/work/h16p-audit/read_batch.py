"""Print bounded complete text chunks; acknowledge only after model inspection.

Usage: python read_batch.py next
       python read_batch.py ack
The caller must NOT acknowledge truncated tool output.
"""
from pathlib import Path
import json,sys

HERE=Path(__file__).resolve().parent
REPO=HERE.parent/'H16P'
LEDGER=HERE/'coverage.json'
PENDING=HERE/'pending_read.json'
d=json.loads(LEDGER.read_text())
if sys.argv[1]=='ack':
    batch=json.loads(PENDING.read_text())
    for item in batch:
        row=next(r for r in d['files'] if r['path']==item['path'])
        row['read_chars']=item['end']
        row['status']='read' if item['end']==item['length'] else 'partial'
    LEDGER.write_text(json.dumps(d,indent=2)+'\n')
    PENDING.unlink()
    print('Acknowledged',len(batch),'chunks.')
else:
    budget=22000
    selected=[]
    # Narrative and source code need actual reading, data receive separate inspection.
    eligible=[r for r in d['files'] if r['status']!='read' and
              Path(r['path']).suffix in ('.md','.py','.cpp','.txt','.log','')]
    eligible.sort(key=lambda r:(0 if r['path'].endswith('.md') else 1,r['path']))
    for row in eligible:
        text=(REPO/row['path']).read_text()
        start=row.get('read_chars',0)
        end=min(len(text),start+budget)
        print('\nFILE:',row['path'],'CHARACTERS:',start,'TO',end,'OF',len(text))
        print(text[start:end])
        selected.append({'path':row['path'],'end':end,'length':len(text)})
        budget-=end-start
        if budget<1000:break
    PENDING.write_text(json.dumps(selected,indent=2))
    print('\nEND OF BATCH. Acknowledge only if the entire output was visible.')
