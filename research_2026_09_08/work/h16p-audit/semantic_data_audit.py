"""Inspect every structured record, preserving schemas, ranges and failures.
Does not label raw files fully read and does not validate ODE integration.
"""
from pathlib import Path
from collections import Counter,defaultdict
from decimal import Decimal,InvalidOperation
import json,re,hashlib,math
ROOT=Path(__file__).resolve().parent
REPO=ROOT.parent/'H16P'
NUM=re.compile(r'^[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?$')
def run(path):
    text=path.read_text();structured=[];unstructured=[];errors=[]
    if path.suffix=='.json':
        try:structured=[json.loads(text)]
        except Exception as e:errors.append(str(e))
    else:
        for line_no,line in enumerate(text.splitlines(),1):
            if not line.strip():continue
            try:structured.append(json.loads(line))
            except json.JSONDecodeError:unstructured.append({'line':line_no,'text':line})
    counts=Counter();schemas=Counter();statuses=Counter();ranges={};texts=defaultdict(Counter);checks=Counter();issues=[]
    def number(v):
        if isinstance(v,bool):return None
        if isinstance(v,(int,float)) or isinstance(v,str) and NUM.fullmatch(v):
            try:return Decimal(str(v))
            except InvalidOperation:return None
    def check_close(d,key,expected,label):
        val=number(d.get(key))
        if val is None:return
        checks[label]+=1
        if not val.is_finite() or not math.isfinite(expected):return
        error=abs(float(val)-expected)/(1+abs(expected))
        if error>1e-7:issues.append({'check':label,'actual':str(val),'expected':expected,'scaled_error':error})
    def walk(v,p='$'):
        counts[type(v).__name__]+=1
        if isinstance(v,dict):
            schemas['|'.join(sorted(v))]+=1
            if isinstance(v.get('status'),str):statuses[v['status']]+=1
            # Only identities declared by labels; ordinary multipliers need
            # section flux and are deliberately not equated off a fixed point.
            g=number(v.get('G'))
            if g is not None and abs(g)<500 and 'multiplier_at_match' in v:
                check_close(v,'multiplier_at_match',math.exp(float(g)),'multiplier_at_match=exp(G)')
            r=number(v.get('r'));rr=number(v.get('return_coordinate'))
            if r is not None and rr is not None and r*rr>0 and 'L' in v:
                check_close(v,'L',math.log(float(rr/r)),'L=log(return_coordinate/r)')
            for k,value in v.items():walk(value,p+'.'+k)
        elif isinstance(v,list):
            for value in v:walk(value,p+'[]')
        else:
            n=number(v)
            if n is not None:
                if not n.is_finite():issues.append({'check':'nonfinite','path':p,'value':str(n)})
                else:
                    if p not in ranges:ranges[p]=[n,n,0]
                    ranges[p][0]=min(ranges[p][0],n);ranges[p][1]=max(ranges[p][1],n);ranges[p][2]+=1
            elif isinstance(v,str):texts[p][v]+=1
    for obj in structured:walk(obj)
    return {'path':str(path.relative_to(REPO)),'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'structured_records':len(structured),'counts':dict(counts),'schemas':dict(schemas),'statuses':dict(statuses),'numeric_ranges':{k:[str(v[0]),str(v[1]),v[2]] for k,v in ranges.items()},'text_values':{k:dict(v) for k,v in texts.items()},'unstructured_lines':unstructured,'parse_errors':errors,'identity_checks':dict(checks),'identity_flags':issues}

def main():
    files=sorted(p for p in REPO.rglob('*') if '.git' not in p.parts and p.suffix in ['.json','.jsonl','.log'])
    rows=[run(p) for p in files]
    (ROOT/'semantic_data_audit.json').write_text(json.dumps(rows,indent=2)+'\n')
    for r in rows:print(json.dumps({k:r[k] for k in ['path','structured_records','statuses','identity_checks','identity_flags']}))
    print('FILES',len(rows),'FLAGS',sum(len(r['identity_flags']) for r in rows),'UNSTRUCTURED_LINES',sum(len(r['unstructured_lines']) for r in rows))
if __name__=='__main__':main()
