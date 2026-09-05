"""Recount, with the fixed compactified counter, every record with total >= 3 in the given jsonl files."""
import sys, json, numpy as np
sys.argv_files = sys.argv[1:]; sys.argv = ['x', 'mvneutral', '/dev/null', '0', '0']
import importlib.util; spec = importlib.util.spec_from_file_location('sl', 'sweep_log.py'); sl = importlib.util.module_from_spec(spec); spec.loader.exec_module(sl)
out = open('data/RECOUNT_fixed_counter.jsonl', 'a'); hist = {}
for f in sys.argv_files:
    for line in open(f):
        r = json.loads(line)
        if r.get('total', 0) < 3: continue
        nests = sl.evaluate(np.array(r['coef'])); tot = sum(len(n['roots']) for n in nests)
        hist[(r['total'], tot)] = hist.get((r['total'], tot), 0)+1
        rec = dict(source=f, old_total=r['total'], new_total=tot, coef=r['coef'], nests=nests); out.write(json.dumps(rec)+"\n"); out.flush()
        flag = " <-- INCREASED" if tot > r['total'] else ""
        print(f"{f}: old {r['total']} -> new {tot} {[(''.join(n['stab']), [f'{x:.3g}' for x in n['roots']]) for n in nests]}{flag}", flush=True)
        if tot >= 5: print("FIVE OR MORE", rec, flush=True)
print("DONE", {f"{k[0]}->{k[1]}": v for k, v in hist.items()})
