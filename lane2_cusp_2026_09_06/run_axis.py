"""Sweep the cusp manifold along a shape axis from seeds in a ledger."""
import json, sys, glob
import mpmath as mp
from engine import Engine
import cusp_axis as ca
mp.mp.dps = 50

def main():
    axis = sys.argv[1]                 # "a20" or "a"
    ledgers = sys.argv[2].split(",")
    x0lo, x0hi = float(sys.argv[3]), float(sys.argv[4])
    maxpts = int(sys.argv[5]); out = sys.argv[6]
    stride = int(sys.argv[7]) if len(sys.argv) > 7 else 40
    seeds = []
    for pat in ledgers:
        for p in sorted(glob.glob(pat)):
            rows = [json.loads(l) for l in open(p)]
            for i, r in enumerate(rows):
                if i % stride: continue
                if not (x0lo < float(r["x0"]) < x0hi): continue
                r["_src"] = p
                seeds.append(r)
    print("%d seeds, axis=%s" % (len(seeds), axis), flush=True)
    e = Engine(quad=True); res = []
    for i, s in enumerate(seeds):
        fixed = ({"a": mp.mpf(s["a"]), "x0": mp.mpf(s["x0"])} if axis == "a20"
                 else {"a20": mp.mpf(s["a20"]), "x0": mp.mpf(s["x0"])})
        w0 = mp.mpf(s["a20"]) if axis == "a20" else mp.mpf(s["a"])
        for d in (+1, -1):
            tag = "%s_%s_i%d_d%+d" % (axis, s["_src"].split("cusp_")[-1][:-6], i, d)
            r = ca.run(e, axis, fixed,
                       [mp.mpf(s["a11"]), mp.mpf(s["a01"]), mp.mpf(s["a10"]), w0],
                       side=s["side"], maxpts=maxpts, direction=d, tag=tag, verbose=False)
            r["tag"] = tag; r["seed_a"] = s["a"]; r["seed_a20"] = s["a20"]; r["seed_x0"] = s["x0"]
            res.append(r); json.dump(res, open(out, "w"), indent=1)
            nsc = len(r.get("Dxxx_sign_changes", []))
            print("%3d %-42s npts=%-4s end=%-14s Dxxx %s -> %s  sc=%d%s"
                  % (i, tag[:42], r["npts"], (r.get("end") or ["?"]*4)[3][:12],
                     (r.get("seed_Dxxx") or "?")[:11], (r.get("end_Dxxx") or "?")[:11], nsc,
                     "   <<<<< SIGN CHANGE" if nsc else ""), flush=True)
    e.close()
    print("total sign changes:", sum(len(r.get("Dxxx_sign_changes", [])) for r in res))

main()
