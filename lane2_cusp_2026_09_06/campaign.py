"""Campaign driver: continue the cusp curve for a list of (a, a20) shapes."""
import json, os, sys, time
import mpmath as mp
import continue_cusp as cc

mp.mp.dps = 50


def main():
    spec = json.load(open(sys.argv[1]))
    out = sys.argv[2]
    res = []
    for s in spec:
        t0 = time.time()
        print("=== %s  a=%s a20=%s side=%d" % (s["tag"], s["a"], s["a20"], s.get("side", 1)), flush=True)
        try:
            r = cc.run(mp.mpf(s["a"]), mp.mpf(s["a20"]), side=s.get("side", 1),
                       maxpts=s.get("maxpts", 400), dsmax=s.get("dsmax", "0.1"),
                       r0_start=s.get("r0", "0.02"), tag=s["tag"],
                       ledger_dir=s.get("ledger_dir", "ledger"), verbose=True)
        except Exception as e:
            r = {"tag": s["tag"], "a": s["a"], "a20": s["a20"], "end_reason": "EXCEPTION:%r" % e}
        r["tag"] = s["tag"]
        r["wall"] = time.time() - t0
        res.append(r)
        print("  -> npts=%s end=%s reason=%s Dxxx_sign_changes=%d  (%.0fs)"
              % (r.get("npts"), (r.get("end") or {}).get("x0"), r.get("end_reason"),
                 len(r.get("Dxxx_sign_changes", [])), r["wall"]), flush=True)
        json.dump(res, open(out, "w"), indent=1)
    json.dump(res, open(out, "w"), indent=1)


if __name__ == "__main__":
    main()
