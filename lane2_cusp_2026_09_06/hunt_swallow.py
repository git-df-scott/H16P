"""Seed the square swallow-tail Newton from every cusp point in a set of ledgers.

Unknowns (a11, a01, a10, a20), equations (D, D_x, D_xx, D_xxx) = 0, at fixed
(a, x0).  a20 -- the second shape parameter -- is the fourth unfolding direction,
so the Newton may leave the fixed-(a,a20) cusp curve in exactly the one direction
that can kill D_xxx.

Seeds must be at NORMAL amplitude: at small amplitude D_xxx = 48 d7 r0^4 and d7
vanishes only on the centre variety, where D vanishes identically, so the Newton
provably diverges (see REPORT_lane2.md II.3).
"""
import glob, json, os, sys, time
import mpmath as mp
from engine import Engine
from swallow import try_from_cusp_point

mp.mp.dps = 50

MIN_R0 = mp.mpf("0.15")          # skip the small-amplitude end


def seeds_from(paths, stride=8):
    out = []
    for p in paths:
        rows = [json.loads(l) for l in open(p)]
        for i, r in enumerate(rows):
            if i % stride:
                continue
            if mp.mpf(r["x0"]) - 1 < MIN_R0:
                continue
            out.append(r)
    return out


def main():
    pats = sys.argv[1:-1] or ["ledger/cusp_*.jsonl"]
    outfile = sys.argv[-1] if len(sys.argv) > 1 else "swallow_hunt.json"
    files = []
    for p in pats:
        files.extend(sorted(glob.glob(p)))
    seeds = seeds_from(files)
    print("%d seeds from %d ledgers" % (len(seeds), len(files)), flush=True)
    eng = Engine(quad=True)
    res, best = [], None
    t0 = time.time()
    for i, s in enumerate(seeds):
        out = try_from_cusp_point(eng, mp.mpf(s["a"]), mp.mpf(s["a20"]),
                                  [mp.mpf(s["a11"]), mp.mpf(s["a01"]), mp.mpf(s["a10"])],
                                  mp.mpf(s["x0"]), side=s["side"])
        out["seed_file"] = os.path.basename(files[0]) if len(files) == 1 else None
        out["seed_x0"] = s["x0"]
        res.append(out)
        flag = ""
        if out["status"] == "OK":
            flag = "  <<< SWALLOW-TAIL CANDIDATE"
            if best is None:
                best = out
        print("%3d/%3d a=%-10s x0=%-10s seedDxxx=%-14s -> %s%s"
              % (i + 1, len(seeds), out["a"][:9], out["x0"][:9],
                 (out["seed_Dxxx"] or "")[:13], out["status"][:52], flag), flush=True)
        json.dump(res, open(outfile, "w"), indent=1)
    eng.close()
    nok = sum(1 for r in res if r["status"] == "OK")
    print("\n%d/%d converged to a swallow-tail   (%.0f s, %d engine calls)"
          % (nok, len(res), time.time() - t0, eng.ncalls))


if __name__ == "__main__":
    main()
