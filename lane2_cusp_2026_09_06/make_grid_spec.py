"""Emit the (a, a20) grid specs for the cusp-manifold survey, split across N workers."""
import json, sys
import mpmath as mp
from grid import a20_centre, admissible, straddle

mp.mp.dps = 40

NW = int(sys.argv[1]) if len(sys.argv) > 1 else 4
MAXPTS = int(sys.argv[2]) if len(sys.argv) > 2 else 160
PREFIX = sys.argv[3] if len(sys.argv) > 3 else "grid"

# a values where the centre curve V7 = 0 is INSIDE the admissible region.
A_CENTRE = ["-2.9", "-2.5", "-2.0", "-1.5", "-1.0", "-0.8", "-0.2", "0.0",
            "0.2", "0.4", "0.5", "0.6", "0.72727272727272727272727272727273",
            "0.8", "0.9", "0.95"]
# offsets on both sides of the centre curve
OFF = ["-1.0", "-0.3", "-0.08", "0.08", "0.3", "1.0"]

# a values where the centre curve is NOT reachable: sweep the admissible band
A_PLAIN = ["-4", "-3", "1.04", "1.5", "2.5", "3", "5", "8"]
BAND = ["0.5", "2.0", "6.0", "20.0"]

jobs = []
for a in A_CENTRE:
    c = a20_centre(a)
    for o in OFF:
        v = c + mp.mpf(o)
        if not admissible(a, v):
            continue
        jobs.append({"tag": "c_a%s_o%s" % (a[:7].replace(".", "p").replace("-", "m"),
                                           o.replace(".", "p").replace("-", "m")),
                     "a": mp.nstr(mp.mpf(a), 34), "a20": mp.nstr(v, 34),
                     "side": 1, "maxpts": MAXPTS, "dsmax": "0.15",
                     "ledger_dir": "ledger_grid"})
for a in A_PLAIN:
    av = mp.mpf(a)
    edge = av - 3
    for b in BAND:
        v = edge - mp.mpf(b) if av > mp.mpf(1) / 3 else edge + mp.mpf(b)
        if not admissible(av, v):
            continue
        jobs.append({"tag": "p_a%s_b%s" % (a.replace(".", "p").replace("-", "m"),
                                           b.replace(".", "p")),
                     "a": mp.nstr(av, 34), "a20": mp.nstr(v, 34),
                     "side": 1, "maxpts": MAXPTS, "dsmax": "0.15",
                     "ledger_dir": "ledger_grid"})

# de-duplicate by tag
seen, uniq = set(), []
for j in jobs:
    if j["tag"] in seen:
        continue
    seen.add(j["tag"]); uniq.append(j)

groups = [[] for _ in range(NW)]
for i, j in enumerate(uniq):
    groups[i % NW].append(j)
for i, g in enumerate(groups):
    json.dump(g, open("%s_g%d.json" % (PREFIX, i), "w"), indent=1)
print("%d jobs -> %s" % (len(uniq), [len(g) for g in groups]))
