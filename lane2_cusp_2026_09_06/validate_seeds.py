"""Validation of the Lane-2 binary128 engine against Cherkas-Artes-Llibre rows 1-8
(PROTOCOL rule 7).  Reports every sign change of D on the section y=-1, x>1."""
import sys, os, json, math, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "engine"))
from eng import Engine, geom, x_of_rho, rho_of_x, engine_hash

ROWS = [
 (1, 3.0,      -12.0, -1.398,   8.4,                15.28,  [1.26,1.98,3.95]),
 (2, 1.5,      -15.0,  0.79993, 3.2,                 9.17,  [1.4,1.9,3.1]),
 (3,-2.0,       12.0, 10.999, -14.0,               -26.1,   [0.32,0.66,0.8]),
 (4,-2.0,       -1.0,  9.49965,-12.5,                6.955, [0.56,0.75,0.87]),
 (5,-4.0,       -1.0, 13.9987, -21.0,               12.4,   [0.63,0.80,0.88]),
 (6, 5.0,      -50.0, -5.49995,16.5,                76.45,  [1.05,1.16,1.5]),
 (7, 8/11,     -12.0,  2.1502,  67/220,            -26.5,   [1.28,2.15,4.43]),
 (8, 1.04,    -120.0,  1.51997, 1.56,              -79.6,   [1.29,2.22,4.63]),
]

def scan(e, mu, r_lo, r_hi, n):
    """Return (list of (rho,D) samples, first rho where the return fails)."""
    out, fail = [], None
    for i in range(n+1):
        r = r_lo + (r_hi-r_lo)*i/n
        if r <= 0: continue
        d = e.D(mu, r)
        if not d["ok"]:
            fail = r; break
        out.append((r, d["D"], d["D1"], d["min_den"]))
    return out, fail

def refine(e, mu, r1, r2, iters=80):
    d1 = e.D(mu, r1)["D"]
    for _ in range(iters):
        rm = 0.5*(r1+r2)
        dm = e.D(mu, rm)
        if not dm["ok"]: return None
        if (dm["D"] > 0) == (d1 > 0): r1, d1 = rm, dm["D"]
        else: r2 = rm
        if abs(r2-r1) < 1e-25*max(1.0, abs(r1)): break
    return 0.5*(r1+r2)

def main():
    e = Engine()
    report = []
    for (rid, a, a20, a11, a01, a10, xpub) in ROWS:
        mu = (a, a20, a11, a01, a10)
        g = geom(mu)
        t0 = time.time()
        # scan out to a generous radius; the second-focus / graphic end will fail
        samples, fail = scan(e, mu, 1e-3, 8.0, 400)
        roots = []
        for i in range(len(samples)-1):
            if samples[i][1] == 0.0: continue
            if (samples[i][1] > 0) != (samples[i+1][1] > 0):
                rr = refine(e, mu, samples[i][0], samples[i+1][0])
                if rr is not None: roots.append(rr)
        xs = [x_of_rho(mu, r) for r in roots]
        rec = dict(row=rid, mu=[repr(v) for v in mu], T=g["T"], L=g["L"], w=g["w"],
                   scale=g["scale"], x_cycles_found=[round(v, 6) for v in xs],
                   x_cycles_published=xpub,
                   rho_max_scanned=(fail if fail else 8.0),
                   fail_at_rho=fail, wall_s=round(time.time()-t0, 1))
        report.append(rec)
        print("row %d: found x = %s  published %s  (fail at rho=%s, %.0fs)" %
              (rid, ["%.4f" % v for v in xs], xpub, fail, time.time()-t0))
        sys.stdout.flush()
    e.close()
    out = dict(engine="cusp128.cpp", engine_sha=engine_hash(), tol="1e-28",
               section="y=-1, x>1 (polar chart round the focus A=(1,-1))",
               rows=report)
    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           "data", "validation_cherkas.json"), "w") as fh:
        json.dump(out, fh, indent=1)

main()
