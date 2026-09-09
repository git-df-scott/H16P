"""PROTOCOL rule 7: reproduce the published Andronov-Hopf function of
Cherkas-Artes-Llibre row 4 qualitatively (2 extrema on x in [0.6,0.9])."""
import sys, os, json
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "engine"))
from eng import Engine, rho_of_x, geom

PUB = [8.89863, 4.39482, -13.5991, 22.9703, -22.4248, 11.9886, -2.72941]
def pub(x): return sum(c*x**i for i, c in enumerate(PUB))

A, A20, A01, A10 = -2.0, -1.0, -12.5, 6.955

def Dval(e, a11, x):
    mu = (A, A20, a11, A01, A10)
    g = geom(mu)
    if g is None: return None
    r = rho_of_x(mu, x, side=1)
    if r <= 0: return None
    d = e.D(mu, r, side=1)
    return d["D"] if d["ok"] else None

def AH(e, x, lo=9.4985, hi=9.5015):
    flo, fhi = Dval(e, lo, x), Dval(e, hi, x)
    if flo is None or fhi is None or (flo > 0) == (fhi > 0): return None
    for _ in range(90):
        m = 0.5*(lo+hi); fm = Dval(e, m, x)
        if fm is None: return None
        if (fm > 0) == (flo > 0): lo, flo = m, fm
        else: hi = m
        if hi-lo < 1e-16: break
    return 0.5*(lo+hi)

e = Engine()
xs = [0.60 + 0.005*i for i in range(61)]
rows = []
for x in xs:
    v = AH(e, x)
    rows.append((x, v, pub(x)))
    print("x=%.2f  AH=%s  published_poly=%.6f" %
          (x, ("%.8f" % v) if v is not None else "None", pub(x)))
e.close()
vals = [(x, v) for x, v in [(r[0], r[1]) for r in rows] if v is not None]
ext = []
for i in range(1, len(vals)-1):
    a, b, c = vals[i-1][1], vals[i][1], vals[i+1][1]
    if (b-a)*(c-b) < 0: ext.append((vals[i][0], b, "max" if b > a else "min"))
pex = []
pv = [(x, pub(x)) for x in xs]
for i in range(1, len(pv)-1):
    a, b, c = pv[i-1][1], pv[i][1], pv[i+1][1]
    if (b-a)*(c-b) < 0: pex.append((pv[i][0], b, "max" if b > a else "min"))
print("computed AH extrema on [0.6,0.9]:", ext)
print("published poly extrema on [0.6,0.9]:", pex)
json.dump(dict(x=[r[0] for r in rows], AH=[r[1] for r in rows],
               pub=[r[2] for r in rows], extrema_computed=ext,
               extrema_published=pex),
          open("data/ah_row4.json", "w"), indent=1)
