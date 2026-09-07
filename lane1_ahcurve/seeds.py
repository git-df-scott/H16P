"""lane1_ahcurve/seeds.py -- the nine fat seeds plus the KKL control.

All coefficient vectors are carried as exact Fractions and only converted to
float at the engine boundary, so a TRIGGER can quote exact rationals.
"""
from fractions import Fraction as F
from engine import local_expand

def _f(x):
    return F(str(x))

# ------------------------------------------------------------- Cherkas 2003 --
# xdot = 1 + x y ; ydot = a00 + a10 x + a20 x^2 + a01 y + a11 x y + a y^2
# a00 = a01 + a11 - a10 - a20 - a ; focus A = (1,-1) ; a11 is a rotating parameter.

CHERKAS_ROWS = [
    # id, a, a20, a11, a01, a10, x_cycles, distribution, second focus x (or None)
    (1, "3",    "-12",  "-1.398",   "8.4",   "15.28", [1.26, 1.98, 3.95], "3",     None),
    (2, "1.5",  "-15",  "0.79993",  "3.2",   "9.17",  [1.4, 1.9, 3.1],    "(3,0)", -0.73),
    (3, "-2",   "12",   "10.999",   "-14",   "-26.1", [0.32, 0.66, 0.8],  "3",     None),
    (4, "-2",   "-1",   "9.49965",  "-12.5", "6.955", [0.56, 0.75, 0.87], "3",     None),
    (5, "-4",   "-1",   "13.9987",  "-21",   "12.4",  [0.63, 0.80, 0.88], "3",     None),
    (6, "5",    "-50",  "-5.49995", "16.5",  "76.45", [1.05, 1.16, 1.5],  "3",     None),
    (7, "8/11", "-12",  "2.1502",   "67/220","-26.5", [1.28, 2.15, 4.43], "(3,1)", -3.2),
    (8, "1.04", "-120", "1.51997",  "1.56",  "-79.6", [1.29, 2.22, 4.63], "(3,1)", -1.79),
]

def cherkas_vec12(a, a20, a11, a01, a10):
    a, a20, a11, a01, a10 = (F(str(v)) for v in (a, a20, a11, a01, a10))
    a00 = a01 + a11 - a10 - a20 - a
    P = [F(1), F(0), F(0), F(0), F(1), F(0)]
    Q = [a00, a10, a01, a20, a11, a]
    return P + Q

# direction of the Cherkas a11 rotating parameter, in LOCAL coordinates at
# A=(1,-1):  d/da11 (P,Q) = (0, 1 + x y) = (0, P), and P|_local = -u + v + u v.
CHERKAS_E11_LOCAL = [F(0)]*6 + [F(0), F(-1), F(1), F(0), F(1), F(0)]

# --------------------------------------------------------------- Perko 1984 --
# P = lam x - y - 10 x^2 + (5+delta) x y + y^2 ; Q = x + x^2 + (8 eps -25 -9 delta) x y
# seed values alpha=-0.0023, lam=-0.005, eps=-0.01, delta=-0.5; focus at (0,0).

def perko_p3_vec12(lam="-0.005", eps="-0.01", delta="-0.5"):
    lam, eps, delta = F(str(lam)), F(str(eps)), F(str(delta))
    P = [F(0), lam, F(-1), F(-10), F(5)+delta, F(1)]
    Q = [F(0), F(1), F(0), F(1), 8*eps - 25 - 9*delta, F(0)]
    return P + Q

PERKO_P3_ALPHA = F("-0.0023")

# ------------------------------------------------------------- KKL control ---
# xdot = y + x^2 + x y ; ydot = a2 x^2 + b2 x y + c2 y^2 + alpha2 x + beta2 y
def kkl_vec12(a2="-10", b2="2.2", c2="0.7", alpha2="-72.7778", beta2="0.0015"):
    a2, b2, c2, alpha2, beta2 = (F(str(v)) for v in (a2, b2, c2, alpha2, beta2))
    P = [F(0), F(0), F(1), F(1), F(1), F(0)]
    Q = [F(0), alpha2, beta2, a2, b2, c2]
    return P + Q

# ------------------------------------------------------------------ seeds ----

def fat_seeds():
    """The nine fat seeds: Cherkas 1-8 and Perko P3."""
    out = []
    for (i, a, a20, a11, a01, a10, xs, dist, sf) in CHERKAS_ROWS:
        v = cherkas_vec12(a, a20, a11, a01, a10)
        x0, y0 = F(1), F(-1)
        loc = local_expand(v, x0, y0)
        outward = 1 if xs[0] > 1 else -1
        out.append(dict(
            name="cherkas%d" % i, family="cherkas", row=i,
            vec12=v, focus=(x0, y0), local=loc,
            direction=(float(outward), 0.0),
            s_cycles=[abs(x - 1.0) for x in xs],
            x_cycles=xs, distribution=dist, base_b=0.0,
            evec_local=CHERKAS_E11_LOCAL, base_t=F(str(a11)),
            second_focus_x=sf,
            params=dict(a=str(a), a20=str(a20), a11=str(a11), a01=str(a01), a10=str(a10)),
        ))
    v = perko_p3_vec12()
    out.append(dict(
        name="perko_p3", family="perko", row=0,
        vec12=v, focus=(F(0), F(0)), local=local_expand(v, F(0), F(0)),
        direction=(0.0, -1.0),
        s_cycles=[0.0425, 0.2160, 1.3838],
        x_cycles=None, distribution="(3,1)", base_b=float(PERKO_P3_ALPHA),
        evec_local=None, base_t=None, second_focus_x=None,
        params=dict(alpha="-0.0023", lam="-0.005", eps="-0.01", delta="-0.5"),
    ))
    return out

def kkl_control():
    v = kkl_vec12()
    return dict(
        name="kkl_control", family="kkl", row=0,
        vec12=v, focus=(F(0), F(0)), local=local_expand(v, F(0), F(0)),
        direction=(1.0, 0.0),
        s_cycles=[0.6832, 2.1837, 15.9628],
        x_cycles=None, distribution="(3,1)", base_b=0.0,
        evec_local=None, base_t=None, second_focus_x=-3711.56,
        params=dict(a2="-10", b2="2.2", c2="0.7", alpha2="-72.7778", beta2="0.0015"),
    )

def exact_strings(vec12):
    return [str(F(v)) for v in vec12]
