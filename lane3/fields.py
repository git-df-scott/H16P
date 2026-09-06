"""Exact rational quadratic vector fields used for Lane 3 validation.

vec12 layout (as in coordination_2026_09_06/SEEDS.json):

    P = c0 + c1 x + c2 y + c3 x^2 + c4 x y + c5 y^2
    Q = c6 + c7 x + c8 y + c9 x^2 + c10 x y + c11 y^2

with xdot = P, ydot = Q.  Every coefficient is a fractions.Fraction, so the
field certified is EXACTLY the field described by the certificate.
"""

from fractions import Fraction as F


def _f(v):
    return v if isinstance(v, F) else F(str(v))


def cherkas(a, a20, a11, a01, a10):
    """Cherkas-Artes-Llibre normal form (paper eq. (2)):
        xdot = 1 + x y
        ydot = a00 + a10 x + a20 x^2 + a01 y + a11 x y + a y^2,
        a00 = a01 + a11 - a10 - a20 - a.
    Antisaddle A = (1, -1); the paper's section is y = -1, x > 1."""
    a, a20, a11, a01, a10 = map(_f, (a, a20, a11, a01, a10))
    a00 = a01 + a11 - a10 - a20 - a
    P = [F(1), F(0), F(0), F(0), F(1), F(0)]
    Q = [a00, a10, a01, a20, a11, a]
    return P + Q


CHERKAS_ROWS = {
    1: dict(a=3,     a20=-12,  a11="-1.398",   a01="8.4",  a10="15.28",
            x_cycles=[1.26, 1.98, 3.95]),
    2: dict(a="1.5", a20=-15,  a11="0.79993",  a01="3.2",  a10="9.17",
            x_cycles=[1.4, 1.9, 3.1]),
    3: dict(a=-2,    a20=12,   a11="10.999",   a01=-14,    a10="-26.1",
            x_cycles=[0.32, 0.66, 0.8]),
    4: dict(a=-2,    a20=-1,   a11="9.49965",  a01="-12.5", a10="6.955",
            x_cycles=[0.56, 0.75, 0.87]),
    5: dict(a=-4,    a20=-1,   a11="13.9987",  a01=-21,    a10="12.4",
            x_cycles=[0.63, 0.80, 0.88]),
    6: dict(a=5,     a20=-50,  a11="-5.49995", a01="16.5", a10="76.45",
            x_cycles=[1.05, 1.16, 1.5]),
    7: dict(a=F(8, 11), a20=-12, a11="2.1502", a01=F(67, 220), a10="-26.5",
            x_cycles=[1.28, 2.15, 4.43]),
    8: dict(a="1.04", a20=-120, a11="1.51997", a01="1.56", a10="-79.6",
            x_cycles=[1.29, 2.22, 4.63]),
}


def cherkas_row(n):
    r = CHERKAS_ROWS[n]
    return cherkas(r["a"], r["a20"], r["a11"], r["a01"], r["a10"])


def kkl(a2="-10", b2="2.2", c2="0.7", alpha2="-72.7778", beta2="0.0015"):
    """Kuznetsov-Kuznetsov-Levitin control from SEEDS.json:
        xdot = y + x^2 + x y
        ydot = a2 x^2 + b2 x y + c2 y^2 + alpha2 x + beta2 y
    (3,1) distribution; origin nest cycles at r ~ 0.6832, 2.1837, 15.9628."""
    a2, b2, c2, alpha2, beta2 = map(_f, (a2, b2, c2, alpha2, beta2))
    P = [F(0), F(0), F(1), F(1), F(1), F(0)]
    Q = [F(0), alpha2, beta2, a2, b2, c2]
    return P + Q


def rotate(vec12, cos_b, sin_b):
    """Uniform rotation X_b = (P cos b - Q sin b, P sin b + Q cos b).
    cos_b, sin_b must be exact rationals with cos^2 + sin^2 = 1 for the rotation
    to be an isometry; the returned field is quadratic and has the same
    equilibria as the base field for any (cos_b, sin_b) != (0,0)."""
    c, s = _f(cos_b), _f(sin_b)
    P, Q = vec12[:6], vec12[6:]
    return [c * P[i] - s * Q[i] for i in range(6)] + \
           [s * P[i] + c * Q[i] for i in range(6)]


def rational_rotation(alpha, denom=10 ** 18):
    """Exact rational (cos, sin) on the unit circle close to the angle alpha,
    via the half-angle parametrisation cos = (1-t^2)/(1+t^2), sin = 2t/(1+t^2)
    with t a rational approximation of tan(alpha/2).  cos^2 + sin^2 = 1 holds
    exactly, so this is a genuine rotation by 2*arctan(t) ~= alpha."""
    import math
    t = F(round(math.tan(alpha / 2.0) * denom), denom)
    d = 1 + t * t
    return (1 - t * t) / d, 2 * t / d


def perko_p3(alpha=-0.0023, lam="-0.005", eps="-0.01", delta="-0.5"):
    """Perko 1984 system P3 (SEEDS.json), rotated by a rational rotation close
    to alpha.  Base field
        P = lam x - y - 10 x^2 + (5+delta) x y + y^2
        Q = x + x^2 + (8 eps - 25 - 9 delta) x y
    Returns (vec12, cos_b, sin_b): the certified field is the rotation by the
    exact rational angle 2*arctan(t), which differs from Perko's alpha by about
    1e-18."""
    lam, eps, delta = map(_f, (lam, eps, delta))
    P = [F(0), lam, F(-1), F(-10), F(5) + delta, F(1)]
    Q = [F(0), F(1), F(0), F(1), 8 * eps - 25 - 9 * delta, F(0)]
    c, s = rational_rotation(alpha)
    return rotate(P + Q, c, s), c, s


def one_cycle_field(mu="1/2"):
    """A constructed control with a single limit cycle: the quadratic Lienard
    system  xdot = y,  ydot = -x + mu y - mu x y.  Its divergence mu(1-x) has a
    single sign change on the line x = 1, so by the Dulac/Bendixson criterion in
    each of the half planes there is no closed orbit lying entirely in one of
    them; the classical analysis of this family gives at most one limit cycle
    around the origin."""
    mu = _f(mu)
    P = [F(0), F(0), F(1), F(0), F(0), F(0)]
    Q = [F(0), F(-1), mu, F(0), -mu, F(0)]
    return P + Q


def poly_str(vec12):
    def side(v, name):
        mons = ["1", "x", "y", "x^2", "x*y", "y^2"]
        parts = [("(%s)*%s" % (c, m) if m != "1" else "(%s)" % c)
                 for c, m in zip(v, mons) if c != 0]
        return "%s = %s" % (name, " + ".join(parts) if parts else "0")
    return side(vec12[:6], "P") + " ; " + side(vec12[6:], "Q")
