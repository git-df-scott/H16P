"""Seed fields for Lane 1, in vec12 form (SEEDS.json convention)."""
import numpy as np
from fractions import Fraction as Fr

# ---------------------------------------------------------------- Cherkas
CHERKAS_ROWS = {
    1: dict(a=Fr(3), a20=Fr(-12), a11=Fr(-1398, 1000), a01=Fr(84, 10), a10=Fr(1528, 100),
            x_cycles=[1.26, 1.98, 3.95], second_focus=None),
    2: dict(a=Fr(15, 10), a20=Fr(-15), a11=Fr(79993, 100000), a01=Fr(32, 10), a10=Fr(917, 100),
            x_cycles=[1.4, 1.9, 3.1], second_focus="x=-0.73"),
    3: dict(a=Fr(-2), a20=Fr(12), a11=Fr(10999, 1000), a01=Fr(-14), a10=Fr(-261, 10),
            x_cycles=[0.32, 0.66, 0.8], second_focus=None),
    4: dict(a=Fr(-2), a20=Fr(-1), a11=Fr(949965, 100000), a01=Fr(-125, 10), a10=Fr(6955, 1000),
            x_cycles=[0.56, 0.75, 0.87], second_focus=None),
    5: dict(a=Fr(-4), a20=Fr(-1), a11=Fr(139987, 10000), a01=Fr(-21), a10=Fr(124, 10),
            x_cycles=[0.63, 0.80, 0.88], second_focus=None),
    6: dict(a=Fr(5), a20=Fr(-50), a11=Fr(-549995, 100000), a01=Fr(165, 10), a10=Fr(7645, 100),
            x_cycles=[1.05, 1.16, 1.5], second_focus=None),
    7: dict(a=Fr(8, 11), a20=Fr(-12), a11=Fr(21502, 10000), a01=Fr(67, 220), a10=Fr(-265, 10),
            x_cycles=[1.28, 2.15, 4.43], second_focus="B=(-3.2,1/3.2)"),
    8: dict(a=Fr(104, 100), a20=Fr(-120), a11=Fr(151997, 100000), a01=Fr(156, 100), a10=Fr(-796, 10),
            x_cycles=[1.29, 2.22, 4.63], second_focus="B=(-1.79,1/1.79)"),
}


def cherkas_vec12(a, a20, a11, a01, a10):
    """xdot = 1 + x y ; ydot = a00 + a10 x + a20 x^2 + a01 y + a11 x y + a y^2."""
    a00 = a01 + a11 - a10 - a20 - a
    return [Fr(0) + 1, Fr(0), Fr(0), Fr(0), Fr(1), Fr(0),
            a00, a10, a01, a20, a11, a]


def cherkas_seed(rid):
    r = CHERKAS_ROWS[rid]
    v = cherkas_vec12(r["a"], r["a20"], r["a11"], r["a01"], r["a10"])
    return dict(name=f"cherkas{rid}", vec12_exact=[str(c) for c in v],
                vec12=np.array([float(c) for c in v]), focus=(1.0, -1.0),
                focus_exact=("1", "-1"), meta=r)


# ------------------------------------------------------------------ Perko P3
def perko_p3_vec12(lam=-0.005, eps=-0.01, delta=-0.5):
    """P = lam x - y - 10x^2 + (5+delta) xy + y^2 ;  Q = x + x^2 + (8eps-25-9delta) xy.
    The rotation alpha of the paper is exactly our uniform rotation b, so the
    base field here is the alpha = 0 member."""
    lam, eps, delta = Fr(lam).limit_denominator(10**9), Fr(eps).limit_denominator(10**9), Fr(delta).limit_denominator(10**9)
    P = [Fr(0), lam, Fr(-1), Fr(-10), Fr(5) + delta, Fr(1)]
    Q = [Fr(0), Fr(1), Fr(0), Fr(1), 8 * eps - 25 - 9 * delta, Fr(0)]
    return P + Q


def perko_p3_seed():
    v = perko_p3_vec12()
    return dict(name="perkoP3", vec12_exact=[str(c) for c in v],
                vec12=np.array([float(c) for c in v]), focus=(0.0, 0.0),
                focus_exact=("0", "0"),
                meta=dict(alpha=-0.0023, cycles_on_neg_y=[-0.0425, -0.2160, -1.3838]))


# ------------------------------------------------------------------- KKL
def kkl_vec12(a2=-10.0, b2=2.2, c2=0.7, alpha2=-72.7778, beta2=0.0015):
    """xdot = y + x^2 + xy ; ydot = a2 x^2 + b2 xy + c2 y^2 + alpha2 x + beta2 y."""
    f = lambda z: Fr(z).limit_denominator(10**9)
    P = [Fr(0), Fr(0), Fr(1), Fr(1), Fr(1), Fr(0)]
    Q = [Fr(0), f(alpha2), f(beta2), f(a2), f(b2), f(c2)]
    return P + Q


def kkl_seed():
    v = kkl_vec12()
    return dict(name="kkl", vec12_exact=[str(c) for c in v],
                vec12=np.array([float(c) for c in v]), focus=(0.0, 0.0),
                focus_exact=("0", "0"),
                meta=dict(origin_cycles_r=[0.6832, 2.1837, 15.9628], remote=-3711.56))


def all_seeds():
    s = [cherkas_seed(i) for i in range(1, 9)]
    s.append(perko_p3_seed())
    return s
