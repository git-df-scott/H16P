"""Land exactly on the cusp manifold instead of searching for it.

PROTOCOL section (c): triple limit cycles are KNOWN to exist at small amplitude
in the Bautin unfolding of a third-order weak focus.  With

        D(r) ~ V1 r + V3 r^3 + V5 r^5 + V7 r^7,

the choice

        D = V7 r (r^2 - r0^2)^3   <=>   V5 = -3 r0^2 V7,
                                        V3 =  3 r0^4 V7,
                                        V1 = -  r0^6 V7

puts a cycle of multiplicity three at radius r0 -- a CUSP of beta*.  Nobody has
continued that cusp manifold out to normal amplitude.

Cherkas-Artes-Llibre give the focal values of their normal form in closed form
(LIT_A section 3.1, their eq. (15)), and the three conditions are almost
triangular in their parameters, so no search is needed at all:

  * V1 = a11 + a01 - 2a - 1                       -- linear in a01
  * V3 = W0(a, a20, a11) - a10 * W(a, a11)        -- linear in a10
  * V5 = (4 - 2a - a11) V / W                     -- depends on a11, a20, a
  * V7 = -(a11 + 2a + 1) U V / W                  -- depends on a11, a20, a

so for a fixed (a, a20) and a target r0 one scalar equation
`V5 + 3 r0^2 V7 = 0` fixes a11, and then a10 and a01 follow by substitution.
The whole cusp manifold is then swept by moving (a, a20, r0).

The r of the Bautin expansion is not the section coordinate s, and the focal
values are given up to the paper's own normalisation, so r0 is a dial rather
than a length: what is checked numerically afterwards is the shape of beta*.
"""
import math
import numpy as np


def focal_values(a, a20, a11, a01, a10):
    """V1, V3, V5, V7 of Cherkas-Artes-Llibre eq. (15) at A = (1,-1)."""
    W0 = (a11 ** 2 * (a + 1) + a11 * (2 * a ** 2 + a - 1)
          - a20 * (a11 * (2 * a - 1) + (2 * a + 1) * (2 * a - 3)))
    W = -1 + 2 * a ** 2 + a * (a11 - 1)
    V = -a11 ** 2 * a * (a + 1) + a20 * (a - 1) * (2 * a + 1) ** 2
    U = ((8 - 2 * a ** 2) * (a11 + 2 * a + 1) ** 2
         - 35 * (2 * a + 1) * (a11 + 2 * a + 1) + 35 * (2 * a + 1) ** 2)
    if W == 0:
        return None
    V1 = a11 + a01 - 2 * a - 1
    V3 = W0 - a10 * W
    V5 = (4 - 2 * a - a11) * V / W
    V7 = -(a11 + 2 * a + 1) * U * V / W
    return V1, V3, V5, V7


def _g(a, a20, a11, r0):
    """V5 + 3 r0^2 V7, the one scalar equation that fixes a11."""
    fv = focal_values(a, a20, a11, 0.0, 0.0)
    if fv is None:
        return None
    _v1, _v3, V5, V7 = fv
    return V5 + 3.0 * r0 * r0 * V7


def triple_cycle_field(a, a20, r0, a11_lo=-30.0, a11_hi=30.0, n=4001):
    """Solve the three cusp conditions.  Returns a list of parameter dicts,
    one per real root a11 of V5 + 3 r0^2 V7 = 0."""
    grid = np.linspace(a11_lo, a11_hi, n)
    vals = np.array([_g(a, a20, x, r0) if _g(a, a20, x, r0) is not None else np.nan
                     for x in grid])
    out = []
    for i in range(len(grid) - 1):
        f1, f2 = vals[i], vals[i + 1]
        if not (np.isfinite(f1) and np.isfinite(f2)) or f1 * f2 >= 0:
            continue
        lo, hi = grid[i], grid[i + 1]
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            fm = _g(a, a20, mid, r0)
            if fm is None:
                break
            if f1 * fm < 0:
                hi = mid
            else:
                lo, f1 = mid, fm
            if hi - lo < 1e-15 * (1 + abs(hi)):
                break
        a11 = 0.5 * (lo + hi)
        fv = focal_values(a, a20, a11, 0.0, 0.0)
        if fv is None:
            continue
        _v1, W0_minus_0, _V5, V7 = fv           # V3 at a10 = 0 is W0
        W = -1 + 2 * a ** 2 + a * (a11 - 1)
        if W == 0 or V7 == 0:
            continue
        # V3 = W0 - a10 W  =  3 r0^4 V7   ->  a10
        a10 = (W0_minus_0 - 3.0 * r0 ** 4 * V7) / W
        # V1 = a11 + a01 - 2a - 1  =  -r0^6 V7  ->  a01
        a01 = -(r0 ** 6) * V7 + 2 * a + 1 - a11
        chk = focal_values(a, a20, a11, a01, a10)
        out.append(dict(a=a, a20=a20, r0=r0, a11=a11, a01=a01, a10=a10,
                        V=list(chk), V7=V7,
                        residual=[abs(chk[0] + r0 ** 6 * V7),
                                  abs(chk[1] - 3 * r0 ** 4 * V7),
                                  abs(chk[2] + 3 * r0 ** 2 * V7)]))
    return out


def vec12(p):
    """The Cherkas normal form as a 12-vector; focus at A = (1,-1)."""
    a, a20, a11, a01, a10 = p["a"], p["a20"], p["a11"], p["a01"], p["a10"]
    a00 = a01 + a11 - a10 - a20 - a
    return np.array([1.0, 0.0, 0.0, 0.0, 1.0, 0.0,
                     a00, a10, a01, a20, a11, a], float)


def unfold(p, e1=0.0, e3=0.0):
    """Move off the cusp inside the Bautin unfolding, exactly and for free.

    a10 enters V3 linearly and a01 enters V1 linearly, with V5 and V7 untouched,
    so the two unfolding directions of a multiplicity-three cycle are available
    in closed form:  V3 -> 3 r0^4 V7 (1 + e3),  V1 -> -r0^6 V7 (1 + e1).
    The (e1, e3) plane at fixed (V5, V7) is the classical Bautin picture; the
    triple root at the origin opens into three simple roots in a cusped wedge.
    """
    a, a20, a11, r0, V7 = p["a"], p["a20"], p["a11"], p["r0"], p["V7"]
    W = -1 + 2 * a ** 2 + a * (a11 - 1)
    fv = focal_values(a, a20, a11, 0.0, 0.0)
    W0 = fv[1]                                  # V3 evaluated at a10 = 0
    a10 = (W0 - 3.0 * r0 ** 4 * V7 * (1.0 + e3)) / W
    a01 = -(r0 ** 6) * V7 * (1.0 + e1) + 2 * a + 1 - a11
    q = dict(p)
    q.update(a10=a10, a01=a01, e1=e1, e3=e3,
             V=list(focal_values(a, a20, a11, a01, a10)))
    return q
