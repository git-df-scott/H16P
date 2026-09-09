"""Build the (a, a20) grid for the cusp-manifold survey.

Two facts shape the grid.

1. Existence of the third-order weak focus (Cherkas eq.(16)) requires
       (a - 3 - a20) / (1 - 3a) < 0,
   i.e.  a20 < a - 3  for a > 1/3,  and  a20 > a - 3  for a < 1/3.
   a = 2 and a = 1/3 are excluded (a11 = 4-2a degenerate / W = 3a-1 = 0).

2. On that stratum,
       V7 = -150 (a-2) [ -4a(a+1)(a-2)^2 + a20 (a-1)(2a+1)^2 ],
   so V7 vanishes on the CENTRE curve
       a20 = a20_c(a) = 4 a (a+1)(a-2)^2 / [ (a-1)(2a+1)^2 ].
   There V1 = V3 = V5 = V7 = 0, i.e. the system has a centre.

At the small-amplitude end of a cusp curve D_xxx = 48 r0^4 d7 with sgn(d7)
constant on each side of the centre curve, so the ENTRY sign of D_xxx flips
across a20 = a20_c(a).  A swallow-tail needs D_xxx = 0 at NONZERO amplitude, so
the interesting places are (i) where a20_c(a) lies inside the admissible region,
and (ii) wherever the far end of a cusp curve has a different D_xxx sign from
its near end.
"""
import json
import mpmath as mp

mp.mp.dps = 40


def a20_centre(a):
    a = mp.mpf(a)
    if a == 1 or a == mp.mpf("-0.5"):
        return None
    return 4 * a * (a + 1) * (a - 2) ** 2 / ((a - 1) * (2 * a + 1) ** 2)


def admissible(a, a20):
    a, a20 = mp.mpf(a), mp.mpf(a20)
    if a == 2 or a == mp.mpf(1) / 3:
        return False
    return (a - 3 - a20) / (1 - 3 * a) < 0


def centre_reachable(a):
    """Is the centre curve inside the admissible region at this a?"""
    c = a20_centre(a)
    return c is not None and admissible(a, c)


def straddle(a, offsets=("-0.30", "-0.10", "-0.03", "0.03", "0.10", "0.30")):
    """Points on both sides of the centre curve at this a (skipping inadmissible)."""
    c = a20_centre(a)
    if c is None:
        return []
    out = []
    for o in offsets:
        v = c + mp.mpf(o)
        if admissible(a, v):
            out.append(v)
    return out


def survey(a_values, n_a20=5, span="3.0"):
    """A broad sweep of a20 inside the admissible region at each a."""
    out = []
    for a in a_values:
        a = mp.mpf(a)
        if a == 2 or a == mp.mpf(1) / 3:
            continue
        edge = a - 3                       # the admissibility boundary
        for i in range(1, n_a20 + 1):
            d = mp.mpf(span) * i / n_a20
            v = edge - d if a > mp.mpf(1) / 3 else edge + d
            out.append((a, v))
    return out


if __name__ == "__main__":
    import sys
    print("a       a20_centre    admissible?  ")
    for a in ("-4", "-3", "-2.5", "-2", "-1.5", "-1", "-0.5", "0", "0.2", "0.4",
              "0.5", "0.6", "0.72727272727272727", "0.9", "1", "1.04", "1.5", "3", "5"):
        c = a20_centre(a)
        print("%-8s %-14s %s" % (a, "n/a" if c is None else mp.nstr(c, 8),
                                 "-" if c is None else ("REACHABLE" if admissible(a, c) else "outside")))
