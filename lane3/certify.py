"""Lane 3 exact certifier: turn numerically located limit cycles of a planar
quadratic field into a machine-checkable Poincare-Bendixson certificate.

Proof skeleton
--------------
For each located cycle we build two closed polygons with RATIONAL vertices, one
just inside and one just outside the cycle, such that

  (T)  along every edge the quadratic  t -> X(A + t(B-A)) . n_out  has one
       strict sign on [0,1]  (decided exactly: the extrema of a quadratic on a
       closed interval are attained at the endpoints or at its vertex);
  (S)  each polygon is simple, counter-clockwise and star-shaped about a
       rational centre O in its interior (decided exactly by cross/dot products
       and a quadrant-advance count);
  (N)  the inner polygon lies strictly inside the outer one (decided exactly);
  (E)  the closed annulus between them contains no equilibrium (decided exactly
       by rational interval arithmetic on the quadrilaterals that tile it);
  (D)  the annuli are pairwise disjoint (decided exactly).

(T) makes each boundary a transversal curve; if the two signs are (inner +1,
outer -1) the annulus is positively invariant, if they are (inner -1, outer +1)
it is negatively invariant.  With (E), the Poincare-Bendixson theorem gives a
periodic orbit inside each annulus (forward time in the first case, reverse time
in the second).  With (D) the periodic orbits are distinct.

Finally, if every finite equilibrium of the field has non-zero divergence, the
field has no centre; the displacement map on a transversal ray from a strong
focus is then analytic and not identically zero, so its zeros are isolated and
each periodic orbit found above is a LIMIT cycle.

Floating point appears only in the proposal stage (module `propose`).  Every
predicate that enters the certificate is evaluated with fractions.Fraction.
"""

import argparse
import hashlib
import json
import math
import os
import sys
import time
from fractions import Fraction as F

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import exactgeom as G
import fields
import propose as pr

HERE = os.path.dirname(os.path.abspath(__file__))
CERTDIR = os.path.join(HERE, "certificates")
LEDGER = os.path.join(HERE, "ledger", "lane3_runs.jsonl")


# ------------------------------------------------------------- exact checks

def check_polygon(vecP, vecQ, O, V):
    ok, why = G.check_star_polygon(O, V)
    if not ok:
        return None, "star: " + why
    sign, bad = G.polygon_transversality(vecP, vecQ, O, V)
    if sign == 0:
        return None, "edge %s is not strictly transversal" % bad
    return sign, "ok"


def failing_edges(vecP, vecQ, V):
    """Indices of edges on which X . n_out is not of the majority strict sign."""
    n = len(V)
    signs = []
    for i in range(n):
        c = G.edge_normal_form(vecP, vecQ, V[i], V[(i + 1) % n])
        signs.append(G.quadratic_sign_on_unit_interval(*c))
    pos = signs.count(1)
    neg = signs.count(-1)
    want = 1 if pos >= neg else -1
    return [i for i, s in enumerate(signs) if s != want], want


def annulus_free_of_equilibria(vecP, vecQ, O, Vin, Vout, eq_boxes,
                               max_split=30):
    """Exact: no equilibrium lies in the closed annulus K_out minus int K_in.

    `eq_boxes` is a superset of the finite equilibria as rational boxes (from
    the resultants).  For each box we either prove it holds no equilibrium at
    all (0 outside the interval enclosure of P or of Q, after splitting), or
    prove it is disjoint from the annulus by placing it strictly inside the
    inner polygon or strictly outside the outer one.
    """
    star_in = G.star_directions(O, Vin)
    star_out = G.star_directions(O, Vout)
    for box0 in eq_boxes:
        stack = [(box0, 0)]
        while stack:
            box, d = stack.pop()
            lo, hi = G.quad_range_on_box(vecP, box)
            if lo > 0 or hi < 0:
                continue
            lo, hi = G.quad_range_on_box(vecQ, box)
            if lo > 0 or hi < 0:
                continue
            if G.box_vs_polygon(O, Vin, star_in, box) == "inside":
                continue
            if G.box_vs_polygon(O, Vout, star_out, box) == "outside":
                continue
            if d >= max_split:
                return False, "equilibrium not excluded from box %s" % (
                    [str(c) for c in box],)
            xlo, xhi, ylo, yhi = box
            xm, ym = (xlo + xhi) / 2, (ylo + yhi) / 2
            stack += [((xlo, xm, ylo, ym), d + 1), ((xm, xhi, ylo, ym), d + 1),
                      ((xlo, xm, ym, yhi), d + 1), ((xm, xhi, ym, yhi), d + 1)]
    return True, "ok"


# --------------------------------------------------- proposal + refinement
#
# The refinement loop runs on integer vertex numerators over a common dyadic
# denominator, using the integer reformulation in exactgeom; that is the same
# exact decision, only faster.  The polygons handed to the certificate are
# converted back to Fractions and re-checked by check_polygon.

def build_polygon_int(orbit, O_num, bits, taus, reverse):
    sol, T1, T2, _ = orbit
    D = 1 << bits
    ox, oy = O_num[0] / float(D), O_num[1] / float(D)
    pts = pr.blend_points(sol, T1, T2, taus)
    V = [(O_num[0] + int(round((x - ox) * D)), O_num[1] + int(round((y - oy) * D)))
         for x, y in pts]
    return V[::-1] if reverse else V


def to_fractions(Vint, bits):
    D = 1 << bits
    return [(F(x, D), F(y, D)) for x, y in Vint]


def bad_edges_int(iP, iQ, D, O_num, V):
    """Edge indices that break either the star conditions or strict
    transversality, plus the majority transversality sign."""
    n = len(V)
    U = [(v[0] - O_num[0], v[1] - O_num[1]) for v in V]
    bad = []
    signs = []
    for i in range(n):
        a, b = U[i], U[(i + 1) % n]
        if a[0] * b[1] - a[1] * b[0] <= 0 or a[0] * b[0] + a[1] * b[1] <= 0:
            bad.append(i)
            signs.append(0)
            continue
        C = G.int_edge_coeffs(iP, iQ, D, V[i], V[(i + 1) % n])
        signs.append(G.int_quadratic_sign_on_unit_interval(*C))
    want = 1 if signs.count(1) >= signs.count(-1) else -1
    bad += [i for i, s in enumerate(signs) if s != want and i not in bad]
    return sorted(set(bad)), want


def refine_polygon(iP, iQ, O_num, orbit, bits, n0, reverse,
                   max_vertices, max_rounds=60, verbose=False, tag=""):
    D = 1 << bits
    taus = [i / float(n0) for i in range(n0)]
    for rnd in range(max_rounds):
        V = build_polygon_int(orbit, O_num, bits, taus, reverse)
        m = len(taus)
        bad, want = bad_edges_int(iP, iQ, D, O_num, V)
        if not bad:
            return V, taus, {"rounds": rnd, "n_vertices": m, "sign": want}
        if reverse:
            bad = [(m - 2 - i) % m for i in bad]
        if m + len(bad) > max_vertices:
            return None, None, "vertex budget exhausted (%d vertices, %d bad)" % (
                m, len(bad))
        new = []
        for j in sorted(set(bad)):
            lo = taus[j]
            hi = taus[j + 1] if j + 1 < m else 1.0
            mid = 0.5 * (lo + hi)
            if lo < mid < hi:
                new.append(mid)
        if not new:
            return None, None, "bisection underflow at %d vertices" % m
        taus = sorted(set(taus) | set(new))
        if verbose:
            print("      %s refine %d: %d bad -> %d vertices"
                  % (tag, rnd, len(bad), len(taus)))
    return None, None, "refinement did not converge"


def propose_annulus(vecx, vf, O_rat, Of, th0, rstar, dsign, delta_in, delta_out,
                    n0=256, bits=40, max_vertices=400000, verbose=True):
    vecP, vecQ = vecx[:6], vecx[6:]
    iP, iQ, _M = G.integer_field(vecP, vecQ)
    D = 1 << bits
    O_num = (int(O_rat[0] * D), int(O_rat[1] * D))
    if F(O_num[0], D) != O_rat[0] or F(O_num[1], D) != O_rat[1]:
        raise ValueError("centre is not dyadic at this precision")
    orb_out = pr.two_turn_orbit(vf, Of, th0, rstar + delta_out, dsign)
    orb_in = pr.two_turn_orbit(vf, Of, th0, rstar - delta_in, dsign)
    if orb_out is None or orb_in is None:
        return None, None, "two-turn integration failed"
    reverse = dsign < 0
    Vout, _t, iout = refine_polygon(iP, iQ, O_num, orb_out, bits, n0, reverse,
                                    max_vertices, verbose=verbose, tag="outer")
    if Vout is None:
        return None, None, "outer: " + str(iout)
    Vin, _t, iin = refine_polygon(iP, iQ, O_num, orb_in, bits, n0, reverse,
                                  max_vertices, verbose=verbose, tag="inner")
    if Vin is None:
        return None, None, "inner: " + str(iin)
    return (to_fractions(Vin, bits), to_fractions(Vout, bits),
            {"outer": iout, "inner": iin})


# ---------------------------------------------------------- no-centre check

def sturm_chain(p):
    """Sturm chain of a squarefree-ified rational polynomial (list of Fractions,
    lowest degree first)."""
    def deriv(a):
        return [a[i] * i for i in range(1, len(a))]

    def trim(a):
        while a and a[-1] == 0:
            a = a[:-1]
        return a

    def divmod_poly(a, b):
        a = list(a)
        q = [F(0)] * max(1, len(a) - len(b) + 1)
        while len(a) >= len(b) and trim(a):
            d = len(a) - len(b)
            c = a[-1] / b[-1]
            q[d] = c
            for i in range(len(b)):
                a[i + d] -= c * b[i]
            a = trim(a)
        return q, a

    p = trim(list(p))
    chain = [p, trim(deriv(p))]
    while len(chain[-1]) > 1:
        _, r = divmod_poly(chain[-2], chain[-1])
        r = trim([-c for c in r])
        if not r:
            break
        chain.append(r)
    return chain


def poly_eval(a, x):
    s = F(0)
    for c in reversed(a):
        s = s * x + c
    return s


def sign_changes(chain, x):
    vals = [poly_eval(a, x) for a in chain]
    vals = [v for v in vals if v != 0]
    return sum(1 for i in range(len(vals) - 1) if (vals[i] > 0) != (vals[i + 1] > 0))


def isolate_roots(p, lo, hi, depth=80):
    """Disjoint rational intervals, each containing exactly one real root of p
    in [lo, hi]; roots of p at the endpoints are reported as degenerate
    intervals."""
    chain = sturm_chain(p)
    out = []
    stack = [(F(lo), F(hi), 0)]
    while stack:
        a, b, d = stack.pop()
        if poly_eval(p, a) == 0:
            out.append((a, a))
        na = sign_changes(chain, a) - sign_changes(chain, b)
        if na <= 0:
            continue
        if na == 1 and (b - a) < F(1, 10 ** 12):
            out.append((a, b))
            continue
        if d > depth:
            out.append((a, b))
            continue
        m = (a + b) / 2
        stack.append((a, m, d + 1))
        stack.append((m, b, d + 1))
    return out


def resultant_y(vecP, vecQ):
    """Res_y(P, Q) as a polynomial in x with rational coefficients (Sylvester
    4x4 with formal degree 2 in y).  A common zero (x0,y0) of P and Q forces
    this polynomial to vanish at x0."""
    p0, p1, p2, p3, p4, p5 = vecP
    q0, q1, q2, q3, q4, q5 = vecQ
    # P = p5 y^2 + (p2 + p4 x) y + (p0 + p1 x + p3 x^2)
    A2 = [p5]
    A1 = [p2, p4]
    A0 = [p0, p1, p3]
    B2 = [q5]
    B1 = [q2, q4]
    B0 = [q0, q1, q3]

    def pm(a, b):
        r = [F(0)] * (len(a) + len(b) - 1)
        for i, ca in enumerate(a):
            if ca == 0:
                continue
            for j, cb in enumerate(b):
                r[i + j] += ca * cb
        return r

    def pa(a, b):
        n = max(len(a), len(b))
        return [(a[i] if i < len(a) else F(0)) + (b[i] if i < len(b) else F(0))
                for i in range(n)]

    def neg(a):
        return [-c for c in a]

    Z = [F(0)]
    M = [[A2, A1, A0, Z], [Z, A2, A1, A0], [B2, B1, B0, Z], [Z, B2, B1, B0]]

    def det(rows, cols):
        if len(cols) == 1:
            return M[rows[0]][cols[0]]
        tot = Z
        for k, c in enumerate(cols):
            sub = det(rows[1:], [cc for cc in cols if cc != c])
            term = pm(M[rows[0]][c], sub)
            tot = pa(tot, term if k % 2 == 0 else neg(term))
        return tot
    r = det([0, 1, 2, 3], [0, 1, 2, 3])
    while r and r[-1] == 0:
        r = r[:-1]
    return r


def finite_equilibria_boxes(vecP, vecQ, span=F(10 ** 6), width=F(1, 10 ** 9)):
    """Rational boxes whose union contains every finite equilibrium in
    [-span, span]^2, together with the exact interval enclosures of P, Q and of
    the divergence on each box."""
    Rx = resultant_y(vecP, vecQ)
    Ry = resultant_y([vecP[0], vecP[2], vecP[1], vecP[5], vecP[4], vecP[3]],
                     [vecQ[0], vecQ[2], vecQ[1], vecQ[5], vecQ[4], vecQ[3]])
    if len(Rx) <= 1 or len(Ry) <= 1:
        return None, "resultant is constant or identically zero (degenerate field)"
    xs = isolate_roots(Rx, -span, span)
    ys = isolate_roots(Ry, -span, span)
    boxes = []
    for (xa, xb) in xs:
        for (ya, yb) in ys:
            xa2, xb2 = xa - width, xb + width
            ya2, yb2 = ya - width, yb + width
            boxes.append((xa2, xb2, ya2, yb2))
    return boxes, "ok"


def divergence_vec(vecP, vecQ):
    """div X = P_x + Q_y, an affine function, as a vec6."""
    return [vecP[1] + vecQ[2], 2 * vecP[3] + vecQ[4], vecP[4] + 2 * vecQ[5],
            F(0), F(0), F(0)]


def no_centre_certificate(vecP, vecQ, max_split=40):
    """Exact: every finite equilibrium has non-zero divergence, hence the field
    has no centre."""
    boxes, why = finite_equilibria_boxes(vecP, vecQ)
    if boxes is None:
        return False, why, []
    dv = divergence_vec(vecP, vecQ)
    kept = []
    for box in boxes:
        stack = [(box, 0)]
        while stack:
            b, d = stack.pop()
            lo, hi = G.quad_range_on_box(vecP, b)
            if lo > 0 or hi < 0:
                continue
            lo, hi = G.quad_range_on_box(vecQ, b)
            if lo > 0 or hi < 0:
                continue
            lo, hi = G.quad_range_on_box(dv, b)
            if lo > 0 or hi < 0:
                kept.append({"box": [str(c) for c in b],
                             "div_range": [str(lo), str(hi)]})
                continue
            if d >= max_split:
                return False, "divergence not separated from 0 on box %s" % (
                    [str(c) for c in b],), kept
            xlo, xhi, ylo, yhi = b
            xm, ym = (xlo + xhi) / 2, (ylo + yhi) / 2
            stack += [((xlo, xm, ylo, ym), d + 1), ((xm, xhi, ylo, ym), d + 1),
                      ((xlo, xm, ym, yhi), d + 1), ((xm, xhi, ym, yhi), d + 1)]
    return True, "ok", kept


# ------------------------------------------------------------------ driver

def certify_field(vec12, section_points=None, name="field", O_hint=None,
                  th0=0.0, scan=(1e-3, 50.0, 160), bits=40, n0=256,
                  delta_frac=0.35, max_vertices=400000, verbose=True,
                  nests=None):
    """Entry point.  `vec12` is a list of 12 exact rationals (or strings).
    `section_points` optionally gives approximate radii on the ray theta = th0
    from the focus at which cycles are expected; when omitted the displacement
    is scanned over `scan`.

    Returns the certificate dictionary."""
    vecx = [c if isinstance(c, F) else F(str(c)) for c in vec12]
    vecP, vecQ = vecx[:6], vecx[6:]
    vf = pr.as_float(vecx)

    eqs = pr.equilibria(vf)
    foci = [e for e in eqs if pr.classify(vf, e)[0] in ("focus",)]
    if O_hint is not None:
        Of = tuple(float(c) for c in O_hint)
    elif foci:
        Of = foci[0]
    else:
        return {"name": name, "status": "no focus found", "equilibria": eqs}
    O_rat = (pr.dyadic(Of[0], 30), pr.dyadic(Of[1], 30))

    eq_boxes, eq_why = finite_equilibria_boxes(vecP, vecQ)
    if eq_boxes is None:
        return {"name": name, "status": "degenerate field: " + eq_why}

    t_start = time.time()
    # --- locate cycles -------------------------------------------------
    if section_points:
        cycles = []
        for s in section_points:
            lo, hi = s * 0.85, s * 1.15
            try:
                d_lo = pr.displacement(vf, Of, th0, lo)
                d_hi = pr.displacement(vf, Of, th0, hi)
            except Exception:
                continue
            if d_lo is None or d_hi is None or (d_lo < 0) == (d_hi < 0):
                continue
            cycles.append(pr.refine_cycle(vf, Of, th0, lo, hi))
    else:
        _, _, br = pr.scan_cycles(vf, Of, th0, scan[0], scan[1], scan[2])
        cycles = [pr.refine_cycle(vf, Of, th0, a, b) for a, b, _, _ in br]
    cycles = sorted(cycles)
    if verbose:
        print("  focus %.10g,%.10g ; cycles at r = %s"
              % (Of[0], Of[1], ["%.10g" % c for c in cycles]))

    annuli = []
    for k, rstar in enumerate(cycles):
        mu = pr.cycle_multiplier(vf, Of, th0, rstar)
        dsign = 1 if (mu is not None and mu < 1.0) else -1
        gap_lo = rstar - (cycles[k - 1] if k > 0 else 0.0)
        gap_hi = ((cycles[k + 1] - rstar) if k + 1 < len(cycles) else rstar)
        rec = None
        for shrink in (1.0, 0.4, 0.15, 0.05, 0.015, 0.005, 0.0015):
            d_in = delta_frac * gap_lo * shrink
            d_out = delta_frac * gap_hi * shrink
            d = (d_in, d_out)
            Vin, Vout, info = propose_annulus(
                vecx, vf, O_rat, Of, th0, rstar, dsign, d_in, d_out,
                n0=n0, bits=bits, max_vertices=max_vertices, verbose=verbose)
            if Vin is not None:
                rec = (Vin, Vout, d, info)
                break
            if verbose:
                print("    cycle %d delta=(%.3g,%.3g) failed: %s"
                      % (k, d_in, d_out, info))
        if rec is None:
            annuli.append({"index": k, "r": rstar, "status": "FAILED",
                           "multiplier": mu})
            continue
        Vin, Vout, d, info = rec
        s_in, why_in = check_polygon(vecP, vecQ, O_rat, Vin)
        s_out, why_out = check_polygon(vecP, vecQ, O_rat, Vout)
        if s_in is None or s_out is None:
            annuli.append({"index": k, "r": rstar, "status": "FAILED",
                           "why": (why_in, why_out), "multiplier": mu})
            continue
        ok_nest, why_nest = G.polygon_strictly_inside(O_rat, Vin, Vout)
        ok_eq, why_eq = annulus_free_of_equilibria(vecP, vecQ, O_rat, Vin, Vout,
                                                   eq_boxes)
        annuli.append({
            "index": k, "r": rstar, "multiplier": mu, "delta": d,
            "status": "OK" if (ok_nest and ok_eq) else "FAILED",
            "why": None if (ok_nest and ok_eq) else (why_nest, why_eq),
            "sign_inner": s_in, "sign_outer": s_out,
            "n_vertices_inner": len(Vin), "n_vertices_outer": len(Vout),
            "invariance": ("positive" if (s_in == 1 and s_out == -1) else
                           "negative" if (s_in == -1 and s_out == 1) else
                           "INCONCLUSIVE"),
            "Vin": Vin, "Vout": Vout,
        })
        if verbose:
            print("    cycle %d r=%.10g mult=%.6f delta=%.3g vertices=%d "
                  "signs=(%+d,%+d) %s" % (k, rstar, mu or float('nan'), d[1],
                                          len(Vin) + len(Vout), s_in, s_out,
                                          annuli[-1]["invariance"]))

    good = [a for a in annuli if a.get("status") == "OK"
            and a.get("invariance") in ("positive", "negative")]
    # --- pairwise disjointness (nested, common centre) -------------------
    disjoint = True
    disjoint_why = []
    for i in range(len(good) - 1):
        ok, why = G.polygon_strictly_inside(O_rat, good[i]["Vout"],
                                            good[i + 1]["Vin"])
        if not ok:
            disjoint = False
            disjoint_why.append("annuli %d,%d: %s" % (i, i + 1, why))

    nc_ok, nc_why, nc_boxes = no_centre_certificate(vecP, vecQ)

    cert = {
        "format": "lane3-polygon-annulus-1",
        "name": name,
        "generated": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "field": {
            "vec12": [str(c) for c in vecx],
            "order": "P: 1,x,y,x^2,xy,y^2 ; Q: 1,x,y,x^2,xy,y^2 ; xdot=P, ydot=Q",
            "printed": fields.poly_str(vecx),
        },
        "centre": [str(O_rat[0]), str(O_rat[1])],
        "annuli": [{
            "index": a["index"],
            "sign_inner": a["sign_inner"], "sign_outer": a["sign_outer"],
            "invariance": a["invariance"],
            "delta_inner": a["delta"][0], "delta_outer": a["delta"][1],
            "n_vertices_inner": a["n_vertices_inner"],
            "n_vertices_outer": a["n_vertices_outer"],
            "approx_radius": a["r"], "approx_multiplier": a["multiplier"],
            "inner": [[str(p[0]), str(p[1])] for p in a["Vin"]],
            "outer": [[str(p[0]), str(p[1])] for p in a["Vout"]],
        } for a in good],
        "failed": [{k: (str(v) if not isinstance(v, (int, float, type(None))) else v)
                    for k, v in a.items() if k not in ("Vin", "Vout")}
                   for a in annuli if a.get("status") != "OK"],
        "nested_disjoint": disjoint,
        "nested_disjoint_why": disjoint_why,
        "no_centre": {"ok": nc_ok, "why": nc_why, "equilibrium_boxes": nc_boxes},
        "claim": {
            "n_periodic_orbits": len(good) if disjoint else 0,
            "n_limit_cycles": (len(good) if (disjoint and nc_ok) else 0),
        },
        "wall_seconds": round(time.time() - t_start, 2),
    }
    return cert


def write_certificate(cert, outdir=CERTDIR):
    os.makedirs(outdir, exist_ok=True)
    path = os.path.join(outdir, "CERTIFICATE_%s.json" % cert["name"])
    with open(path, "w") as fh:
        json.dump(cert, fh, indent=1, sort_keys=True)
    h = hashlib.sha256(open(path, "rb").read()).hexdigest()[:16]
    return path, h


def append_ledger(rec):
    os.makedirs(os.path.dirname(LEDGER), exist_ok=True)
    with open(LEDGER, "a") as fh:
        fh.write(json.dumps(rec, sort_keys=True) + "\n")
