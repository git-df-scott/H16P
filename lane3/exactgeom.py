"""Exact rational primitives for the Lane 3 polygon certifier.

Every function here works with fractions.Fraction only.  No floating point.
Used by certify.py; verify.py deliberately re-implements the same predicates
by different routes so that the two programs are independent checks.

Conventions
-----------
A point is a 2-tuple of Fractions.  A *star polygon* is a pair (O, V) where
O is a rational centre and V = [v_0, ..., v_{n-1}] a rational vertex list; the
polygon is the closed curve v_0 -> v_1 -> ... -> v_{n-1} -> v_0.  It is
required to be star-shaped about O with strictly increasing polar angle and
total winding number 1 about O (predicate `check_star_polygon`).  Under those
conditions the polygon is a simple closed curve, positively (counter-clockwise)
oriented, and O lies in its interior.

For a counter-clockwise simple polygon the outward normal of the edge with
direction d = (dx, dy) is n = (dy, -dx).
"""

from fractions import Fraction as F

ZERO = F(0)


# ---------------------------------------------------------------- vectors

def sub(a, b):
    return (a[0] - b[0], a[1] - b[1])


def cross(a, b):
    return a[0] * b[1] - a[1] * b[0]


def dot(a, b):
    return a[0] * b[0] + a[1] * b[1]


def quadrant(u):
    """Quadrant index 0..3 of a non-zero rational vector.

    0: x>0, y>=0    1: x<=0, y>0    2: x<0, y<=0    3: x>=0, y<0
    The four cases partition R^2 \\ {0}.
    """
    x, y = u
    if x > 0 and y >= 0:
        return 0
    if x <= 0 and y > 0:
        return 1
    if x < 0 and y <= 0:
        return 2
    if x >= 0 and y < 0:
        return 3
    raise ValueError("zero vector has no quadrant")


def check_star_polygon(O, V):
    """Certify that (O, V) is a simple, CCW, star-shaped-about-O closed polygon.

    Returns (ok, reason).  The three conditions checked are

      (a) every direction u_i = v_i - O is non-zero;
      (b) consecutive directions satisfy cross(u_i, u_{i+1}) > 0 and
          dot(u_i, u_{i+1}) > 0, i.e. each angular step lies strictly in
          (0, pi/2);
      (c) the quadrant index advances by a total of exactly 4 around the loop.

    (b) forces the polar angle to increase strictly by less than pi/2 per step,
    so the quadrant index advances by 0 or 1 at each step and (c) pins the total
    turning at exactly 2*pi.  Hence theta is a strictly increasing continuous
    parametrisation of the vertex directions covering the circle once, the
    polygon is star-shaped about O (each open wedge contains exactly one edge,
    crossed once), simple, and positively oriented.
    """
    n = len(V)
    if n < 3:
        return False, "fewer than 3 vertices"
    U = []
    for i, v in enumerate(V):
        u = sub(v, O)
        if u[0] == 0 and u[1] == 0:
            return False, "vertex %d coincides with the centre" % i
        U.append(u)
    adv = 0
    for i in range(n):
        a, b = U[i], U[(i + 1) % n]
        if cross(a, b) <= 0:
            return False, "directions %d,%d not strictly CCW" % (i, (i + 1) % n)
        if dot(a, b) <= 0:
            return False, "angular step %d,%d not below pi/2" % (i, (i + 1) % n)
        qa, qb = quadrant(a), quadrant(b)
        step = (qb - qa) % 4
        if step not in (0, 1):
            return False, "quadrant jump at step %d" % i
        adv += step
    if adv != 4:
        return False, "total winding is %d, not 1" % (adv // 4 if adv % 4 == 0 else -1)
    return True, "ok"


# ------------------------------------------------- edge transversality

def _quad_along_edge(v6, A, D):
    """Coefficients (c0,c1,c2) of  g(t) = G(A + t D)  for the quadratic
    G(x,y) = v6[0] + v6[1] x + v6[2] y + v6[3] x^2 + v6[4] x y + v6[5] y^2."""
    g0, gx, gy, gxx, gxy, gyy = v6
    ax, ay = A
    dx, dy = D
    c0 = g0 + gx * ax + gy * ay + gxx * ax * ax + gxy * ax * ay + gyy * ay * ay
    c1 = gx * dx + gy * dy + 2 * gxx * ax * dx + gxy * (ax * dy + ay * dx) \
        + 2 * gyy * ay * dy
    c2 = gxx * dx * dx + gxy * dx * dy + gyy * dy * dy
    return c0, c1, c2


def edge_normal_form(vecP, vecQ, A, B):
    """Return (c0,c1,c2) with  f(t) = X(A + t(B-A)) . n  for the outward normal
    n = (dy, -dx) of the CCW edge A->B, as a quadratic in t."""
    D = sub(B, A)
    nx, ny = D[1], -D[0]
    p0, p1, p2 = _quad_along_edge(vecP, A, D)
    q0, q1, q2 = _quad_along_edge(vecQ, A, D)
    return (p0 * nx + q0 * ny, p1 * nx + q1 * ny, p2 * nx + q2 * ny)


def quadratic_sign_on_unit_interval(c0, c1, c2):
    """Exact strict sign of c0 + c1 t + c2 t^2 on the closed interval [0,1].

    Returns +1 if the quadratic is > 0 throughout, -1 if < 0 throughout, and
    0 if it vanishes somewhere on [0,1].

    A quadratic attains its extrema on [0,1] at t = 0, t = 1 and, when c2 != 0
    and the vertex t* = -c1/(2 c2) lies in (0,1), at t*.  So the minimum and
    maximum over [0,1] are the min and max of that finite set of exact values.
    """
    vals = [c0, c0 + c1 + c2]
    if c2 != 0:
        ts = F(-c1, 1) / (2 * c2)
        if 0 < ts < 1:
            vals.append(c0 + c1 * ts + c2 * ts * ts)
    lo, hi = min(vals), max(vals)
    if lo > 0:
        return 1
    if hi < 0:
        return -1
    return 0


def polygon_transversality(vecP, vecQ, O, V):
    """Strict sign of X . n_out over the whole boundary, or 0 with the index of
    the first offending edge.  Returns (sign, bad_index_or_None)."""
    n = len(V)
    sign = 0
    for i in range(n):
        A, B = V[i], V[(i + 1) % n]
        c = edge_normal_form(vecP, vecQ, A, B)
        s = quadratic_sign_on_unit_interval(*c)
        if s == 0:
            return 0, i
        if sign == 0:
            sign = s
        elif s != sign:
            return 0, i
    return sign, None


# --------------------------------------------------- radial containment

def star_directions(O, V):
    """Directions, quadrants and a rotation that makes the quadrant blocks
    contiguous.  Requires `check_star_polygon` to have passed, which forces
    every quadrant to be non-empty and the sequence to wrap 3 -> 0 exactly once.
    Returns (U, quads, rot, blocks) with blocks[q] = (lo, hi) in rotated indices."""
    n = len(V)
    U = [sub(v, O) for v in V]
    quads = [quadrant(u) for u in U]
    rot = None
    for i in range(n):
        if quads[i] == 0 and quads[i - 1] == 3:
            rot = i
            break
    if rot is None:
        rot = 0
    rq = [quads[(rot + i) % n] for i in range(n)]
    blocks = [None] * 4
    for q in range(4):
        idx = [i for i in range(n) if rq[i] == q]
        blocks[q] = (idx[0], idx[-1]) if idx else None
    return U, quads, rot, blocks


def _wedge_index(U, quads, rot, blocks, u):
    """Index i (in the ORIGINAL vertex numbering) with u in the closed wedge
    spanned by U[i], U[i+1].  Bisection inside the quadrant block of u."""
    n = len(U)
    q = quadrant(u)
    blk = blocks[q]
    if blk is None:
        for i in range(n):
            if cross(U[i], u) >= 0 and cross(u, U[(i + 1) % n]) >= 0:
                return i
        raise ValueError("direction not covered by the polygon wedges")
    lo, hi = blk
    if cross(U[(rot + lo) % n], u) < 0:
        return (rot + lo - 1) % n
    a, b = lo, hi
    while a < b:
        m = (a + b + 1) // 2
        if cross(U[(rot + m) % n], u) >= 0:
            a = m
        else:
            b = m - 1
    return (rot + a) % n


def _side_of_wedge_chord(O, V, U, quads, rot, blocks, p):
    """Signed value dot(n_out, p - A) for the chord of the wedge containing p-O.
    Negative: strictly inside.  Zero: on the boundary.  Positive: strictly outside."""
    u = sub(p, O)
    if u[0] == 0 and u[1] == 0:
        return F(-1)
    i = _wedge_index(U, quads, rot, blocks, u)
    A, B = V[i], V[(i + 1) % len(V)]
    d = sub(B, A)
    n = (d[1], -d[0])
    return dot(n, sub(p, A))


def strictly_inside(O, V, U, quads, rot, blocks, p):
    return _side_of_wedge_chord(O, V, U, quads, rot, blocks, p) < 0


def strictly_outside(O, V, U, quads, rot, blocks, p):
    return _side_of_wedge_chord(O, V, U, quads, rot, blocks, p) > 0


def polygon_strictly_inside(O, Vin, Vout):
    """Closed region of star polygon (O,Vin) strictly inside the interior of
    (O,Vout).  Both must be certified star-shaped about the same centre O.

    In a wedge of the merged direction set both boundaries are single straight
    chords; a chord lies strictly in the open half-plane bounded by the other
    chord and containing O as soon as both of its endpoints do, by convexity of
    the half-plane.  The two vertex loops below supply exactly those endpoint
    comparisons at every direction of the merged set.
    """
    Uin, Qin, Rin, Bin = star_directions(O, Vin)
    Uout, Qout, Rout, Bout = star_directions(O, Vout)
    for k, p in enumerate(Vin):
        if not strictly_inside(O, Vout, Uout, Qout, Rout, Bout, p):
            return False, "inner vertex %d is not strictly inside the outer polygon" % k
    for k, p in enumerate(Vout):
        if not strictly_outside(O, Vin, Uin, Qin, Rin, Bin, p):
            return False, "outer vertex %d is not strictly outside the inner polygon" % k
    return True, "ok"


# ------------------------------------------------------------- boxes

def bbox(V):
    xs = [v[0] for v in V]
    ys = [v[1] for v in V]
    return (min(xs), max(xs), min(ys), max(ys))


def boxes_disjoint(b1, b2):
    return b1[1] < b2[0] or b2[1] < b1[0] or b1[3] < b2[2] or b2[3] < b1[2]


def quad_range_on_box(v6, box):
    """Exact interval enclosure of the quadratic v6 over the rational box
    (xlo,xhi,ylo,yhi), by interval arithmetic on the monomials."""
    xlo, xhi, ylo, yhi = box
    g0, gx, gy, gxx, gxy, gyy = v6

    def imul(a, b):
        p = (a[0] * b[0], a[0] * b[1], a[1] * b[0], a[1] * b[1])
        return (min(p), max(p))

    def iscale(c, a):
        return (c * a[0], c * a[1]) if c >= 0 else (c * a[1], c * a[0])

    X = (xlo, xhi)
    Y = (ylo, yhi)
    XX = imul(X, X)
    XY = imul(X, Y)
    YY = imul(Y, Y)
    lo = g0
    hi = g0
    for c, iv in ((gx, X), (gy, Y), (gxx, XX), (gxy, XY), (gyy, YY)):
        if c == 0:
            continue
        a, b = iscale(c, iv)
        lo += a
        hi += b
    return lo, hi


def segment_meets_box(A, B, box):
    """Exact: does the closed segment AB meet the closed rational box?
    Separating-axis test on the two box normals and the segment normal."""
    xlo, xhi, ylo, yhi = box
    if max(A[0], B[0]) < xlo or min(A[0], B[0]) > xhi:
        return False
    if max(A[1], B[1]) < ylo or min(A[1], B[1]) > yhi:
        return False
    d = sub(B, A)
    n = (d[1], -d[0])
    if n[0] == 0 and n[1] == 0:
        return xlo <= A[0] <= xhi and ylo <= A[1] <= yhi
    vals = [n[0] * (cx - A[0]) + n[1] * (cy - A[1])
            for cx in (xlo, xhi) for cy in (ylo, yhi)]
    if min(vals) > 0 or max(vals) < 0:
        return False
    return True


def box_vs_polygon(O, V, star, box):
    """Exact position of a closed rational box relative to a certified star
    polygon: "inside", "outside", or "crosses".  `star` is star_directions(O,V)."""
    U, quads, rot, blocks = star
    n = len(V)
    bb = bbox(V)
    if boxes_disjoint(bb, box):
        return "outside"
    for i in range(n):
        if segment_meets_box(V[i], V[(i + 1) % n], box):
            return "crosses"
    corner = (box[0], box[2])
    return "inside" if strictly_inside(O, V, U, quads, rot, blocks, corner) \
        else "outside"


# ---------------------------------------------- integer fast path (exact)
#
# All proposed vertices are dyadic with one common denominator D = 2**bits, and
# the field may be scaled by a positive integer M without changing any sign.
# Writing a vertex as (ax, ay)/D, the edge quadratic
#     f(t) = X(A + t(B-A)) . n_out,      n_out = (dy, -dx)/D
# becomes an integer quadratic divided by the positive constant M * D^3, so the
# integer sign test below is exactly the rational one, only faster.  It is used
# inside the proposal loop; the certificate itself is re-checked with Fractions
# by check_polygon and, independently, by verify.py.

def _gcd(a, b):
    while b:
        a, b = b, a % b
    return a


def integer_field(vecP, vecQ):
    """Scale (P, Q) by a positive integer so that all 12 coefficients are ints.
    Returns (iP, iQ, M)."""
    den = 1
    for c in list(vecP) + list(vecQ):
        d = c.denominator
        den = den * d // _gcd(den, d)
    return ([int(c * den) for c in vecP], [int(c * den) for c in vecQ], den)


def int_edge_coeffs(iP, iQ, D, A, B):
    """Integer (C0, C1, C2) proportional (by M*D^3 > 0) to the coefficients of
    f(t) = X(A + t(B-A)) . n_out for integer-numerator vertices over D."""
    ax, ay = A
    bx, by = B
    dx, dy = bx - ax, by - ay
    p0, p1, p2, p3, p4, p5 = iP
    q0, q1, q2, q3, q4, q5 = iQ
    axD, ayD = ax * D, ay * D
    axax, axay, ayay = ax * ax, ax * ay, ay * ay
    dxdx, dxdy, dydy = dx * dx, dx * dy, dy * dy
    axdx, aydy = ax * dx, ay * dy
    mix = ax * dy + ay * dx
    P0 = p0 * D * D + p1 * axD + p2 * ayD + p3 * axax + p4 * axay + p5 * ayay
    P1 = p1 * D * dx + p2 * D * dy + 2 * p3 * axdx + p4 * mix + 2 * p5 * aydy
    P2 = p3 * dxdx + p4 * dxdy + p5 * dydy
    Q0 = q0 * D * D + q1 * axD + q2 * ayD + q3 * axax + q4 * axay + q5 * ayay
    Q1 = q1 * D * dx + q2 * D * dy + 2 * q3 * axdx + q4 * mix + 2 * q5 * aydy
    Q2 = q3 * dxdx + q4 * dxdy + q5 * dydy
    return (P0 * dy - Q0 * dx, P1 * dy - Q1 * dx, P2 * dy - Q2 * dx)


def int_quadratic_sign_on_unit_interval(C0, C1, C2):
    """Strict sign of C0 + C1 t + C2 t^2 on [0,1]; 0 if it vanishes there."""
    v1 = C0 + C1 + C2
    lo = C0 if C0 < v1 else v1
    hi = C0 if C0 > v1 else v1
    if C2 != 0:
        num, den = -C1, 2 * C2
        inside = (0 < num < den) if den > 0 else (den < num < 0)
        if inside:
            disc = 4 * C0 * C2 - C1 * C1          # f(t*) = disc / (4 C2)
            if disc == 0:
                return 0
            ext_positive = (disc > 0) == (C2 > 0)
            if C2 > 0:                            # vertex is the minimum
                return 1 if ext_positive else (-1 if hi < 0 else 0)
            return -1 if not ext_positive else (1 if lo > 0 else 0)
    if lo > 0:
        return 1
    if hi < 0:
        return -1
    return 0
