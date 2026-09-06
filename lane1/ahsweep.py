"""Lane 1: the Andronov-Hopf curve sweep.

Search space.  A base field is carried as the LOCAL 10-vector
    L = (p1,p2,p3,p4,p5, q1,q2,q3,q4,q5)
of the field re-expanded about its focus, so the focus sits at the origin of
the integration chart by construction.  Two directions in this 10-space are
inert for our purposes and are projected out of every perturbation:

  * L itself  -- an overall time rescaling, which changes nothing;
  * V_rot = (-q, +p)  -- the uniform rotation, which is exactly the parameter
    beta that beta*(s) already integrates out.

That leaves 8 live directions.  ||L|| is held fixed (the "one overall scale").

Objective.
  primary   : number of interior extrema of beta*(s)  (3 cycles <=> 2 extrema,
              4 cycles <=> 3 extrema with overlapping height windows)
  secondary : `fold_margin` -- the smallest |d beta*/d log s|, normalised by the
              curve's height range, at an interior local minimum of that
              derivative lying OUTSIDE the interval spanned by the two known
              extrema.  A third pair of extrema is born exactly when this
              reaches zero, so minimising it walks towards a cusp of beta*.
"""
import json, math, os, time
import numpy as np
import engine as E

HERE = os.path.dirname(os.path.abspath(__file__))


# ------------------------------------------------------------------ geometry
def rot_dir(L):
    """d/db of the rotated local vector at b = 0."""
    p, q = L[:5], L[5:]
    return np.concatenate([-q, p])


def project_live(L, v):
    """Remove the scale and rotation components from a perturbation."""
    for d in (np.asarray(L, float), rot_dir(L)):
        nd = np.dot(d, d)
        if nd > 0:
            v = v - d * (np.dot(v, d) / nd)
    return v


def is_focus(L):
    A = np.array([[L[0], L[1]], [L[5], L[6]]])
    tr, dt = A[0, 0] + A[1, 1], A[0, 0] * A[1, 1] - A[0, 1] * A[1, 0]
    return dt > 0 and tr * tr < 4 * dt, tr, dt


def length_scale(L):
    """Radius at which the quadratic part becomes comparable to the linear part."""
    lin = math.hypot(math.hypot(L[0], L[1]), math.hypot(L[5], L[6]))
    quad = math.sqrt(sum(x * x for x in (L[2], L[3], L[4], L[7], L[8], L[9])))
    if quad <= 0:
        return 1.0
    return lin / quad


# ---------------------------------------------------------------- objective
def longest_run(st, b):
    ok = (st == 0) & np.isfinite(b)
    idx = np.where(ok)[0]
    if idx.size == 0:
        return np.array([], int)
    runs, cur = [], [idx[0]]
    for i in idx[1:]:
        if i == cur[-1] + 1:
            cur.append(i)
        else:
            runs.append(cur); cur = [i]
    runs.append(cur)
    return np.array(max(runs, key=len), int)


def curve_features(s, b, st, min_prom_rel=2e-4):
    """Extrema of beta* and the fold margin, on the longest resolved run."""
    run = longest_run(st, b)
    out = dict(n_run=int(run.size), n_extrema=0, extrema=[], height_range=0.0,
               fold_margin=float("inf"), fold_margin_s=None, s_lo=None, s_hi=None)
    if run.size < 15:
        return out
    ss, bb = s[run], b[run]
    out["s_lo"], out["s_hi"] = float(ss[0]), float(ss[-1])
    rng = float(bb.max() - bb.min())
    out["height_range"] = rng
    if rng <= 0:
        return out
    t = np.log(ss)
    B = (bb - bb.min()) / rng
    # derivative on a possibly non-uniform log grid
    G = np.gradient(B, t)
    # turning points of B  (sign changes of G)
    ext = []
    for j in range(1, len(bb) - 1):
        if (bb[j] - bb[j - 1]) * (bb[j + 1] - bb[j]) < 0:
            ext.append([j, "max" if bb[j] > bb[j - 1] else "min", 0.0])
    marks = [0] + [e[0] for e in ext] + [len(bb) - 1]
    for k, e in enumerate(ext):
        j = e[0]
        e[2] = min(abs(B[j] - B[marks[k]]), abs(B[j] - B[marks[k + 2]]))
    ext = [e for e in ext if e[2] > min_prom_rel]
    out["n_extrema"] = len(ext)
    out["extrema"] = [dict(s=float(ss[e[0]]), b=float(bb[e[0]]), kind=e[1],
                           prominence=float(e[2])) for e in ext]

    # fold margin: interior local minima of |G| strictly outside [first,last] extremum
    if len(ext) >= 1:
        lo, hi = ext[0][0], ext[-1][0]
    else:
        lo = hi = len(bb) // 2
    A = np.abs(G)
    best, best_s = float("inf"), None
    guard = 3
    for j in range(guard, len(A) - guard):
        if lo - guard <= j <= hi + guard:
            continue
        if A[j] < A[j - 1] and A[j] < A[j + 1]:
            if A[j] < best:
                best, best_s = float(A[j]), float(ss[j])
    out["fold_margin"] = best
    out["fold_margin_s"] = best_s
    return out


# ------------------------------------------------------------------ evaluate
def ah_domain(L, phi, s_ref, hi_mult=300.0, lo_mult=1e-3, **kw):
    """Bracket the outer end of the domain where the return map exists."""
    lo = lo_mult * s_ref
    hi = hi_mult * s_ref
    D, st, _ = E.d_curve(L, phi, np.array([lo]), **kw)
    if st[0] != 0:
        # walk outwards until the return works
        for m in (1e-2, 1e-1, 1.0):
            lo = m * s_ref
            D, st, _ = E.d_curve(L, phi, np.array([lo]), **kw)
            if st[0] == 0:
                break
        else:
            return None, None
    D, st, _ = E.d_curve(L, phi, np.array([hi]), **kw)
    if st[0] == 0:
        return lo, hi
    a, c = lo, hi
    for _ in range(45):
        m = math.sqrt(a * c)
        D, st, _ = E.d_curve(L, phi, np.array([m]), **kw)
        if st[0] == 0:
            a = m
        else:
            c = m
        if c / a < 1 + 1e-8:
            break
    return lo, a


def evaluate(L, phi=0.0, n=200, refine=True, bmax=1.5, btol=1e-10, **kw):
    """beta*(s) for one base field.  Returns a feature dict."""
    L = np.ascontiguousarray(np.asarray(L, float))
    okf, tr, dt = is_focus(L)
    if not okf:
        return dict(status="not_a_focus", trace=float(tr), det=float(dt),
                    n_extrema=-1, fold_margin=float("inf"))
    s_ref = length_scale(L)
    lo, hi = ah_domain(L, phi, s_ref, **kw)
    if lo is None or not (hi > lo):
        return dict(status="no_domain", n_extrema=-1, fold_margin=float("inf"),
                    s_ref=float(s_ref))
    s = np.geomspace(lo, hi * (1 - 1e-9), n)
    dirh = E.rotation_direction(L, phi, math.sqrt(lo * hi), **kw)
    b, st, d0, nf = E.betastar(L, phi, s, dirhint=dirh, bmax=bmax, btol=btol, **kw)
    feat = curve_features(s, b, st)

    if refine and feat["n_run"] >= 15:
        # densify around every turning point and around the best fold margin
        marks = [e["s"] for e in feat["extrema"]]
        if feat["fold_margin_s"]:
            marks.append(feat["fold_margin_s"])
        extra = []
        span = (math.log(hi) - math.log(lo)) / n
        for m in marks:
            extra.append(np.exp(np.log(m) + np.linspace(-6 * span, 6 * span, 25)))
        if extra:
            s2 = np.unique(np.concatenate([s] + extra))
            s2 = np.ascontiguousarray(s2[(s2 >= lo) & (s2 <= hi)])
            b2, st2, d02, nf2 = E.betastar(L, phi, s2, dirhint=dirh,
                                           bmax=bmax, btol=btol, **kw)
            f2 = curve_features(s2, b2, st2)
            if f2["n_run"] >= feat["n_run"]:
                s, b, st, feat = s2, b2, st2, f2
                nf = np.concatenate([nf, nf2])

    feat.update(status="ok", s_ref=float(s_ref), phi=float(phi),
                trace=float(tr), det=float(dt), n=int(len(s)),
                n_returns=float(np.nansum(nf)), dirhint=int(dirh))
    return feat, s, b, st


def evaluate_only(L, **kw):
    r = evaluate(L, **kw)
    return r[0] if isinstance(r, tuple) else r


# ---------------------------------------------------------------- confirmation
def count_at_level(L, phi, b_level, s_lo, s_hi, n=400, **kw):
    """Sign changes of D(., b_level) in one nest, with the two-tolerance rule."""
    s = np.geomspace(s_lo, s_hi, n)
    D, st, noise, T = E.d_curve_noisy(L, phi, s, b=b_level, **kw)
    br = E.count_sign_changes(s, D, st, noise)
    return br, s, D, st, noise


def overlap_window(extrema):
    """Height window on which a horizontal line meets the curve 4 times,
    given three interior extrema listed in s-order."""
    if len(extrema) < 3:
        return None
    h = [e["b"] for e in extrema]
    lo = max(min(h[0], h[1]), min(h[1], h[2]))
    hi = min(max(h[0], h[1]), max(h[1], h[2]))
    return (lo, hi) if hi > lo else None


# -------------------------------------------------------------------- ledger
class Ledger:
    def __init__(self, path):
        self.path = path
        self.n = 0
        if os.path.exists(path):
            with open(path) as f:
                self.n = sum(1 for _ in f)

    def write(self, rec):
        rec = dict(rec)
        rec.setdefault("engine", E.ENGINE_NAME)
        rec.setdefault("engine_sha256", E.ENGINE_HASH)
        rec.setdefault("t", time.time())
        with open(self.path, "a") as f:
            f.write(json.dumps(rec, default=_jd) + "\n")
        self.n += 1


def _jd(o):
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, np.ndarray):
        return o.tolist()
    return str(o)


def coeff_strings(L):
    """Exact decimal strings (repr round-trips a double exactly) + hex."""
    return dict(dec=[repr(float(x)) for x in L],
                hex=[float(x).hex() for x in L])
