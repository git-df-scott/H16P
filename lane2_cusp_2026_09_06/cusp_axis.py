"""Continuation of the cusp manifold along a chosen axis.

The cusp set {D = D_x = D_xx = 0} is 3 equations in the 6 coordinates
(a, a20, a11, a01, a10, x0), hence 3-dimensional.  continue_cusp.py sweeps it in
the AMPLITUDE direction (unknowns (a11,a01,a10,x0) at fixed shape (a,a20)).
This module sweeps the two SHAPE directions instead:

    axis "a20" : unknowns (a11, a01, a10, a20) at fixed (a, x0)
    axis "a"   : unknowns (a11, a01, a10, a)   at fixed (a20, x0)
    axis "x0"  : the amplitude sweep, for comparison

Why this matters.  On the third-order stratum V7 vanishes on the centre curve
a20_c(a) = 4a(a+1)(a-2)^2/[(a-1)(2a+1)^2], and the SMALL-amplitude sign of D_xxx
flips across it.  The amplitude sweeps at the Cherkas shapes never cross that
curve, so it is unsurprising that they show no sign change.  Sweeping a20 (or a)
at NORMAL amplitude crosses the shape directions directly, and it is a sign change
of D_xxx on THAT sweep, at nonzero amplitude, that is a swallow-tail.
"""
import json, os, time
import mpmath as mp
from engine import Engine, third_order, V7_of, L_of, V1_of
from cusp import solve4, nullvec, wres

mp.mp.dps = 50

FD = mp.mpf("1e-13")
# slot order in z: 0,1,2 are always (a11, a01, a10); slot 3 is the axis
AXES = ("x0", "a20", "a")


def js(v, n=34):
    return mp.nstr(mp.mpf(v), n, strip_zeros=False)


class CuspAxis:
    def __init__(self, eng, axis, fixed, side=1):
        assert axis in AXES
        self.eng = eng
        self.axis = axis
        self.fixed = {k: mp.mpf(v) for k, v in fixed.items()}
        self.side = side

    def unpack(self, z):
        a11, a01, a10, w = z
        a = self.fixed.get("a"); a20 = self.fixed.get("a20"); x0 = self.fixed.get("x0")
        if self.axis == "x0":
            x0 = w
        elif self.axis == "a20":
            a20 = w
        else:
            a = w
        return a, a20, a11, a01, a10, x0

    def val(self, z):
        a, a20, a11, a01, a10, x0 = self.unpack(z)
        return self.eng.D(a, a20, a11, a01, a10, x0, side=self.side)

    def F(self, z):
        r = self.val(z)
        if r["status"] != "OK":
            return None, r
        return [r["D"], r["Dx"], r["Dxx"]], r

    def jac(self, z, r=None):
        """3x4 Jacobian.  For the x0 axis the last column is exact (D_x, D_xx, D_xxx)."""
        J = [[mp.mpf(0)] * 4 for _ in range(3)]
        ncol = 3 if self.axis == "x0" else 4
        for j in range(ncol):
            h = FD * max(mp.mpf(1), abs(z[j]))
            zp = list(z); zp[j] += h
            zm = list(z); zm[j] -= h
            rp = self.val(zp); rm = self.val(zm)
            if rp["status"] != "OK" or rm["status"] != "OK":
                return None
            for i, k in enumerate(("D", "Dx", "Dxx")):
                J[i][j] = (rp[k] - rm[k]) / (2 * h)
        if self.axis == "x0":
            if r is None:
                r = self.val(z)
                if r["status"] != "OK":
                    return None
            for i, k in enumerate(("Dx", "Dxx", "Dxxx")):
                J[i][3] = r[k]
        return J

    def ok(self, z, parmax=mp.mpf("5e3")):
        a, a20, a11, a01, a10, x0 = self.unpack(z)
        if max(abs(t) for t in (a, a20, a11, a01, a10)) > parmax:
            return False
        # NOTE: a = 2 and a = 1/3 are excluded only from the third-order
        # WEAK-FOCUS parametrisation (a11 = 4-2a is degenerate at a=2; W = 3a-1
        # vanishes at a=1/3).  Nothing about the vector field or the cusp manifold
        # is singular there, so the sweep is free to cross them.
        if x0 <= 1 and self.side > 0:
            return False
        if x0 >= 1 and self.side < 0:
            return False
        return L_of(a, a20, a11, a01, a10) > 0            # A must stay an antisaddle

    def newton(self, z, tang, zpred, tol="1e-26", itmax=16, refresh=4):
        z = list(z)
        J = None
        for it in range(itmax):
            F, r = self.F(z)
            if F is None:
                return None, r["status"], None
            if J is None or it % refresh == refresh - 1:
                J = self.jac(z, r)
                if J is None:
                    return None, "jac-fail", None
            elif self.axis == "x0":
                J = [J[i][:3] + [[r["Dx"], r["Dxx"], r["Dxxx"]][i]] for i in range(3)]
            A = J + [list(tang)]
            b = F + [sum(tang[k] * (z[k] - zpred[k]) for k in range(4))]
            s = solve4(A, b)
            if s is None:
                return None, "singular", None
            z = [z[k] - s[k] for k in range(4)]
            if not self.ok(z):
                return None, "out-of-region", None
            sc = max(abs(s[k]) / max(mp.mpf(1), abs(z[k])) for k in range(4))
            if sc < mp.mpf(tol):
                F, r = self.F(z)
                if F is None:
                    return None, r["status"], None
                x0 = self.unpack(z)[5]
                if wres(F, x0) > mp.mpf(tol) * 1000:
                    return None, "residual-too-large", None
                return z, r, J
        return None, "no-converge", None

    def tangent(self, z, prev=None, J=None):
        if J is None:
            J = self.jac(z)
        if J is None:
            return None, None
        t = nullvec(J)
        if t is None:
            return None, None
        if prev is not None and sum(t[k] * prev[k] for k in range(4)) < 0:
            t = [-v for v in t]
        return t, J


def run(eng, axis, fixed, seed_z, side=1, ds0="0.01", dsmax="0.12", dsmin="1e-7",
        maxpts=200, direction=+1, tag="", ledger_dir="ledger_axis", verbose=True):
    """seed_z = (a11, a01, a10, w) already ON the cusp manifold."""
    os.makedirs(ledger_dir, exist_ok=True)
    path = os.path.join(ledger_dir, "axis_%s.jsonl" % (tag or axis))
    f = open(path, "a")
    c = CuspAxis(eng, axis, fixed, side=side)
    z = [mp.mpf(v) for v in seed_z]
    F, r = c.F(z)
    if F is None:
        return {"end_reason": "SEED_EVAL_FAIL:%s" % r["status"], "npts": 0, "ledger": path}
    tang, J = c.tangent(z)
    if tang is None:
        return {"end_reason": "NO_TANGENT", "npts": 0, "ledger": path}
    if (tang[3] > 0) != (direction > 0):
        tang = [-v for v in tang]

    summ = {"axis": axis, "fixed": {k: js(v) for k, v in c.fixed.items()}, "side": side,
            "direction": direction, "seed": [js(v) for v in z],
            "seed_Dxxx": js(r["Dxxx"], 12), "ledger": path,
            "Dxxx_sign_changes": [], "npts": 0, "end_reason": None}
    ds = mp.mpf(ds0); dsmax_, dsmin_ = mp.mpf(dsmax), mp.mpf(dsmin)
    prevD3 = r["Dxxx"]
    npts = 1
    t0 = time.time()
    minabs, minat = abs(prevD3), list(z)

    def rec(z, r, J=None, extra=None):
        a, a20, a11, a01, a10, x0 = c.unpack(z)
        d = {"axis": axis, "side": side, "a": js(a), "a20": js(a20), "a11": js(a11),
             "a01": js(a01), "a10": js(a10), "x0": js(x0),
             "D": js(r["D"], 12), "Dx": js(r["Dx"], 12), "Dxx": js(r["Dxx"], 12),
             "Dxxx": js(r["Dxxx"], 20),
             "Dxxxx": js(r["Dxxxx"], 16) if "Dxxxx" in r else None,
             "res": js(wres([r["D"], r["Dx"], r["Dxx"]], x0), 6),
             "T": js(r["T"], 16), "transv": js(r["transv"], 8),
             "V1": js(V1_of(a, a11, a01), 10), "L": js(L_of(a, a20, a11, a01, a10), 10)}
        if extra:
            d.update(extra)
        return d

    f.write(json.dumps(rec(z, r, extra={"kind": "seed"})) + "\n"); f.flush()

    while npts < maxpts:
        zpred = [z[k] + ds * tang[k] for k in range(4)]
        if not c.ok(zpred):
            summ["end_reason"] = "PRED_OUT_OF_REGION"
            break
        zn, rn, Jn = c.newton(zpred, tang, zpred)
        if zn is None:
            ds /= 3
            if ds < dsmin_:
                summ["end_reason"] = "NEWTON_FAIL:%s" % rn
                break
            continue
        tn, _ = c.tangent(zn, prev=tang, J=Jn)
        if tn is None:
            summ["end_reason"] = "TANGENT_LOST"
            break
        npts += 1
        extra = None
        if (rn["Dxxx"] > 0) != (prevD3 > 0):
            extra = {"EVENT": "DXXX_SIGN_CHANGE"}
            summ["Dxxx_sign_changes"].append(
                {"between": [js(z[3], 20), js(zn[3], 20)],
                 "Dxxx": [js(prevD3, 12), js(rn["Dxxx"], 12)],
                 "npt": npts,
                 "point": [js(v, 30) for v in zn]})
            if verbose:
                print("  *** D_xxx SIGN CHANGE  %s: %s -> %s  at %s = %s"
                      % (axis, js(prevD3, 6), js(rn["Dxxx"], 6), axis, js(zn[3], 14)),
                      flush=True)
        if abs(rn["Dxxx"]) < minabs:
            minabs, minat = abs(rn["Dxxx"]), list(zn)
        f.write(json.dumps(rec(zn, rn, extra=extra)) + "\n"); f.flush()
        prevD3 = rn["Dxxx"]
        z, tang = zn, tn
        ds = min(dsmax_, ds * mp.mpf("1.5"))
        if verbose and npts % 20 == 0:
            print("  n=%-4d %s=%-20s Dxxx=%-20s L=%-10s %.0fs"
                  % (npts, axis, js(z[3], 12), js(prevD3, 10),
                     js(L_of(*[c.unpack(z)[i] for i in (0, 1, 2, 3, 4)]), 6),
                     time.time() - t0), flush=True)
    else:
        summ["end_reason"] = "MAXPTS"
    summ["npts"] = npts
    summ["end"] = [js(v) for v in z]
    summ["end_Dxxx"] = js(prevD3, 12)
    summ["Dxxx_min_abs"] = {"value": js(minabs, 12), "z": [js(v) for v in minat]}
    summ["wall_s"] = time.time() - t0
    f.close()
    return summ
