"""TASK 3 — pseudo-arclength continuation of the cusp curve (triple limit cycles)
in amplitude at fixed shape (a, a20), logging D_xxx and Perko's nondegeneracy
Jacobians at every accepted point, and watching for a SIGN CHANGE of D_xxx
(= a multiplicity-four limit cycle = Perko's swallow-tail C4).

Ledger: append-only JSONL, one record per accepted point (PROTOCOL rule 5).
"""
import json, os, sys, time, argparse
import mpmath as mp
from engine import Engine, third_order, V7_of, L_of, V1_of
from cusp import Cusp, perko_data, wres, nullvec, solve4, det

mp.mp.dps = 50


def js(v, n=34):
    return mp.nstr(mp.mpf(v), n, strip_zeros=False)


class Ledger:
    def __init__(self, path):
        self.path = path
        self.f = open(path, "a")
        self.n = 0

    def write(self, rec):
        self.f.write(json.dumps(rec) + "\n")
        self.f.flush()
        self.n += 1


def record(c, z, r, tang=None, extra=None):
    mu, x0 = z[:3], z[3]
    M = c.jac_mu(mu, x0)
    pk = perko_data(M) if M else None
    rec = {
        "a": js(c.a), "a20": js(c.a20), "side": c.side,
        "a11": js(mu[0]), "a01": js(mu[1]), "a10": js(mu[2]), "x0": js(x0),
        "r0": js(mp.mpf(x0) - 1),
        "D": js(r["D"], 12), "Dx": js(r["Dx"], 12), "Dxx": js(r["Dxx"], 12),
        "Dxxx": js(r["Dxxx"], 20),
        "res": js(wres([r["D"], r["Dx"], r["Dxx"]], x0), 6),
        "T": js(r["T"], 18), "transv": js(r["transv"], 8), "nsteps": r["nsteps"],
        "V1": js(V1_of(c.a, mu[0], mu[1]), 12),
        "L": js(L_of(c.a, c.a20, mu[0], mu[1], mu[2]), 12),
        "perko": {k: js(v, 12) for k, v in pk.items()} if pk else None,
    }
    if tang is not None:
        rec["tangent"] = [js(t, 12) for t in tang]
    if extra:
        rec.update(extra)
    return rec


def tangent_at(c, z, prev=None):
    J = c.jac_full(z[:3], z[3])
    if J is None:
        return None
    t = nullvec(J)
    if t is None:
        return None
    if prev is not None:
        if sum(t[k] * prev[k] for k in range(4)) < 0:
            t = [-v for v in t]
    elif t[3] < 0:          # start by going OUTWARD in amplitude
        t = [-v for v in t]
    return t


def run(a, a20, side=1, r0_start="0.02", ds0="0.004", dsmax="0.25", dsmin="1e-7",
        maxpts=400, x0max=200.0, tag="", ledger_dir="ledger", verbose=True):
    os.makedirs(ledger_dir, exist_ok=True)
    name = tag or ("a%s_a20%s_s%d" % (js(a, 8), js(a20, 8), side))
    name = name.replace(" ", "").replace("/", "_")
    led = Ledger(os.path.join(ledger_dir, "cusp_%s.jsonl" % name))
    eng = Engine(quad=True)
    c = Cusp(eng, a, a20, side=side)

    summary = {"a": js(a), "a20": js(a20), "side": side,
               "V7_stratum": js(V7_of(a, a20), 16),
               "engine": eng.banner, "start": None, "end": None,
               "end_reason": None, "npts": 0,
               "Dxxx_sign_changes": [], "Dxxx_min_abs": None,
               "ledger": led.path}

    # ---- enter the cusp manifold at small amplitude (Bautin) -------------
    mu0 = third_order(a, a20)
    r0 = mp.mpf(r0_start)
    mu = None
    for trial in (r0, r0 / 2, r0 * 2, r0 / 4):
        mu, r = c.newton_mu(mu0, 1 + trial)
        if mu is not None:
            r0 = trial
            break
    if mu is None:
        summary["end_reason"] = "ENTRY_FAILED:%s" % r
        eng.close()
        return summary

    z = list(mu) + [1 + r0]
    summary["start"] = {"x0": js(z[3]), "mu": [js(v) for v in mu],
                        "Dxxx": js(r["Dxxx"], 16), "res": js(wres([r["D"], r["Dx"], r["Dxx"]], z[3]), 6)}
    led.write(record(c, z, r, extra={"kind": "entry"}))

    tang = tangent_at(c, z)
    if tang is None:
        summary["end_reason"] = "NO_TANGENT_AT_ENTRY"
        eng.close()
        return summary

    ds = mp.mpf(ds0)
    dsmax_, dsmin_ = mp.mpf(dsmax), mp.mpf(dsmin)
    prevDxxx = r["Dxxx"]
    prevz = list(z)
    minabs = abs(r["Dxxx"])
    minat = dict(x0=js(z[3]), mu=[js(v) for v in mu])
    npts = 1
    t0 = time.time()
    fails = 0

    while npts < maxpts:
        zpred = [z[k] + ds * tang[k] for k in range(4)]
        zn, rn = c.newton_arc(zpred, tang, zpred)
        if zn is None:
            fails += 1
            ds /= 3
            if ds < dsmin_:
                summary["end_reason"] = "NEWTON_FAIL:%s (ds<%s)" % (rn, dsmin)
                break
            continue
        # sanity
        if not (abs(zn[3]) < x0max and all(abs(v) < mp.mpf("1e7") for v in zn[:3])):
            summary["end_reason"] = "PARAM_BLOWUP x0=%s" % js(zn[3], 8)
            break
        tn = tangent_at(c, zn, prev=tang)
        if tn is None:
            summary["end_reason"] = "TANGENT_LOST"
            break
        npts += 1
        rec = record(c, zn, rn, tang=tn)
        # ---- D_xxx watch -------------------------------------------------
        if (rn["Dxxx"] > 0) != (prevDxxx > 0):
            rec["EVENT"] = "DXXX_SIGN_CHANGE"
            summary["Dxxx_sign_changes"].append(
                {"between_x0": [js(z[3], 20), js(zn[3], 20)],
                 "Dxxx": [js(prevDxxx, 12), js(rn["Dxxx"], 12)],
                 "npt": npts})
            if verbose:
                print("  *** D_xxx SIGN CHANGE between x0=%s and %s : %s -> %s"
                      % (js(z[3], 12), js(zn[3], 12), js(prevDxxx, 6), js(rn["Dxxx"], 6)))
        if abs(rn["Dxxx"]) < minabs:
            minabs = abs(rn["Dxxx"])
            minat = {"x0": js(zn[3]), "mu": [js(v) for v in zn[:3]],
                     "Dxxx": js(rn["Dxxx"], 16)}
        led.write(rec)
        prevDxxx = rn["Dxxx"]
        prevz, z, tang = list(z), zn, tn
        ds = min(dsmax_, ds * mp.mpf("1.5"))
        if verbose and npts % 10 == 0:
            print("  n=%-4d x0=%-22s a11=%-22s Dxxx=%-22s ds=%-10s %.0fs"
                  % (npts, js(z[3], 14), js(z[0], 14), js(prevDxxx, 10), js(ds, 6), time.time() - t0))
    else:
        summary["end_reason"] = "MAXPTS"

    summary["npts"] = npts
    summary["end"] = {"x0": js(z[3]), "mu": [js(v) for v in z[:3]], "Dxxx": js(prevDxxx, 16)}
    summary["Dxxx_min_abs"] = {"value": js(minabs, 16), **minat}
    summary["wall_s"] = time.time() - t0
    summary["calls"] = eng.ncalls
    eng.close()
    return summary


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--a", required=True)
    ap.add_argument("--a20", required=True)
    ap.add_argument("--side", type=int, default=1)
    ap.add_argument("--maxpts", type=int, default=400)
    ap.add_argument("--tag", default="")
    ap.add_argument("--out", default="")
    a = ap.parse_args()
    s = run(mp.mpf(a.a), mp.mpf(a.a20), side=a.side, maxpts=a.maxpts, tag=a.tag)
    print(json.dumps(s, indent=1))
    if a.out:
        json.dump(s, open(a.out, "w"), indent=1)
