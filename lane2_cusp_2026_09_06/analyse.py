"""Summarise cusp-curve ledgers: endpoints, D_xxx sign along each curve,
closest approach of D_xxx to zero, and Perko nondegeneracy margins."""
import glob, json, os, sys
import mpmath as mp

mp.mp.dps = 40


def summarise(path):
    rows = [json.loads(l) for l in open(path)]
    if not rows:
        return None
    d = [mp.mpf(r["Dxxx"]) for r in rows]
    x = [mp.mpf(r["x0"]) for r in rows]
    signs = [1 if v > 0 else -1 for v in d]
    changes = [i for i in range(1, len(signs)) if signs[i] != signs[i - 1]]
    imin = min(range(len(d)), key=lambda i: abs(d[i]))

    # ---- the CORRECT swallow-tail detector -------------------------------
    # {D = D_x = D_xx = 0} contains the CENTRE VARIETY as a spurious component
    # of the same dimension: where the field has a centre, D vanishes
    # identically and the three cusp equations hold trivially.  Crossing it
    # flips the sign of D_xxx -- and of D_xxxx, and of D itself -- all at once,
    # so a bare D_xxx sign watch fires on centres, not on swallow-tails.
    # A genuine multiplicity-four cycle needs D_xxx = 0 with D_xxxx != 0, i.e.
    # a sign change of  nu = D_xxx / (D_xxxx * r0).
    nu = [mp.mpf(r["nu"]) if r.get("nu") else None for r in rows]
    nus = [None if v is None else (1 if v > 0 else -1) for v in nu]
    nuch = [i for i in range(1, len(nus))
            if nus[i] is not None and nus[i - 1] is not None and nus[i] != nus[i - 1]]
    have_nu = any(v is not None for v in nu)
    numin = (min((i for i in range(len(nu)) if nu[i] is not None),
                 key=lambda i: abs(nu[i])) if have_nu else None)
    return {
        "file": os.path.basename(path),
        "a": rows[0]["a"], "a20": rows[0]["a20"], "side": rows[0]["side"],
        "npts": len(rows),
        "x0_start": mp.nstr(x[0], 12), "x0_end": mp.nstr(x[-1], 12),
        "x0_max": mp.nstr(max(x), 12), "x0_min": mp.nstr(min(x), 12),
        "Dxxx_start": mp.nstr(d[0], 8), "Dxxx_end": mp.nstr(d[-1], 8),
        "sign_start": signs[0], "sign_end": signs[-1],
        "n_sign_changes": len(changes),
        "sign_change_at": [{"x0": mp.nstr(x[i], 16),
                            "Dxxx_before": mp.nstr(d[i - 1], 8),
                            "Dxxx_after": mp.nstr(d[i], 8)} for i in changes],
        "Dxxx_closest_to_zero": {
            "abs": mp.nstr(abs(d[imin]), 10), "x0": mp.nstr(x[imin], 20),
            "a11": rows[imin]["a11"], "a01": rows[imin]["a01"],
            "a10": rows[imin]["a10"], "index": imin,
            "ratio_to_start": mp.nstr(abs(d[imin]) / abs(d[0]), 6) if d[0] != 0 else None},
        "n_nu_sign_changes": len(nuch),
        "nu_sign_change_at": [{"x0": mp.nstr(x[i], 16),
                               "nu_before": mp.nstr(nu[i - 1], 8),
                               "nu_after": mp.nstr(nu[i], 8),
                               "Dxxxx": rows[i].get("Dxxxx")} for i in nuch],
        "nu_start": (mp.nstr(nu[0], 8) if have_nu and nu[0] is not None else None),
        "nu_end": (mp.nstr(nu[-1], 8) if have_nu and nu[-1] is not None else None),
        "nu_closest_to_zero": ({"abs": mp.nstr(abs(nu[numin]), 8),
                                "x0": mp.nstr(x[numin], 20),
                                "a11": rows[numin]["a11"], "a01": rows[numin]["a01"],
                                "a10": rows[numin]["a10"], "index": numin}
                               if numin is not None else None),
        "V1_end": rows[-1]["V1"], "L_end": rows[-1]["L"],
        "L_min": mp.nstr(min(mp.mpf(r["L"]) for r in rows), 8),
        "res_max": mp.nstr(max(mp.mpf(r["res"]) for r in rows), 6),
        "perko_min_abs_min": (mp.nstr(min(mp.mpf(r["perko"]["min_abs"])
                                          for r in rows if r.get("perko")), 6)
                              if rows[0].get("perko") else None),
    }


def main():
    pats = sys.argv[1:] or ["ledger/cusp_*.jsonl"]
    files = []
    for p in pats:
        files.extend(sorted(glob.glob(p)))
    out = []
    for f in files:
        s = summarise(f)
        if s:
            out.append(s)
    hdr = "%-30s %5s %-10s %-10s %-11s %-11s %4s %-11s %4s" % (
        "curve", "npts", "x0 start", "x0 end", "nu start", "nu end", "D3sc",
        "min|nu|", "NUsc")
    print(hdr); print("-" * len(hdr))
    for s in out:
        print("%-30s %5d %-10.5g %-10.5g %-11s %-11s %4d %-11s %4d"
              % (s["file"][:30], s["npts"], float(s["x0_start"]), float(s["x0_end"]),
                 (s["nu_start"] or "-")[:10], (s["nu_end"] or "-")[:10],
                 s["n_sign_changes"],
                 (s["nu_closest_to_zero"] or {}).get("abs", "-")[:10],
                 s["n_nu_sign_changes"]))
        for ch in s["nu_sign_change_at"]:
            print("      *** nu SIGN CHANGE (SWALLOW-TAIL CANDIDATE) at x0 ~ %s : %s -> %s  D_xxxx=%s"
                  % (ch["x0"], ch["nu_before"], ch["nu_after"], (ch["Dxxxx"] or "?")[:14]))
    json.dump(out, open("analysis.json", "w"), indent=1)
    nsc = sum(s["n_sign_changes"] for s in out)
    nnu = sum(s["n_nu_sign_changes"] for s in out)
    print("\n%d curves. %d bare D_xxx sign changes (these include CENTRE-VARIETY"
          " crossings and are NOT swallow-tails);  %d nu sign changes"
          " (= genuine multiplicity-four candidates)." % (len(out), nsc, nnu))
    return out


if __name__ == "__main__":
    main()
