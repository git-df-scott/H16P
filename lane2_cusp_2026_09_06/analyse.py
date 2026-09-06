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
    hdr = "%-28s %5s %-11s %-11s %-13s %-13s %3s %-11s" % (
        "curve", "npts", "x0 start", "x0 end", "Dxxx start", "Dxxx end", "sc", "min|Dxxx|")
    print(hdr); print("-" * len(hdr))
    for s in out:
        print("%-28s %5d %-11.5g %-11.5g %-13.4g %-13.4g %3d %-11.4g"
              % (s["file"][:28], s["npts"], float(s["x0_start"]), float(s["x0_end"]),
                 float(s["Dxxx_start"]), float(s["Dxxx_end"]), s["n_sign_changes"],
                 float(s["Dxxx_closest_to_zero"]["abs"])))
        for ch in s["sign_change_at"]:
            print("      *** D_xxx SIGN CHANGE at x0 ~ %s : %s -> %s"
                  % (ch["x0"], ch["Dxxx_before"], ch["Dxxx_after"]))
    json.dump(out, open("analysis.json", "w"), indent=1)
    nsc = sum(s["n_sign_changes"] for s in out)
    print("\n%d curves, %d D_xxx sign changes in total" % (len(out), nsc))
    return out


if __name__ == "__main__":
    main()
