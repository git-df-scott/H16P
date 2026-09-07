"""Summarise the line-scan and cusp ledgers."""
import glob, json, os, sys
import numpy as np
from collections import Counter, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))


def load(pat):
    for p in sorted(glob.glob(os.path.join(HERE, pat))):
        with open(p) as f:
            for line in f:
                line = line.strip()
                if line:
                    try:
                        yield json.loads(line)
                    except json.JSONDecodeError:
                        pass


def lines():
    rows = [r for r in load("ledger/line_*.jsonl") if r.get("kind") == "line"]
    if not rows:
        print("no line-scan rows yet")
        return
    by = defaultdict(list)
    for r in rows:
        by[r["seed"]].append(r)
    nev = sum(len(r["counts"]) for r in rows)
    n3 = sum(1 for r in rows if r["max_extrema"] >= 3)
    print(f"LINE SCANS: {len(rows)} directions, {nev} field evaluations, "
          f"{n3} directions reaching 3 interior extrema")
    print()
    print(f"{'seed':10s} {'dirs':>5s} {'evals':>7s} {'max ext':>8s} "
          f"{'survival + (median)':>20s} {'survival - (median)':>20s} {'widest window':>14s}")
    for s, rs in sorted(by.items()):
        pos = np.array([r["lam_break_pos"] for r in rs])
        neg = np.array([abs(r["lam_break_neg"]) for r in rs])
        width = pos + neg
        print(f"{s:10s} {len(rs):5d} {sum(len(r['counts']) for r in rs):7d} "
              f"{max(r['max_extrema'] for r in rs):8d} "
              f"{np.median(pos):20.4e} {np.median(neg):20.4e} {width.max():14.4e}")
    print()
    allw = np.array([r["lam_break_pos"] + abs(r["lam_break_neg"]) for r in rows])
    print(f"survival width of the three-cycle configuration along a random line,")
    print(f"as a relative displacement of the coefficient vector:")
    for q in (10, 25, 50, 75, 90, 100):
        print(f"   {q:3d}th percentile  {np.percentile(allw, q):.4e}")
    cnt = Counter()
    for r in rows:
        for c in r["counts"]:
            cnt[c[1]] += 1
    print()
    print("interior extremum count over every field evaluated:", dict(sorted(cnt.items())))


def cusps(pat="ledger/cusp_c1.jsonl"):
    rows = [r for r in load(pat) if r.get("kind") == "cusp_try"]
    if not rows:
        print("no cusp rows yet")
        return
    print()
    print(f"CUSP CONTINUATION ({pat}): {len(rows)} Newton runs, "
          f"{sum(1 for r in rows if r['ok'])} converged")
    agg = defaultdict(list); why = defaultdict(Counter)
    for r in rows:
        k = (r["seed"], r.get("region"))
        why[k][r.get("why") or "CONVERGED"] += 1
        if r["ok"]:
            agg[k].append((abs(r["lam"]), r.get("max_extrema_near_cusp", -1)))
    print()
    print(f"{'seed':10s} {'region':9s} {'runs':>5s} {'conv':>5s} "
          f"{'med |lam*|':>12s} {'min |lam*|':>12s} {'max ext near':>13s}")
    for k in sorted(why):
        v = agg.get(k, [])
        lam = np.array([x[0] for x in v]) if v else np.array([])
        m = max([x[1] for x in v] or [-1])
        print(f"{k[0]:10s} {str(k[1]):9s} {sum(why[k].values()):5d} {len(v):5d} "
              f"{(f'{np.median(lam):.3e}' if lam.size else '-'):>12s} "
              f"{(f'{lam.min():.3e}' if lam.size else '-'):>12s} {m:13d}")
    print()
    for reg in ("inner", "outer", "between", "whole"):
        c = Counter()
        for k in why:
            if k[1] == reg:
                c.update(why[k])
        if c:
            print(f"  {reg:8s} {dict(c)}")


if __name__ == "__main__":
    lines()
    cusps()
