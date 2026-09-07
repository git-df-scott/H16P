"""Summarise the Lane 1 ledgers."""
import glob, json, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))


def load(pattern="ledger/*.jsonl"):
    rows = []
    for p in sorted(glob.glob(os.path.join(HERE, pattern))):
        with open(p) as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                try:
                    r = json.loads(line)
                except json.JSONDecodeError:
                    continue
                r["_file"] = os.path.basename(p)
                rows.append(r)
    return rows


def summary(rows):
    out = {}
    by_seed = {}
    for r in rows:
        if r.get("kind") not in ("seed", "perturb", "climb", "climb_start", "climb_end"):
            continue
        s = r.get("seed")
        by_seed.setdefault(s, []).append(r)
    for s, rs in sorted(by_seed.items()):
        ok = [r for r in rs if r.get("status") == "ok"]
        ne = np.array([r.get("n_extrema", -1) for r in ok]) if ok else np.array([])
        fm = np.array([r.get("score", np.inf) if r.get("score") is not None
                       else np.inf for r in ok]) if ok else np.array([])
        two = fm[ne == 2] if ne.size else np.array([])
        best_i = int(np.argmin(two)) if two.size else None
        out[s] = dict(
            n_rows=len(rs), n_ok=len(ok),
            max_extrema=int(ne.max()) if ne.size else None,
            n_with_3plus=int((ne >= 3).sum()) if ne.size else 0,
            best_fold_margin=float(two.min()) if two.size else None,
            median_fold_margin=float(np.median(two)) if two.size else None,
            seed_fold_margin=next((r.get("score") for r in rs
                                   if r.get("kind") in ("seed", "climb_start")), None),
        )
    return out


def best_records(rows, k=10):
    cand = [r for r in rows
            if r.get("status") == "ok" and r.get("n_extrema") == 2
            and r.get("score") is not None and np.isfinite(r["score"])]
    cand.sort(key=lambda r: r["score"])
    return cand[:k]


def three_plus(rows):
    return [r for r in rows if r.get("status") == "ok" and r.get("n_extrema", 0) >= 3]


if __name__ == "__main__":
    pat = sys.argv[1] if len(sys.argv) > 1 else "ledger/*.jsonl"
    rows = load(pat)
    print(f"{len(rows)} ledger rows from {pat}")
    tp = three_plus(rows)
    print(f"rows with >=3 interior extrema of beta*: {len(tp)}")
    for r in tp[:20]:
        print("  ", r.get("seed"), r.get("kind"), r.get("n_extrema"),
              [round(e["prominence"], 5) for e in (r.get("extrema") or [])])
    print()
    print(f"{'seed':10s} {'rows':>6s} {'ok':>6s} {'maxext':>7s} {'#>=3':>5s} "
          f"{'seed score':>11s} {'best score':>11s} {'median':>10s}")
    for s, v in summary(rows).items():
        sf, bf, mf = (v["seed_fold_margin"], v["best_fold_margin"],
                      v["median_fold_margin"])
        fmt = lambda x, d=6: ("-" if x is None else f"{x:.{d}f}")
        print(f"{s:10s} {v['n_rows']:6d} {v['n_ok']:6d} {str(v['max_extrema']):>7s} "
              f"{v['n_with_3plus']:5d} {fmt(sf):>11s} {fmt(bf):>11s} {fmt(mf, 4):>10s}")
    print()
    print("best 10 fields by fold margin:")
    for r in best_records(rows):
        print(f"  {r.get('seed'):10s} {r.get('kind'):12s} eps={r.get('eps')} "
              f"score={r['score']:.6f} fold={r.get('fold_margin')} ({r.get('fold_margin_kind')}) "
              f"rng={r.get('height_range'):.3e} run={r.get('n_run')}")
