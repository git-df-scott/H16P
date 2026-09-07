"""Append-only JSONL ledgers (PROTOCOL rule 5)."""
import json, os, time, threading
from fractions import Fraction

_HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(_HERE, "data")
os.makedirs(DATA, exist_ok=True)
_lock = threading.Lock()

def _plain(o):
    if isinstance(o, Fraction):
        return str(o)
    if isinstance(o, (list, tuple)):
        return [_plain(v) for v in o]
    if isinstance(o, dict):
        return {k: _plain(v) for k, v in o.items()}
    try:
        import numpy as np
        if isinstance(o, np.generic):
            return o.item()
        if isinstance(o, np.ndarray):
            return [_plain(v) for v in o.tolist()]
    except Exception:
        pass
    if isinstance(o, float):
        if o != o or o in (float("inf"), float("-inf")):
            return str(o)
    return o

def append(name, rec):
    rec = dict(rec)
    rec.setdefault("ts", time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()))
    path = os.path.join(DATA, name + ".jsonl")
    line = json.dumps(_plain(rec), allow_nan=False, sort_keys=True)
    with _lock:
        with open(path, "a") as f:
            f.write(line + "\n")
    return path

def sizes():
    out = {}
    for fn in sorted(os.listdir(DATA)):
        if fn.endswith(".jsonl"):
            p = os.path.join(DATA, fn)
            with open(p) as f:
                n = sum(1 for _ in f)
            out[fn] = (n, os.path.getsize(p))
    return out
