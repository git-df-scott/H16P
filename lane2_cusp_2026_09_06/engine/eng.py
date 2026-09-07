"""Python driver for the binary128 Lane-2 engine (persistent subprocess).

All numbers cross the pipe as decimal strings with >= 36 significant digits, so
no precision is lost either way.  Values come back as mpmath mpf (dps 50)."""
import subprocess, os, math, hashlib
from mpmath import mp, mpf, nstr, sqrt as msqrt, atan2 as matan2

mp.dps = 50
HERE = os.path.dirname(os.path.abspath(__file__))
BIN = os.path.join(HERE, "cusp128")

def fmt(v):
    if isinstance(v, str): return v
    return nstr(mpf(v), 40, strip_zeros=False)

class Engine:
    def __init__(self, tol="1e-28"):
        self.tol = tol
        self.p = subprocess.Popen([BIN, "--tol", tol], stdin=subprocess.PIPE,
                                  stdout=subprocess.PIPE, text=True, bufsize=1)
        self.calls = 0
    def _cmd(self, line):
        self.p.stdin.write(line + "\n"); self.p.stdin.flush()
        return self.p.stdout.readline().strip()
    def info(self, mu):
        return self._cmd("INFO " + " ".join(fmt(x) for x in mu))
    def D(self, mu, rho, side=0):
        """mu=(a,a20,a11,a01,a10); side 0 = ray {y=-1,x>1}, 1 = {y=-1,x<1}."""
        self.calls += 1
        r = self._cmd("D " + " ".join(fmt(x) for x in mu) + " " + fmt(rho)
                      + " " + str(int(side)))
        f = r.split()
        if f[0] != "OK":
            return {"ok": False, "why": f[1] if len(f) > 1 else "?",
                    "th": float(f[2]) if len(f) > 2 else 0.0}
        return {"ok": True, "D": mpf(f[1]), "D1": mpf(f[2]), "D2": mpf(f[3]),
                "D3": mpf(f[4]), "min_den": mpf(f[5]), "min_rho": mpf(f[6]),
                "x": mpf(f[7]), "w": mpf(f[8]), "steps": int(f[9]),
                "raw": f[1:5]}
    def close(self):
        try:
            self.p.stdin.write("QUIT\n"); self.p.stdin.flush(); self.p.wait(5)
        except Exception:
            self.p.kill()

def geom(mu):
    a, a20, a11, a01, a10 = [mpf(v) if not isinstance(v, str) else mpf(v) for v in mu]
    A = a10 + 2*a20 - a11
    B = a01 + a11 - 2*a
    T = B - 1
    L = -B - A
    disc = L - T*T/4
    if L <= 0 or disc <= 0: return None
    w = msqrt(disc); k1 = 1 + T/2
    nrm = msqrt(w*w + k1*k1)
    return dict(T=T, L=L, w=w, th0=matan2(k1, w), nrm=nrm, scale=w/nrm)

def x_of_rho(mu, rho, side=0):
    g = geom(mu)
    return 1 + mpf(rho)*g["scale"] if side == 0 else 1 - mpf(rho)*g["scale"]
def rho_of_x(mu, x, side=0):
    g = geom(mu)
    return (mpf(x) - 1)/g["scale"] if side == 0 else (1 - mpf(x))/g["scale"]

def engine_hash():
    h = hashlib.sha256()
    with open(os.path.join(HERE, "cusp128.cpp"), "rb") as fh: h.update(fh.read())
    return h.hexdigest()[:16]
