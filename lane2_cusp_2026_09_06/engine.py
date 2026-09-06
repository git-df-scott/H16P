"""Python driver for the Lane-2 cusp engine (binary128 / long double Taylor-jet).

All parameters and returned values are carried as mpmath mpf at dps=50.  This
matters: the cusp Newton takes steps of relative size ~1e-25 in (a11,a01,a10),
which a Python float cannot represent, and with float parameters the Newton
stagnates at a residual of ~1e-17 for purely representational reasons.
"""
import subprocess, os, json, time
import mpmath as mp

mp.mp.dps = 50

HERE = os.path.dirname(os.path.abspath(__file__))


def s(v, digits=42):
    """Decimal string at `digits` significant figures, safe for strtoflt128."""
    if isinstance(v, str):
        return v
    return mp.nstr(mp.mpf(v), digits, strip_zeros=False)


class Engine:
    """Persistent process wrapper.  D() returns mpf values."""

    def __init__(self, quad=True, exe=None, env=None):
        self.quad = quad
        self.exe = exe or os.path.join(HERE, "cusp_engine")
        args = [self.exe] + ([] if quad else ["--ld"])
        e = dict(os.environ)
        if env:
            e.update(env)
        self.p = subprocess.Popen(args, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                  text=True, bufsize=1, env=e)
        self.p.stdin.write("V\n"); self.p.stdin.flush()
        self.banner = self.p.stdout.readline().strip()
        self.ncalls = 0

    def D(self, a, a20, a11, a01, a10, x0, side=1):
        line = ("D" if side > 0 else "L") + " %s %s %s %s %s %s\n" % (s(a), s(a20), s(a11), s(a01), s(a10), s(x0))
        self.p.stdin.write(line); self.p.stdin.flush()
        out = self.p.stdout.readline().strip()
        self.ncalls += 1
        f = out.split()
        if not f:
            raise RuntimeError("engine died")
        if f[0] != "OK":
            return {"status": f[0]}
        return {"status": "OK",
                "D": mp.mpf(f[1]), "Dx": mp.mpf(f[2]), "Dxx": mp.mpf(f[3]),
                "Dxxx": mp.mpf(f[4]), "T": mp.mpf(f[5]), "transv": mp.mpf(f[6]),
                "nsteps": int(f[7])}

    def close(self):
        try:
            self.p.stdin.write("Q\n"); self.p.stdin.flush(); self.p.wait(timeout=5)
        except Exception:
            self.p.kill()


# --------------------------------------------------------------------------
# Cherkas third-order weak focus family (SEEDS.json / LIT_A eq.(16))
# --------------------------------------------------------------------------
def third_order(a, a20):
    """(a11, a01, a10) making A=(1,-1) a weak focus of order exactly 3."""
    a, a20 = mp.mpf(a), mp.mpf(a20)
    a11 = 4 - 2 * a
    a01 = 2 * a + 1 - a11
    a10 = (6 * (a * a - a - 2) + a20 * (6 * a - 7)) / (1 - 3 * a)
    return [a11, a01, a10]


def V7_of(a, a20):
    """V7 restricted to the third-order stratum (derived from LIT_A eq.(15)):
       V7 = -150 (a-2) [ -4 a (a+1)(a-2)^2 + a20 (a-1)(2a+1)^2 ]."""
    a, a20 = mp.mpf(a), mp.mpf(a20)
    return -150 * (a - 2) * (-4 * a * (a + 1) * (a - 2) ** 2 + a20 * (a - 1) * (2 * a + 1) ** 2)


def L_of(a, a20, a11, a01, a10):
    """det J at A = (1,-1);  L > 0 <=> antisaddle."""
    return 2 * mp.mpf(a) - mp.mpf(a01) - mp.mpf(a10) - 2 * mp.mpf(a20)


def V1_of(a, a11, a01):
    """trace J at A = (1,-1)."""
    return mp.mpf(a11) + mp.mpf(a01) - 2 * mp.mpf(a) - 1


def a00_of(a, a20, a11, a01, a10):
    return mp.mpf(a01) + mp.mpf(a11) - mp.mpf(a10) - mp.mpf(a20) - mp.mpf(a)
