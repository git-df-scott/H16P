#!/usr/bin/env python3
"""Lane C: independent numerical check of the weak-focus order at Shi's seed
x'=-y-10x^2+5xy+y^2, y'=x+x^2-25xy, via the displacement of the first
return map d(r) = r1 - r for small r at 30 digits. A weak focus of order k
gives d(r) ~ c r^(2k+1); a center gives d(r)=0 to precision."""
import mpmath as mp
mp.mp.dps = 30
def return_disp(r, l=-10, m=5, a=1, b=-25):
    f = lambda t, u: [-u[1]+l*u[0]**2+m*u[0]*u[1]+u[1]**2, u[0]+a*u[0]**2+b*u[0]*u[1]]
    sol = mp.odefun(f, 0, [mp.mpf(r), mp.mpf(0)], tol=mp.mpf('1e-28'))
    # find return to positive x-axis after one revolution: y crosses 0 from below near t≈2pi
    T = mp.findroot(lambda t: sol(t)[1], 2*mp.pi)
    return sol(T)[0]-r, T
vals = []
for r in ("0.02", "0.04", "0.08"):
    d, T = return_disp(mp.mpf(r)); vals.append(d)
    print(f"r={r}: d(r)={mp.nstr(d,10)}  T={mp.nstr(T,8)}")
print("ratios d(2r)/d(r):", mp.nstr(vals[1]/vals[0], 6), mp.nstr(vals[2]/vals[1], 6), " (2^7=128 for order 3, 2^5=32 order 2, 2^9=512 order 4)")
