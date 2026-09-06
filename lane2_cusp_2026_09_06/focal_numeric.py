#!/usr/bin/env python3
"""Focal values at the Newton limit, computed at 40 digits by the formal first
integral, and compared with a deliberately perturbed nearby point."""
import mpmath as mp
mp.mp.dps = 40

def focal(a, a20, a11, a10, nmax=8):
    a01 = -a11 - 3 if a == -2 else None
    a01 = -a11 - 2*a - 1 + 0           # V1 = 0  =>  a01 = 2a+1-a11
    a01 = 2*a + 1 - a11
    a00 = a01 + a11 - a10 - a20 - a
    # translate to A=(1,-1):  x=u+1, y=v-1
    # P = 1+xy = uv - u + v ;  Q = ... (expanded below)
    # quadratic system in (u,v)
    Pc = {(1,1): mp.mpf(1), (1,0): mp.mpf(-1), (0,1): mp.mpf(1)}
    Qc = {(1,0): a10 + 2*a20 - a11, (0,1): a01 + a11 - 2*a,
          (2,0): a20, (1,1): a11, (0,2): a}
    A11 = Pc.get((1,0), 0); A12 = Pc.get((0,1), 0)
    A21 = Qc.get((1,0), 0); A22 = Qc.get((0,1), 0)
    det = A11*A22 - A12*A21; tr = A11 + A22
    if abs(tr) > mp.mpf(10)**(-mp.mp.dps+8): raise ValueError("V1 != 0: %s" % mp.nstr(tr,6))
    w = mp.sqrt(det)
    # u = U, v = (-A11 U + w V)/A12
    def ev(c, U, Vv):
        uu = U; vv = (-A11*U + w*Vv)/A12
        return sum(co*uu**i*vv**j for (i, j), co in c.items())
    # build polynomial coefficients of Ud, Vd in (U,V) up to degree 2 by sampling
    import itertools
    mons = [(2,0),(1,1),(0,2),(1,0),(0,1),(0,0)]
    def poly_of(fun):
        pts = [(mp.mpf(i+1)/7, mp.mpf(j+2)/5) for i in range(3) for j in range(3)][:6]
        A = mp.matrix(6,6); b = mp.matrix(6,1)
        for r,(Uu,Vv) in enumerate(pts):
            for cidx,(i,j) in enumerate(mons): A[r,cidx] = Uu**i*Vv**j
            b[r] = fun(Uu,Vv)
        s = mp.lu_solve(A,b)
        return {mons[k]: s[k] for k in range(6)}
    Ud = poly_of(lambda Uu,Vv: ev(Pc,Uu,Vv)/w)
    Vd = poly_of(lambda Uu,Vv: (A12*ev(Qc,Uu,Vv) + A11*ev(Pc,Uu,Vv))/(w*w))
    # strip the linear part (should be V and -U)
    f2 = {k:v for k,v in Ud.items() if sum(k)==2}
    g2 = {k:v for k,v in Vd.items() if sum(k)==2}
    # formal integral
    Fc = {2: {(2,0): mp.mpf(1)/2, (0,2): mp.mpf(1)/2}}
    out = {}
    for k in range(3, nmax+1):
        src = {}
        prev = Fc[k-1]
        for (i,j),c in prev.items():
            if i>0:
                for (p,q),d in f2.items():
                    key=(i-1+p, j+q); src[key]=src.get(key,0) - c*i*d
            if j>0:
                for (p,q),d in g2.items():
                    key=(i+p, j-1+q); src[key]=src.get(key,0) - c*j*d
        n = k+1
        A = mp.matrix(n+1, n+1 if k%2 else n+2) if False else None
        cols = n + (1 if k%2==0 else 0)
        M = mp.matrix(n+1, cols); rhs = mp.matrix(n+1,1)
        # L(U^p V^q) = p U^(p-1) V^(q+1) - q U^(p+1) V^(q-1)
        basis = [(k-i, i) for i in range(k+1)]
        rows = [(k-i, i) for i in range(k+1)]
        ridx = {m:i for i,m in enumerate(rows)}
        for ci,(p,q) in enumerate(basis):
            if p>0: M[ridx[(p-1,q+1)], ci] += p
            if q>0: M[ridx[(p+1,q-1)], ci] -= q
        if k%2==0:
            # column for -v*(U^2+V^2)^(k/2)
            m = k//2
            from math import comb
            for t in range(m+1):
                M[ridx[(2*t, k-2*t)], n] = -mp.binomial(m, t)
        for m0, c in src.items():
            if m0 in ridx: rhs[ridx[m0]] += c
        sol = mp.lu_solve(M, rhs) if M.rows == M.cols else mp.qr_solve(M, rhs)[0]
        Fc[k] = {basis[i]: sol[i] for i in range(n)}
        if k%2==0: out[k] = sol[n]
    return out

pts = [("Newton limit", mp.mpf(-2), -mp.mpf(3254)/675,
        mp.mpf("8.067217612039481856440047"), mp.mpf("16.53363941409600843877444")),
       ("a11 shifted 1e-6", mp.mpf(-2), -mp.mpf(3254)/675,
        mp.mpf("8.067218612039481856440047"), mp.mpf("16.53363941409600843877444")),
       ("a10 shifted 1e-6", mp.mpf(-2), -mp.mpf(3254)/675,
        mp.mpf("8.067217612039481856440047"), mp.mpf("16.53363841409600843877444"))]
print(" point                V3 (=V4 coeff)        V5 (=V6 coeff)        V7 (=V8 coeff)")
for nm, a, a20, a11, a10 in pts:
    try:
        o = focal(a, a20, a11, a10, nmax=8)
        print("  %-19s %-21s %-21s %s"
              % (nm, mp.nstr(o.get(4,mp.mpf('nan')),8), mp.nstr(o.get(6,mp.mpf('nan')),8),
                 mp.nstr(o.get(8,mp.mpf('nan')),8)))
    except Exception as e:
        print("  %-19s ERROR %s" % (nm, e))
