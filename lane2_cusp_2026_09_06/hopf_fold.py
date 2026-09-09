#!/usr/bin/env python3
"""Task 5 alternative organizer:  T(mu)=0, D(s_f,mu)=0, D_s(s_f,mu)=0.

An ordinary Hopf at the equilibrium A simultaneous with a FOLD of a separate
periodic orbit.  Row 7's continued fold already sits at T = -4.7e-5, so solve
the three equations exactly and then test every nondegeneracy the brief lists:
  det J(A) > 0 ; l1 != 0 ; D_ss(s_f) != 0 ;
  rank[grad_mu T ; d_mu D(s_f)] = 2 ;
  and at that SAME mu: one extra hyperbolic local cycle + one remote cycle.
"""
import json
import mpmath as mp
from engine import Engine
from cusp import Cusp
import foldcont as FC
mp.mp.dps = 40

def Tval(mu): return mu["a11"] + mu["a01"] - 2*mu["a"] - 1
def detJ(mu):  return 2*mu["a"] - mu["a01"] - mu["a10"] - 2*mu["a20"]

START = dict(a=mp.mpf(8)/11, a20=mp.mpf(-12), a11=mp.mpf("2.149952807637873858354431"),
             a01=mp.mpf(67)/220, a10=mp.mpf("-26.5"))
S_F = mp.mpf("3.468994183028599845460731")

eng = Engine(); print("engine:", eng.banner)
mu = dict(START); s = S_F
c = Cusp(eng, mu["a"], mu["a20"], side=1)
# Unknowns (a11, a01, s): T=0 is linear in (a11,a01); D,D_s by Newton.
print("\nsolving  T=0, D(s_f)=0, D_s(s_f)=0  in (a11, a01, s_f)")
for it in range(30):
    r = FC.val(c, mu, s)
    if r["status"] != "OK": print("   return failed:", r["status"]); break
    F = [Tval(mu), r["D"], r["Dx"]]
    m = max(abs(F[0]), abs(F[1]), abs(F[2])*abs(s-1))
    print("   it%-2d |T|=%.3e |D|=%.3e |D_s|=%.3e" % (it, float(abs(F[0])), float(abs(F[1])), float(abs(F[2]))))
    if m < mp.mpf("1e-26"): break
    g = FC.grad_mu(c, mu, s, which=("D", "Dx"))
    if g is None: print("   grad failed"); break
    sc = FC.scales(mu)
    # columns: a11 (idx2), a01 (idx3), s
    J = [[mp.mpf(1)*sc[2], mp.mpf(1)*sc[3], mp.mpf(0)],
         [g["D"][2],       g["D"][3],       r["Dx"]],
         [g["Dx"][2],      g["Dx"][3],      r["Dxx"]]]
    M = mp.matrix(3, 3)
    for i in range(3):
        for j in range(3): M[i, j] = J[i][j]
    try: d = mp.lu_solve(M, mp.matrix(F))
    except Exception as e: print("   singular:", e); break
    mu["a11"] -= d[0]*sc[2]/sc[2]; mu["a01"] -= d[1]*sc[3]/sc[3]; s -= d[2]
    mu["a11"] = mu["a11"]; mu["a01"] = mu["a01"]
print("\nresult:")
for p in FC.PARAMS: print("   %-4s = %s" % (p, mp.nstr(mu[p], 25)))
print("   s_f  = %s" % mp.nstr(s, 20))
r = FC.val(c, mu, s)
print("   T = %s   det J(A) = %s" % (mp.nstr(Tval(mu), 8), mp.nstr(detJ(mu), 10)))
print("   D = %s  D_s = %s  D_ss = %s"
      % (mp.nstr(r["D"],4), mp.nstr(r["Dx"],4), mp.nstr(r["Dxx"],10)))
inv = FC.inventory(c, mu)
print("   local cycles (D=0): %s" % [mp.nstr(v,10) for v in inv["cycles"]])
print("   stationary        : %s" % [mp.nstr(v,10) for v in inv["stationary"]])
rem = FC.remote_cycles(mu)
print("   remote: focus %s  cycles=%s"
      % ("none" if rem["second_focus"] is None else "(%.4f,%.4f) tr %.4f" % rem["second_focus"],
         rem["cycles"]))
json.dump(dict(**{p: mp.nstr(mu[p],25) for p in FC.PARAMS}, s_f=mp.nstr(s,25),
               T=mp.nstr(Tval(mu),20), detJ=mp.nstr(detJ(mu),20),
               Dss=mp.nstr(r["Dxx"],20), local_cycles=len(inv["cycles"]),
               remote_cycles=rem["cycles"]),
          open("ledger_opus/hopf_fold.json", "w"), indent=1)
print("\ncalls:", eng.ncalls); eng.close()
