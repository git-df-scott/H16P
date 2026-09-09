#!/usr/bin/env python3
"""Counting rules for the raw displacement D on a section from the focus,
with D(1) = 0 (the focus itself is a root).

LEMMA 1.  If D has three distinct cycle roots 1 < c1 < c2 < c3, then D has at
least THREE stationary points: one in (1,c1), one in (c1,c2), one in (c2,c3).
  Proof: Rolle on each of the four consecutive roots 1 < c1 < c2 < c3.
  => three stationary points is what three cycles ALREADY force.  It is NOT an
     entry certificate for anything further.

LEMMA 2.  If D has two separated double roots 1 < s1 < s2 (each a double cycle),
then D has at least FOUR distinct stationary points: s1 itself, s2 itself, one
in (1,s1) and one in (s1,s2).
  Proof: s1,s2 are stationary as double roots.  Rolle on the roots 1 < s1 gives
  a stationary point in (1,s1); Rolle on s1 < s2 gives one in (s1,s2).  These
  four are distinct because the two interior ones lie in disjoint open intervals
  neither containing s1 nor s2.

CONSEQUENCE.  Two separated finite-amplitude double cycles need >= 4 stationary
points of raw D.  Rows 2, 7, 8 have exactly 3 -- the Rolle minimum for their
three cycles -- so the earlier 'encouraging census' claim was wrong.

QUALIFICATION of the opposite-type rejection.  Adjacent max/min stationary
points on a nonconstant MONOTONE segment cannot share a value.  But opposite-type
double zeros coexist freely when an intervening root is present."""
import sympy as S

r = S.symbols('r')
d = r*(r-1)**2*(r-2)*(r-3)**2
print("REGRESSION for the rejection logic:  d(r) = r(r-1)^2(r-2)(r-3)^2")
print("   roots:", S.roots(S.Poly(d, r)))
dp = S.expand(S.diff(d, r))
st = sorted([x for x in S.real_roots(S.Poly(dp, r))], key=lambda z: float(z))
print("   stationary points of d:", [S.N(x, 8) for x in st])
for x in st:
    print("      s=%-12s d=%-14s d''=%s"
          % (S.N(x, 8), S.N(d.subs(r, x), 6), S.N(S.diff(d, r, 2).subs(r, x), 6)))
print("""
   d has double zeros at r=1 and r=3 of OPPOSITE type (max at 1, d''=-8; min at 3, d''=+24),
   with a simple zero at r=2 between them, and d is not identically zero.
   => 'two opposite-type folds cannot both sit at zero' is FALSE in general.
      It holds only across a monotone segment with no intervening root.
   The CENTRE_AND_DOUBLEFOLD_REPORT claim is corrected accordingly.""")
n_double = 2
print("   Lemma 2 check: two double roots (1,3) + simple root 2 =>",
      len(st), "stationary points found; Lemma 2 requires >= 4 when both double")
