#!/usr/bin/env python3
"""Reproduce Fable's reported neutral loop at a=-2.5 (b=1.28125, l=-1.13719,
sigma=-1.3e-9, eta_2=+2.0e-8) as a cross-check of this machinery."""
import numpy as np
from d2core import *
a, b, l = -2.5, 1.28125, -1.13719
print("a=%g b=%g l=%g   m=%.10g" % (a, b, l, mval(a, b, l)))
print("eta_2 = %+.6e   C3 = %+.6e" % (eta2(a, b, l), C3(a, b, l)))
for (x, y, sg) in saddles(a, b, l):
    sp = splitting(a, b, l, (x, y))
    print("  saddle (%+.8f, %+.8f)  sigma=%+.6e  splitting=%s"
          % (x, y, sg, "None" if sp is None else "%+.3e" % sp))
