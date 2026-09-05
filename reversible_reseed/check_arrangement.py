"""Controls for the finite-curve line-arrangement algorithm.

The artificial curves are software controls, not Abelian integrals or fields.
"""
import json
import numpy as np
from moment_search import ROOT,arrangement

t=np.arange(7,dtype=float)
wave=np.column_stack([t,t+1,(-1.)**t])
line=np.column_stack([t,t+1,t-3])
positive_upper=arrangement(wave,line)
positive_lower=arrangement(line,wave)
assert positive_upper['five_candidate'] is not None
assert positive_lower['five_candidate'] is not None
parabola=np.column_stack([t,t+1,(t+1)**2])
negative=arrangement(parabola,parabola)
assert negative['maxima']['upper']<=2
assert negative['five_candidate'] is None
out={'synthetic_upper_five':positive_upper,'synthetic_lower_five':positive_lower,'parabola_negative':negative}
(ROOT/'data'/'arrangement_controls.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
