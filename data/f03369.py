#!/usr/bin/env python3
"""Numerical interior stationary polish; not a global norm proof."""
import json,numpy as np
from pathlib import Path
from scipy.optimize import root
import screen_odd_tail as S
HERE=Path(__file__).resolve().parent
q3=1/np.sqrt(6)
def fun(x,z=S.Z,w=S.W):
 rec=S.witness([x[0],x[1],q3],z,w)
 return [rec['r']+rec['v'],rec['M']-rec['U']]
x=root(fun,[.055,2.5077],tol=1e-12);assert x.success
r=S.witness([*x.x,q3]);b=r['u']/(r['r']*(x.x[1]+3/np.sqrt(6)));t=x.x[0]*b*r['r'];norm=S.scalar(r,b,t)/(1+t)
z,w=S.setup(192);fine=S.witness([*x.x,q3],z,w);fnorm=S.scalar(fine,b,t)/(1+t)
out={'status':'NUMERICAL_STATIONARY_SCREEN_ONLY','b':float(b),'t':float(t),'q0':float(x.x[0]),'q1':float(x.x[1]),'q3':float(q3),'residual':list(map(float,fun(x.x))),'normalized_witness_norm':float(norm),'heldout_normalized_witness_norm':float(fnorm),'quadrature_difference':float(abs(norm-fnorm)),'ceiling':float(1/fnorm),'target_margin':float(1/1.775-fnorm),'witness':r}
(HERE/'stationary_polish.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out),flush=True)
