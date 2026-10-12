"""Discovery for genuine all-law sufficient criterion; local numerical maxima NOT certificates."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import numpy as np,json,argparse
from scipy.special import ndtr
from scipy.optimize import minimize
from pathlib import Path
EP=.72;AL=1.0;DE=0.0;SQ2=np.sqrt(2);SQ6=np.sqrt(6)


def realroots(c):
 c=np.array(c,dtype=float)
 while len(c)>1 and abs(c[-1])<1e-14:c=c[:-1]
 if len(c)<2:return []
 rr=np.polynomial.polynomial.polyroots(c)
 return rr.real[abs(rr.imag)<1e-7].tolist()
def pol(c):return np.array([c[0]-c[2]/SQ2,c[1]-3*c[3]/SQ6,c[2]/SQ2,c[3]/SQ6])
def moments(ell,beta,v):
 A=pol([0,ell,0,beta]);B=pol(v)
 rr=[]
 for sa in [-1,1]:
  for sb in [-1,1]:
   co=sa*A+sb*B;co[0]-=DE;rr+=realroots(co)
 rr+=realroots(A)+realroots(B)
 rr=np.sort(np.array([r for r in rr if abs(r)<39]));rr=np.r_[-np.inf,rr,np.inf]
 lo,hi=rr[:-1],rr[1:];mid=(lo+hi)/2
 mid[np.isneginf(lo)]=np.minimum(-40,hi[np.isneginf(lo)]-1)
 mid[np.isposinf(hi)]=np.maximum(40,lo[np.isposinf(hi)]+1)
 av=np.polynomial.polynomial.polyval(mid,A);bv=np.polynomial.polynomial.polyval(mid,B)
 ish=abs(av)-DE>=abs(bv); hs=np.sign(av)*ish;ks=np.sign(bv)*(~ish)
 # integrals H0,H1,H2,H3 over intervals
 def edge(t):
  z=np.where(np.isfinite(t),t,0);ph=np.exp(-z*z/2)/np.sqrt(2*np.pi);ph[~np.isfinite(t)]=0
  return np.array([ph,z*ph/SQ2,(z*z-1)*ph/SQ6])
 ww=np.vstack((ndtr(hi)-ndtr(lo),edge(lo)-edge(hi)))
 hm=ww@hs;km=ww@ks;mass=ww[0]@ish
 return np.r_[hm[1],hm[3],km,mass]