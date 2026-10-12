"""Discover reflected-sign CLT witnesses across coefficients; no norm upper claims."""
from pathlib import Path
from fractions import Fraction as F
import json, numpy as np
from scipy.special import roots_legendre
from scipy.optimize import minimize
import sys
x,w=roots_legendre(160);z=(x+1)*2.5;w=w*2.5

def moments(a,c):
 q=a*z+c*z**3;dp=a+3*c*z*z;e=np.exp(-(z*z+q*q)/2)
 return np.array([2/np.pi*np.dot(w,e),2/np.pi*np.dot(w,dp*e),2/(np.pi*np.sqrt(6))*np.dot(w,(z*z-1)*dp*e)])
def norm(p,b):
 y,u,v=moments(*p);return u*u-b*(v*v-y*y)
rows=[];start=np.array([2.6412631,.4084047])
for bs in ['.699','.703','.705','.707','.710','.713','.715','.716','.7171']:
 b=float(bs);r=minimize(lambda p:-norm(p,b),start,method='Nelder-Mead',options={'xatol':1e-11,'fatol':1e-14,'maxiter':500})
 assert np.min(r.x)>0
 start=r.x
 rows.append({'b':str(F(bs)),'a':str(F(round(r.x[0]*10**8),10**8)),'c':str(F(round(r.x[1]*10**8),10**8)), 'sampled_norm':norm(r.x,b),'sampled_method_ceiling':1/norm(r.x,b),'status':'DISCOVERY_ONLY'})
p=Path(__file__).with_suffix('.json');p.write_text(json.dumps(rows,indent=2)+'\n')
print(json.dumps(rows,indent=2))
