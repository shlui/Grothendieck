#!/usr/bin/env python3
"""Numerical discovery only: cubic Stein caps for J_a(x)=E max(a|Z|,|x|)."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
from pathlib import Path
import json,time,argparse
import numpy as np
from numpy.polynomial import chebyshev as ch
from scipy.special import erfcx
from scipy.optimize import minimize,Bounds
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
C=np.sqrt(2/np.pi);K=194/np.sqrt(10)
def basis(x,n,L=2):
 t=(x-L)/(x+L);dt=2*L/(x+L)**2;T=ch.chebvander(t,n)
 DT=np.column_stack([ch.chebval(t,ch.chebder(np.eye(n+1)[j]))*dt for j in range(n+1)])
 p=x/(L+x);dp=L/(L+x)**2
 return np.column_stack([x,x*x/2,x**3/3,x**4/4,p[:,None]*T]),np.column_stack([np.ones_like(x),x,x*x,x**3,dp[:,None]*T+p[:,None]*DT])
def fprime(x,a):
 s=np.sqrt(1+a*a);b=s/a
 return s*C*(x*np.sqrt(np.pi/2)*erfcx(x/np.sqrt(2))-1)-s*x*np.exp(-x*x/(2*a*a))*erfcx(b*x/np.sqrt(2))+C*a*np.exp(-x*x/(2*a*a))
def solve(a,k,old=None,n=12):
 start=time.time();x=np.unique(np.r_[np.linspace(0,3,1001),np.linspace(3,12,601),np.geomspace(12,1e3,251)])
 xx=np.unique(np.r_[np.linspace(0,64,64001),np.geomspace(64,1e6,4001)])
 m6=15+15*k+np.array([K,K-633/15])*k**1.5;hv=k*np.array([7/90,1/9]);hc=C*np.sqrt(1+a*a)
 if old is None:
  old=json.loads((ROOT/'reduced_search/energy/stein_mean_discovery.json').read_text())['rows'][2]
  current=np.array(old['phi_coefficients'])*.05;q=np.array(old['majorant'])*.05;eta=old['eta']*.05
 else:current=np.array(old['phi_coefficients']);q=np.array(old['majorant']);eta=old['eta']
 hist=[]
 for it in range(5):
  w=1/(1+x**6);sw=np.sqrt(w);V,D=basis(x,n);A=D-x[:,None]*V
  _,sing,Vt=np.linalg.svd(np.vstack([A*w[:,None],D*sw[:,None]]),full_matrices=False);keep=sing>sing[0]*1e-8;Tr=Vt[keep].T/sing[keep];rank=Tr.shape[1]
  Av=w[:,None]*(A@Tr);Dv=sw[:,None]*(D@Tr);ss=sw*fprime(x,a);W=np.column_stack([w,w*x*x,w*x**4,w*x**6]);mom=np.column_stack([np.ones(2),np.ones(2),np.full(2,3+k),m6])
  def con(z):
   v=z[:rank];qq=z[rank:rank+4];e=z[-2];obj=z[-1];r=ss-Dv@v;tt=-(Tr@v)[3]
   return np.r_[W@qq-Av@v-r*r/(2*e),obj-hc-mom@qq-hv*e,qq[3]-tt*tt/(2*e)]
  def jac(z):
   v=z[:rank];e=z[-2];r=ss-Dv@v;tt=-(Tr@v)[3];J=np.zeros((len(x)+3,rank+6));dR=r/e
   J[:len(x),:rank]=-Av+dR[:,None]*Dv;J[:len(x),rank:rank+4]=W;J[:len(x),-2]=dR*dR/2
   J[len(x):-1,rank:rank+4]=-mom;J[len(x):-1,-2]=-hv;J[len(x):-1,-1]=1
   J[-1,:rank]=tt/e*Tr[3];J[-1,rank+3]=1;J[-1,-2]=tt*tt/(2*e*e)
   return J
  v0=np.linalg.lstsq(Tr,current,rcond=None)[0];z0=np.r_[v0,q,eta,hc+np.max(mom@q+hv*eta)+1e-5];g=np.zeros(rank+6);g[-1]=1
  lo=np.full(rank+6,-np.inf);hi=-lo;lo[-2]=1e-9;hi[-2]=2;lo[rank+3]=0;lo[:rank]=v0-10;hi[:rank]=v0+10
  res=minimize(lambda z:z[-1],z0,jac=lambda z:g,method='SLSQP',bounds=Bounds(lo,hi),constraints=[{'type':'ineq','fun':con,'jac':jac}],options={'maxiter':350,'ftol':5e-13})
  current=Tr@res.x[:rank];q=res.x[rank:rank+4];eta=res.x[-2];VV,DD=basis(xx,n);rr=fprime(xx,a)-DD@current
  err=((DD-xx[:,None]*VV)@current+rr*rr/(2*eta)-np.polynomial.polynomial.polyval(xx*xx,q))/(1+xx**6)
  inds=np.argsort(err)[-80:];vi=float(max(err));hist.append(dict(iteration=it,status=res.message,iterations=res.nit,rank=rank,cap=float(res.x[-1]),constraint_min=float(min(con(res.x))),violation=vi,seconds=time.time()-start));print(a,k,hist[-1],flush=True)
  if vi<1e-10:break
  x=np.unique(np.r_[x,xx[inds]])
 return dict(status='NUMERICAL_DISCOVERY_ONLY',degree=n,scale=2.,kappa=k,a=a,eta=float(eta),phi_coefficients=current.tolist(),majorant=q.tolist(),J_upper_sampled=float(res.x[-1]),J_Gaussian=hc,max_scaled_violation_dense=vi,history=hist)
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--a',type=float,nargs='+',default=[2.,2.4,2.8,3.2]);ap.add_argument('--kappa',type=float,default=.53);ap.add_argument('--out',default='Ja_source_discovery.json');args=ap.parse_args()
 rows=[];old=None
 for a in args.a:
  old=solve(a,args.kappa,old);rows.append(old);(HERE/args.out).write_text(json.dumps({'status':'NUMERICAL_DISCOVERY_ONLY','rows':rows},indent=2))
