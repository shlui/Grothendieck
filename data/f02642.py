#!/usr/bin/env python3
"""Numerical discovery for source retaining longitudinal h3. Not a certificate."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
from pathlib import Path
import json,argparse,numpy as np
from scipy.special import roots_hermitenorm,roots_laguerre,erf,erfcx
import discover_Ja_source as M
HERE=Path(__file__).resolve().parent
class Target:
 def __init__(self,a,d,laguerre=96,hermite=256):
  self.C=np.sqrt(2/np.pi);self.a=a;self.d=d;self.B=d/np.sqrt(6);self.A=a-3*self.B
  if self.A<=0 or d<0:raise ValueError('Require monotone odd cubic, positive d')
  self.nodes,self.weights=roots_laguerre(laguerre)
  z,w=roots_hermitenorm(hermite);y=np.abs(self.A*z+self.B*z**3)
  self.mean=float(np.dot(w,(self.C*np.exp(-y*y/2)+y*erf(y/np.sqrt(2))))/np.sqrt(2*np.pi))
 def inverse(self,y):
  if self.B==0:return y/self.A
  return 2*np.sqrt(self.A/(3*self.B))*np.sinh(np.arcsinh((3*y/(2*self.A))*np.sqrt(3*self.B/self.A))/3)
 def deriv(self,x):
  x=np.asarray(x,dtype=float);z=self.inverse(x);EYtail=self.C*np.exp(-z*z/2)*(self.a+self.B*(z*z-1))
  G=1-x*np.sqrt(np.pi/2)*erfcx(x/np.sqrt(2))
  ss=self.nodes[None,:];xx=x[:,None];t=ss/(xx+1);y=xx+t;zz=self.inverse(y)
  Gy=1-y*np.sqrt(np.pi/2)*erfcx(y/np.sqrt(2))
  integ=Gy*np.exp(ss/(xx+1)-t*t/2-zz*zz/2)/(self.A+3*self.B*zz*zz)
  term=self.C*x/(x+1)*(integ@self.weights)
  return -G*(self.mean-EYtail)-term
 def value(self,x):
  # Not needed by Stein fitting; Gaussian quadrature used only for diagnostics.
  z,w=roots_hermitenorm(256);y=np.abs(self.A*z+self.B*z**3)
  return np.maximum(np.asarray(x)[:,None],y[None,:])@w/np.sqrt(2*np.pi)
def solve(a,d,k,n=16,old=None):
 target=Target(a,d)
 orig_f=M.fprime;orig_C=M.C
 # Reuse discovery SLSQP exactly; its Gaussian mean gets replaced below.
 def deriv(x,aa):return target.deriv(x)
 M.fprime=deriv;M.C=target.mean/np.sqrt(1+a*a)
 try:r=M.solve(a,k,old,n)
 finally:M.fprime=orig_f;M.C=orig_C
 r['status']='NUMERICAL_CUBIC_HEAD_SOURCE_ONLY';r['a']=a;r['d']=d;r['Gaussian_mean']=target.mean;r['source_formula']='a*u+d*v+r <= A_ad(kappa,V)'
 # Recompute mean with twice quadrature order as numerical convergence evidence.
 high=Target(a,d,192,512);r['mean_quadrature_difference']=high.mean-target.mean
 grid=np.r_[np.linspace(0,8,401),np.linspace(8,30,100),np.geomspace(30,1e3,30)]
 r['derivative_quadrature_difference']=float(max(abs(high.deriv(grid)-target.deriv(grid))))
 q=r['majorant'];eta=r['eta'];K=194/np.sqrt(10);base=target.mean+q[0]+q[1]+(3+k)*q[2]+(15+15*k+K*k**1.5)*q[3]+7*eta*k/90;slope=eta*k/30-633*q[3]*k**1.5/15
 r['shape_endpoint_prices']=[base,base+slope];r['zeta_slope']=slope
 print('CUBIC HEAD',a,d,'mean',target.mean,'shape_prices',r['shape_endpoint_prices'],'quad_convergence',r['mean_quadrature_difference'],r['derivative_quadrature_difference'],flush=True)
 return r
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--a',type=float,default=2.925);ap.add_argument('--d',type=float,nargs='+',default=[0.,.1,.2,.4]);ap.add_argument('--kappa',type=float,default=.53);ap.add_argument('--out',default='cubic_head_source_discovery.json');args=ap.parse_args()
 # The inherited warm start has degree12; pad rational Chebyshev block.
 old=None;rows=[]
 for d in args.d:
  r=solve(args.a,d,args.kappa,n=12,old=old);old=r;rows.append(r);(HERE/args.out).write_text(json.dumps({'status':'NUMERICAL_DISCOVERY_ONLY','rows':rows},indent=2))
