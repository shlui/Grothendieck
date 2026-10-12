#!/usr/bin/env python3
"""Floating sign-witness obstruction screen. No Grothendieck lower bound."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
from pathlib import Path
import json,time,argparse
import numpy as np
from scipy.special import erf
from scipy.optimize import minimize,linprog,differential_evolution
from numpy.polynomial.legendre import leggauss
HERE=Path(__file__).resolve().parent;C=np.sqrt(2/np.pi);TARGET=1.775

def setup(n=64):
 xx,ww=leggauss(n);z=[];w=[]
 for lo,hi in [(0,1),(1,2),(2,4),(4,8),(8,12)]:
  p=lo+(xx+1)*(hi-lo)/2;z.extend(p);w.extend(2*ww*(hi-lo)/2*np.exp(-p*p/2)/np.sqrt(2*np.pi))
 return np.array(z),np.array(w)
Z,W=setup();H3=(Z**3-3*Z)/np.sqrt(6);FEATURE=np.column_stack([Z,Z**3,Z**5])

def witness(ac,z=Z,w=W,jac=False):
 feat=FEATURE if z is Z else np.column_stack([z,z**3,z**5]);h3=H3 if z is Z else(z**3-3*z)/np.sqrt(6)
 p=feat@np.asarray(ac);ep=np.exp(-p*p/2);s=erf(p/np.sqrt(2));u=w@(z*s);v=w@(h3*s);r=C*(w@ep);m=2*w@np.abs(s)-1
 out=np.array([u*u,r*r-v*v,m])
 if jac:
  du=C*((w*z*ep)@feat);dv=C*((w*h3*ep)@feat);dr=-C*((w*p*ep)@feat);dm=2*C*((w*np.sign(p)*ep)@feat)
  return out,np.array([2*u*du,2*r*dr-2*v*dv,dm]),[u,v,r]
 return out,[u,v,r]

def optimize(b,t,seed=None,degree=5,global_search=False):
 bounds=[(0,8),(-2,2),(-.25,.25) if degree==5 else(0,0)]
 def fun(ac):
  row,J,_=witness(ac,jac=True);return -(row[0]+b*row[1]+t*row[2]),-(J[0]+b*J[1]+t*J[2])
 initials=[[2.55,.39,0],[.6,.6,0],[4.,0.,0],[1.,-.05,0],[.1,1.,0]]
 if seed is not None:initials.insert(0,seed)
 best=None
 for start in initials:
  opt=minimize(fun,start,jac=True,method='L-BFGS-B',bounds=bounds,options={'ftol':1e-15,'gtol':1e-11,'maxiter':500})
  if best is None or opt.fun<best.fun:best=opt
 if global_search:
  glob=differential_evolution(lambda x:fun(x)[0],bounds,seed=20261001,popsize=12,maxiter=350,tol=2e-10,polish=False)
  opt=minimize(fun,glob.x,jac=True,method='L-BFGS-B',bounds=bounds,options={'ftol':1e-15,'gtol':1e-11,'maxiter':500})
  if opt.fun<best.fun:best=opt
 row,mom=witness(best.x);return dict(coefficients=best.x.tolist(),features=row.tolist(),moments=mom,scalar=row[0]+b*row[1]+t*row[2],optimizer_success=bool(best.success))

HEMI=np.array([C*C,-C*C/6,1.])
# Opposite normalizedcubicGaussianCLT signfields and oppositesignH3(Z).
ANTI_G=np.array([0,C*C,-1.])
u3=C*(2*np.exp(-1.5)-1)
v3=(C+4*C*np.exp(-1.5))/np.sqrt(6)
ANTI_H3=np.array([-u3*u3,v3*v3,-1.])

def evaluated(b,t,q,verify_n=192):
 z,w=setup(verify_n);row,mom=witness(q['coefficients'],z,w);vals={'reflected_quintic':row[0]+b*row[1]+t*row[2], 'hemisphere':HEMI[0]+b*HEMI[1]+t,'opposite_cubic_CLT':b*C*C-t,'opposite_H3':ANTI_H3[0]+b*ANTI_H3[1]-t}
 lower=max(vals.values());norm=lower/(1+t)
 return dict(b=b,t=t,vector_numerator_limit=1+t,sampled_scalar_norm_lower=lower,normalized_norm_lower=norm,ratio_ceiling=(1+t)/lower,target_norm_headroom=1/TARGET-norm,unnormalized_target_headroom=(1+t)/TARGET-lower,witness_values=vals,reflected_coefficients=q['coefficients'],reflected_moments=mom,quadrature_discrepancy=abs(vals['reflected_quintic']-q['scalar']))

def grid(blo=.72,bhi=.80,tlo=-.02,thi=.02,nb=9,nt=17):
 records=[];seed=None
 for b in np.linspace(blo,bhi,nb):
  for t in np.linspace(tlo,thi,nt):
   q=optimize(float(b),float(t),seed);seed=q['coefficients'];rec=evaluated(float(b),float(t),q,96);records.append(rec)
 return records

def minimax(blo=.72,bhi=.80,tlo=-.02,thi=.02,seed_rows=None):
 # Normalize the operator first: x=b/(1+t),y=t/(1+t).
 # Every fixedsignwitness is affine: c0+x*c1+y*(c2-c0).
 rows=[HEMI,ANTI_G,ANTI_H3]
 if seed_rows:
  for rec in seed_rows:
   row,_=witness(rec['reflected_coefficients']);rows.append(row)
 history=[];latest=None
 for i in range(80):
  cuts=np.array([[r[1],r[2]-r[0],-1]for r in rows]);rhs=-np.array([r[0]for r in rows])
  cuts=np.vstack([cuts,[1,blo,0],[-1,-bhi,0]]);rhs=np.r_[rhs,bhi,-blo] # corrected below
  # x>=blo*(1-y): -x-blo*y<=-blo; x<=bhi*(1-y):x+bhi*y<=bhi.
  cuts[-2]=[1,bhi,0];rhs[-2]=bhi;cuts[-1]=[-1,-blo,0];rhs[-1]=-blo
  lp=linprog([0,0,1],A_ub=cuts,b_ub=rhs,bounds=[(None,None),(tlo/(1+tlo),thi/(1+thi)),(0,None)],method='highs')
  if not lp.success:raise RuntimeError(lp.message)
  x,y,lower=lp.x;b=x/(1-y);t=y/(1-y);q=optimize(b,t,latest['coefficients'] if latest else None,global_search=(i%12==0));upper=max(q['scalar']/(1+t),HEMI[0]*(1-y)+HEMI[1]*x+y,ANTI_G[1]*x-y,ANTI_H3[0]*(1-y)+ANTI_H3[1]*x-y)
  latest=q;gap=upper-lower;history.append(dict(iteration=i,b=b,t=t,LP_lower_norm=lower,found_witness_norm=upper,separation_gap=gap))
  print(history[-1],flush=True)
  rows.append(np.array(q['features']))
  if gap<3e-10:break
 return dict(status='FLOATING_MINIMAX_SIGN_WITNESS_OBSTRUCTION_ONLY',history=history,best=evaluated(b,t,q),cut_count=len(rows),scope='Witness norms are lower bounds on true scalar norms. Ratio ceilings are NOT Grothendieck lower bounds; positive targetheadroom does NOT prove globalcap.')

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--bmin',type=float,default=.72);ap.add_argument('--bmax',type=float,default=.80);ap.add_argument('--tmin',type=float,default=-.02);ap.add_argument('--tmax',type=float,default=.02);args=ap.parse_args();start=time.time()
 records=grid(args.bmin,args.bmax,args.tmin,args.tmax);out={'status':'FLOATING_SIGN_WITNESS_SCREEN_ONLY','records':records,'grid_best':max(records,key=lambda x:x['ratio_ceiling'])};(HERE/'grid_screen.json').write_text(json.dumps(out,indent=2)+'\n');print('GRIDBEST',out['grid_best'],flush=True)
 best=minimax(args.bmin,args.bmax,args.tmin,args.tmax,records);best['seconds']=time.time()-start;(HERE/'minimax_screen.json').write_text(json.dumps(best,indent=2)+'\n');print('MINIMAXBEST',best['best'],flush=True)
