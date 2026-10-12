#!/usr/bin/env python3
"""Floating witness minimax screen for P1-bP3+tPodd. No lower G endpoint.

The diffuse cubic source is a Gaussian CLT limit of normalized third chaoses.
For odd q(Z), use f=sign(q(Z)+Y), g=sign(q(Z)-Y), Y the cubic CLT.
All reported scalar values are adversarial *lower* tests for operator norm.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
from pathlib import Path
import json,time,numpy as np
from scipy.special import erf
from scipy.optimize import minimize,linprog,differential_evolution
from numpy.polynomial.legendre import leggauss
HERE=Path(__file__).resolve().parent
C=np.sqrt(2/np.pi);C2=2/np.pi

def setup(n=80):
 x,w=leggauss(n);zz=[];ww=[]
 for lo,hi in [(0,.25),(.25,.5),(.5,1),(1,2),(2,4),(4,8),(8,12)]:
  z=lo+(x+1)*(hi-lo)/2;zz.extend(z);ww.extend(2*w*(hi-lo)/2*np.exp(-z*z/2)/np.sqrt(2*np.pi))
 return np.array(zz),np.array(ww)
Z,W=setup();H=(Z**3-3*Z)/np.sqrt(6)

def basis(z,quintic=False):
 return np.column_stack([np.ones_like(z),z,z**3]+([z**5]if quintic else[]))
def witness(qcoef,z=Z,w=W):
 V=basis(z,len(qcoef)==4);q=V@np.asarray(qcoef);h=(z**3-3*z)/np.sqrt(6);e=np.exp(-q*q/2);s=erf(q/np.sqrt(2))
 u=w@(z*s);v=w@(h*s);r=C*(w@e);m=2*(w@np.abs(s))-1
 return {'q_coefficients_signz_z_z3_z5':list(map(float,qcoef)),'u':float(u),'v':float(v),'r':float(r),'identity_pairing':float(m),'U':float(u*u),'D':float(r*r-v*v),'M':float(m)}
def scalar(rec,b,t):return rec['U']+b*rec['D']+t*rec['M']
def optimize(b,t,old=None,quintic=False):
 V=basis(Z,quintic);nn=V.shape[1]
 bounds=[(-.75,.75),(-2,8),(-.4,2)]+([(-.15,.15)]if quintic else[])
 def fun(x):
  q=V@x;e=np.exp(-q*q/2);sf=erf(q/np.sqrt(2));u=W@(Z*sf);v=W@(H*sf);r=C*(W@e);m=2*(W@abs(sf))-1
  val=u*u+b*(r*r-v*v)+t*m
  fac=2*C*e*(u*Z-b*v*H-b*r*q+t*np.sign(q))
  grad=V.T@(W*fac)
  return -val,-grad
 starts=[np.array([0.,2.55,.39]),np.array([.03,2.5,.4]),np.array([-.03,2.5,.4]),np.array([0.,1.,.1])]
 if old is not None:starts.insert(0,np.array(old['q_coefficients_signz_z_z3_z5'][:3]))
 best=None
 for s in starts:
  if quintic:s=np.r_[s,0.]
  opt=minimize(fun,s,jac=True,bounds=bounds,method='L-BFGS-B',options={'ftol':1e-15,'gtol':1e-10,'maxiter':500})
  if best is None or opt.fun<best.fun:best=opt
 rec=witness(best.x);rec.update(b=float(b),t=float(t),scalar_test=scalar(rec,b,t),optimizer_status=str(best.message),iterations=int(best.nit))
 return rec

def simple_witnesses():
 # Exact hemisphere and two pure cubic anti-sign CLT/fixed-coordinate limits.
 e=np.exp(-1.5);u=C*(2*e-1);v=C*(1+4*e)/np.sqrt(6)
 return [{'label':'hemisphere_equal_signs','U':C2,'D':-C2/6,'M':1.},
 {'label':'pure_diffuse_cubic_anti_signs','U':0.,'D':C2,'M':-1.},
 {'label':'pure_coordinate_H3_anti_signs','U':-u*u,'D':v*v,'M':-1.}]

def minimax_screen(blo=.70,bhi=.84,tlo=-.02,thi=.04):
 records=simple_witnesses();history=[];old=None
 # normalized variables x=b/(1+t), y=t/(1+t); normalized C=U+D*x+(M-U)*y.
 for it in range(40):
  mat=[];rhs=[]
  for r in records:mat.append([r['D'],r['M']-r['U'],-1.]);rhs.append(-r['U'])
  mat.extend([[-1.,-blo,0.],[1.,bhi,0.] ]);rhs.extend([-blo,bhi])
  lp=linprog([0.,0.,1.],A_ub=mat,b_ub=rhs,bounds=[(None,None),(tlo/(1+tlo),thi/(1+thi)),(0.,None)],method='highs')
  assert lp.success,lp.message
  x,y,cap=lp.x;t=y/(1-y);b=x/(1-y);old=optimize(b,t,old)
  actual=max(scalar(r,b,t) for r in records+[old]);normalized=actual/(1+t);gap=normalized-cap
  history.append(dict(iteration=it,b=float(b),t=float(t),LP_normalized_norm_lower=float(cap),new_normalized_witness=float(normalized),separation_gap=float(gap),screen_ceiling=float((1+t)/actual),new_witness=old))
  records.append(old)
  print(json.dumps({k:history[-1][k]for k in ['iteration','b','t','separation_gap','screen_ceiling']}),flush=True)
  if gap<2e-11:break
 return {'status':'FLOATING_SELECTED_WITNESS_MINIMAX_CEILING_ONLY','b_range':[blo,bhi],'t_range':[tlo,thi],'history':history,'records':records,'best_parameters':{'b':float(b),'t':float(t)},'normalized_witness_norm_lower':float(normalized),'screen_method_ceiling':float(1/normalized)}

def asymmetric_probe(b,t,quintic=False,seed=20261002):
 V=basis(Z,quintic);nn=V.shape[1];bounds=[(-.3,.3),(-2,8),(-.3,1.5)]+([(-.05,.05)]if quintic else[])
 def fun(xx):
  qa=V@xx[:nn];qb=V@xx[nn:];sa=erf(qa/np.sqrt(2));sb=erf(qb/np.sqrt(2))
  ua=W@(Z*sa);ub=W@(Z*sb);va=W@(H*sa);vb=W@(H*sb);ra=C*(W@np.exp(-qa*qa/2));rb=C*(W@np.exp(-qb*qb/2));mm=W@(abs(sa+sb)-1)
  return -(ua*ub-b*va*vb+b*ra*rb+t*mm)
 opt=differential_evolution(fun,bounds*2,seed=seed,popsize=12,maxiter=500,tol=1e-10,polish=True)
 za,wa=setup(192);VV=basis(za,quintic);Ha=(za**3-3*za)/np.sqrt(6);qa=VV@opt.x[:nn];qb=VV@opt.x[nn:];sa=erf(qa/np.sqrt(2));sb=erf(qb/np.sqrt(2));val=(wa@(za*sa))*(wa@(za*sb))-b*(wa@(Ha*sa))*(wa@(Ha*sb))+b*C2*(wa@np.exp(-qa*qa/2))*(wa@np.exp(-qb*qb/2))+t*(wa@(abs(sa+sb)-1))
 return {'label':'asymmetric_quintic'if quintic else'asymmetric_cubic_with_odd_constant','coefficients_f':opt.x[:nn].tolist(),'coefficients_g':opt.x[nn:].tolist(),'scalar_test':float(-opt.fun),'heldout_scalar_test':float(val),'quadrature_difference':float(abs(val+opt.fun)),'method_ceiling':float((1+t)/val),'evaluations':int(opt.nfev)}

if __name__=='__main__':
 start=time.time();out=minimax_screen();b=out['best_parameters']['b'];t=out['best_parameters']['t'];out['symmetric_quintic_check']=optimize(b,t,out['history'][-1]['new_witness'],True)
 z,w=setup(192);fine=witness(out['symmetric_quintic_check']['q_coefficients_signz_z_z3_z5'],z,w);out['heldout_symmetric_quintic']=dict(fine,scalar_test=scalar(fine,b,t))
 out['asymmetric_checks']=[asymmetric_probe(b,t,False),asymmetric_probe(b,t,True)]
 allvals=[out['normalized_witness_norm_lower']*(1+t),out['heldout_symmetric_quintic']['scalar_test']]+[r['heldout_scalar_test']for r in out['asymmetric_checks']]
 out['final_normalized_witness_norm_lower']=max(allvals)/(1+t);out['final_method_ceiling']=(1+t)/max(allvals)
 last=out['history'][-1]['new_witness'];out['stationary_odd_constant']={'optimized_q0':last['q_coefficients_signz_z_z3_z5'][0],'expected_t_over_b_r':t/(b*last['r'])};out['target_1775_normalized_margin']=1/1.775-out['final_normalized_witness_norm_lower'];out['seconds']=time.time()-start;out['scope']='Numerical adversarial LOWER tests for scalar norm, hence CEILINGS for this matrixfamily. Does not upper-bound global all-sign norm or prove any Grothendieck lowerendpoint.'
 (HERE/'odd_tail_screen.json').write_text(json.dumps(out,indent=2)+'\n')
 print(json.dumps({k:out[k]for k in ['best_parameters','final_normalized_witness_norm_lower','final_method_ceiling','target_1775_normalized_margin','stationary_odd_constant','seconds']}),flush=True)
