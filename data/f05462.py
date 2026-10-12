"""Expanded numerical local/source screen. This does not certify a supremum."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
import json,time,sys,numpy as np
from pathlib import Path
from scipy.optimize import minimize
from operator_gaussian_moments import moments,EP,AL,DE
from prior_joint_source_fit import features,matrix

ROOT=Path(__file__).resolve().parents[1]
file=ROOT/'results/positive_cov_residual_n16_clip1_res1_u0.755_r0.2578_k0.1_pub1_peNone.json'
if len(sys.argv)>1:file=Path(sys.argv[1])
suffix='_refined' if file.name.startswith('refined_') else ''
d=json.load(open(file));x=np.array(d['history'][-1]['x']);n=d['order'];r=d['r'];kap=d['kappa'];p=np.array(d['p']);pp=np.outer(p,p);theta=kap/(3*(kap+2))
ia=14+n;ib=ia+1;ie=ia+2;iz=ia+3;ic=ia+4;id=ia+5;ig=ia+6;it=ia+7;igst=ia+8;ng=8;iag=igst+ng;ieg=iag+1
M=matrix(x);MM=M.copy();MM[np.ix_([0,2],[0,2])]+=x[iz]**2/x[ig]*pp
qs=np.unique(np.r_[np.arange(0,8.0001,.025),np.arange(8.25,30.01,.25),np.geomspace(32,1000,61)])
records=[];start=time.time();last=np.zeros(4)
for count,q in enumerate(qs):
 g=x[ic]*q+x[id]*(q**3-3*q);v0=q*np.array([x[8],0,x[9]-x[3]*r,0]);v0[[0,2]]+=x[iz]*(-(q**3-3*q)+g/(2*x[ig]))*p
 def fg(v):
  B=v0+MM@v;row=moments(x[0],x[1],B);k=row[2:6];val=2*(x[0]*row[0]+x[1]*row[1]+B@k-DE*row[6])-v@MM@v
  return -val,2*MM@(v-k)
 rng=np.random.default_rng(114+count);seeds=[last,np.zeros(4)]+[s*np.eye(4)[j]for j in range(4)for s in [.8,-.8]]
 best=(-np.inf,None,None)
 for v in seeds:
  sol=minimize(fg,v,jac=True,method='L-BFGS-B',bounds=[(-1,1)]*4,options={'maxiter':180,'ftol':5e-14,'gtol':3e-9});row=moments(x[0],x[1],v0+MM@sol.x);k=row[2:6];K=p@k[[0,2]]
  val=2*(x[0]*row[0]+x[1]*row[1]+q*x[8]*k[0]+q*(x[9]-x[3]*r)*k[2]-DE*row[6])+k@M@k-2*x[iz]*(q**3-3*q)*K+(x[iz]*K+g/2)**2/x[ig]
  if val>best[0]:best=(val,row,sol.x)
 last=best[2];ph,dp=features(q,n);hh=q*q-1;e=x[14:14+n]@dp-x[ia]-x[ib]*hh;lo=-1-theta*hh;eta=x[ie];pen=e*e/(2*eta)if e>=eta*lo else lo*e-eta*lo*lo/2
 D=q*(x[14:14+n]@ph)-(1+theta*hh)*(x[14:14+n]@dp)
 major=x[10:14]@np.array([1,q*q,q**4,q**6])+D-pen;gap=best[0]-major
 records.append({'q':float(q),'local_value':float(best[0]),'majorant':float(major),'gap':float(gap),'scaled_gap':float(gap/(1+q**6)),'row':best[1].tolist()})
 if count%50==0:
  (file.with_name(file.stem+'_expanded_partial.json')).write_text(json.dumps({'status':'INCOMPLETE_NUMERICAL_SCREEN','records':records},indent=2))
  print(count,len(qs),'worst',max(records,key=lambda t:t['gap']),'seconds',time.time()-start,flush=True)
out={'status':'NUMERICAL_SCREEN_ONLY_NOT_A_CERTIFICATE','source_file':file.name,'q_nodes':len(qs),'starts_per_node':10,'seconds':time.time()-start,'worst_raw':max(records,key=lambda t:t['gap']),'worst_scaled':max(records,key=lambda t:t['scaled_gap']),'matrix_eigenvalues':np.linalg.eigvalsh(M).tolist(),'records':records}
(file.with_name(file.stem+'_expanded.json')).write_text(json.dumps(out,indent=2));print(json.dumps({k:v for k,v in out.items()if k!='records'},indent=2),flush=True)
