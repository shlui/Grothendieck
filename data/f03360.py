#!/usr/bin/env python3
"""Small numerical feasibility refit; never a scalar certificate."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
from pathlib import Path
import json,time
import numpy as np
import cvxpy as cp
from numpy.polynomial import chebyshev as ch
from scipy.special import erfcx

HERE=Path(__file__).resolve().parent
C=np.sqrt(2/np.pi);K=194/np.sqrt(10)


def basis(x,n=12):
    z=(x-2)/(x+2);dz=4/(x+2)**2;T=ch.chebvander(z,n)
    DT=np.column_stack([ch.chebval(z,ch.chebder(np.eye(n+1)[j]))*dz for j in range(n+1)])
    p=x/(x+2);dp=2/(x+2)**2
    V=np.column_stack([x,x*x/2,x**3/3,x**4/4,p[:,None]*T])
    D=np.column_stack([np.ones_like(x),x,x*x,x**3,dp[:,None]*T+p[:,None]*DT])
    return V,D


def solve(lo,hi,iterations=4):
    start=time.time()
    x=np.unique(np.r_[np.linspace(0,4,801),np.linspace(4,16,401),np.geomspace(16,1000,160)])
    dense=np.unique(np.r_[np.linspace(0,64,64001),np.geomspace(64,1e5,2001)])
    history=[]
    for it in range(iterations):
        w=1/(1+x**6);sw=np.sqrt(w);V,D=basis(x);A=D-x[:,None]*V
        _,sig,vt=np.linalg.svd(np.vstack([A*w[:,None],D*sw[:,None]]),full_matrices=False)
        keep=sig>sig[0]*1e-8;Tr=vt[keep].T/sig[keep]
        Av=w[:,None]*(A@Tr);Dv=sw[:,None]*(D@Tr)
        gg=1-x*np.sqrt(np.pi/2)*erfcx(x/np.sqrt(2));base=-C*sw*gg
        W=np.column_stack([w,w*x*x,w*x**4,w*x**6])
        v=cp.Variable(Tr.shape[1]);q=cp.Variable(4);eta=cp.Variable();cap=cp.Variable()
        z=W@q-Av@v;res=base-Dv@v
        cons=[eta>=1e-6,eta<=20,q[3]>=1e-9,z>=0,
              cp.SOC(eta+z,cp.vstack([np.sqrt(2)*res,eta-z]),axis=0),
              cp.SOC(eta+q[3]-1e-9,cp.hstack([np.sqrt(2)*(Tr@v)[3],eta-q[3]+1e-9]))]
        for k in [lo,hi]:
            for s in [0,1]:
                mom=np.array([1,1,3+k,15+15*k+(K-633*s/15)*k**1.5])
                cons.append(C+mom@q+eta*k*(7/90+s/30)<=cap)
        prob=cp.Problem(cp.Minimize(cap),cons)
        prob.solve(solver='CLARABEL',tol_gap_abs=1e-10,tol_gap_rel=1e-10,tol_feas=1e-10,max_iter=350,max_threads=1)
        if q.value is None:raise RuntimeError(prob.status)
        phi=Tr@v.value;qq=q.value;e=float(eta.value)
        VV,DD=basis(dense);rr=-C*(1-dense*np.sqrt(np.pi/2)*erfcx(dense/np.sqrt(2)))-DD@phi
        gap=((DD-dense[:,None]*VV)@phi+rr*rr/(2*e)-np.polynomial.polynomial.polyval(dense*dense,qq))/(1+dense**6)
        index=np.argsort(gap)[-80:];vi=float(gap.max())
        rec=dict(iteration=it,kappa_interval=[lo,hi],sampled_cap=float(cap.value),scaled_violation=vi,rank=int(sum(keep)),seconds=time.time()-start,status=prob.status)
        history.append(rec);print(json.dumps(rec),flush=True)
        out=dict(status='NUMERICAL_COUPLED_633_CUBIC_L1_REFIT',degree=12,scale=2,
                 kappa=(lo+hi)/2,sigma=1,odd_weights=[1,0,0],sixth_offset=0,combined_sixth=False,
                 eta=e,phi_coefficients=phi.tolist(),majorant=qq.tolist(),kappa_interval=[lo,hi],history=history)
        (HERE/f'refit_{str(lo).replace(".","_")}_{str(hi).replace(".","_")}.json').write_text(json.dumps(out,indent=2)+'\n')
        if vi<1e-10:break
        x=np.unique(np.r_[x,dense[index]])
    return out


if __name__=='__main__':
    for lo,hi in [(6.25,7.3),(8.85,10.3)]:solve(lo,hi)
