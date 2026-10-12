"""Exploratory joint, centered Gaussian Stein constraints for a fixed operator.
This file does NOT produce a complete norm certificate. Local searches are lower
witnesses to a local supremum; whole-line and source-parameter covers are absent.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import sys,json,time,argparse
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
ROOT=Path(__file__).resolve().parents[1]
from operator_gaussian_moments import moments,EP,AL,DE

def matrix(x):
    M=np.zeros((4,4))
    M[np.ix_([0,2],[0,2])]=[[x[2],x[3]],[x[3],x[4]]]
    M[np.ix_([1,3],[1,3])]=[[x[5],x[6]],[x[6],x[7]]]
    return M

def repair(z):
    z=z.copy()
    blocks=[np.array([[z[2],z[3]],[z[3],z[4]]]),np.array([[z[5],z[6]],[z[6],z[7]]]),np.array([[z[2]+AL,z[3]],[z[3],z[4]-EP]]),np.array([[z[5]+AL,z[6]],[z[6],z[7]-EP]]),np.array([[z[5]-EP,z[6]],[z[6],z[7]]])]
    mineig=min(np.linalg.eigvalsh(b)[0] for b in blocks)
    shift=max(0.,2e-8-mineig)
    z[[2,4,5,7]]+=shift
    return z,shift

def source_oracle(x,q,r,starts=5,seed=0):
    # The physical-moment completion avoids inverting a singular matrix.
    x,_=repair(x)
    M=matrix(x);v0=q*np.array([x[8],0,x[9]-x[3]*r,0])
    def fg(v):
        B=v0+M@v;row=moments(x[0],x[1],B);k=row[2:6]
        sup=2*(x[0]*row[0]+x[1]*row[1]+B@k-DE*row[6])
        return -sup+v@M@v,2*M@(v-k)
    rng=np.random.default_rng(seed)
    seeds=[np.zeros(4)]+[s*np.eye(4)[j] for j in range(4) for s in [.8,-.8]]+[rng.uniform(-1,1,4) for _ in range(starts)]
    best=(-np.inf,None,None)
    for v in seeds:
        sol=minimize(fg,v,jac=True,method='L-BFGS-B',bounds=[(-1,1)]*4,options={'maxiter':150,'ftol':2e-14,'gtol':2e-9})
        row=moments(x[0],x[1],v0+M@sol.x);k=row[2:6]
        val=2*(x[0]*row[0]+x[1]*row[1]+v0@k-DE*row[6])+k@M@k
        if val>best[0]:best=(val,row,sol.x)
    return best

def features(q,order):
    x=q*q;t=(x-1)/(x+1)
    base=q/(1+x);der=(1-x)/(1+x)**2; dt=4*q/(1+x)**2
    T=np.polynomial.chebyshev.chebvander(t,order-1).reshape(-1)
    Tp=np.array([np.polynomial.chebyshev.chebval(t,np.polynomial.chebyshev.chebder([0]*j+[1])) if j else 0. for j in range(order)])
    return base*T,der*T+base*Tp*dt

def run(kap=.1,u=.755,r=.258,order=8,rounds=16,ceiling=8.,project=True,center=True,start_file=None):
    vbar=2*kap/9-kap*kap/(9*(kap+2)) if project else 2*kap/9
    theta=kap/(3*(kap+2)) if project else 0.
    m6=15+15*kap+(194/np.sqrt(10))*kap**1.5
    ew=np.array([1,1,3+kap,m6])
    prev=json.load(open(ROOT/'inputs/previous_numerical_seed.json'))['history'][-1]['parameters']
    x=np.r_[np.array(prev[:14]),np.zeros(order),0.,0.,.01]
    if start_file:
        old=json.load(open(start_file))['history'][-1]
        a=np.array(old['parameters']); oo=len(a)-17
        x[:14]=a[:14]; x[14:14+min(order,oo)]=a[14:14+min(order,oo)];x[-3:]=a[-3:]
    def objective(z):
        return AL*u*u-2*z[0]*u+DE+z[1]**2/EP+(EP-z[2])*r*r-2*z[8]*r+z[9]**2/z[4]+z[10:14]@ew+z[-1]*vbar/2
    def gradient(z):
        out=np.zeros_like(z);out[:14]=np.r_[-2*u,2*z[1]/EP,-r*r,0,-z[9]**2/z[4]**2,0,0,0,-2*r,2*z[9]/z[4],ew];out[-1]=vbar/2;return out
    def psd(z):
        a,b,c,d,e,f=z[2:8]
        blocks=[(a-1e-8,b,c-1e-8),(d-1e-8,e,f-1e-8),(a+AL,b,c-EP),(d+AL,e,f-EP),(d-EP,e,f)]
        return np.array([v for a,b,c in blocks for v in [a,c,a*c-b*b]])
    def psd_jac(z):
        J=np.zeros((15,len(z)))
        blocks=[(2,3,4,z[2]-1e-8,z[3],z[4]-1e-8),(5,6,7,z[5]-1e-8,z[6],z[7]-1e-8),(2,3,4,z[2]+AL,z[3],z[4]-EP),(5,6,7,z[5]+AL,z[6],z[7]-EP),(5,6,7,z[5]-EP,z[6],z[7])]
        for i,(ai,bi,ci,a,b,c) in enumerate(blocks):
            J[3*i,ai]=1;J[3*i+1,ci]=1;J[3*i+2,ai]=c;J[3*i+2,ci]=a;J[3*i+2,bi]=-2*b
        return J
    qs=np.r_[np.arange(0,8.001,.1),10.,14.,20.,30.]
    bank=[]; rng=np.random.default_rng(5121)
    for q in qs:
        bank.append((q,source_oracle(x,q,r,1)[1]))
        for j in range(3):bank.append((q,moments(.4,.2,q*np.array([.08,0,0,0])+rng.normal(size=4)*.3)))
    history=[]
    if start_file and Path(start_file).with_suffix('.npz').exists():
        old_bank=np.load(Path(start_file).with_suffix('.npz'));bank.extend(zip(old_bank['q'],old_bank['row']))
        history=json.load(open(start_file))['history']
    offset=len(history)
    name=f'joint_u{u}_r{r}_k{kap}_o{order}_p{int(project)}_c{int(center)}_L{ceiling}.json'
    for step in range(rounds):
        it=offset+step
        rows=[]; rhs=[]; H=[]
        for q,row in bank:
            h1,h3=row[:2];a,b,c,d=row[2:6];phi,dp=features(q,order); hh=q*q-1
            A=np.zeros(len(x));A[:14]=[2*h1,2*h3,a*a,2*a*c-2*q*r*c,c*c,b*b,2*b*d,d*d,2*q*a,2*q*c,-1,-q*q,-q**4,-q**6]
            A[14:14+order]=-q*phi+(1+theta*hh)*dp
            h=np.r_[np.zeros(14),dp,-1 if center else 0.,-hh if project else 0.,0.]
            sc=max(1,np.linalg.norm(A[:14]))
            rows.append(A/sc);H.append(h/np.sqrt(sc));rhs.append(2*DE*row[6]/sc)
        A=np.array(rows);H=np.array(H);rhs=np.array(rhs)
        def constraint(z):
            v=H@z;return rhs-A@z-v*v/(2*z[-1])
        def jac(z):
            v=H@z;J=-A-v[:,None]*H/z[-1];J[:,-1]+=v*v/(2*z[-1]**2);return J
        bounds=[(0,ceiling),(-5,5),(1e-7,15),(-15,15),(1e-7,15),(1e-7,15),(-15,15),(1e-7,15),(-10,10),(-10,10),(-100,100),(-100,100),(-10,10),(0,10)]+[(-5,5)]*order+[(-5,5) if center else (0,0),(-5,5) if project else (0,0),(1e-8,10)]
        start=time.time()
        sol=minimize(objective,x,jac=gradient,constraints=[{'type':'ineq','fun':constraint,'jac':jac},{'type':'ineq','fun':psd,'jac':psd_jac}],bounds=bounds,method='SLSQP',options={'maxiter':1000,'ftol':2e-12})
        x,shift=repair(sol.x);gaps=[]
        for q in qs:
            val,row,v=source_oracle(x,q,r,2 if it<7 else 4,seed=it)
            phi,dp=features(q,order);der=x[14:14+order]@dp;ph=x[14:14+order]@phi;hh=q*q-1
            d=q*ph-(1+theta*hh)*der
            sq=der-(x[-3] if center else 0)-(x[-2]*hh if project else 0)
            major=x[10:14]@np.array([1,q*q,q**4,q**6])+d-sq*sq/(2*x[-1])
            gaps.append(val-major);bank.append((q,row))
        rec={'iteration':it,'objective':float(objective(x)),'gap':float(max(gaps)),'q_worst':float(qs[np.argmax(gaps)]),'min_psd':float(min(psd(x))),'sample_constraint_min':float(min(constraint(x))),'matrix_diagonal_repair':shift,'success':bool(sol.success),'message':str(sol.message),'time':time.time()-start,'parameters':x.tolist()}
        history.append(rec);print(name,json.dumps(rec),flush=True)
        (ROOT/'results'/name).write_text(json.dumps({'status':'NUMERICAL ONLY: no whole-line, local-max, or continuum-parameter certificate.','kappa':kap,'u':u,'r':r,'order':order,'project':project,'center':center,'ceiling':ceiling,'vbar':vbar,'m6':m6,'q_grid':qs.tolist(),'history':history},indent=2))
        np.savez_compressed((ROOT/'results'/name).with_suffix('.npz'),q=np.array([q for q,row in bank]),row=np.array([row for q,row in bank]))
        if max(gaps)<2e-8 and sol.success:break
    return name
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--u',type=float,default=.755);p.add_argument('--r',type=float,default=.258);p.add_argument('--order',type=int,default=8);p.add_argument('--rounds',type=int,default=12);p.add_argument('--kappa',type=float,default=.1);p.add_argument('--ceiling',type=float,default=8);p.add_argument('--no-project',action='store_true');p.add_argument('--no-center',action='store_true');p.add_argument('--start');a=p.parse_args()
    run(kap=a.kappa,u=a.u,r=a.r,order=a.order,rounds=a.rounds,ceiling=a.ceiling,project=not a.no_project,center=not a.no_center,start_file=a.start)
