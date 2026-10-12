"""Discovery only: full cubic projection residual and clipped centered Stein price.
No numerical result here is a uniform norm upper certificate.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
import numpy as np,json,time,argparse
from pathlib import Path
from scipy.optimize import minimize
from operator_gaussian_moments import moments,EP,AL,DE
from prior_joint_source_fit import features,matrix
ROOT=Path(__file__).resolve().parents[1]

def run(order=8,rounds=6,u=.755,r=.2578,kap=.1,clipped=True,residual=True,startfile=None,source68=True,pe=None,klo=0.,khi=.1,ulo=None,uhi=None,rlo=None,rhi=None):
 seed=ROOT/'inputs/joint_u0.755_r0.2578_k0.1_o8_p1_c1_L30.0.json'
 prev=json.load(open(seed))['history'][-1]['parameters']; oldbank=np.load(seed.with_suffix('.npz'))
 n=order;ia=14+n;ib=ia+1;ie=ia+2;iz=ia+3;ic=ia+4;id=ia+5;ig=ia+6;it=ia+7;iu=ia+8;sz=ia+9
 x=np.zeros(sz);x[:14]=prev[:14];x[14:14+min(order,8)]=prev[14:14+min(order,8)];x[ia:ie+1]=prev[-3:];x[ig]=.002
 # Frozen rank-one form that can be subtracted from the even matrix deficit.
 e=float(prev[2]);p=np.array([np.sqrt(e),np.sqrt(EP*(1+e/AL))]); pp=np.outer(p,p)
 theta=kap/(3*(kap+2));nu=kap*(kap+4)/(9*(kap+2));k6=(68 if source68 else 194/np.sqrt(10))*kap**1.5
 m6=15+15*kap+k6;ew=np.array([1,1,3+kap,m6]);Ds=np.array([kap*kap,kap*kap+2.25*kap**1.5]);Bc=6+9*kap+k6
 kactual=np.repeat([klo,khi],4);shz=np.tile(np.repeat([.4,1.],2),2);Dend=kactual**2+np.tile([0.,1.],4)*2.25*kactual**1.5*shz**1.5;k6end=(194/np.sqrt(10)-(633/9)*(shz-.4))*kactual**1.5;M6end=15+15*kactual+k6end;Bend=6+9*kactual+k6end-Dend;Vare=kactual/9*(1+shz)-2*theta*kactual/3+theta**2*(kactual+2);cov=kactual/3-theta*(kactual+2);ew[-1]=0; x[it]=float(max(x[12]*(kactual-kap)+x[ib]*cov+x[13]*M6end+x[ie]*Vare/2+x[iz]**2*Bend+(x[ig]-x[id])*Dend-x[ic]*kactual))
 if startfile:
  dat=json.load(open(startfile));arr=dat['history'][-1]['x'];oo=dat['order'];x[:14]=arr[:14];x[14:14+min(n,oo)]=arr[14:14+min(n,oo)];x[ia:it+1]=arr[14+oo:14+oo+8];p=np.array(dat['p']);pp=np.outer(p,p)
 if pe is not None:
  p=np.array([np.sqrt(pe),np.sqrt(EP*(1+pe/AL))]);pp=np.outer(p,p);x[2]=pp[0,0];x[3]=pp[0,1];x[4]=pp[1,1]
 def allblocks(z):
  a,b,c,d,e,f=z[2:8]
  return [(a,b,c),(d,e,f),(a+AL,b,c-EP),(d+AL,e,f-EP),(d-EP,e,f),(a-pp[0,0],b-pp[0,1],c-pp[1,1])]
 def repair(z):
  z=z.copy();mineig=min(np.linalg.eigvalsh([[a,b],[b,c]])[0] for a,b,c in allblocks(z));shift=max(0,2e-9-mineig);z[[2,4,5,7]]+=shift;return z,shift
 x,_=repair(x)
 def psd(z):return np.array([v for a,b,c in allblocks(z) for v in [a,c,a*c-b*b]])
 def psdj(z):
  out=np.zeros((18,sz))
  for i,((a,b,c),(ai,bi,ci)) in enumerate(zip(allblocks(z),[(2,3,4),(5,6,7),(2,3,4),(5,6,7),(5,6,7),(2,3,4)])):
   out[3*i,ai]=1;out[3*i+1,ci]=1;out[3*i+2,ai]=c;out[3*i+2,bi]=-2*b;out[3*i+2,ci]=a
  return out
 def outer(z,U,R):
  y=z[9]+z[3]*(R-r)
  return AL*U*U-2*z[0]*U+(EP-z[2])*R*R-2*z[8]*R+y*y/z[4]
 us=[u if ulo is None else ulo,u if uhi is None else uhi];rs=[r if rlo is None else rlo,r if rhi is None else rhi]
 x[iu]=max(outer(x,U,R)for U in us for R in rs)
 def objective(z):return DE+z[1]**2/EP+z[10:14]@ew+z[it]+z[iu]
 def objj(z):
  out=np.zeros(sz);out[:14]=[0,2*z[1]/EP,0,0,0,0,0,0,0,0,*ew];out[ie]=0;out[it]=1;out[iu]=1;return out
 def outercon(z):return np.array([z[iu]-outer(z,U,R)for U in us for R in rs])
 def outerjac(z):
  J=[]
  for U in us:
   for R in rs:
    y=z[9]+z[3]*(R-r);v=np.zeros(sz);v[0]=2*U;v[2]=R*R;v[3]=-2*y*(R-r)/z[4];v[4]=y*y/z[4]**2;v[8]=2*R;v[9]=-2*y/z[4];v[iu]=1;J.append(v)
  return np.array(J)
 def rcon(z):return z[it]-z[12]*(kactual-kap)-z[ib]*cov-z[13]*M6end-z[ie]*Vare/2-z[iz]**2*Bend-(z[ig]-z[id])*Dend+z[ic]*kactual
 def rjac(z):
  J=np.zeros((len(kactual),sz));J[:,it]=1;J[:,12]=-(kactual-kap);J[:,ib]=-cov;J[:,13]=-M6end;J[:,ie]=-Vare/2;J[:,iz]=-2*z[iz]*Bend;J[:,ig]=-Dend;J[:,id]=Dend;J[:,ic]=kactual;return J
 def oracle(z,q,extra=0):
  M=matrix(z);g=z[ic]*q+z[id]*(q**3-3*q);coef=z[iz]*(-(q**3-3*q)+g/(2*z[ig]));v0=q*np.array([z[8],0,z[9]-z[3]*r,0]);v0[[0,2]]+=coef*p
  MM=M.copy();MM[np.ix_([0,2],[0,2])]+=z[iz]**2/z[ig]*pp
  def fg(v):
   B=v0+MM@v;row=moments(z[0],z[1],B);k=row[2:6];val=2*(z[0]*row[0]+z[1]*row[1]+B@k-DE*row[6])-v@MM@v
   return -val,2*MM@(v-k)
  seeds=[np.zeros(4)]+[s*np.eye(4)[j] for j in range(4) for s in [.75,-.75]]
  if extra: seeds.extend(np.random.default_rng(34+int(100*q)).uniform(-1,1,(extra,4)))
  best=(-np.inf,None)
  for v in seeds:
   sol=minimize(fg,v,jac=True,method='L-BFGS-B',bounds=[(-1,1)]*4,options={'maxiter':120,'ftol':1e-13,'gtol':5e-9});B=v0+MM@sol.x;row=moments(z[0],z[1],B);k=row[2:6];K=p@k[[0,2]]
   val=2*(z[0]*row[0]+z[1]*row[1]+q*z[8]*k[0]+q*(z[9]-z[3]*r)*k[2]-DE*row[6])+k@M@k-2*z[iz]*(q**3-3*q)*K+(z[iz]*K+g/2)**2/z[ig]
   if val>best[0]:best=(val,row)
  return best
 qs=np.unique(np.r_[np.arange(0,8.001,.1),np.arange(3.8,4.601,.05),np.arange(8,16.001,.25),np.arange(16,64.001,4.)])
 bank=list(zip(oldbank['q'],oldbank['row']))
 if startfile and Path(startfile).with_suffix('.npz').exists():
  saved=np.load(Path(startfile).with_suffix('.npz'));bank=list(zip(saved['q'],saved['row']))
 name=f'robust_interval_shape_e{EP}_n{n}_u{u}_r{r}_k{klo}_{khi}_pe{pe}'
 history=[]
 for step in range(rounds):
  tstart=time.time();rows=[];R=[];H=[];lows=[];scales=[];rhs=[]
  for q,row in bank:
   h1,h3=row[:2];a,b,c,d=row[2:6];K=p@np.array([a,c]);S=q**3-3*q;ph,dp=features(q,n);hh=q*q-1
   A=np.zeros(sz);A[:14]=[2*h1,2*h3,a*a,2*a*c-2*q*r*c,c*c,b*b,2*b*d,d*d,2*q*a,2*q*c,-1,-q*q,-q**4,-q**6];A[14:14+n]=-q*ph+(1+theta*hh)*dp;A[iz]=-2*S*K
   h=np.zeros(sz);h[14:14+n]=dp;h[ia]=-1;h[ib]=-hh
   rr=np.zeros(sz);rr[iz]=K;rr[ic]=q/2;rr[id]=S/2
   scale=max(1,np.linalg.norm(A[:14]));rows.append(A);R.append(rr);H.append(h);lows.append(-1-theta*hh);scales.append(scale);rhs.append(2*DE*row[6])
  A=np.array(rows);R=np.array(R);H=np.array(H);low=np.array(lows);scale=np.array(scales);rhs=np.array(rhs)
  def cons(z):
   v=H@z;rr=R@z;eta=z[ie];pen=v*v/(2*eta)
   if clipped:pen=np.where(v>=eta*low,pen,low*v-eta*low*low/2)
   return (rhs-A@z-rr*rr/z[ig]-pen)/scale
  def cjac(z):
   v=H@z;rr=R@z;eta=z[ie];dv=v/eta;de=-v*v/(2*eta**2)
   if clipped:inside=v>=eta*low;dv=np.where(inside,dv,low);de=np.where(inside,de,-low*low/2)
   J=-A-2*rr[:,None]*R/z[ig]-dv[:,None]*H;J[:,ig]+=rr*rr/z[ig]**2;J[:,ie]-=de;return J/scale[:,None]
  bounds=[(0,30),(-5,5),(1e-7,15),(-15,15),(1e-7,15),(1e-7,15),(-15,15),(1e-7,15),(-10,10),(-10,10),(-100,100),(-100,100),(-10,10),(0,10)]+[(-5,5)]*n+[(-5,5),(-5,5),(1e-8,10),(-1,1),(-5,5),(-5,5),(1e-8,10),(-10,10),(-100,100)]
  if not residual:
   for idx in [iz,ic,id,it]:bounds[idx]=(0,0)
   bounds[ig]=(1e-8,1e-8);x[iz]=x[ic]=x[id]=x[it]=0.;x[ig]=1e-8
  sol=minimize(objective,x,jac=objj,method='SLSQP',bounds=bounds,constraints=[{'type':'ineq','fun':psd,'jac':psdj},{'type':'ineq','fun':rcon,'jac':rjac},{'type':'ineq','fun':outercon,'jac':outerjac},{'type':'ineq','fun':cons,'jac':cjac}],options={'ftol':2e-12,'maxiter':450})
  x,shift=repair(sol.x);x[it]=max(x[it],max(x[12]*(kactual-kap)+x[ib]*cov+x[13]*M6end+x[ie]*Vare/2+x[iz]**2*Bend+(x[ig]-x[id])*Dend-x[ic]*kactual))
  gaps=[]
  for q in qs:
   val,row=oracle(x,q,extra=1);ph,dp=features(q,n);hh=q*q-1;e=x[14:14+n]@dp-x[ia]-x[ib]*hh;lo=-1-theta*hh;eta=x[ie];pen=e*e/(2*eta) if not clipped or e>=eta*lo else lo*e-eta*lo*lo/2;D=q*(x[14:14+n]@ph)-(1+theta*hh)*(x[14:14+n]@dp);major=x[10:14]@np.array([1,q*q,q**4,q**6])+D-pen;gaps.append(val-major);bank.append((q,row))
  rec={'step':step,'objective':float(objective(x)),'max_grid_violation':float(max(gaps)),'q_worst':float(qs[np.argmax(gaps)]),'min_psd':float(min(psd(x))),'min_constraints':float(min(cons(x))),'rcon':rcon(x).tolist(),'repair':shift,'solver_success':bool(sol.success),'message':str(sol.message),'seconds':time.time()-tstart,'x':x.tolist()};history.append(rec)
  (ROOT/'results'/f'{name}.json').write_text(json.dumps({'status':'NUMERICAL DISCOVERY ONLY. No full local/source/parameter certificate.','u':u,'r':r,'kappa':kap,'kappa_interval':[klo,khi],'u_interval':us,'r_interval':rs,'order':n,'clipped':clipped,'residual':residual,'published_moment_constant':False,'uses_inherited_633_coupling':True,'p':p.tolist(),'m6':m6,'D_interval':Ds.tolist(),'history':history},indent=2))
  np.savez_compressed(ROOT/'results'/f'{name}.npz',q=np.array([q for q,row in bank]),row=np.array([row for q,row in bank]));print(json.dumps({k:v for k,v in rec.items() if k!='x'}),flush=True)
  if max(gaps)<1e-7 and sol.success:break
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--u-low',type=float);p.add_argument('--u-high',type=float);p.add_argument('--r-low',type=float);p.add_argument('--r-high',type=float);p.add_argument('--epsilon',type=float,default=.72);p.add_argument('--kappa-low',type=float,default=0.);p.add_argument('--kappa-high',type=float,default=.1);p.add_argument('--order',type=int,default=8);p.add_argument('--rounds',type=int,default=4);p.add_argument('--u',type=float,default=.755);p.add_argument('--r',type=float,default=.2578);p.add_argument('--kappa',type=float,default=.1);p.add_argument('--start');p.add_argument('--no-clip',action='store_true');p.add_argument('--no-residual',action='store_true');p.add_argument('--inherited-moment',action='store_true');p.add_argument('--pe',type=float);a=p.parse_args();EP=a.epsilon;run(a.order,a.rounds,a.u,a.r,a.kappa,not a.no_clip,not a.no_residual,a.start,not a.inherited_moment,a.pe,a.kappa_low,a.kappa_high,a.u_low,a.u_high,a.r_low,a.r_high)
