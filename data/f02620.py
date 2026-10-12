"""Directed whole-line Stein majorants for smoothed absolute-value source ceilings."""
from pathlib import Path
from fractions import Fraction as F
import json,hashlib,time
from functools import lru_cache
from flint import arb,arb_series,ctx,fmpq,fmpq_poly
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"reduced_search/energy"))
from halfline_stein_verifier import aa,qq,pp,horner,reduced_coefficients,ceilq
HERE=Path(__file__).resolve().parent;ctx.prec=512
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()

def Rval(x):return (arb.pi()/2).sqrt()*(x*x/2).exp()*(x/arb(2).sqrt()).erfc()
@lru_cache(None)
def M(q):return Rval(aa(q))
def jetR(pos,r,scale=arb(1)):
 co=[r,pos*r-1]
 for j in range(1,4):co.append((pos*co[-1]+co[-2])/(j+1))
 return arb_series([v*scale**j for j,v in enumerate(co)])
def rawjet(row,pos,r0,r1):
 ctx.cap=5;X=arb_series([pos,1]);d=X+2;n=row['degree'];N,B=reduced_coefficients(n,'2',tuple(row['phi_coefficients']))
 phi=horner(N,X)/d**(n+1);der=horner(B,X)/d**(n+2)
 a=aa(row['a']);bb=(1+a*a).sqrt();c=bb/a;C=(2/arb.pi()).sqrt()
 G=1-X*jetR(pos,r0);G1=1-c*X*jetR(c*pos,r1,c)
 fh=-C*bb*G+C*a*(-X*X/(2*a*a)).exp()*G1
 Q=sum((aa(v)*X**(2*j) for j,v in enumerate(row['majorant'])),arb_series([0]))
 return Q-der+X*phi,fh-der

def exact_tail(row,R):
 X=pp([0,1]);d=X+2;z=X-2;n=row['degree'];v=list(map(F,row['phi_coefficients']));A=[pp([1]),z]
 for j in range(2,n+1):A.append(2*z*A[-1]-d*d*A[-2])
 N=sum((qq(v[j]/(j+1))*X**(j+1) for j in range(4)),pp([0]))*d**(n+1)
 for j,c in enumerate(v[4:]):N+=qq(c)*X*A[j]*d**(n-j)
 B=N.derivative()*d-(n+1)*N;D=d**(n+2);Q=sum((qq(t)*X**(2*j) for j,t in enumerate(map(F,row['majorant']))),pp([0]));eta=F(row['eta'])
 a=aa(row['a']);cu=ceilq((2/arb.pi()).sqrt()*(1+a*a).sqrt());checks=[]
 if v[3]<0:
  # f'_J lies in[-cu/x²,0]. Prove residual>=eta, then use the clipped branch.
  polys=[X**2*(-B-qq(eta)*D)-qq(cu)*D,(Q+qq(eta/2))*d**(n+1)+X*N]
  mode='POSITIVE_RESIDUAL_CLIPPED_LINEAR'
 else:
  polys=[2*qq(eta)*X**4*D*(Q*D-B+X*N*d)-(-X**2*B-qq(e*cu)*D)**2 for e in [0,1]]
  mode='QUADRATIC_YOUNG_ENDPOINTS'
 for i,P in enumerate(polys):
  co=P(pp([R,1])).coeffs();assert co and co[0]>0 and all(t>=0 for t in co),('tailfail',i,str(R),sum(t<0 for t in co))
  checks.append(dict(index=i,degree=P.degree(),minimum_translated_coefficient=str(min(co)),polynomial_sha256=hashlib.sha256('|'.join(map(str,P.coeffs())).encode()).hexdigest()))
 return dict(domain=[str(R),'infinity'],mode=mode,c_upper=str(cu),stein_derivative_bounds=['-c_upper/x²','0'],checks=checks)

def certify(raw,index):
 start=time.time();margin=F(1,10**9)
 row=dict(a=str(F(str(raw['a']))),u=str(F(str(raw['u']))),kappa_discovery=str(F(str(raw['kappa']))),degree=raw['degree'],scale='2',eta=str(F(format(raw['eta'],'.16f'))),phi_coefficients=[str(F(format(v,'.18f'))) for v in raw['phi_coefficients']],majorant=[str(F(format(v,'.18f'))+(margin if j in[0,3]else 0)) for j,v in enumerate(raw['majorant'])],added_q0_q6_margin=str(margin))
 R=F(64)
 while True:
  try:tail=exact_tail(row,R);break
  except AssertionError:R*=2;assert R<=16384
 a=aa(row['a']);c=(1+a*a).sqrt()/a;eta=aa(row['eta']);stack=[(F(0),R,0)];leaves=[];minimum=None
 while stack:
  l,r,depth=stack.pop();mid=(l+r)/2;rad=(r-l)/2
  try:
   whole0,wholeR=rawjet(row,aa(l).union(aa(r)),M(r).union(M(l)),Rval(c*aa(r)).union(Rval(c*aa(l))))
   branch='linear' if wholeR[0]>eta else 'quadratic'
   center0,centerR=rawjet(row,aa(mid),M(mid),Rval(c*aa(mid)))
   if branch=='linear':center=center0-centerR+eta/2;whole=whole0-wholeR+eta/2
   else:center=center0-centerR*centerR/(2*eta);whole=whole0-wholeR*wholeR/(2*eta)
   gap=center[0]-sum((abs(center[j])*aa(rad**j)for j in range(1,4)),arb(0))-abs(whole[4])*aa(rad**4)
  except ValueError:gap=arb(-1);branch="quadratic"
  if gap>0:
   low=F(str(gap.lower().fmpq()));minimum=low if minimum is None else min(minimum,low);leaves.append([depth,branch])
  else:
   assert depth<48,('compactfail',index,str(l),str(r),str(gap));stack.extend([(mid,r,depth+1),(l,mid,depth+1)])
 k=F(row['kappa_discovery']);q=list(map(F,row['majorant']));et=F(row['eta']);C=(2/arb.pi()).sqrt();B=(1+a*a).sqrt();K=194/arb(10).sqrt();prices=[]
 for c6,hv in [(K,F(7,90)),(K-aa(F(633,15)),F(1,9))]:
  val=C*B+aa(q[0]+q[1]+q[2]*(3+k))+aa(q[3])*(15+15*aa(k)+c6*aa(k)*aa(k).sqrt())+aa(et*k*hv)
  prices.append(F(str(val.upper().fmpq())))
 upper=max(prices)+F(1,10**20);source=upper-F(row['a'])*F(row['u'])
 row.update(status='PASS_DIRECTED_GENERIC_J_WHOLE_LINE',precision_bits=512,compact_domain=['0',str(R)],leaf_records=leaves,compact_leaves=len(leaves),minimum_compact_gap=str(minimum),tail=tail,J_upper=str(upper),J_upper_decimal=str(float(upper)),source_r_upper=str(source),source_r_upper_decimal=str(float(source)),seconds=time.time()-start)
 print(index,row['compact_leaves'],row['J_upper_decimal'],row['source_r_upper_decimal'],row['seconds'],flush=True);return row

if __name__=='__main__':
 src=HERE/'archived_Ja_discovery.json';raw=json.loads(src.read_text());out=dict(status='PASS_DIRECTED_GENERIC_J_SOURCE_CEILINGS',source_sha256=sha(src),producer_sha256=sha(Path(__file__)),scope='Whole-line scalar Stein majorants plus633moments; physicalsource supportlemma inJ_source_ceiling.tex. No fulloperator certificate.',rows=[])
 for i in [0,2]:
  r=raw['rows'][i]
  out['rows'].append(json.loads((HERE/'J_source_fresh_directed.json').read_text())['rows'][0] if i==0 else certify(r,i));(HERE/'J_source_k053_fresh_directed.json').write_text(json.dumps(out,indent=2)+'\n')
