"""Independent 768-bit whole-domain order-eight Taylor proof of J ceilings."""
from pathlib import Path
from fractions import Fraction as F
from functools import lru_cache
from math import factorial
from flint import arb,arb_series,ctx,fmpq,fmpq_poly
import json,hashlib,time
P=Path(__file__).resolve().parent;ctx.prec=768;ctx.cap=9
A=lambda x:arb(F(x).numerator)/F(x).denominator
Q=lambda x:fmpq(F(x).numerator,F(x).denominator)
poly=lambda cs:fmpq_poly([Q(c)for c in cs])
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def R(x):return (arb.pi()/2).sqrt()*(x*x/2).exp()*(x/arb(2).sqrt()).erfc()
@lru_cache(None)
def RM(x):return R(A(x))
def horner(cs,x):
 z=arb_series([0])
 for c in reversed(cs):z=z*x+A(str(c))
 return z

def audit(row):
 n=row['degree'];a=F(row['a']);eta=F(row['eta']);co=list(map(F,row['phi_coefficients']));q=list(map(F,row['majorant']));assert len(co)==n+5 and q[3]>=0 and eta>0
 x=poly([0,1]);den=x+2;diff=x-2
 # Explicit Chebyshev coefficients, not the producer's recurrence.
 T=[poly([1])]
 for j in range(1,n+1):
  v=poly([0])
  for m in range(j//2+1):
   coeff=F((-1)**m*j*factorial(j-m-1)*2**(j-2*m),2*factorial(m)*factorial(j-2*m))
   v+=Q(coeff)*diff**(j-2*m)*den**(2*m)
  T.append(v)
 N=sum((Q(co[j]/(j+1))*x**(j+1)for j in range(4)),poly([0]))*den**(n+1)
 N+=sum((Q(co[4+j])*x*T[j]*den**(n-j)for j in range(n+1)),poly([0]))
 D=den**(n+2);B=den*N.derivative()-(n+1)*N
 C=(2/arb.pi()).sqrt();beta=(1+A(a)**2).sqrt();scale=beta/A(a)
 def seriesR(pos,r,mul):
  cc=[r,mul*(pos*r-1)]
  for j in range(1,8):cc.append((mul*pos*cc[j]+mul*mul*cc[j-1])/(j+1))
  return arb_series(cc)
 def jet(pos,r0,r1,linear):
  X=arb_series([pos,1]);dd=X+2
  ph=horner(N.coeffs(),X)/dd**(n+1);dp=horner(B.coeffs(),X)/dd**(n+2)
  m0=seriesR(pos,r0,arb(1));m1=seriesR(scale*pos,r1,scale);e=(-X**2/(2*A(a*a))).exp()
  target=C*(-beta+A(a)*e+beta*X*m0-A(a)*scale*X*e*m1)
  res=target-dp
  maj=sum((A(q[j])*X**(2*j)for j in range(4)),arb_series([0]))
  return maj-dp+X*ph-(res-A(eta)/2 if linear else res*res/(2*A(eta))),res
 last=F(row['compact_domain'][1]);fresh=0;minimum=None;stack=[(F(0),last,0)]
 while stack:
  l,r,d=stack.pop();mid=(l+r)/2;rad=(r-l)/2
  try:
   # Independently choose the clipped branch over this whole interval.
   _,res=jet(A(l).union(A(r)),RM(r).union(RM(l)),R(scale*A(r)).union(R(scale*A(l))),False)
   linear=res[0]>A(eta)
   center,_=jet(A(mid),RM(mid),R(scale*A(mid)),linear)
   whole,_=jet(A(l).union(A(r)),RM(r).union(RM(l)),R(scale*A(r)).union(R(scale*A(l))),linear)
   gap=center[0]-sum((abs(center[j])*A(rad**j)for j in range(1,8)),arb(0))-abs(whole[8])*A(rad**8)
  except ValueError:gap=arb(-1)
  if gap>0:
   low=F(str(gap.lower().fmpq()));minimum=low if minimum is None else min(minimum,low);fresh+=1
  else:
   assert d<48,('compact',str(l),str(r),str(gap));stack.extend([(mid,r,d+1),(l,mid,d+1)])
 assert len(row['leaf_records'])==row['compact_leaves']
 cu=F(row['tail']['c_upper']);assert A(cu)>C*beta
 even=sum((Q(q[j])*x**(2*j)for j in range(4)),poly([0]))
 if row['tail']['mode']=='POSITIVE_RESIDUAL_CLIPPED_LINEAR':
  assert co[3]<0
  tail=[-x*x*B-Q(eta)*x*x*D-Q(cu)*D,(even+Q(eta/2))*den**(n+1)+x*N]
 else:
  tail=[2*Q(eta)*x**4*D*(even*D-B+x*N*den)-(x*x*B+Q(j*cu)*D)**2 for j in [0,1]]
 for j,v in enumerate(tail):
  tr=v(poly([last,1])).coeffs();assert tr[0]>0 and all(z>=0 for z in tr)
  digest=hashlib.sha256('|'.join(map(str,v.coeffs())).encode()).hexdigest();assert digest==row['tail']['checks'][j]['polynomial_sha256']
 k=F(row['kappa_discovery']);K=194/arb(10).sqrt()
 for j in [0,1]:
  coeff=K-A(F(633*j,15));hv=F(7,90)if j==0 else F(1,9)
  val=C*beta+A(q[0]+q[1]+3*q[2]+15*q[3])+A(k*(q[2]+15*q[3]+eta*hv))+A(q[3])*coeff*(A(k)**3).sqrt()
  assert val<A(row['J_upper'])
 assert F(row['source_r_upper'])==F(row['J_upper'])-a*F(row['u'])
 return dict(a=row['a'],kappa=row['kappa_discovery'],source_upper=row['source_r_upper'],producer_leaves=row['compact_leaves'],fresh_leaves=fresh,minimum_fresh_gap=str(minimum),tail='PASS_EXACT_RECONSTRUCTION',price='PASS_768')

if __name__=='__main__':
 source=P/'J_source_k053_fresh_directed.json';data=json.loads(source.read_text());start=time.time();out=dict(status='PASS_INDEPENDENT_GENERIC_J_SOURCE_CEILINGS',precision_bits=768,compact_method='FRESH_ADAPTIVE_TAYLOR_ORDER_8',source_sha256=sha(source),auditor_sha256=sha(Path(__file__)),results=[])
 for r in data['rows']:
  out['results'].append(audit(r));print(out['results'][-1]['a'],out['results'][-1]['kappa'],'PASS',time.time()-start,flush=True)
 out['seconds']=time.time()-start;(P/'J_source_fresh_independent.json').write_text(json.dumps(out,indent=2)+'\n')
