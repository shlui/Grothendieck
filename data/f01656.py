"""Directed scalar norm LOWER witness for the fresh cubic coefficient b=.72.
This is NOT a uniform scalar norm upper certificate or Grothendieck lower bound.
The integration uses Taylor coefficients (derivatives divided by factorial),
order-N whole-cell remainders, and explicit Gaussian tail recurrences. Only
Python's standard library and the system MPFR shared library are required.
"""
from fractions import Fraction as F
from math import comb
from pathlib import Path
import argparse, hashlib, json, time
from mpfr_interval import I, selftest
HERE=Path(__file__).resolve().parent
A=F(511,200); B=F(393,1000); EPS=F(18,25)
def add(a,b): return [x+y for x,y in zip(a,b)]
def scale(a,b): return [x*b for x in a]
def mul(a,b):
 return [sum((a[j]*b[k-j] for j in range(k+1)),I(0)) for k in range(len(a))]
def expjet(a):
 r=[a[0].exp()]
 for k in range(1,len(a)):
  r.append(sum((j*a[j]*r[k-j] for j in range(1,k+1)),I(0))/k)
 return r

def jets(x,N):
 q=[I(B)*comb(3,j)*x**(3-j)+(I(A)*x if j==0 else I(A) if j==1 else I(0)) if j<=3 else I(0) for j in range(N+1)]
 xx=[x.square(),2*x,I(1)]+[I(0)]*(N-2)
 ee=expjet(scale(add(xx,mul(q,q)),F(-1,2)))
 dp=[I(A)+3*I(B)*x.square(),6*I(B)*x,3*I(B)]+[I(0)]*(N-2)
 u1=mul(dp,ee); hm=xx.copy();hm[0]-=1
 return [ee,u1,mul(hm,u1)]

def compute(cells=320,N=10,L=5,check_fixed_witness=True):
 assert N>=4 and N%2==0 and L>1
 start=time.time();h=F(L,2*cells);sums=[I(0) for _ in range(3)]
 for k in range(cells):
  mid=F((2*k+1)*L,2*cells)
  center=jets(I(mid),N);whole=jets(I(I(mid-h).lo,I(mid+h).hi),N)
  for j in range(3):
   sums[j]+=sum((center[j][i]*2*h**(i+1)/F(i+1) for i in range(0,N,2)),I(0))+whole[j][N]*2*h**(N+1)/F(N+1)
 pi=I.pi();c=I(1+A*A);el=(-c*L*L/2).exp();J0=el/(c*L);J2=L*el/c+J0/c;J4=L**3*el/c+3*J2/c
 tails=[2/pi*J0,2/pi*(A*J0+3*B*J2),2/(pi*I(6).sqrt())*(A*(J2+J0)+3*B*(J4+J2))]
 moments=[2/pi*sums[0],2/pi*sums[1],2/(pi*I(6).sqrt())*sums[2]]
 moments=[v+I(I(0).lo,t.hi) for v,t in zip(moments,tails)]
 mY,m1,m3=moments;E1=m1.square();E3=m3.square()-mY.square();value=E1-EPS*E3
 if check_fixed_witness:
  assert value.above('0.56324') and value.below(F(40,71))
  assert (1/value).above('1.77543') and (1/value).below('1.77545')
 return {'status':'PASS_DIRECTED_CLT_SCALAR_NORM_LOWER_WITNESS','engine':selftest(),'cells':cells,'Taylor_order':N,'L':L,'parameters':{'a':str(A),'c':str(B),'epsilon':str(EPS)},'moments':{k:v.record() for k,v in zip(['diffuse_cubic','first','ordinary_third'],moments)},'P1_pairing':E1.record(),'P3_pairing':E3.record(),'scalar_norm_lower_witness':value.record(),'reciprocal_norm_ceiling':(1/value).record(),'norm_margin_below_40_71':(I(F(40,71))-value).record(),'tail_bounds':[t.record() for t in tails],'seconds':time.time()-start,'definition':'Y_n=n^(-1/2) sum_j H3(X_j); f_n=sign(Y_n+aZ+cZ^3), g_n=sign(-Y_n+aZ+cZ^3), all Gaussians independent. The moments are the n-to-infinity limits.','scope':'A rigorous limiting scalar norm lower test for P1-(18/25)P3. Its reciprocal is an upper ceiling for this reciprocal-norm sequence method, not an achieved Grothendieck lower bound. This witness is below40/71 and does not exclude target1.775.','source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--cells',type=int,default=320);ap.add_argument('--order',type=int,default=10);ap.add_argument('--output',default='cubic720_norm_witness.json');a=ap.parse_args();r=compute(a.cells,a.order)
 p=Path(a.output);p=p if p.is_absolute() else HERE/p;p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(r,indent=2)+'\n')
 print(json.dumps({k:r[k] for k in ['status','scalar_norm_lower_witness','reciprocal_norm_ceiling','norm_margin_below_40_71','seconds']},indent=2))
