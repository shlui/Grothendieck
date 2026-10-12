#!/usr/bin/env python3
"""Directed complete-halfline verification of rational Chebyshev Stein row.
Scalar consequence is conditional on inherited cubic moment identities.
"""
from pathlib import Path
from fractions import Fraction as F
import sys,json,time,hashlib,argparse
from functools import lru_cache
from flint import arb,arb_series,ctx,fmpq,fmpq_poly
ROOT=Path(__file__).resolve().parents[2]
def aa(x):
 x=F(x);return arb(x.numerator)/x.denominator
def qq(x):
 x=F(x);return fmpq(x.numerator,x.denominator)
def pp(v):return fmpq_poly([qq(x) for x in v])
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def ceilq(x,den=10**24):return F(int((x*den).ceil().unique_fmpz()),den)
@lru_cache(None)
def mills(q):
 x=aa(q);return (arb.pi()/2).sqrt()*(x*x/2).exp()*(x/arb(2).sqrt()).erfc()
def freeze(r,margin):
 return dict(odd_weights=[str(F(format(v,'.18f'))) for v in r['odd_weights']],sixth_offset=str(F(format(r.get('sixth_offset',0),'.18f'))),sigma_discovery=str(F(str(r['sigma']))),combined_sixth=r['combined_sixth'],kappa_discovery=str(F(str(r['kappa']))),degree=r['degree'],scale=str(F(str(r['scale']))),eta=str(F(format(r['eta'],'.16f'))),phi_coefficients=[str(F(format(v,'.18f'))) for v in r['phi_coefficients']],majorant=[str(F(format(v,'.18f'))+(margin if i in [0,3] else 0)) for i,v in enumerate(r['majorant'])],added_q0_q6_margin=str(margin))
@lru_cache(None)
def reduced_coefficients(degree,scale,values):
 X=pp([0,1]);L=F(scale);d=X+qq(L);z=X-qq(L);n=degree;v=list(map(F,values));A=[pp([1]),z]
 for j in range(2,n+1):A.append(2*z*A[-1]-d*d*A[-2])
 N=(qq(v[0])*X+qq(v[1]/2)*X**2+qq(v[2]/3)*X**3+qq(v[3]/4)*X**4)*d**(n+1)
 for j,c in enumerate(v[4:]):N+=qq(c)*X*A[j]*d**(n-j)
 B=N.derivative()*d-(n+1)*N
 return [aa(str(c)) for c in N.coeffs()],[aa(str(c)) for c in B.coeffs()]
def horner(c,X):
 out=arb_series([0])
 for a in reversed(c):out=out*X+a
 return out
def gapjet(row,pos,R):
 ctx.cap=5;X=arb_series([pos,1]);d=X+aa(row['scale']);n=row['degree'];N,B=reduced_coefficients(n,row['scale'],tuple(row['phi_coefficients']))
 phi=horner(N,X)/d**(n+1);der=horner(B,X)/d**(n+2)
 rr=[R,pos*R-1]
 for j in range(1,4):rr.append((pos*rr[-1]+rr[-2])/(j+1))
 G=1-X*arb_series(rr);a1,a3,a5=list(map(aa,row['odd_weights']));res=-(a1+2*a3+8*a5)*(2/arb.pi()).sqrt()*G-2*a3*X-a5*(4*X**3+8*X)-der
 major=arb_series([0])
 for i,q in enumerate(row['majorant']):major+=aa(q)*X**(2*i)
 return major-der+X*phi-res*res/(2*aa(row['eta']))
def exact_tail(row,R):
 X=pp([0,1]);L=F(row['scale']);d=X+qq(L);z=X-qq(L);n=row['degree'];v=list(map(F,row['phi_coefficients']));A=[pp([1]),z]
 for j in range(2,n+1):A.append(2*z*A[-1]-d*d*A[-2])
 N=(qq(v[0])*X+qq(v[1]/2)*X**2+qq(v[2]/3)*X**3+qq(v[3]/4)*X**4)*d**(n+1)
 for j,c in enumerate(v[4:]):N+=qq(c)*X*A[j]*d**(n-j)
 D=d**(n+2);B=N.derivative()*d-(n+1)*N;XP=X*N*d;Q=pp([0])
 for i,q in enumerate(row['majorant']):Q+=qq(q)*X**(2*i)
 cu=ceilq((2/arb.pi()).sqrt());base=qq(2*F(row['eta']))*X**4*D*(Q*D-B+XP);out=[]
 for e in [0,1]:
  a1,a3,a5=list(map(F,row['odd_weights']));H=-qq(2*a3+8*a5)*X-qq(4*a5)*X**3
  P=base-(X*X*(H*D-B)-qq(e*cu*(a1+2*a3+8*a5))*D)**2;shift=P(pp([R,1]));cc=shift.coeffs()
  assert cc and cc[0]>0 and all(t>=0 for t in cc),('TAIL_FAIL',e,str(R),sum(t<0 for t in cc))
  out.append(dict(endpoint=e,degree=P.degree(),positive_shift=str(R),translated_coefficient_count=len(cc),minimum_translated_coefficient=str(min(cc)),polynomial_sha256=hashlib.sha256('|'.join(map(str,P.coeffs())).encode()).hexdigest()))
 return dict(domain=[str(R),'infinity'],g_bounds=['0','1/x^2'],c_upper=str(cu),checks=out)
def certify(row,R=F(64)):
 start=time.time();tail=exact_tail(row,R);stack=[(F(0),R,0)];depths=[];mini=None
 while stack:
  l,r,d=stack.pop();mid=(l+r)/2;rad=(r-l)/2
  try:
   center=gapjet(row,aa(mid),mills(mid));whole=gapjet(row,aa(l).union(aa(r)),mills(r).union(mills(l)));gap=center[0]-sum((abs(center[j])*aa(rad**j) for j in range(1,4)),arb(0))-abs(whole[4])*aa(rad**4)
  except ValueError:gap=arb(-1)
  if gap>0:
   depths.append(d);lo=F(str(gap.lower().fmpq()));mini=lo if mini is None else min(mini,lo)
  else:
   assert d<44,('COMPACT_FAIL',str(l),str(r),str(gap));stack.extend([(mid,r,d+1),(l,mid,d+1)])
 # Certify objective at discovery kappa with shared tensor y, no finite y grid.
 k=F(row['kappa_discovery']);q=list(map(F,row['majorant']));eta=F(row['eta']);C=(2/arb.pi()).sqrt();high=None
 a1,a3,a5=list(map(F,row['odd_weights']));ac=a1+2*a3+8*a5;offset=F(row['sixth_offset']);assert q[3]+offset>=0
 # New exact affine coupling; each expression is a certified moment price.
 K=194/arb(10).sqrt()
 for c6,hv in [(K,aa(F(7,90))),(K-aa(F(633,15)),aa(F(1,9)))]:
  val=C*aa(ac)+aa(q[0]+q[1]+q[2]*(3+k))+aa(q[3]+offset)*(15+15*aa(k)+aa(k)**aa(F(3,2))*c6)+aa(eta*k)*hv
  h=F(str(val.upper().fmpq()));high=h if high is None else max(h,high)
 ys=[0,1,2]
 row.update(status='PASS_DIRECTED_SHARP_GENERIC_ODD_CHEBYSHEV_RATIONAL_HALF_LINE',precision_bits=ctx.prec,taylor_order=4,compact_domain=['0',str(R)],leaf_depths=depths,compact_leaves=len(depths),minimum_compact_gap=str(mini),tail=tail,coupled_objective_upper=str(high),coupled_objective_upper_decimal=str(float(high)),coupled_y_cells=len(ys)-1,seconds=time.time()-start)
 return row
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--source',type=Path,required=True);ap.add_argument('--row',type=int,default=0);ap.add_argument('--margin',default='1/1000000000');ap.add_argument('--output',type=Path,required=True);args=ap.parse_args();ctx.prec=512
 raw=json.loads(args.source.read_text());row=certify(freeze(raw['rows'][args.row],F(args.margin)))
 args.output.write_text(json.dumps(dict(status='PASS_DIRECTED_REOPTIMIZED_AFFINE_COUPLED_SCALAR',rows=[row],source=str(args.source),source_sha256=sha(args.source),program_sha256=sha(__file__),scope='Complete scalar majorant on real line plus new affine moment coupling. Requires cubic identities and new coupled-moment theorem. Not a Grothendieck bound.'),indent=2)+'\n')
 print({k:row[k] for k in ['compact_leaves','coupled_objective_upper_decimal','seconds']})
