#!/usr/bin/env python3
"""Independent Arb512 Taylor12 witnesses at three exact coefficients.

No producer modules are imported. Compact integrals and Gaussian tails
are recomputed directly; the producer's intervals are only overlap checks.
"""
from pathlib import Path
from fractions import Fraction as F
from flint import arb,arb_series,ctx
import hashlib,json,time
ROOT=Path(__file__).resolve().parents[2];HERE=Path(__file__).resolve().parent
SOURCE=ROOT/'tonight/selected_witness_coefficient_certificates.json'
ctx.prec=512;ORDER=12;CELLS=400;L=F(5);ctx.cap=ORDER+1
aa=lambda x:arb(F(x).numerator)/F(x).denominator
def jets(pos,A,C):
 z=arb_series([pos,1]);q=aa(A)*z+aa(C)*z**3
 e=(-(z*z+q*q)/2).exp();qp=aa(A)+3*aa(C)*z*z
 return e,qp*e,(z*z-1)*qp*e
def interval(v):return dict(lower=str(F(str(v.lower().fmpq()))),upper=str(F(str(v.upper().fmpq()))))
def compute(A,C,b):
 h=L/(2*CELLS);sums=[arb(0)]*3
 for k in range(CELLS):
  mid=F(2*k+1)*h;center=jets(aa(mid),A,C);whole=jets(aa(mid-h).union(aa(mid+h)),A,C)
  for j in range(3):
   sums[j]+=sum((center[j][i]*aa(2*h**(i+1)/F(i+1))for i in range(0,ORDER,2)),arb(0))
   sums[j]+=whole[j][ORDER]*aa(2*h**(ORDER+1)/F(ORDER+1))
 c=aa(1+A*A);decay=(-c*aa(L*L)/2).exp();J0=decay/(c*aa(L));J2=aa(L)*decay/c+J0/c;J4=aa(L**3)*decay/c+3*J2/c
 factors=[2/arb.pi(),2/arb.pi(),2/(arb.pi()*arb(6).sqrt())]
 tails=[factors[0]*J0,factors[1]*(aa(A)*J0+3*aa(C)*J2),factors[2]*(aa(A)*(J2+J0)+3*aa(C)*(J4+J2))]
 moments=[factors[j]*sums[j]+arb(0).union(tails[j].upper())for j in range(3)]
 my,m1,m3=moments;P1=m1*m1;P3=m3*m3-my*my;W=P1-aa(b)*P3
 return dict(parameters=dict(a=str(A),c=str(C),b=str(b)),moments={k:interval(v)for k,v in zip(('diffuse_cubic','first','ordinary_third'),moments)},tails=list(map(interval,tails)),P1_pairing=interval(P1),P3_pairing=interval(P3),scalar_norm_lower_witness=interval(W),reciprocal_norm_ceiling=interval(1/W),excludes_target1775_for_this_coefficient=W>aa(F(40,71)),display=dict(witness=float(W),ceiling=float(1/W))),moments,W
def main():
 start=time.time();digest=hashlib.sha256(SOURCE.read_bytes()).hexdigest();raw=json.loads(SOURCE.read_text());expected=[(F(31091537,10**7),F(49105809,10**8),F(141,200)),(F(277966397,10**8),F(1353131,3125000),F(713,1000)),(F(264126401,10**8),F(8168077,2*10**7),F(7171,10000))];assert len(raw)==len(expected)==3
 rows=[]
 for item,(A,C,b)in zip(raw,expected):
  assert item['status']=='PASS_DIRECTED_CLT_SCALAR_NORM_LOWER_WITNESS'
  assert {k:F(v)for k,v in item['parameters'].items()}==dict(a=A,c=C,epsilon=b)
  r,moments,W=compute(A,C,b)
  for key,new in zip(('diffuse_cubic','first','ordinary_third'),moments):
   q=item['moments'][key];assert new.overlaps(aa(F(q['lower'])).union(aa(F(q['upper']))))
  for key,new in [('scalar_norm_lower_witness',W),('reciprocal_norm_ceiling',1/W)]:
   q=item[key];assert new.overlaps(aa(F(q['lower'])).union(aa(F(q['upper']))))
  if b==F(141,200):assert W>aa(F(40,71)) and 1/W<aa(F(1774988761828,10**12))
  else:assert W<aa(F(40,71))
  rows.append(r)
 assert hashlib.sha256(SOURCE.read_bytes()).hexdigest()==digest
 out=dict(status='PASS_INDEPENDENT_ARB_SELECTED_CUBIC_COEFFICIENT_WITNESSES',bits=ctx.prec,Taylor_order=ORDER,cells=CELLS,domain=['0','5'],source=dict(path=str(SOURCE.relative_to(ROOT)),sha256=digest),producer_modules_imported=False,rows=rows,seconds=time.time()-start,analytic_review='The reflected-sign CLT witness uses the included original-degree projection proof. q_prime>=a>0 gives its uniform small-strip control. Here all three univariate compact integrals and rigorous Gaussian tails are freshly recomputed. For fixed b the limiting scalar pairing is m1^2-b*(m3^2-my^2), so it is a norm LOWER witness; its reciprocal is a construction-family CEILING.',scope='At b=.705 the pure-cubic reciprocal-norm sequence cannot reach1.775. This neither bounds G_R from above nor excludes1.775 for another coefficient or matrix family. No Grothendieck lower bound is inferred from these witnesses.')
 (HERE/'selected_cubic_witnesses_independent_replay.json').write_text(json.dumps(out,indent=2)+'\n')
 print(json.dumps(dict(status=out['status'],rows=[dict(b=r['parameters']['b'],**r['display'],excludes=r['excludes_target1775_for_this_coefficient'])for r in rows],seconds=out['seconds']),indent=2))
if __name__=='__main__':main()
