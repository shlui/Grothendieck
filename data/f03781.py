#!/usr/bin/env python3
"""Independent B=0 contact coefficients and all-source domain exclusion."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json
from mpfr_interval import I,BITS
ROOT=Path(__file__).resolve().parents[2];src=ROOT/'tonight/moments/cubic_head_contact/linear_contact_exclusion.json';d=json.loads(src.read_text())
assert d['status']=='PASS_EXACT_GLOBAL_LINEAR_CONTACT_PHYSICAL_SOURCE_EXCLUSION'
A,B,W=map(F,(d['parameters']['A'],d['parameters']['B'],d['parameters']['w']))
assert (A,B,W)==(F(147,50),F(0),F(11,2))
family=ROOT/'tonight/moments/cubic_head_contact/certify_direct_contact_majorant.py';assert hashlib.sha256(family.read_bytes()).hexdigest()==d['analytic_family_source_sha256']
accepted=json.loads((ROOT/'tonight/audit/cubic_head_contact_independent.json').read_text());assert accepted['status']=='PASS_INDEPENDENT_ALL_LINE_CUBIC_HEAD_CONTACT_GATE'
nu=(I(2)/I.pi()).sqrt();s=W*A*A;x=I(W).sqrt()*A;e=(I(W)/2).sqrt().erf()
J=x*e+nu*(-I(W)/2).exp()*A;Jp=e/(2*x);p0=nu*A;p2=nu/(2*A);delta=J-p0-p2*s;de=Jp-p2
p4=(3*delta-s*de)/(s*s);p6=(s*de-2*delta)/(s**3);assert p6.above(0)
row=d['contact_row'];assert F(row['contact_s'])==s
q=list(map(F,row['polynomial_coefficients_q0_q2_q4_q6']))
for p,Q in zip((p0,p2,p4,p6),q):assert p.below(Q)
assert q[3]>0
K=I(194)/I(10).sqrt();assert(K-F(633,15)).above(0)
ulo=F(d['u_lower']);assert ulo==F(151,200)
assert list(map(F,d['kappa_interval']))==[F(2),F(4)]and list(map(F,d['shape_interval']))==[F(0),F(1)]
seen=set();corners=[]
for c in d['moment_price_corners']:
 k,z=F(c['kappa']),F(c['shape']);seen.add((k,z))
 v=I(q[0]+q[1]+(3+k)*q[2])+I(q[3])*(15+15*k+(K-F(633,15)*z)*I(k)*I(k).sqrt())
 assert v.below(F(c['price_upper']))and F(c['source_upper_at_ulo'])==F(c['price_upper'])-A*ulo
 corners.append(F(c['price_upper']))
assert seen=={(k,z)for k in(F(2),F(4))for z in(F(0),F(1))}
price=max(corners);r=price-A*ulo;assert price==F(d['uniform_price_upper'])and r==F(d['source_r_uniform_upper'])
clean,excluded=F(d['clean_source_r_upper']),F(d['excluded_r_lower']);assert r<clean<excluded==F(129,500)
out={'status':'PASS_INDEPENDENT_COMPLETE_PHYSICAL_SOURCE_EXCLUSION','input_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'precision_bits':BITS,'parameters':{'A':str(A),'B':'0','w':str(W)},'u_lower':str(ulo),'kappa_interval':['2','4'],'all_native_shapes':['0','1'],'uniform_source_upper':str(r),'clean_source_upper':str(clean),'excluded_r_lower':str(excluded),'directed_coefficients_checked':4,'directed_corner_prices_checked':4,'producer_code_imported':False,'analytic_review':'The accepted positive-polynomial derivative proof includes B=0. The fixed contact Hermite polynomial is an all-line majorant; q6>=0 makes its source moment price convex in kappa and affine in shape. The exact physical gate A*u+r<=price and A>0 imply the complete exclusion for u>=.755, kappa[2,4], all shapes. No operator coefficient or law assumption is used.','scope':'Physical domain exclusion only; clipping remaining source cells and a complete global bound require a separate envelope replay.'}
path=Path(__file__).parent/'linear_contact_exclusion_independent.json';path.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
