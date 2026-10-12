#!/usr/bin/env python3
"""Independent global contacts, homogeneity, prices, and frozen box caps."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json
from mpfr_interval import I,BITS
ROOT=Path(__file__).resolve().parents[2];src=ROOT/'tonight/moments/linear_contact_bank/source_bank.json';d=json.loads(src.read_text());assert d['status']=='PASS_EXACT_CONTACT_PHYSICAL_SOURCE_BANK_AND_FROZEN_BOX_EXCLUSIONS'
A=F(d['base_A']);assert A==F(147,50)and F(d['base_B'])==0
assert hashlib.sha256((ROOT/'tonight/moments/cubic_head_contact/certify_direct_contact_majorant.py').read_bytes()).hexdigest()==d['analytic_family_sha256']
proof=json.loads((ROOT/'tonight/audit/cubic_head_contact_independent.json').read_text());assert proof['status']=='PASS_INDEPENDENT_ALL_LINE_CUBIC_HEAD_CONTACT_GATE'
nu=(I(2)/I.pi()).sqrt();K=I(194)/I(10).sqrt();assert(K-F(633,15)).above(0)
base=[]
for row in d['base_rows']:
 W=F(row['w']);s=W*A*A;x=I(W).sqrt()*A;e=(I(W)/2).sqrt().erf();J=x*e+nu*(-I(W)/2).exp()*A;Jp=e/(2*x)
 p0=nu*A;p2=nu/(2*A);delta=J-p0-p2*s;de=Jp-p2;p4=(3*delta-s*de)/(s*s);p6=(s*de-2*delta)/(s**3)
 q=list(map(F,row['polynomial_coefficients_q0_q2_q4_q6']));assert p6.above(0)
 for p,Q in zip((p0,p2,p4,p6),q):assert p.below(Q)
 assert F(row['contact_s'])==s;assert q[3]>0;base.append(q)
assert len(base)==6
rows=[]
for r in d['rows']:
 a,scale=F(r['a']),F(r['scale']);i=r['base_row'];assert scale==a/A>0
 qs=[q*scale**(1-2*j)for j,q in enumerate(base[i])];assert qs==list(map(F,r['polynomial_coefficients_q0_q2_q4_q6']))
 assert F(r['w'])==F(d['base_rows'][i]['w'])and F(r['contact_s'])==F(r['w'])*a*a
 rows.append((a,qs))
assert len(rows)==60
cache={}
for p in d['all_endpoint_prices']:
 i,k,z=p['row'],F(p['kappa']),F(p['shape']);assert 0<=k<=90 and z in(F(0),F(1));q=rows[i][1]
 v=I(q[0]+q[1]+(3+k)*q[2])+I(q[3])*(15+15*k+(K-F(633,15)*z)*I(k)*I(k).sqrt())
 upper=F(p['upper']);assert v.below(upper);key=(i,k,z);assert key not in cache;cache[key]=upper
assert len(cache)==1800
bindings={x['path']:x['sha256']for x in d['source_ledger_bindings']};boxes=0;reports=[]
for ledger in d['ledgers']:
 path=ROOT/ledger['source'];assert hashlib.sha256(path.read_bytes()).hexdigest()==bindings[ledger['source']]
 source=json.loads(path.read_text());assert len(source['unclosed_boxes'])==len(ledger['boxes'])==ledger['box_count']
 reduced=excluded=0
 for raw,b in zip(source['unclosed_boxes'],ledger['boxes']):
  for k,v in raw.items():assert b[k]==v
  ulo,uhi=map(F,b['u_interval']);rlo,rhi=map(F,b['r_interval']);klo,khi=map(F,b['kappa_interval']);candidates=[]
  for i,(a,_)in enumerate(rows):
   corners=[cache[i,k,z]for k in(klo,khi)for z in(F(0),F(1))];candidates.append((max(corners)-a*ulo,i,corners))
  cap,i,corners=min(candidates,key=lambda x:x[0]);assert i==b['selected_row']and cap==F(b['uniform_source_upper'])and corners==list(map(F,b['moment_price_corners']))
  assert b['excluded']==(cap<rlo) and F(b['clipped_r_upper'])==min(rhi,cap) and F(b['source_clip_gain'])==max(F(0),rhi-cap)
  reduced+=min(rhi,cap)<rhi;excluded+=cap<rlo;boxes+=1
 assert reduced==ledger['source_reduced_count']and excluded==ledger['excluded_count']
 reports.append({'source':ledger['source'],'boxes':ledger['box_count'],'reduced':reduced,'excluded':excluded})
assert boxes==227
out={'status':'PASS_INDEPENDENT_COMPLETE_LINEAR_CONTACT_SOURCE_BANK','input_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'precision_bits':BITS,'base_rows':6,'base_coefficients_checked':24,'exact_scaled_rows':60,'directed_endpoint_prices':1800,'frozen_boxes_checked':boxes,'ledger_reports':reports,'producer_code_imported':False,'analytic_review':'Positive rational homogeneity preserves the accepted all-line majorants, including negative q4. Every selected row is fixed over its entire source box. q6>0 makes price convex in cumulant and affine in shape; subtracting positive a*u makes u_lo its uniform upper source cap. Only strict cap<r_lo yields a complete physical exclusion.','scope':'Additional all-source physical inequalities and frozen-box clipping only, not a complete global norm certificate.'}
p=Path(__file__).parent/'linear_contact_bank_independent.json';p.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
