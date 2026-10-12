"""Independent MPFR768 replay of all physical-curved B0 source planes."""
import os,sys
sys.dont_write_bytecode=True
os.environ['GR_MPFR_BITS']='768'
from pathlib import Path
from fractions import Fraction as F
from itertools import product
from functools import lru_cache
import argparse,json,gzip,hashlib
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
BASE=ROOT/'verified_2026-10-01/tonight/operator/preferred_global'
sys.path.insert(0,str(BASE/'tonight/audit'))
from mpfr_interval import I,UP,DOWN,BITS,selftest
assert BITS==768
ap=argparse.ArgumentParser();ap.add_argument('source',type=Path);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
assert args.output.resolve().is_relative_to(HERE)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):
 with(gzip.open(p,'rt')if p.suffix=='.gz'else p.open())as h:return json.load(h)
def rec(p):return {'path':str(p.relative_to(ROOT)),'sha256':sha(p)}
source=args.source.resolve();initial=sha(source);d=read(source)
assert d['physical_constraint']=='4*kappa*zeta-kappa<=54'
assert d['physical_kappa_shape_constraint']is True and F(d['operator_b'])==F(711,1000)
bankpath=BASE/'tonight/moments/linear_contact_bank/source_bank.json';bank=read(bankpath)
acceptpath=BASE/'tonight/audit/linear_contact_bank_independent.json';accepted=read(acceptpath)
assert accepted['status']=='PASS_INDEPENDENT_COMPLETE_LINEAR_CONTACT_SOURCE_BANK'and accepted['input_sha256']==sha(bankpath)
gatepath=BASE/'tonight/audit/cubic_head_contact_independent.json';gate=read(gatepath)
assert gate['status']=='PASS_INDEPENDENT_ALL_LINE_CUBIC_HEAD_CONTACT_GATE'
family=BASE/'tonight/moments/cubic_head_contact/certify_direct_contact_majorant.py'
assert bank['analytic_family_sha256']==sha(family)
A=F(bank['base_A']);assert A==F(147,50)and F(bank['base_B'])==0
nu=(I(2)/I.pi()).sqrt();K=I(194)/I(10).sqrt();assert (K-F(633,15)).above(0)
base=[];basechecks=[]
for item in bank['base_rows']:
 W=F(item['w']);s=W*A*A;x=I(W).sqrt()*A;e=(I(W)/2).sqrt().erf()
 J=x*e+nu*(-I(W)/2).exp()*A;Jp=e/(2*x);p0=nu*A;p2=nu/(2*A)
 delta=J-p0-p2*s;de=Jp-p2;coeffs=[p0,p2,(3*delta-s*de)/(s*s),(s*de-2*delta)/(s**3)]
 q=list(map(F,item['polynomial_coefficients_q0_q2_q4_q6']));assert F(item['contact_s'])==s and coeffs[3].above(0)and q[3]>0
 for value,Q in zip(coeffs,q):assert value.below(Q)
 base.append(q);basechecks.append({'w':str(W),'contact_s':str(s),'directed_coefficients':[v.record()for v in coeffs],'accepted_majorant':list(map(str,q))})
assert len(base)==6
rows=[]
for item in bank['rows']:
 a,scale=F(item['a']),F(item['scale']);i=item['base_row'];assert scale==a/A>0
 q=[Q*scale**(1-2*j)for j,Q in enumerate(base[i])]
 assert q==list(map(F,item['polynomial_coefficients_q0_q2_q4_q6']))and q[3]>0
 assert F(item['w'])==F(bank['base_rows'][i]['w'])and F(item['contact_s'])==F(item['w'])*a*a
 rows.append((a,q))
assert len(rows)==60
def candidates(ka,kb,za,zb):
 assert 0<=ka<kb<=90 and F(2,5)<=za<zb<=1
 points={(k,z)for k,z in product((ka,kb),(za,zb))if k*(4*z-1)<=54}
 for k in(ka,kb):
  if k>0 and za<=(k+54)/(4*k)<=zb:points.add((k,(k+54)/(4*k)))
 for z in(za,zb):
  if ka<=F(54)/(4*z-1)<=kb:points.add((F(54)/(4*z-1),z))
 assert points
 return sorted(points)
@lru_cache(None)
def price(i,k,z):
 q=rows[i][1];sixth=(K-F(633,9)*(z-F(2,5)))*I(k)*I(k).sqrt()
 return I(q[0]+q[1]+(3+k)*q[2]+(15+15*k)*q[3])+I(q[3])*sixth
tables=[];seen=set();selected_evaluations=0;all_evaluations=0
for table in d['source_plane_cache']:
 ka,kb=map(F,table['kappa_interval']);za,zb=map(F,table['shape_interval']);key=(ka,kb,za,zb)
 assert key not in seen;seen.add(key)
 assert table['physical_constraint']=='4*kappa*zeta-kappa<=54'
 ends=candidates(*key);assert table['evaluation_points']==[list(map(str,x))for x in ends]
 # Evaluate every accepted fixed row, independently of producer selections.
 for i in range(60):
  for k,z in ends:price(i,k,z);all_evaluations+=1
 planes=[];slopes=set()
 for plane in table['planes']:
  i=plane['source_row'];a,q=rows[i];assert F(plane['a'])==a and a not in slopes;slopes.add(a)
  saved=list(map(F,plane['corner_uppers']));cap=F(plane['cap']);assert len(saved)==len(ends)and max(saved)==cap
  values=[]
  for(k,z),upper in zip(ends,saved):
   value=price(i,k,z);assert value.below(upper)
   values.append({'kappa':str(k),'shape':str(z),'lower':value.lo.text(DOWN,90),'upper':value.hi.text(UP,90),'saved_upper':str(upper)})
   selected_evaluations+=1
  planes.append({'a':str(a),'source_row':i,'cap':str(cap),'coefficients':list(map(str,q)),'corners':values})
 assert slopes=={a for a,q in rows}
 tables.append({'kappa_interval':[str(ka),str(kb)],'shape_interval':[str(za),str(zb)],'planes':planes})
cacheboxes={tuple(map(F,x['kappa_interval']+x['shape_interval']))for x in d['cache']}
assert seen==cacheboxes and len(tables)==70 and sha(source)==initial
out={'status':'PASS_INDEPENDENT_CURVED_B0_UNIFORM_SOURCE_PLANES','source':rec(source),'input_sha256':initial,'precision_bits':BITS,'operator_b':d['operator_b'],
 'source_bank':rec(bankpath),'accepted_all_line_bank_report':rec(acceptpath),'accepted_analytic_contact_gate':rec(gatepath),'analytic_family':rec(family),
 'base_coefficients_freshly_verified':24,'exact_homogeneous_rows':60,'uniform_source_boxes':len(tables),'selected_uniform_planes':sum(len(t['planes'])for t in tables),
 'selected_directed_prices':selected_evaluations,'all60_row_directed_prices':all_evaluations,'distinct_MPFR_prices':price.cache_info().currsize,'base_coefficient_checks':basechecks,'tables':tables,
 'boundary_proof':rec(ROOT/'new_round_2026-10-02/audit/LINEAR_TENSOR_SHAPE_CONSTRAINT.md'),'physical_tensor_proof':rec(ROOT/'new_round_2026-10-02/audit/TENSOR_SHAPE_NECESSARY_CONSTRAINT.md'),
 'backend':selftest(),'script_sha256':sha(Path(__file__)),'producer_modules_imported':False,
 'scope':'All70 declared necessary physical kappa/shape boxes and700 selected fixed B0 planes. Every selected cap bounds one already accepted global source row over the entire curved physical subdomain. Exact homogeneity and fresh Gaussian contact coefficients are verified; global clipping, rank minimum, partition and norm composition remain separate.'}
args.output.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:out[k]for k in ('status','uniform_source_boxes','selected_uniform_planes','selected_directed_prices','all60_row_directed_prices')},indent=2))
