#!/usr/bin/env python3
"""Independent directed contact price/normalization bindings on selected cells."""
from pathlib import Path
from fractions import Fraction as F
from functools import lru_cache
import argparse,hashlib,json,sys,gzip
sys.path.insert(0,str(Path(__file__).resolve().parent.parent))
from mpfr_interval import I,BITS
ROOT=Path(__file__).resolve().parents[3]
ap=argparse.ArgumentParser();ap.add_argument('ledger',type=Path);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--r3-report',type=Path);args=ap.parse_args();path=args.ledger.resolve();d=json.loads(gzip.decompress(path.read_bytes()) if path.suffix=='.gz' else path.read_text());b=F(d['operator_b']);assert 0<b<=F(7171,10000)
contact=ROOT/d['contact_gate_binding']['path'];assert hashlib.sha256(contact.read_bytes()).hexdigest()==d['contact_gate_binding']['sha256']
C=json.loads(contact.read_text());a=F(C['parameters']['a']);B=F(C['parameters']['B']);d2=6*B*B
accepted=json.loads((ROOT/'tonight/audit/cubic_head_contact_independent.json').read_text());assert accepted['input_sha256']==d['contact_gate_binding']['sha256']and accepted['status'].startswith('PASS')
K=I(194)/I(10).sqrt();native={v['mask']:v for v in d['native_endpoint_inputs']}
r3=json.loads((ROOT/'tonight/audit/global_b7171/ninebank_R3_corner_replay.json').read_text())['complete_whole_q_bindings']
r3.update(json.loads((ROOT/'tonight/audit/global_b7171/elevenbank_R3_corner_replay.json').read_text())['complete_new_whole_q_bindings'])
if args.r3_report:
 report=json.loads(args.r3_report.read_text());assert report['status'].startswith('PASS')
 assert report['source']['sha256']==hashlib.sha256(path.read_bytes()).hexdigest()
 r3.update(report['complete_whole_q_bindings'])
@lru_cache(None)
def point(index,k,z):
 q=list(map(F,C['rows'][index]['polynomial_coefficients_q0_q2_q4_q6']));assert q[3]>=0
 return I(q[0]+q[1]+(3+k)*q[2])+I(q[3])*(15+15*k+(K-F(633,15)*z)*I(k)*I(k).sqrt())
seen=set();count=0;inactive=0;nativebeta={};rawbeta={}
for cell in d['cells']:
 if 'contact_credit'not in cell:continue
 if not cell['winner'].startswith(('native-contact:','R3-contact:')):
  inactive+=1;continue
 cc=cell['contact_credit'];index=cc['contact_row'];assert cc['contact_w']==C['rows'][index]['w'];upper=F(cc['gate_price_upper']);klo,khi=map(F,cell['kappa_interval'])
 for k in(klo,khi):
  for z in(F(0),F(1)):
   assert point(index,k,z).below(upper);seen.add((index,k,z))
 winner=cell['winner'];beta_upper=F(cc['beta_d_over_b_upper'])
 if winner.startswith('native-contact:'):
  name=winner[len('native-contact:'):];binding=native[name];local=ROOT/binding['package']/'data/local_summary.json'
  if name not in nativebeta:nativebeta[name]=F(json.loads(local.read_text())['profile']['beta'])
  value=I(nativebeta[name])*I(6).sqrt()*B/b;assert value.below(beta_upper)
 else:
  assert winner.startswith('R3-contact:');name=winner[len('R3-contact:'):]
  if name not in rawbeta:
   prof=ROOT/r3[name]['profile']['path'];assert hashlib.sha256(prof.read_bytes()).hexdigest()==r3[name]['profile']['sha256'];rawbeta[name]=F(json.loads(prof.read_text())['profile']['beta'])
  assert beta_upper==6*rawbeta[name]*B/b
 ulo=F(cell['u_interval'][0]);rlo=F(cell['r_interval'][0]);S=upper-a*ulo-rlo+beta_upper
 assert S==F(cc['S_upper'])<0
 lam=F(cc['lambda_value']);assert lam>=0
 # This identity is an implementation check only; validity holds for EVERY fixed lambda>=0.
 star=-2*b*S/d2;scale=10**10;floor=F(star.numerator*scale//star.denominator,scale);assert lam==floor
 count+=1
out={'status':'PASS_INDEPENDENT_SELECTED_CONTACT_CREDIT_NUMERIC_BINDINGS','input_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'precision_bits':BITS,'selected_contact_cells':count,'inactive_contact_metadata_cells':inactive,'distinct_contact_corner_prices':len(seen),'native_normalized_beta_terms_checked':len(nativebeta),'R3_raw_beta_terms_checked':len(rawbeta),'producer_code_imported':False,'analytic_review':'Every contact row and lambda is fixed over its whole cell. All four coupled cumulant/shape gate prices are bounded; normalized native beta*d/b uses directed sqrt6, raw R3 beta*d/b equals6*beta_raw*B/b. The recorded S has PLUS beta*d/b. The exact cell outer/chord/quadratic pricing is separately audited.','scope':'Numeric contact-credit premises on all selected cells; complete physical partition and all options are separate.'}
parentpath=ROOT/'tonight/audit/global_b7171/full_upper57_b713_22_contact_credits.json'
parent=json.loads(parentpath.read_text());assert parent['status']=='PASS_INDEPENDENT_SELECTED_CONTACT_CREDIT_NUMERIC_BINDINGS'
assert F(d['operator_b'])==F(713,1000)
out['accepted_parent_contact_premise']={'path':str(parentpath.relative_to(ROOT)),'sha256':hashlib.sha256(parentpath.read_bytes()).hexdigest(),'input_sha256':parent['input_sha256']}
out['inactive_metadata_review']='Only active native-contact/R3-contact winners are repriced here. Contact fields on inherited or newly replaced plain-profile winners are inactive. Exact parent/child geometry proves their chosen bound without using those fields; inherited parent contact bounds are additionally bound to the accepted parent numeric proof.'
args.output.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
