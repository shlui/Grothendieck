#!/usr/bin/env python3
"""Independent binding of the strengthened nine-profile/contact cover."""
from pathlib import Path
from fractions import Fraction as F
import json,hashlib
ROOT=Path(__file__).resolve().parents[3];HERE=Path(__file__).resolve().parent
paths={
 'first_transfers':'tonight/audit/native_b7171_union/first_four_transfers_independent.json',
 'lower_base':'tonight/audit/native_b7171_union/u600_682_independent_replay.json',
 'lower_R3_union':'tonight/audit/native_b7171_union/u682_6905_independent_replay.json',
 'lower_R3_price':'tonight/audit/diffuse_b7171_R3/independent_extended_remainder_box.json',
 'upper_target':'tonight/audit/native_b7171_union/u775_c_independent_replay.json',
 'middle_native_lower':'tonight/audit/global_b7171/native_middle_lower_input_replay.json',
 'middle_native_upper':'tonight/audit/global_b7171/native_middle_upper_input_replay.json',
 'critical_native_endpoints':'tonight/audit/global_b7171/critical_native_endpoint_replay.json',
 'nine_R3_corners':'tonight/audit/global_b7171/ninebank_R3_corner_replay.json',
 'contact_global_majorant':'tonight/audit/cubic_head_contact_independent.json',
 'linear_exclusion':'tonight/audit/linear_contact_exclusion_independent.json',
 'inverse_bounds':'tonight/audit/global_b7171/critical_ninebank_inverse_bounds.json',
 'selected_contact_prices':'tonight/audit/global_b7171/critical_ninebank_contact_prices.json',
 'exact_complete_clipped_cover':'tonight/theory/clipped_contact_envelope_audit.json',
 'cubic_L1_prices':'tonight/audit/refit847/independent_complete_price_cover.json'}
def read(p):return json.loads((ROOT/p).read_text())
def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
r={k:read(p)for k,p in paths.items()};assert all(x['status'].startswith('PASS')for x in r.values())
b=F(7171,10000);B=F(56513883,10**8);E=r['exact_complete_clipped_cover']
assert F(E['operator_b'])==b and F(E['strict_rational_cap'])==B
assert E['original_cells']==39662==E['accepted_price_cells']+E['exact_physically_empty_cells']
assert E['complete_groups']==224 and E['u_bands']==16 and E['kappa_bands']==14
assert F(E['maximum_cell_cap'])<B and F(40,71)<B
source=E['source'];eh=sha(source);assert eh==E['source_sha256']
for k in['inverse_bounds','selected_contact_prices']:assert r[k]['input_sha256']==eh
assert r['critical_native_endpoints']['envelope_input_sha256']==eh
assert r['nine_R3_corners']['source']['sha256']==eh
assert r['nine_R3_corners']['profiles']==9 and r['nine_R3_corners']['directed_corner_prices']==504
assert r['critical_native_endpoints']['endpoint_prices_recomputed']==672
assert r['selected_contact_prices']['selected_contact_cells']==226
assert r['inverse_bounds']['bounds_checked']==238
for k in['lower_base','lower_R3_union','upper_target']:assert F(r[k]['maximum_rounded_cap'])<F(40,71)and F(r[k]['b'])==b
assert r['first_transfers']['coverage']==['0','319/500']
assert list(map(F,r['lower_base']['interval']))==[F(3,5),F(341,500)]
assert list(map(F,r['lower_R3_union']['interval']))==[F(341,500),F(1381,2000)]
assert list(map(F,E['amplitude_interval']))==[F(1381,2000),F(31,40)]
assert list(map(F,r['upper_target']['interval']))==[F(31,40),F(159577,200000)]
out={'status':'PASS_INDEPENDENT_STRENGTHENED_COMPLETE_GLOBAL_B7171_SCALAR_BOUND','operator_b':str(b),'norm_upper_bound':str(B),'constructive_G_R_lower_bound':str(1/B),'constructive_G_R_lower_bound_display':float(1/B),'global_amplitude_coverage':['0','sqrt(2/pi)'],'source_rho':'847/1000','middle_price_cells':39657,'proved_physically_empty_middle_cells':5,'independent_native_endpoint_prices':9120,'independent_R3_corner_prices':504,'independent_physical_inverse_caps':238,'selected_contact_credit_cells':226,'envelope':{'path':source,'sha256':eh},'bindings':{k:{'path':p,'sha256':sha(p),'status':r[k]['status']}for k,p in paths.items()},'historical_1773_or_1780498_premises_used':False,'scope':'Fresh same-coefficient global norm and explicit matrix lower bound, independent of inherited global premises. The target1.775 and three-significant-digit value of G_R remain unproved. Newer middle-profile or joint-credit components are not used.'}
(HERE/'STRENGTHENED_GLOBAL_B7171_ACCEPTANCE.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items()if k!='bindings'},indent=2))
