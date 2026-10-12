#!/usr/bin/env python3
"""Bind independently accepted arithmetic/gates/partitions for ONE operator."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json
ROOT=Path(__file__).resolve().parents[3]
D=ROOT/'tonight/audit/global_b7171'
def read(rel):return json.loads((ROOT/rel).read_text())
def sha(rel):return hashlib.sha256((ROOT/rel).read_bytes()).hexdigest()
paths={
 'first_transfers':'tonight/audit/native_b7171_union/first_four_transfers_independent.json',
 'lower_base':'tonight/audit/native_b7171_union/u600_682_independent_replay.json',
 'lower_R3_union':'tonight/audit/native_b7171_union/u682_6905_independent_replay.json',
 'lower_R3_price':'tonight/audit/diffuse_b7171_R3/independent_extended_remainder_box.json',
 'upper_target':'tonight/audit/native_b7171_union/u775_c_independent_replay.json',
 'middle_native_lower':'tonight/audit/global_b7171/native_middle_lower_input_replay.json',
 'middle_native_upper':'tonight/audit/global_b7171/native_middle_upper_input_replay.json',
 'middle_R3_corners':'tonight/audit/global_b7171/R3_corner_replay.json',
 'middle_exact_partition':'tonight/theory/global_native_R3_envelope_algebra_audit.json',
 'cubic_L1_prices':'tonight/audit/refit847/independent_complete_price_cover.json'}
records={k:read(p)for k,p in paths.items()}
assert all(x['status'].startswith('PASS')for x in records.values())
b=F(7171,10000);B=F(70756627,125000000);small=F(40,71)
for key in['lower_base','lower_R3_union','upper_target','middle_native_lower','middle_native_upper']:
 assert F(records[key]['b'])==b
assert F(records['first_transfers']['operator_b'])==b
assert F(records['middle_R3_corners']['operator_b'])==b
assert F(records['middle_exact_partition']['operator_b'])==b
assert records['middle_R3_corners']['all_used_R3_profiles_whole_q_bound']
assert records['middle_R3_corners']['profiles']==7 and records['middle_R3_corners']['directed_corner_prices']==392
assert records['middle_exact_partition']['cell_count']==39662
assert F(records['middle_exact_partition']['strict_rational_cap'])==B
assert F(records['middle_exact_partition']['maximum_cell_cap'])<B
assert small<B
for key in ['lower_base','lower_R3_union','upper_target']:assert F(records[key]['maximum_rounded_cap'])<small
assert records['first_transfers']['coverage']==['0','319/500']
assert list(map(F,records['lower_base']['interval']))==[F(3,5),F(341,500)]
assert list(map(F,records['lower_R3_union']['interval']))==[F(341,500),F(1381,2000)]
assert list(map(F,records['middle_exact_partition']['amplitude_interval']))==[F(1381,2000),F(31,40)]
assert list(map(F,records['upper_target']['interval']))==[F(31,40),F(159577,200000)]
envelope=ROOT/'tonight/recovery/profile_cover/global_middle_b7171_native_R3_envelope.json'
eh=hashlib.sha256(envelope.read_bytes()).hexdigest()
assert records['middle_R3_corners']['source_sha256']==eh==records['middle_exact_partition']['envelope_sha256']
for name in['middle_native_lower','middle_native_upper']:
 assert records[name]['missing_native_price_cells']==0
 assert records[name]['input_sha256']in[x['sha256']for x in records['middle_exact_partition']['native_inputs']]
# This fixed rational upper endpoint .797885 exceeds sqrt(2/pi); verified below.
import sys
sys.path.insert(0,str(ROOT/'tonight/audit'))
from mpfr_interval import I,BITS
nu=(I(2)/I.pi()).sqrt();assert nu.below(F(159577,200000))
result={'status':'PASS_INDEPENDENT_COMPLETE_GLOBAL_B7171_SCALAR_BOUND','operator_b':str(b),
 'norm_upper_bound':str(B),'strictly_smaller_cell_cap':records['middle_exact_partition']['maximum_cell_cap'],
 'constructive_G_R_lower_bound':str(1/B),'constructive_G_R_lower_bound_display':float(1/B),
 'global_amplitude_coverage':['0','sqrt(2/pi)'],'global_source_rho':'847/1000',
 'scope':'Fresh global bound for P1-(7171/10000)P3 and its explicit Gaussian-cell/reduced/rational matrix sequence. This does not prove G_R>=1.775, the retained1.773 premise, the retained upper1.780498 premise, or a three-digit value of G_R.',
 'historical_1773_or_1780498_premises_used':False,'new_contact_credit_used':False,
 'independent_endpoint_prices':sum(records[k]['selected_endpoint_prices_recomputed']for k in['lower_base','lower_R3_union','upper_target','middle_native_lower','middle_native_upper'])+records['first_transfers']['endpoint_prices_recomputed'],
 'independent_R3_corner_prices':392,'global_middle_cells':39662,
 'bindings':{k:{'path':p,'sha256':sha(p),'status':records[k]['status']}for k,p in paths.items()},
 'frozen_envelope':{'path':str(envelope.relative_to(ROOT)),'sha256':eh}}
(D/'GLOBAL_B7171_ACCEPTANCE.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items()if k!='bindings'},indent=2))
