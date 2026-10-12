#!/usr/bin/env python3
"""Exact consequences of frozen aggregate caps; pure-P3 replacement stays conditional."""
from fractions import Fraction as F
from math import isqrt
from pathlib import Path
import json,hashlib
HERE=Path(__file__).resolve().parent
A=F(2403891349687,2000000000000)
B=F(12706580435619,10000000000000)
E=F(699,1000)
U=F(1000,1773)
TARGET=F(40,71)
RHO=F(437,500)
ETA_UP=(1-A*U)/(B-A*E)
ETA_DOWN=(1-A*U)/(B+A*E)
assert ETA_UP>0 and ETA_UP<1/B
assert ETA_DOWN>0 and ETA_DOWN<1/B
def floor_sqrt(q,den=10**9):
 if q<0:return None
 k=isqrt(q.numerator*den**2//q.denominator)
 out=F(k,den)
 assert out*out<=q
 assert F(k+1,den)**2>q
 return out
perturb=[]
for eps in [F(698999,1000000),F(699001,1000000),F(7,10),F(707,1000),F(18,25)]:
 if eps>=E:
  delta=eps-E;eta=ETA_UP
  assert delta<=(B-A*E)/A
 else:
  delta=E-eps;eta=ETA_DOWN
  assert delta<=(B+A*E)/A
 cap=U+delta*eta
 assert cap>=U
 perturb.append({'epsilon':str(eps),'norm_upper_from_baseline_and_energy':str(cap),
  'norm_upper_decimal':float(cap),'reciprocal_lower':str(1/cap),
  'reciprocal_decimal':float(1/cap),'improves_baseline':cap<U})
conditional=[]
for eps in [F(699,1000),F(707,1000),F(18,25)]:
 pure_cut2=TARGET-eps*RHO**2
 mixed_cut2=(TARGET-eps/B)/(1-eps*A/B)
 pure=floor_sqrt(pure_cut2);mixed=floor_sqrt(mixed_cut2)
 best=max(x for x in [pure,mixed] if x is not None)
 conditional.append({'epsilon':str(eps),'target':str(TARGET),
  'pure_cutoff_squared':str(pure_cut2),'pure_safe_rational_right':str(pure),
  'mixed_cutoff_squared':str(mixed_cut2),'mixed_safe_rational_right':str(mixed),
  'combined_safe_right':str(best),'combined_safe_right_decimal':float(best)})
# Necessary condition for the pure estimate alone to cover the old first three intervals.
old_energy_right=F(414137,2048000)
assert float(old_energy_right)==0.20221533203125
rho_needed2=(TARGET-old_energy_right**2)/F(18,25)
out={'status':'PASS_EXACT_AGGREGATE_PERTURBATION_ARITHMETIC',
 'frozen_accepted_inputs':{'baseline_norm_cap':str(U),'baseline_epsilon':str(E),
  'mixed_energy_A':str(A),'mixed_energy_B':str(B)},
 'perturbation_penalty_positive_epsilon':str(ETA_UP),
 'perturbation_penalty_positive_epsilon_decimal':float(ETA_UP),
 'perturbation_penalty_negative_epsilon':str(ETA_DOWN),
 'perturbation_penalty_negative_epsilon_decimal':float(ETA_DOWN),
 'triangle_penalty_from_same_energy':str(1/B),'rows':perturb,
 'aggregate_allowed_saturating_point_positive':{'E1_h':str(U-E*ETA_UP),
  'E1_k':'0','E3_h':'0','E3_k':str(ETA_UP),
  'warning':'Feasible for these necessary aggregate constraints only; not claimed to be realizable by signs.'},
 'new_pure_P3_section_status':'CONDITIONAL_UNACCEPTED_DISCOVERY_INPUT',
 'assumed_uniform_cubic_L1_cap':str(RHO),
 'pure_mixed_switch_amplitude_squared':str((1-B*RHO**2)/A),
 'pure_mixed_switch_safe_rational':str(floor_sqrt((1-B*RHO**2)/A)),
 'conditional_low_u_rows':conditional,
 'rho_squared_threshold_to_improve_aggregate_positive_penalty':str(ETA_UP),
 'rho_threshold_to_improve_aggregate_positive_penalty_decimal':float(ETA_UP)**.5,
 'rho_squared_needed_for_pure_low_u_old_cover_at_720':str(rho_needed2),
 'rho_needed_for_pure_low_u_old_cover_at_720_decimal':float(rho_needed2)**.5,
 'norm_decrease_needed_for_1_77301':str(U-F(100000,177301)),
 'scope':'Rigorous perturbation upper bounds from frozen accepted premises; they yield no endpoint improvement. New pure-P3 consequences are conditional until its global certificate passes.',
 'producer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
(HERE/'aggregate_perturbation_audit.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'status':out['status'],'positive_penalty':float(ETA_UP),
 'negative_penalty':float(ETA_DOWN),'conditional_low_u_rows':conditional,
 'needed_rho_for_old_low_u_720':out['rho_needed_for_pure_low_u_old_cover_at_720_decimal']},indent=2))
