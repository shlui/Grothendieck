#!/usr/bin/env python3
"""Bind the independently accepted eleven-bank polygon global theorem."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json
ROOT=Path(__file__).resolve().parents[3];HERE=Path(__file__).resolve().parent
paths={
 'accepted_previous_global':'tonight/audit/global_b7171/STRENGTHENED_GLOBAL_B7171_ACCEPTANCE.json',
 'extended_native_endpoints':'tonight/audit/global_b7171/extended_critical_native_endpoint_replay.json',
 'eleven_R3_corners':'tonight/audit/global_b7171/elevenbank_R3_corner_replay.json',
 'all_line_source_bank':'tonight/audit/linear_contact_bank_independent.json',
 'uniform_planes_and_inverse_caps':'tonight/audit/global_b7171/elevenbank_polygon_plane_prices.json',
 'selected_contact_credits':'tonight/audit/global_b7171/elevenbank_selected_contact_prices.json',
 'exact_complete_polygon_cover':'tonight/theory/polygon_envelope_audit.json',
 'first_transfers':'tonight/audit/native_b7171_union/first_four_transfers_independent.json',
 'lower_base':'tonight/audit/native_b7171_union/u600_682_independent_replay.json',
 'lower_R3_union':'tonight/audit/native_b7171_union/u682_6905_independent_replay.json',
 'upper_target':'tonight/audit/native_b7171_union/u775_c_independent_replay.json'}
sha=lambda p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
read=lambda p:json.loads((ROOT/p).read_text())
r={k:read(p)for k,p in paths.items()}
assert all(v['status'].startswith('PASS')for v in r.values())
old=r['accepted_previous_global'];assert old['historical_1773_or_1780498_premises_used'] is False
b=F(7171,10000);E=r['exact_complete_polygon_cover']
assert F(E['operator_b'])==b
assert E['original_cells']==58282==E['accepted_price_cells']+E['exact_physically_empty_cells']
assert E['accepted_price_cells']==57857 and E['exact_physically_empty_cells']==425
assert E['complete_groups']==294 and E['u_bands']==21 and E['kappa_bands']==14
assert E['recorded_complete_polygons']==57857 and E['native_knot_pieces']==14228
assert E['checked_edge_stationary_points']==3190 and E['uniform_source_planes']==140
source=E['source'];eh=sha(source);assert eh==E['source_sha256'];envelope=read(source)
assert F(envelope['operator_b'])==b and F(envelope['maximum_cell_cap'])==F(E['maximum_cell_cap'])
assert len(envelope['cells'])==57857 and len(envelope['physically_empty_cells'])==425
for key in('uniform_planes_and_inverse_caps','selected_contact_credits'):
 assert r[key]['input_sha256']==eh
R=r['eleven_R3_corners'];assert R['source']['sha256']==eh and R['profiles']==11 and R['fresh_all_corner_prices']==616
P=r['uniform_planes_and_inverse_caps'];assert P['uniform_polygon_planes']==140 and P['fresh_directed_corner_prices']==560 and P['inverse_caps']==250
C=r['selected_contact_credits'];assert C['selected_contact_cells']==26 and C['distinct_contact_corner_prices']==32
assert r['extended_native_endpoints']['endpoint_prices_recomputed']==224
oldenv=read(old['envelope']['path']);oldframes={v['mask']:v for v in oldenv['native_endpoint_inputs']}
extension=r['extended_native_endpoints']['frame'];newframes=envelope['native_endpoint_inputs']
assert len(newframes)==7 and len(oldframes)==6
for frame in newframes:
 assert frame==oldframes.get(frame['mask'],extension)
 assert sha(frame['mask'])==frame['mask_sha256']
assert extension in newframes
for key in('lower_base','lower_R3_union','upper_target'):
 assert F(r[key]['b'])==b and F(r[key]['maximum_rounded_cap'])<F(40,71)
assert r['first_transfers']['coverage']==['0','319/500']
assert list(map(F,r['lower_base']['interval']))==[F(3,5),F(341,500)]
assert list(map(F,r['lower_R3_union']['interval']))==[F(341,500),F(1381,2000)]
assert list(map(F,E['amplitude_interval']))==[F(1381,2000),F(31,40)]
assert list(map(F,r['upper_target']['interval']))==[F(31,40),F(159577,200000)]
maximum=F(E['maximum_cell_cap']);scale=10**12
strict=F((maximum.numerator*scale)//maximum.denominator+1,scale)
assert maximum<strict and F(40,71)<strict<F(old['norm_upper_bound'])
saved=F(envelope['strict_rational_cap']);assert maximum<strict<=saved
out=dict(status='PASS_INDEPENDENT_COMPLETE_POLYGON_GLOBAL_B7171_SCALAR_BOUND',
 operator_b=str(b),norm_upper_bound=str(saved),constructive_G_R_lower_bound=str(1/saved),
 constructive_G_R_lower_bound_display=float(1/saved),proved_middle_maximum=str(maximum),
 saved_envelope_strict_cap=str(saved),saved_envelope_reciprocal=str(1/saved),
 optional_twelve_decimal_strict_cap=str(strict),optional_tighter_reciprocal=str(1/strict),
 global_amplitude_coverage=['0','sqrt(2/pi)'],source_rho='847/1000',
 middle_price_cells=57857,proved_physically_empty_middle_cells=425,
 complete_middle_groups=294,independent_R3_corner_prices=616,
 uniform_physical_source_planes=140,independent_plane_corner_prices=560,
 independent_physical_inverse_caps=250,selected_contact_credit_cells=26,
 native_knot_pieces=14228,checked_edge_stationary_points=3190,
 envelope=dict(path=source,sha256=eh),
 bindings={k:dict(path=p,sha256=sha(p),status=r[k]['status'])for k,p in paths.items()},
 historical_1773_or_1780498_premises_used=False,
 scope='Complete fresh same-coefficient global norm bound and explicit matrix lower bound. The primary norm cap is the simple strict rational saved in the frozen envelope; the strict twelve-decimal upward successor is also recorded as an optional stronger cap. Fine39-kappa and robust-profile refinements are not used. Target1.775 and a three-significant-digit determination of G_R remain unproved.')
(HERE/'POLYGON_GLOBAL_B7171_ACCEPTANCE.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:v for k,v in out.items()if k!='bindings'},indent=2))
