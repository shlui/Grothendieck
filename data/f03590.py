#!/usr/bin/env python3
"""Complete preferred b=.713 theorem from fresh component acceptances."""
from pathlib import Path
from fractions import Fraction as F
import gzip,hashlib,json
ROOT=Path(__file__).resolve().parents[3];HERE=Path(__file__).resolve().parent
paths={
 'baseline_directed_native_prices':'tonight/audit/global_b7171/final46_kappa_native_endpoint_replay.json',
 'exact_target_b_native_transport':'tonight/audit/global_b7171/b713_final46_native_transport.json',
 'target_b_R3_corners_and_whole_q':'tonight/audit/global_b7171/conditional_b713_fine46_R3_corner_replay.json',
 'uniform_planes_and_inverse':'tonight/audit/global_b7171/b713_final46_polygon_plane_prices.json',
 'fixed_contact_credits':'tonight/audit/global_b7171/b713_final46_selected_contact_prices.json',
 'complete_exact_middle_geometry':'tonight/theory/final46_b713_polygon_envelope_audit.json',
 'complete_exact_same_b_outer_transport':'tonight/audit/global_b7171/b713_complete_outer_transport.json',
 'all_line_source_bank':'tonight/audit/linear_contact_bank_independent.json',
 'all_line_contact_rows':'tonight/audit/cubic_head_contact_independent.json',
 'cubic_L1_moment_prices':'tonight/audit/refit847/independent_complete_price_cover.json'}
sha=lambda path:hashlib.sha256((ROOT/path).read_bytes()).hexdigest()
read=lambda path:json.loads((ROOT/path).read_text())
r={key:read(path)for key,path in paths.items()};assert all(v['status'].startswith('PASS')for v in r.values())
b=F(713,1000);E=r['complete_exact_middle_geometry'];assert F(E['operator_b'])==b
source=E['source'];digest=sha(source);assert digest==E['source_sha256']
envelope=json.loads(gzip.decompress((ROOT/source).read_bytes()))
assert F(envelope['operator_b'])==b
assert E['original_cells']==187335==E['accepted_price_cells']+E['exact_physically_empty_cells']
assert E['accepted_price_cells']==185761 and E['exact_physically_empty_cells']==1574
assert E['complete_groups']==945 and E['u_bands']==21 and E['kappa_bands']==45
assert len(envelope['cells'])==185761 and len(envelope['physically_empty_cells'])==1574
for key in('uniform_planes_and_inverse','fixed_contact_credits','exact_target_b_native_transport'):
 assert r[key]['input_sha256']==digest
R=r['target_b_R3_corners_and_whole_q'];assert R['source']['sha256']==digest
assert R['profiles']==14 and R['corner_tables']==630 and R['directed_corner_prices']==2520
assert F(R['operator_b'])==b
N=r['exact_target_b_native_transport'];assert F(N['operator_b'])==b and N['accepted_baseline_directed_prices']==9360
assert N['baseline_numeric_audit']['sha256']==sha(paths['baseline_directed_native_prices'])
assert len(N['frames'])==7 and N['kappa_bands']==45
P=r['uniform_planes_and_inverse'];assert P['uniform_polygon_planes']==450 and P['fresh_directed_corner_prices']==1800 and P['inverse_caps']==250
C=r['fixed_contact_credits'];assert C['selected_contact_cells']==97 and C['distinct_contact_corner_prices']==80
O=r['complete_exact_same_b_outer_transport'];assert F(O['operator_b'])==b
assert O['complete_outer_amplitude_coverage']==[['0','1381/2000'],['31/40','sqrt(2/pi)']]
assert O['accepted_native_cells_transported']==45025 and O['accepted_R3_cells_transported']==41
assert list(map(F,E['amplitude_interval']))==[F(1381,2000),F(31,40)]
B=F(envelope['strict_rational_cap']);maximum=F(E['maximum_cell_cap'])
assert maximum==F(envelope['maximum_cell_cap'])<B==F(564479481,10**9)
assert F(O['outer_norm_upper_bound'])<B
out=dict(status='PASS_INDEPENDENT_COMPLETE_PREFERRED_GLOBAL_B713_SCALAR_BOUND',
 operator_b=str(b),norm_upper_bound=str(B),constructive_G_R_lower_bound=str(1/B),
 constructive_G_R_lower_bound_display=float(1/B),proved_middle_maximum=str(maximum),
 proved_outer_cap=O['outer_norm_upper_bound'],global_amplitude_coverage=['0','sqrt(2/pi)'],
 source_rho='847/1000',middle_price_cells=185761,proved_physically_empty_middle_cells=1574,
 complete_middle_groups=945,cumulant_grid_points=46,
 baseline_directed_native_endpoint_prices=9360,exact_target_b_native_frames=7,
 independent_R3_corner_prices=2520,complete_whole_q_profiles=14,
 uniform_physical_source_planes=450,independent_plane_corner_prices=1800,
 independently_certified_physical_inverse_caps=250,fixed_contact_credit_cells=97,
 envelope=dict(path=source,sha256=digest),
 bindings={key:dict(path=path,sha256=sha(path),status=r[key]['status'])for key,path in paths.items()},
 historical_1773_or_1780498_premises_used=False,
 actual_sequence_limit='1/N(713/1000), bounded below by1000000000/564479481; the cap is not identified with N',
 scope='Complete same-coefficient global scalar norm proof and explicit reduced rational matrix lower bound. The target1.775 and a certified three-significant-digit value of G_R remain unproved. Higher-profile and alternative-coefficient discovery data are not used.')
(HERE/'PREFERRED_GLOBAL_ACCEPTANCE.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({key:value for key,value in out.items()if key!='bindings'},indent=2))
