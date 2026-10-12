#!/usr/bin/env python3
"""Bind fresh geometry/prices, exact uncapping, and strict sequence arithmetic."""
from pathlib import Path
from fractions import Fraction as F
import argparse,json,hashlib,subprocess
R=Path('/workspace/scratch/8242d848994f');O=Path(__file__).resolve().parent
EXPECTED_ARCHIVE_SHA='35f1119b00800862571051990e2e89eb97f6d32cbcf28d69842bccb4374b4d33'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text())
def binding(p):return {'path':str(Path(p).resolve()),'sha256':sha(p)}
def fraction(q):return F(int(q['numerator']),int(q['denominator']))
ap=argparse.ArgumentParser();ap.add_argument('--archive',type=Path);a=ap.parse_args()
archive={'expected_sha256':EXPECTED_ARCHIVE_SHA,'archive_bytes_checked_in_this_run':False}
if a.archive:
 assert sha(a.archive)==EXPECTED_ARCHIVE_SHA
 archive.update(binding(a.archive));archive['archive_bytes_checked_in_this_run']=True
manifest_path=R/'new_round_2026-10-02/audit/portable_replay/FINAL_GOAL_INPUT_MANIFEST.json';manifest=read(manifest_path)
for n,h in manifest['sha256'].items():assert sha(R/n)==h
assert manifest['extension_filter_used']is False
patch=read(O/'GEOMETRY_COMPATIBILITY_PATCH.json');assert patch['numeric_and_geometry_arithmetic_changes']is False and patch['immutable_input_changes']is False
assert sha(patch['original_auditor'])==patch['original_sha256'] and sha(patch['patched_auditor'])==patch['patched_sha256']
s=Path(patch['original_auditor']).read_text();assert s.count(patch['exact_removed_line'])==1
assert s.replace(patch['exact_removed_line'],patch['exact_added_line'])==Path(patch['patched_auditor']).read_text()
helper=read(O/'NEW_INDEPENDENT_AFFINE_HELPER_BINDING.json');assert sha(helper['path'])==helper['sha256'] and helper['original_byte_identity_claimed']is False
numeric=read(O/'ALL_DEFAULT_NUMERIC_REPLAYS.json');assert numeric['status']=='PASS_FRESH_ALL_DEFAULT_NUMERIC_PRICE_REPLAYS' and numeric['jobs']==31 and numeric['bits']==768
for job in numeric['results']:
 p=Path(job['fresh_report']);assert sha(p)==job['fresh_sha256'] and read(p)['status']==job['status']
geometry=[]
for name in ('strong43_geometry','u_refinement_geometry','layer_00_geometry','layer_01_geometry','layer_02_geometry','final_geometry'):
 p=O/(name+'.json');g=read(p);assert g['status'].startswith('PASS_PRODUCER_FREE_')
 script=R/'new_round_2026-10-02/sequence_theory'/({'strong43_geometry':'audit_strong_shape_geometry.py','u_refinement_geometry':'audit_u_refinement_geometry.py','layer_00_geometry':'audit_monotone_rank_children_v2.py'}[name])if name in {'strong43_geometry','u_refinement_geometry','layer_00_geometry'}else O/'portable_geometry'/({'layer_01_geometry':'audit_rank_child_moment_refinement_v2.py','layer_02_geometry':'audit_affine_hybrid_delta_v2.py','final_geometry':'audit_affine_moment_subrefinement.py'}[name])
 assert sha(script)==g['entrypoint_sha256']
 if 'parent'in g:
  parent=g['parent'];assert sha(parent['accepted_geometry_report'])==parent['accepted_geometry_report_sha256']
  assert read(parent['accepted_geometry_report'])['source']['sha256']==parent['sha256']
 if name in {'layer_02_geometry','final_geometry'}:assert g['affine_replay_sha256']==helper['sha256']
 geometry.append({'name':name,**binding(p),'actual_checker':binding(script)})
# Exact source assembly and finite arithmetic are rerun, not merely read.
for name in ('check_uncapped.py','finite_sequence_arithmetic.py','test_affine_geometry_independent.py'):
 with(O/(name+'.verification.log')).open('w')as log:
  result=subprocess.run(['python',str(O/name)],stdout=log,stderr=subprocess.STDOUT)
 assert result.returncode==0,name
assembly=read(O/'STRICT_UNCAPPED_ASSEMBLY.json');finite=read(O/'finite_sequence_arithmetic.json');fixtures=read(O/'AFFINE_HELPER_MANUAL_FIXTURES.json')
raw=F(assembly['exact_upper_maximum']);assert raw==F(281690137769807,500000000000000)
assert assembly['counts']['final_leaf_count']==305039 and assembly['strict_all_leaves']is True and assembly['threshold_equal_leaf_count']==0
B0=F(563380276,10**9);b=F(711000001,10**9);B=F(563380277,10**9)
assert raw<=B0 and B==B0+F(1,10**9)<F(40,71)
assert F(563307766601,10**12)<raw
assert fraction(finite['b'])==b and fraction(finite['perturbed_norm_cap'])==B
assert fraction(finite['finite_lower_bound_limit'])==1/B and 1/B>F(71,40)
for n,row in enumerate(finite['first_three'],1):
 d=int(row['gaussian_dimension']);assert d==10**9*2**(n-1)
 value=(1-(1+b)/(2*d)-F(10,d*d))/(B+F(1,d*d));assert fraction(row['lower_bound'])==value>F(71,40)
assert len(finite['first_three'])==3
# Norm perturbation uses only the probability-space contraction
# ||P3 f||1 <= ||P3 f||2 <= ||f||2 <= ||f||infinity.
answer={'status':'PASS_STRICT_MODIFIED_SEQUENCE_CERTIFICATE','unperturbed_coefficient':'711/1000','fresh_uncapped_operator_norm_cap':str(raw),'unperturbed_simple_cap':str(B0),'modified_coefficient':str(b),'modified_norm_cap':str(B),'modified_cap_margin_below_40_over_71':str(F(40,71)-B),'actual_matrix_ratio_limit_lower_bound':str(1/B),'actual_ratio_limit_lower_bound_decimal':str(float(1/B)),'strict_ratio_margin_above_71_over_40':str(1/B-F(71,40)),'sequence':'H_n=H^*_{10^9*2^(n-1)}(711000001/10^9), n>=1','first_three':finite['first_three'],'fresh_directed_price_jobs':31,'fresh_closed_geometry_stages':6,'independently_checked_final_leaves':305039,'immutable_manifest_files_rechecked':len(manifest['sha256']),'archive':archive,'bindings':{'numeric_replay':binding(O/'ALL_DEFAULT_NUMERIC_REPLAYS.json'),'geometry_stages':geometry,'uncapped_assembly':binding(O/'STRICT_UNCAPPED_ASSEMBLY.json'),'uncapped_checker':binding(O/'check_uncapped.py'),'finite_arithmetic':binding(O/'finite_sequence_arithmetic.json'),'finite_checker':binding(O/'finite_sequence_arithmetic.py'),'new_helper':binding(helper['path']),'manual_helper_fixtures':binding(O/'AFFINE_HELPER_MANUAL_FIXTURES.json'),'manual_helper_fixture_checker':binding(O/'test_affine_geometry_independent.py'),'compatibility_patch':binding(O/'GEOMETRY_COMPATIBILITY_PATCH.json'),'input_manifest':binding(manifest_path)},'scope':'All current directed price jobs and six exact current geometry stages freshly replayed. Whole-line local inequalities and the protected original all-local norm proof remain inherited accepted premises; the entire 219-stage local replay was not rerun. New affine helper independently reconstructed and reviewed; no claim to absent original byte identity. No frozen proof input changed.'}
(O/'STRICT_SEQUENCE_ACCEPTANCE.json').write_text(json.dumps(answer,indent=2)+'\n')
print(json.dumps({k:v for k,v in answer.items()if k not in('bindings','first_three')},indent=2))
