#!/usr/bin/env python3
"""Transparent geometry replay with new helper and optional-legacy guard.
No existing immutable input is edited. Copied final checkers use the new
helper's actual digest and the actual patched predecessor's digest.
"""
from pathlib import Path
import subprocess,os,json,hashlib
R=Path('/workspace/scratch/8242d848994f');O=Path(__file__).resolve().parent
S=R/'new_round_2026-10-02/sequence_theory';D=O/'portable_geometry';D.mkdir(exist_ok=True)
helper=D/'audit_affine_hybrid_geometry.py';assert helper.is_file()
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
for p in S.glob('*.py'):(D/p.name).write_bytes(p.read_bytes())
p=D/'audit_rank_child_moment_refinement_v2.py';s=p.read_text()
old="allowed_parent_entries={H.sha(Path(__file__).with_name(n))for n in ('audit_monotone_rank_children.py','audit_monotone_rank_children_v2.py')}"
new="allowed_parent_entries={H.sha(Path(__file__).with_name(n))for n in ('audit_monotone_rank_children.py','audit_monotone_rank_children_v2.py')if Path(__file__).with_name(n).is_file()}"
assert s.count(old)==1;p.write_text(s.replace(old,new))
(O/'GEOMETRY_COMPATIBILITY_PATCH.json').write_text(json.dumps({'missing_unused_legacy_filename':'audit_monotone_rank_children.py','original_auditor':str(S/p.name),'original_sha256':sha(S/p.name),'patched_auditor':str(p),'patched_sha256':sha(p),'exact_removed_line':old,'exact_added_line':new,'actual_active_parent_hash_still_checked':True,'numeric_and_geometry_arithmetic_changes':False,'immutable_input_changes':False},indent=2)+'\n')
(O/'NEW_INDEPENDENT_AFFINE_HELPER_BINDING.json').write_text(json.dumps({'kind':'NEW_INDEPENDENT_RECONSTRUCTION','path':str(helper),'sha256':sha(helper),'original_helper_present_in_archive':False,'original_byte_identity_claimed':False,'specification':'Exact source-knot splitting; uniform native chords; quadratic secant/tangent majorants; closed affine dominance polygons; rational vertices and maxima.'},indent=2)+'\n')
c=json.loads((R/'new_round_2026-10-02/audit/portable_replay/FINAL_GOAL_REPLAY_CONFIG.json').read_text())
env=dict(os.environ,GR_FROZEN_BASELINE=str(R/'verified_2026-10-01/tonight/operator/preferred_global'),GR_MPFR_BITS='768',PYTHONDONTWRITEBYTECODE='1')
original={'strong43_geometry','u_refinement_geometry','layer_00_geometry'}
names=original|{'layer_01_geometry','layer_02_geometry','final_geometry'}
for j in c['jobs']:
 if j['name']not in names:continue
 def arg(x):return x.replace('{output}',str(O))if'{output}'in x else str(R/x)if x.startswith(('new_round','verified_'))else x
 script=R/j['script']if j['name']in original else D/Path(j['script']).name
 print('START',j['name'],flush=True)
 with(O/(j['name']+'.log')).open('w')as f:p=subprocess.run(['python',str(script),*map(arg,j['args'])],env=env,cwd=R,stdout=f,stderr=subprocess.STDOUT)
 assert p.returncode==0,(j['name'],p.returncode)
 report=json.loads((O/(j['name']+'.json')).read_text());assert report['status'].startswith('PASS')
 assert report['entrypoint_sha256']==sha(script)
 print('PASS',j['name'],flush=True)
