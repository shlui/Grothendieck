#!/usr/bin/env python3
"""Independent selected prices on all seven final46-cumulant endpoint tables."""
from pathlib import Path
from fractions import Fraction as F
import sys,json,hashlib,time
ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT/'tonight/audit'))
import replay_native_R3_unions as N
src=ROOT/'tonight/recovery/profile_cover/kappa_refinement_2_4/complete.json';d=json.loads(src.read_text());assert d['status']=='PASS_ALL_FRESH_FINE_KAPPA_ENDPOINT_LEDGERS'
points=list(map(F,d['kappa_points']));assert len(points)==46 and points[0]==0 and points[-1]==90 and all(a<c for a,c in zip(points,points[1:]))
rows=[];start=time.time()
for bind in d['frames']:
    binding={**bind,'mask_sha256':N.sha(ROOT/bind['mask'])};before=N.POINTS
    f=N.frame(binding,F(7171,10000));assert f['keys']==list(zip(points,points[1:]))
    rows.append({'mask':bind['mask'],'mask_sha256':binding['mask_sha256'],'package':bind['package'],'source_sha256':bind['source_sha256'],'selected_corner_prices_recomputed':N.POINTS-before,'r_grid_points':len(f['rg'])})
    print(json.dumps({'frame':bind['index'],'prices':N.POINTS-before,'status':'PASS'}),flush=True)
assert len(rows)==7
out={'status':'PASS_INDEPENDENT_ALL_FINAL46_KAPPA_NATIVE_SELECTED_ENDPOINT_PRICES','input_sha256':N.sha(src),'operator_b':'7171/10000','precision_bits':N.BITS,'cumulant_points':list(map(str,points)),'frames':rows,'selected_endpoint_prices_recomputed':N.POINTS,'accepted_package_bindings':[{'package':name,'manifest_sha256':p['manifest_sha256'],'replay_sha256':p['replay_sha256'],'source_sha256':p['source_sha256']}for name,p in N.PACKAGES.items()],'producer_modules_imported':False,'seconds':time.time()-start,'scope':'Numeric native endpoint premises only. Fine-cumulant physical partition, source constraints, all profile-bank caches, and final global envelope require their own audits.'}
p=Path(__file__).parent/'final46_kappa_native_endpoint_replay.json';p.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items()if k not in['frames','accepted_package_bindings','cumulant_points']},indent=2))
