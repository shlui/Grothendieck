#!/usr/bin/env python3
"""Bind fresh accepted endpoint tables to exact smaller-b native prices."""
from pathlib import Path
from fractions import Fraction as F
import argparse,json,gzip,sys
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/'tonight/audit'))
import replay_native_R3_unions as N
ap=argparse.ArgumentParser();ap.add_argument('ledger',type=Path);ap.add_argument('--output',type=Path,required=True)
args=ap.parse_args();src=args.ledger.resolve();d=json.loads(gzip.decompress(src.read_bytes())if src.suffix=='.gz'else src.read_text())
b=F(d['operator_b']);b0=F(7171,10000);assert 0<b<=b0
oldpath=ROOT/'tonight/audit/global_b7171/final57_kappa_native_endpoint_replay.json'
old=json.loads(oldpath.read_text());assert old['status']=='PASS_INDEPENDENT_ALL_FINAL57_KAPPA_NATIVE_SELECTED_ENDPOINT_PRICES'
assert F(old['operator_b'])==b0 and old['selected_endpoint_prices_recomputed']==11648
accepted={r['mask']:r for r in old['frames']};rows=[];reference=None
for bind in d['native_endpoint_inputs']:
 prior=accepted[bind['mask']]
 for key in('mask_sha256','package','source_sha256'):assert bind[key]==prior[key]
 maskpath=ROOT/bind['mask'];assert N.sha(maskpath)==bind['mask_sha256']
 raw=json.loads(maskpath.read_text());assert F(raw['operator_b'])==b0
 profile=N.package(bind['package']);p=profile['profile'];M=profile['M'];beta=F(p['beta'])
 assert profile['source_sha256']==bind['source_sha256'] and profile['replay_sha256']==bind['replay_sha256']
 for a,o,c in ((M[0][0]+1,M[0][2],M[2][2]-b),(M[1][1]+1,M[1][3],M[3][3]-b),(M[1][1]-b,M[1][3],M[3][3])):N.psd(a,o,c)
 oldcost=(F(p['nu'])-M[0][2])**2/M[2][2]-M[0][0]+b0-2*F(p['delta'])
 assert oldcost==F(raw['exact_gram_cost'])
 oldshift=beta**2*(1/b0-1/F(p['epsilon']))
 assert oldshift==F(raw['local_constant_transport'])
 delta=beta**2*(1/b-1/b0)
 targetcost=oldcost+b-b0;targetshift=oldshift+delta
 assert targetshift==beta**2*(1/b-1/F(p['epsilon']))
 keys=sorted(tuple(map(F,r['kappa_interval']))for r in raw['endpoint_prices'])
 if reference is None:reference=keys
 else:assert keys==reference
 assert keys[0][0]==0 and keys[-1][1]==90 and all(a[1]==c[0]for a,c in zip(keys,keys[1:]))
 for endpoint in raw['endpoint_prices']:
  assert len(endpoint['prices'])==profile['rg'].__len__()
  for corner in endpoint['prices']:assert len(corner)==4
 rows.append(dict(**bind,target_gram_cost=str(targetcost),baseline_gram_cost=str(oldcost),
                  exact_endpoint_constant_addition=str(delta),target_local_constant_transport=str(targetshift),
                  normalized_beta_squared=str(beta*beta),r_grid_points=len(profile['rg'])))
assert len(rows)==7 and len(reference)==56
assert list(map(F,d['kappa_grid']))==[reference[0][0]]+[k[1]for k in reference]
out=dict(status='PASS_INDEPENDENT_EXACT_REFINED_NATIVE_SMALLER_B_TRANSPORT',
 source=dict(path=str(src.relative_to(ROOT)),sha256=N.sha(src)),input_sha256=N.sha(src),
 operator_b=str(b),baseline_b=str(b0),baseline_numeric_audit=dict(path=str(oldpath.relative_to(ROOT)),sha256=N.sha(oldpath)),
 accepted_baseline_directed_prices=11648,kappa_bands=56,frames=rows,
 target_spectral_gates_checked=True,producer_code_imported=False,
 analytic_review='The frozen target, matrix and scalar rows are unchanged. Every saved four-corner endpoint bound receives exactly beta_norm^2*(1/b-1/b0), while the radial quadratic coefficient receives b-b0. All target-b PSD gates are independently checked. The physical r/kappa/shape domain is unchanged. No new scalar integration is required because the only b-dependence is this exact completion and spectral transport.',
 scope='Exact native numeric premises at the target coefficient only. Complete middle geometry and SAME-b outer amplitude coverage require their separate accepted audits.')
args.output.write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:v for k,v in out.items()if k!='frames'},indent=2))
