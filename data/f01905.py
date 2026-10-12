#!/usr/bin/env python3
"""Replay every selected native endpoint in any frozen same-b envelope."""
from pathlib import Path
from fractions import Fraction as F
import argparse,sys,json,time,gzip
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/'tonight/audit'))
import replay_native_R3_unions as N
ap=argparse.ArgumentParser();ap.add_argument('ledger',type=Path);ap.add_argument('--output',type=Path,required=True)
args=ap.parse_args();src=args.ledger.resolve();d=json.loads(gzip.decompress(src.read_bytes()) if src.suffix=='.gz' else src.read_text());digest=N.sha(src)
assert F(d['operator_b'])==F(7171,10000)
rows=[];reference=None;start=time.time()
for binding in d['native_endpoint_inputs']:
 before=N.POINTS;f=N.frame(binding,F(7171,10000))
 if reference is None:reference=f['keys']
 else:assert f['keys']==reference
 rows.append(dict(mask=binding['mask'],mask_sha256=binding['mask_sha256'],
                  package=binding['package'],source_sha256=binding['source_sha256'],
                  replay_sha256=binding['replay_sha256'],
                  selected_corner_prices_recomputed=N.POINTS-before,
                  r_grid_points=len(f['rg'])))
 print(json.dumps(dict(frame=len(rows)-1,prices=N.POINTS-before,status='PASS')),flush=True)
assert reference[0][0]==0 and reference[-1][1]==90
assert all(a[1]==b[0]for a,b in zip(reference,reference[1:]))
assert N.sha(src)==digest
points=[reference[0][0]]+[key[1]for key in reference]
out=dict(status='PASS_INDEPENDENT_ALL_REFINED_NATIVE_SELECTED_ENDPOINT_PRICES',
 source=dict(path=str(src.relative_to(ROOT)),sha256=digest),input_sha256=digest,
 operator_b='7171/10000',precision_bits=N.BITS,
 cumulant_points=list(map(str,points)),kappa_bands=len(reference),frames=rows,
 selected_endpoint_prices_recomputed=N.POINTS,
 accepted_package_bindings=[dict(package=name,manifest_sha256=p['manifest_sha256'],
                                replay_sha256=p['replay_sha256'],source_sha256=p['source_sha256'])
                            for name,p in N.PACKAGES.items()],
 producer_modules_imported=False,seconds=time.time()-start,
 scope='Numeric native endpoint premises only. Physical source geometry, profile caches, and global amplitude glue remain separately checked.')
args.output.write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:v for k,v in out.items()if k not in('frames','accepted_package_bindings','cumulant_points')},indent=2))
