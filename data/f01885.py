#!/usr/bin/env python3
"""Path-independent complete preferred-global replay, without optimizers."""
from pathlib import Path
from fractions import Fraction as F
import argparse
import gzip
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time

ROOT=Path(__file__).resolve().parent;OUT=ROOT/'replay_outputs'
def read(path):
    with(gzip.open(path,'rt')if path.suffix=='.gz'else path.open())as stream:return json.load(stream)
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def run(work,tag,args,bits=768):
    env=dict(os.environ);env['GR_MPFR_BITS']=str(bits);env['PYTHONINTMAXSTRDIGITS']='0'
    env['PYTHONPATH']=os.pathsep.join(str(work/p)for p in('tonight/audit','tonight/theory','tonight/operator/projection/certificate'))
    with(OUT/(tag+'.log')).open('w')as stream:
        result=subprocess.run([sys.executable,*map(str,args)],cwd=work,env=env,stdout=stream,stderr=subprocess.STDOUT)
    if result.returncode:raise RuntimeError(f'{tag} failed; inspect replay_outputs/{tag}.log')
    print(json.dumps(dict(stage='PASS',job=tag)),flush=True)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--full',action='store_true')
    parser.add_argument('--workers',type=int,default=4);args=parser.parse_args();start=time.time()
    OUT.mkdir(exist_ok=True);manifest=read(ROOT/'MANIFEST.json')
    assert all(sha(ROOT/p)==h for p,h in manifest['sha256'].items())
    provenance=read(ROOT/'PROVENANCE.json');preferred=provenance['preferred']
    original_archive=ROOT/preferred['original_source_archive'];assert sha(original_archive)==preferred['original_source_sha256']
    original_bytes=original_archive.read_bytes()
    if original_bytes[:2]==b'\x1f\x8b':original_bytes=gzip.decompress(original_bytes)
    original=json.loads(original_bytes);portable=read(ROOT/preferred['source'])
    for field in('cells','physically_empty_cells','cached_corner_constants','profile_coefficients',
                  'maximum_cell_cap','strict_rational_cap','source_polygon_planes','physical_u_bounds'):
        assert original[field]==portable[field]
    # The strengthened component entrypoint includes all earlier same-b
    # branches and scalar/L1 source-price prerequisites. Data remain immutable.
    base=[ROOT/'verify_strengthened_components.py']
    if args.full:base+=['--full','--workers',str(args.workers)]
    run(ROOT,'accepted_strengthened_global_components',base)
    with tempfile.TemporaryDirectory(prefix='preferred_global_replay_')as temporary:
        work=Path(temporary)/'proof';shutil.copytree(ROOT,work,ignore=shutil.ignore_patterns('replay_outputs','replay_outputs_targeted','__pycache__'))
        config=read(work/'REPLAY_CONFIG.json')
        if args.full:
            run(work,'all_promoted_R3_q_prerequisites',[work/'verify_targeted_prerequisites.py','--workers',str(args.workers)])
        for job in config['jobs']:
            arguments=[work/arg if arg.startswith(('tonight/','reduced_search/'))else arg for arg in job['args']]
            run(work,job['name'],arguments,job.get('bits',768))
        run(work,'preferred_complete_global_assembly',[work/preferred['assembly']])
        result=read(work/preferred['acceptance'])
        assert result['status'].startswith('PASS_INDEPENDENT_COMPLETE_')
        assert F(result['norm_upper_bound'])==F(preferred['norm_cap'])
        assert F(result['constructive_G_R_lower_bound'])==F(preferred['GR_lower_bound'])
        for name,binding in result['bindings'].items():
            shutil.copy2(work/binding['path'],OUT/('preferred_'+name+'.json'))
        summary=dict(status='PASS_PORTABLE_COMPLETE_PREFERRED_GLOBAL_B7171',
            norm_upper_bound=result['norm_upper_bound'],constructive_G_R_lower_bound=result['constructive_G_R_lower_bound'],
            original_source_sha256=preferred['original_source_sha256'],
            mathematical_source_fields_exactly_equal=True,immutable_evidence_manifest_pass=True,
            mode='complete_scalar_local_primitives_and_prices'if args.full else'accepted_fresh_primitive_reports_plus_fresh_prices_and_complete_cover',
            scalar_local_primitives_freshly_replayed_this_run=args.full,seconds=time.time()-start,
            extension_filter_used=False,
            scope='Complete preferred same-coefficient global proof only. All claimed local/scalar inputs and exact domain data are present; the manifest includes every prerequisite file without extension filtering.')
        for key in('middle_price_cells','proved_physically_empty_middle_cells','complete_middle_groups',
            'independent_R3_corner_prices','uniform_physical_source_planes','independent_plane_corner_prices',
            'independent_physical_inverse_caps','selected_contact_credit_cells','native_knot_pieces','checked_edge_stationary_points'):
            if key in result:summary[key]=result[key]
        (OUT/'result.json').write_text(json.dumps(summary,indent=2)+'\n')
        print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
