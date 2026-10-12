#!/usr/bin/env python3
"""Complete path-independent full-upper/adaptive preferred certificate replay."""
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
    with(gzip.open(path,'rt')if path.suffix=='.gz'else path.open())as f:return json.load(f)
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def run(work,tag,args,bits=768):
    env=dict(os.environ);env['GR_MPFR_BITS']=str(bits);env['PYTHONINTMAXSTRDIGITS']='0'
    env['PYTHONPATH']=os.pathsep.join(str(work/p)for p in('tonight/audit','tonight/theory','tonight/operator/projection/certificate'))
    with(OUT/(tag+'.log')).open('w')as stream:
        r=subprocess.run([sys.executable,*map(str,args)],cwd=work,env=env,stdout=stream,stderr=subprocess.STDOUT)
    if r.returncode:raise RuntimeError(f'{tag} failed; inspect replay_outputs/{tag}.log')
    print(json.dumps(dict(stage='PASS',job=tag)),flush=True)
def arguments(work,args):return[work/v if v.startswith(('tonight/','reduced_search/'))else v for v in args]

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--full',action='store_true');ap.add_argument('--workers',type=int,default=4)
    ap.add_argument('--resume-base',action='store_true',help='Reuse a byte-bound complete base replay checkpoint; all preferred stages are fresh.')
    ap.add_argument('--manifest-only',action='store_true',help='Extraction smoke test only; does not replay mathematical prices or partitions.')
    x=ap.parse_args();start=time.time();OUT.mkdir(exist_ok=True)
    manifest=read(ROOT/'MANIFEST.json');assert manifest['extension_filter_used']is False
    assert all(sha(ROOT/rel)==h for rel,h in manifest['sha256'].items())
    provenance=read(ROOT/'PROVENANCE.json');preferred=provenance['preferred']
    assert sha(ROOT/preferred['original_source_archive'])==preferred['original_source_sha256']
    assert sha(ROOT/preferred['source'])==preferred['portable_source_sha256']
    if x.manifest_only:
        print(json.dumps(dict(status='PASS_PORTABLE_EXTRACTION_MANIFEST_AND_ORIGINAL_SOURCE_SMOKE',
            manifest_files=len(manifest['sha256']),operator_b=preferred['operator_b'],
            norm_cap=preferred['norm_cap'],scope='Immutable file identities and relative package lookup only; complete mathematical replay is a separate mode.'),indent=2))
        return
    for item in preferred['archived_operational_json']:
        raw=ROOT/item['archive'];assert sha(raw)==item['original_sha256']
        data=raw.read_bytes();data=gzip.decompress(data)if data[:2]==b'\x1f\x8b'else data
        # All portable metadata changes are the explicitly bound deterministic
        # path-prefix and metadata-hash substitutions. No optimizer is invoked.
        for before,after in item['ordered_replacements']:data=data.replace(before.encode(),after.encode())
        assert json.loads(data)==read(ROOT/item['path'])
    baseline=[ROOT/'verify_strengthened_components.py']
    if x.full:baseline+=['--full','--workers',str(x.workers)]
    if x.resume_base:
        assert not x.full,'Full scalar mode requires a fresh base replay'
        checkpoint=read(ROOT/'BASE_REPLAY_CHECKPOINT.json')
        assert all(sha(ROOT/rel)==h for rel,h in checkpoint['unchanged_base_inputs'].items())
        if 'baseline_provenance'in checkpoint:
            assert {k:v for k,v in provenance.items()if k not in('preferred','proof_scope')}==checkpoint['baseline_provenance']
        assert sha(ROOT/'BASE_REPLAY_RESULT.json')==checkpoint['base_result_sha256']
        assert read(ROOT/'BASE_REPLAY_RESULT.json')['status']=='PASS_PORTABLE_STRENGTHENED_GLOBAL_B7171'
        print(json.dumps(dict(stage='PASS',job='accepted_original_and_strengthened_components_byte_bound_checkpoint')),flush=True)
    else:
        with tempfile.TemporaryDirectory(prefix='preferred_immutable_base_replay_')as temporary:
            base_work=Path(temporary)/'proof'
            shutil.copytree(ROOT,base_work,ignore=shutil.ignore_patterns('replay_outputs','replay_outputs_targeted','replay_outputs_promoted','__pycache__'))
            isolated_baseline=[base_work/'verify_strengthened_components.py']
            if x.full:isolated_baseline+=['--full','--workers',str(x.workers)]
            run(base_work,'accepted_original_and_strengthened_components',isolated_baseline)
        assert all(sha(ROOT/rel)==h for rel,h in manifest['sha256'].items())
    with tempfile.TemporaryDirectory(prefix='preferred_fullupper_replay_')as temporary:
        work=Path(temporary)/'proof'
        shutil.copytree(ROOT,work,ignore=shutil.ignore_patterns('replay_outputs','replay_outputs_targeted','replay_outputs_promoted','__pycache__'))
        config=read(work/'REPLAY_CONFIG.json')
        if x.full:
            run(work,'original_fourteen_profile_R3_prerequisites',[work/'verify_targeted_prerequisites.py','--workers',str(x.workers)])
            run(work,'all_additional_promoted_rank_prerequisites',[work/'verify_promoted_prerequisites.py',
                '--cache-report',work/config['accepted_cache_report'],'--workers',str(x.workers)])
            for number,report in enumerate(config.get('additional_gate_binding_reports',[])):
                run(work,f'final_delta_promoted_rank_prerequisites_{number}',[work/'verify_promoted_prerequisites.py',
                    '--cache-report',work/report,'--workers',str(x.workers)])
        for job in config['jobs']:run(work,job['name'],arguments(work,job['args']),job.get('bits',768))
        run(work,'complete_fullupper_same_b_global_assembly',arguments(work,config['assembly_args']))
        result=read(work/preferred['acceptance'])
        assert result['status'].startswith('PASS_INDEPENDENT_COMPLETE_')
        assert F(result['operator_b'])==F(preferred['operator_b'])
        assert F(result['norm_upper_bound'])==F(preferred['norm_cap'])
        assert F(result['constructive_G_R_lower_bound'])==F(preferred['GR_lower_bound'])
        for name,binding in result['bindings'].items():shutil.copy2(work/binding['path'],OUT/('preferred_'+name+'.json'))
        for name,binding in result.get('refinement_bindings',{}).items():shutil.copy2(work/binding['path'],OUT/('preferred_refinement_'+name+'.json'))
        for name,binding in result.get('monotone_bindings',{}).items():
            if isinstance(binding,dict)and 'path'in binding:shutil.copy2(work/binding['path'],OUT/('preferred_monotone_'+name+'.json'))
        shutil.copy2(work/preferred['acceptance'],OUT/'preferred_complete_acceptance.json')
        summary=dict(status='PASS_PORTABLE_COMPLETE_FULLUPPER_PREFERRED_GLOBAL',
            operator_b=result['operator_b'],norm_upper_bound=result['norm_upper_bound'],
            constructive_G_R_lower_bound=result['constructive_G_R_lower_bound'],
            original_source_sha256=preferred['original_source_sha256'],
            mathematical_and_operational_rewrite_exact=True,immutable_evidence_manifest_pass=True,
            mode='all_local_scalar_gates_and_complete_global_replay'if x.full else'accepted_fresh_primitive_reports_plus_fresh_prices_and_complete_cover',
            scalar_local_primitives_freshly_replayed_this_run=x.full,seconds=time.time()-start,
            original_and_strengthened_base_reused_from_complete_fresh_checkpoint=x.resume_base,
            extension_filter_used=False,full_amplitude_coverage=result['full_amplitude_coverage'],
            scope='Complete same-coefficient full-amplitude norm and constructive lower bound. The actual matrix sequence limit is 1/N(b); a three-digit enclosure of unrestricted G_R is not asserted.')
        for key in('cumulant_points','complete_whole_q_profiles','full_upper_original_cells','full_upper_priced_cells',
            'proved_physically_empty_full_upper_cells','full_upper_complete_groups','baseline_native_endpoint_prices',
            'independent_R3_corner_prices','independent_plane_corner_prices','uniform_physical_source_planes',
            'forward_verified_inverse_source_caps','fixed_contact_credit_cells','new_native_endpoint_prices','new_uniform_source_planes'):
            if key in result:summary[key]=result[key]
        assert all(sha(ROOT/rel)==h for rel,h in manifest['sha256'].items())
        (OUT/'result.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
