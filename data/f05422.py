#!/usr/bin/env python3
"""Portable replay of the frozen original b7171 global scalar certificate.

Default: hashes, exact633, accepted primitive bindings, all directed endpoint
and corner prices, exact source partitions and complete global assembly.
--full additionally recomputes every scalar/local primitive and whole-q gate.
All mutable replay work takes place in a temporary copy of the bundle.
"""
from concurrent.futures import ThreadPoolExecutor
from fractions import Fraction as F
from pathlib import Path
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

ROOT=Path(__file__).resolve().parent
OUTPUT=ROOT/'replay_outputs'


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def read(path):return json.loads(path.read_text())


def verify_inputs():
    manifest=read(ROOT/'MANIFEST.json')
    for rel,digest in manifest['sha256'].items():assert sha(ROOT/rel)==digest,rel
    prov=read(ROOT/'PROVENANCE.json')
    with gzip.open(ROOT/prov['archival_copy'],'rb')as stream:
        original=stream.read()
    assert hashlib.sha256(original).hexdigest()==prov['original_envelope_sha256']
    old=json.loads(original);new=read(ROOT/'tonight/recovery/profile_cover/global_middle_b7171_native_R3_envelope.json')
    assert {k:v for k,v in old.items()if k!='accepted_R3_binding'}=={
            k:v for k,v in new.items()if k!='accepted_R3_binding'}
    assert all(a['sha256']==b['sha256']for a,b in zip(old['accepted_R3_binding'],new['accepted_R3_binding']))
    assert len(old['accepted_R3_binding'])==len(new['accepted_R3_binding'])
    for row in new['accepted_R3_binding']:
        assert not Path(row['path']).is_absolute()
        assert sha(ROOT/row['path'])==row['sha256']
    for name in prov['primitive_packages']:
        package=ROOT/name;hashes=read(package/'SHA256.json')
        for rel,digest in hashes.items():assert sha(package/rel)==digest,(name,rel)
        assert read(package/'replay.json')['status']in {
            'PASS_COMPLETE_LOCAL_AND_SCALAR_PRIMITIVES','PASS_COMPLETE_REGIONAL_BOUND'}
    accepted=read(ROOT/'reduced_search/energy/audit_stein_L1_independent.json')
    assert accepted['status']=='PASS_INDEPENDENT_COMPLETE_HALF_LINE_CUBIC_L1'
    assert len(accepted['rows'])==8
    for row in accepted['rows']:
        assert sha(ROOT/f"reduced_search/energy/stein_L1_row_{row['row']}.json")==row['row_sha256']
        assert row['compact']['domain']==['0','80']and row['tail']['domain']==['80','infinity']
    return manifest,prov


def execute(work,tag,args,bits=768):
    env=dict(os.environ);env['GR_MPFR_BITS']=str(bits)
    env['PYTHONPATH']=os.pathsep.join(str(work/rel)for rel in (
        'tonight/audit','tonight/theory','tonight/operator/projection/certificate'))
    log=OUTPUT/(tag+'.log')
    with log.open('w')as stream:
        result=subprocess.run([sys.executable,*map(str,args)],cwd=work,env=env,
                              stdout=stream,stderr=subprocess.STDOUT)
    if result.returncode:raise RuntimeError(f'Replay {tag} failed; inspect replay_outputs/{tag}.log')
    print(json.dumps(dict(stage='PASS',job=tag)),flush=True)


def full_primitives(work,prov,workers):
    jobs=[]
    for name in prov['primitive_packages']:
        script=work/name/'code/verify_region.py'
        args=[script,'--output',OUTPUT/('full_'+Path(name).name+'.json')]
        if '--primitives-only'in script.read_text():args.append('--primitives-only')
        jobs.append(('primitive_'+Path(name).name,args))
    with ThreadPoolExecutor(workers)as pool:
        for result in pool.map(lambda job:execute(work,*job),jobs):pass
    execute(work,'L1_old_eight_whole_line',[work/'reduced_search/energy/audit_stein_L1_independent.py'])
    for label in ('68','96'):
        execute(work,'L1_new_'+label+'_whole_line',[
            work/'tonight/audit/refit847/certify_refit_rows_independently.py','--row',label])
    # Critical755 compact replays also independently recheck the exact tail.
    for script in ('replay_r3_local_cover.py','replay_r3_bernstein_center_cover.py'):
        execute(work,'critical_'+script[:-3],[work/'tonight/theory'/script],384)
    for label,folder in (('68','lower_amp_certificate'),('70','lower_amp_certificate'),
                         ('interval_u7478_k0_01','diffuse_certificate'),
                         ('interval_b7171_newworst','diffuse_certificate')):
        audit_folder='lower_amp_R3'if label in('68','70')else(
            'diffuse_R3'if label=='interval_u7478_k0_01'else'diffuse_b7171_R3')
        for part in ('body','exterior','tail'):
            execute(work,'R3_'+label+'_'+part,[
                work/'tonight/audit/replay_R3_physical_q_trees.py','--profile',label,
                '--branch','finer','--part',part,'--data-dir',
                work/'tonight/operator/projection'/folder,'--output-dir',
                work/'tonight/audit'/audit_folder],384)
    for label in ('68','70'):
        execute(work,'fullE_'+label+'_whole_q',[
            work/'tonight/theory/replay_general_r3_cover.py','--label',label,
            '--workers',str(workers)],384)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--full',action='store_true')
    parser.add_argument('--workers',type=int,default=4);args=parser.parse_args()
    assert args.workers>=1
    OUTPUT.mkdir(exist_ok=True);start=time.time();manifest,prov=verify_inputs()
    print(json.dumps(dict(stage='MANIFEST_AND_ORIGINAL_ARCHIVE_PASS',files=len(manifest['sha256']))),flush=True)
    with tempfile.TemporaryDirectory(prefix='grothendieck_b7171_replay_')as tmp:
        work=Path(tmp)/'proof';shutil.copytree(ROOT,work,ignore=shutil.ignore_patterns('replay_outputs','__pycache__'))
        if args.full:full_primitives(work,prov,args.workers)
        for name in ('independent_replay.py','audit_cover_independent_fraction.py','certify_wick_identity.py'):
            execute(work,'633_'+name[:-3],[work/'reduced_search/energy/633'/name])
        execute(work,'L1_complete_ten_row_prices',[work/'tonight/audit/refit847/replay_complete_refit_cover.py'])
        execute(work,'first_four_transfers',[work/'tonight/audit/check_b7171_first_transfers.py'])
        native=work/'tonight/audit/replay_native_R3_unions.py'
        cover=work/'tonight/recovery/profile_cover'
        r3=work/'tonight/operator/projection/lower_amp_certificate/b7171_native_remainder_extended_box.json'
        jobs=[
            ('lower_base','native_lower_base_b7171_target_union.json',['.6','.682'],None,False,
             work/'tonight/audit/native_b7171_union/u600_682_independent_replay.json'),
            ('lower_R3','native_union_b7171_target.json',['.682','.6905'],r3,False,
             work/'tonight/audit/native_b7171_union/u682_6905_independent_replay.json'),
            ('upper','native_upper_b7171_target_union.json',['.775','.797885'],None,False,
             work/'tonight/audit/native_b7171_union/u775_c_independent_replay.json'),
            ('middle_native_lower','native_union_b7171_target.json',['.6905','.720'],None,True,
             work/'tonight/audit/global_b7171/native_middle_lower_input_replay.json'),
            ('middle_native_upper','native_upper_b7171_target_union.json',['.720','.775'],None,True,
             work/'tonight/audit/global_b7171/native_middle_upper_input_replay.json')]
        for tag,ledger,interval,box,envelope,output in jobs:
            command=[native,'--ledger',cover/ledger,'--u',*interval,'--output',output]
            if box is not None:command+=['--r3',box]
            if envelope:command+=['--envelope']
            execute(work,tag,command)
        execute(work,'extended_remainder_source_price',[work/'replay_extended_remainder.py'])
        execute(work,'R3_all_392_corners_and_gate_bindings',[
            work/'tonight/audit/global_b7171/replay_R3_corners.py'],384)
        execute(work,'exact_lower_amplitude_glue',[work/'tonight/theory/audit_same_b_native_R3_glue.py'])
        execute(work,'exact_global_middle_partition',[work/'tonight/theory/audit_global_native_R3_envelope.py'])
        execute(work,'global_acceptance_assembly',[work/'tonight/audit/global_b7171/assemble_independent_acceptance.py'])
        accepted=read(work/'tonight/audit/global_b7171/GLOBAL_B7171_ACCEPTANCE.json')
        assert accepted['status']=='PASS_INDEPENDENT_COMPLETE_GLOBAL_B7171_SCALAR_BOUND'
        assert F(accepted['norm_upper_bound'])==F(70756627,125000000)
        assert F(accepted['constructive_G_R_lower_bound'])==F(125000000,70756627)
        for name,entry in accepted['bindings'].items():
            source=work/entry['path'];target=OUTPUT/(name+'.json');shutil.copy2(source,target)
        result=dict(status='PASS_PORTABLE_FROZEN_ORIGINAL_GLOBAL_B7171',
            mode='complete_fresh_primitives_prices_and_partition'if args.full else
                 'accepted_fresh_primitive_reports_plus_fresh_prices_and_partition',
            norm_upper_bound=accepted['norm_upper_bound'],
            constructive_G_R_lower_bound=accepted['constructive_G_R_lower_bound'],
            original_envelope_sha256=prov['original_envelope_sha256'],
            immutable_bundle_manifest_pass=True,archival_mathematical_content_equal=True,
            scalar_local_primitives_freshly_replayed_this_run=args.full,
            independent_endpoint_prices=accepted['independent_endpoint_prices'],
            independent_R3_corner_prices=accepted['independent_R3_corner_prices'],
            global_middle_cells=accepted['global_middle_cells'],seconds=time.time()-start,
            scope='Frozen original global norm cap and constructive lower bound only. It does not prove a three-digit value of the Grothendieck constant or use later research refinements.')
        (OUTPUT/'result.json').write_text(json.dumps(result,indent=2)+'\n')
        print(json.dumps(result,indent=2))


if __name__=='__main__':main()
