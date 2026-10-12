#!/usr/bin/env python3
"""Optional fresh replay of every additional accepted rank-profile gate."""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import argparse
import hashlib
import json
import os
import subprocess
import sys
import time

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'replay_outputs_promoted'
OLD_FOLDERS={'targeted_middle_R3_recheck','targeted_middlehigh_R3_recheck',
    'robust_middlehigh_R3_recheck','robust_currentworst_R3_recheck','robust_broadcurrent_R3_recheck'}

def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def resolved(value):
    p=Path(value);return p if p.is_absolute()else ROOT/p

def replay(job):
    tag,adapter,which,branch,part,data_dir,radius,expected,profile=job
    destination=OUT/tag;destination.mkdir(parents=True,exist_ok=True)
    args=[sys.executable,str(adapter),'--profile',which,'--branch',branch,'--part',part,
        '--data-dir',str(data_dir),'--tail-R',str(radius),'--output-dir',str(destination)]
    env=dict(os.environ);env['PYTHONINTMAXSTRDIGITS']='0'
    env['PYTHONPATH']=os.pathsep.join(str(ROOT/p)for p in('tonight/audit','tonight/theory','tonight/operator/projection/certificate'))
    with(destination/'replay.log').open('w')as log:
        done=subprocess.run(args,cwd=ROOT,env=env,stdout=log,stderr=subprocess.STDOUT)
    assert done.returncode==0,(tag,'see replay_outputs_promoted/'+tag+'/replay.log')
    file=destination/f'independent_{which}_{branch}_{part}.json';fresh=read(file)
    assert fresh['status']==expected['status']
    if part=='tail':
        for key in('R','difference_polynomial','shifted_coefficients'):assert fresh[key]==expected[key]
    else:
        for key in('profile_sha256','tree_sha256','accepted','excluded','nodes'):assert fresh[key]==expected[key]
        assert fresh['profile_sha256']==sha(profile)
    assert fresh['cubic_integral_backend_sha256']==expected['cubic_integral_backend_sha256']
    print(json.dumps(dict(stage='PASS',job=tag)),flush=True)
    return dict(job=tag,fresh_report=str(file.relative_to(ROOT)),sha256=sha(file),status=fresh['status'])

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--cache-report',type=Path,required=True)
    ap.add_argument('--workers',type=int,default=4);ap.add_argument('--quick',action='store_true');x=ap.parse_args()
    start=time.time();OUT.mkdir(exist_ok=True)
    manifest=read(ROOT/'MANIFEST.json')
    assert all(sha(ROOT/rel)==h for rel,h in manifest['sha256'].items())
    cache=read(resolved(x.cache_report));assert cache['status'].startswith('PASS_INDEPENDENT_GENERIC_COMPLETE_')
    jobs=[];profiles=[]
    for name,binding in cache['complete_whole_q_bindings'].items():
        if 'accepted_adapter'not in binding:continue
        gate=resolved(binding['accepted_whole_q_gate']['path'])
        if gate.parent.name in OLD_FOLDERS:continue
        ledger=read(gate);assert ledger['status'].startswith('PASS_INDEPENDENT_COMPLETE_')
        assert sha(gate)==binding['accepted_whole_q_gate']['sha256']
        adapter=resolved(binding['accepted_adapter']['path'])
        assert sha(adapter)==binding['accepted_adapter']['sha256']==ledger['adapter_sha256']
        profile=resolved(binding['profile']['path']);assert sha(profile)==ledger['profile_sha256']
        body=read(resolved(ledger['reports']['body']['path']))
        which=body['profile'];branch=body['branch'];radius=binding['fresh_exact_tail']['R']
        data_dirs={resolved(region['tree']['path']).parent for region in binding['q_regions']}
        assert len(data_dirs)==1;data_dir=next(iter(data_dirs))
        assert sha(data_dir/f'profile_{which}.json')==sha(profile)
        for part in(['tail']if x.quick else['body','exterior','tail']):
            premise=resolved(ledger['reports'][part]['path']);expected=read(premise)
            assert sha(premise)==ledger['reports'][part]['sha256']
            assert expected['profile']==which and expected['branch']==branch
            jobs.append((name+'_'+part,adapter,which,branch,part,data_dir,radius,expected,profile))
        profiles.append(name)
    with ThreadPoolExecutor(max_workers=x.workers)as pool:reports=list(pool.map(replay,jobs))
    assert all(sha(ROOT/rel)==h for rel,h in manifest['sha256'].items())
    out=dict(status='PASS_PORTABLE_ALL_ADDITIONAL_PROMOTED_RANK_GATES',
        mode='tails_only'if x.quick else'complete_body_exterior_exact_tail',
        seconds=time.time()-start,profiles=profiles,fresh_reports=reports,
        scope='Only the accepted additional rank-profile whole-q prerequisites. The separate global numeric and exact coverage entrypoints assemble the same-coefficient operator bound.')
    (OUT/'result.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))

if __name__=='__main__':main()
