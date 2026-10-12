#!/usr/bin/env python3
"""Read-only goal replay, including the exact universal upper and three-digit rounding."""
from pathlib import Path
from fractions import Fraction as F
from concurrent.futures import ThreadPoolExecutor
import argparse,hashlib,json,os,shutil,subprocess,sys,tempfile,time,uuid
sys.dont_write_bytecode=True
EXPECTED=Path('/workspace/scratch/8242d848994f')
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
ap=argparse.ArgumentParser();ap.add_argument('--manifest',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--materialize',type=Path);ap.add_argument('--manifest-only',action='store_true');ap.add_argument('--all-local',action='store_true');ap.add_argument('--workers',type=int,default=4);a=ap.parse_args()
m=read(a.manifest.resolve());assert m['extension_filter_used']is False
if a.materialize:
    for rr,h in m['sha256'].items():
        src=a.materialize.resolve()/rr;assert sha(src)==h
        dst=EXPECTED/rr
        if dst.exists():assert sha(dst)==h
        else:dst.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(src,dst)
R=EXPECTED;assert R.is_dir(),'Replay requires the recorded original root; use the documented Docker materialization.'
assert all(sha(R/rr)==h for rr,h in m['sha256'].items())
c=read(R/m['config']);assert c['original_root']==str(R)
out=a.output.resolve();out.mkdir(parents=True,exist_ok=True);start=time.time();results=[]
if a.manifest_only:
    print(json.dumps({'status':'PASS_PORTABLE_NEW_ROUND_IMMUTABLE_MANIFEST','files':len(m['sha256']),'operator_b':c['operator_b'],'norm_cap':c['norm_cap'],'mathematical_replay_performed':False},indent=2));sys.exit(0)
env=dict(os.environ,GR_FROZEN_BASELINE=str(R/c['baseline']),GR_MPFR_BITS='768',PYTHONDONTWRITEBYTECODE='1',PYTHONINTMAXSTRDIGITS='0')
numeric=R/'new_round_2026-10-02/audit/numeric'/('_fresh_replay_'+uuid.uuid4().hex);numeric.mkdir()
def arg(v):
    if '{output}'in v:return v.replace('{output}',str(out))
    if '{numeric}'in v:return v.replace('{numeric}',str(numeric))
    return str(R/v)if v.startswith(('new_round_2026-10-02/','verified_2026-10-01/'))else v
def run(job):
    log=out/(job['name']+'.log')
    with log.open('w')as f:r=subprocess.run([sys.executable,str(R/job['script']),*map(arg,job['args'])],cwd=R,env=env,stdout=f,stderr=subprocess.STDOUT)
    assert r.returncode==0,('FAILED_REPLAY',str(log))
    p=Path(arg(job['result']));d=read(p);assert d['status'].startswith('PASS'),(job['name'],d['status'])
    if 'expected_status'in job:assert d['status']==job['expected_status'] and d.get('accepted')==job['accepted']
    dst=out/(job['name']+'.json')
    if p!=dst:shutil.copy2(p,dst)
    q={'job':job['name'],'status':d['status'],'report':str(dst),'sha256':sha(dst)};print(json.dumps(q),flush=True);return q
# The original portable baseline is executed in a disposable copy, so even its
# mutable replay_outputs remain outside the frozen preferred proof package.
with tempfile.TemporaryDirectory(prefix='gr_frozen_baseline_replay_')as tmp:
    work=Path(tmp)/'proof';shutil.copytree(R/c['baseline'],work,ignore=shutil.ignore_patterns('replay_outputs*','__pycache__'))
    cmd=[sys.executable,str(work/'verify_bundle.py'),'--workers',str(a.workers)]
    if a.all_local:cmd.append('--full')
    with(out/'baseline.log').open('w')as f:r=subprocess.run(cmd,cwd=work,env=env,stdout=f,stderr=subprocess.STDOUT)
    assert r.returncode==0,('BASELINE_FAILED',str(out/'baseline.log'))
    results.append({'job':'protected_baseline','status':'PASS_FRESH_PORTABLE_BASELINE','all_local':a.all_local})
if a.all_local:
    with ThreadPoolExecutor(a.workers)as pool:results.extend(pool.map(run,c['all_local_new_gates']))
for job in c['jobs']:results.append(run(job))
g=read(out/'final_geometry.json');assert F(g['strict_upper_cap'])==F(c['norm_cap'])
if c['norm_cap_relation']=='LEQ_EXACT_TARGET_40_OVER_71':
    assert F(c['norm_cap'])==F(40,71) and F(g['maximum'])<=F(c['norm_cap']) and g['cap_relation']==c['norm_cap_relation']
else:assert F(g['maximum'])<F(c['norm_cap'])
accepted=read(R/c['accepted_global']);assert accepted['operator_b']==c['operator_b']and accepted['norm_upper_bound']==c['norm_cap']
upper=read(out/'three_digit_rounding.json');assert upper['three_significant_digit_rounding_proved']is True and upper['rounding_result']==c['rounding_goal']=='1.78'
assert F(c['norm_cap'])<=F(40,71) and F(upper['universal_strict_upper'])==F(c['universal_strict_upper'])<F(357,200)
assert F(accepted['lower_amplitude_upper_bound'])<=F(c['norm_cap'])
assert all(sha(R/rr)==h for rr,h in m['sha256'].items()),'Immutable input changed during replay'
summary={'status':'PASS_PORTABLE_COMPLETE_NEW_ROUND_GLOBAL_REPLAY','operator_b':c['operator_b'],'norm_upper_bound':c['norm_cap'],
 'constructive_G_R_lower_bound':c['constructive_G_R_lower_bound'],'norm_cap_relation':c['norm_cap_relation'],'immutable_input_files':len(m['sha256']),
 'all_local_replayed':a.all_local,'three_significant_digit_rounding':'1.78','universal_strict_upper':upper['universal_strict_upper'],'certificate_bytes_rewritten':False,'original_root_preserved':str(R),
 'seconds':time.time()-start,'jobs':results,'scope':c['all_local_mode']if a.all_local else c['mode_default']}
(out/'REPLAY_RESULT.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
shutil.rmtree(numeric)
