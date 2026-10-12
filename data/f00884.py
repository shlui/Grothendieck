#!/usr/bin/env python3
"""Freeze the goal lower certificate and exact universal rounding replay graph."""
from pathlib import Path
from fractions import Fraction as F
import argparse,gzip,hashlib,json
R=Path(__file__).resolve().parents[3];N=R/'new_round_2026-10-02';H=N/'audit';P=Path(__file__).resolve().parent;O=R/'verified_2026-10-01/tonight/operator/preferred_global'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):
    with (gzip.open(p,'rt') if p.suffix=='.gz' else p.open()) as f:return json.load(f)
def rel(p):return str(p.relative_to(R))
def resolve(p):
    p=Path(p);p=p if p.is_absolute() else R/p
    if not p.exists() and p.is_relative_to(R/'tonight'):p=O/p.relative_to(R)
    return p
def pairs(v):
    if isinstance(v,dict):
        if isinstance(v.get('path'),str) and isinstance(v.get('sha256'),str):yield resolve(v['path']),v['sha256']
        for x in v.values():yield from pairs(x)
    elif isinstance(v,list):
        for x in v:yield from pairs(x)
ap=argparse.ArgumentParser();ap.add_argument('--acceptance',type=Path,required=True);ap.add_argument('--tag',required=True);a=ap.parse_args()
accept=a.acceptance.resolve();final=read(accept);assert final['status'].startswith('PASS_INDEPENDENT_COMPLETE_') and final['operator_b']=='711/1000'
assert F(final['norm_upper_bound'])<=F(40,71),'Goal builder requires an independently accepted complete lower bound at least71/40'
bank=final['complete_whole_q_profiles'];tag=a.tag;assert tag.replace('_','').isalnum()
cp=P/(tag+'_REPLAY_CONFIG.json');mp=P/(tag+'_INPUT_MANIFEST.json');assert not cp.exists() and not mp.exists()
layers=[];p=accept
while True:
    d=read(p);assert d['operator_b']==final['operator_b']
    if 'accepted43_u_global' in d['bindings']:
        layers.append((p,d));break
    layers.append((p,d))
    predecessors=[v for k,v in d['bindings'].items() if k.startswith('accepted') and (k.endswith('_global') or k.endswith('_full_peak_global'))]
    assert len(predecessors)==1,(p,list(d['bindings']))
    p=resolve(predecessors[0]['path']);assert sha(p)==predecessors[0]['sha256']
layers.reverse()
corepath,core=layers[0]
corecache=resolve(core['bindings']['fresh_R3_prices_and_local_gates']['path'])
basecache=H/'numeric/BANK43_LINEAR_TENSOR_R3_CACHE_REPLAY.json'
strong=N/'cubic_search/b711_bank43_strong_shape_rank_envelope_proposal.json.gz';urefine=N/'cubic_search/b711_bank43_u0005_refinement_proposal.json.gz'
scripts={}
for q in (H/'numeric').glob('*.py'):scripts.setdefault(sha(q),[]).append(q)
jobs=[];visited=set()
def numeric(p):
    p=resolve(p)
    if p in visited:return
    d=read(p);assert d['status'] in {'PASS_INDEPENDENT_NEW_ROUND_COMPLETE_R3_CACHE_PRICES','PASS_INDEPENDENT_CURVED_B0_UNIFORM_SOURCE_PLANES'}
    visited.add(p)
    for dep,digest in pairs(d):
        assert sha(dep)==digest
        if dep.is_relative_to(H/'numeric') and dep.suffix=='.json':
            dd=read(dep)
            if dd.get('status') in {'PASS_INDEPENDENT_NEW_ROUND_COMPLETE_R3_CACHE_PRICES','PASS_INDEPENDENT_CURVED_B0_UNIFORM_SOURCE_PLANES'}:numeric(dep)
    candidates=sorted(scripts[d['script_sha256']]);assert candidates
    result='{numeric}/'+p.stem+'.json'
    jobs.append({'name':p.stem,'script':rel(candidates[0]),'args':[rel(resolve(d['source']['path'])),'--output',result],'result':result})
numeric(basecache);numeric(H/'numeric/BANK43_LINEAR_TENSOR_B0_PLANE_REPLAY.json')
local_binding_reports=[basecache]
for _,layer in layers:
    for k,rec in layer['bindings'].items():
        if k in ('fresh_R3_prices_and_local_gates','shared_moment_R3_prices_and_local_gates','shared_moment_B0_planes'):
            q=resolve(rec['path']);numeric(q)
            if k!='shared_moment_B0_planes':local_binding_reports.append(q)
jobs.append({'name':'native_curved','script':rel(H/'numeric/replay_curved_native_prices.py'),'args':['--source',rel(N/'native_curved_endpoint_prices.json'),'--output','{numeric}/native_curved.json'],'result':'{numeric}/native_curved.json'})
jobs.extend([
 {'name':'strong43_geometry','script':rel(N/'sequence_theory/audit_strong_shape_geometry.py'),'args':[rel(strong),'{output}/strong43_geometry.json'],'result':'{output}/strong43_geometry.json'},
 {'name':'u_refinement_geometry','script':rel(N/'sequence_theory/audit_u_refinement_geometry.py'),'args':[rel(urefine),'{output}/u_refinement_geometry.json','--parent-report','{output}/strong43_geometry.json'],'result':'{output}/u_refinement_geometry.json'}])
def geometry(record,name,parent):
    p=resolve(record['path']);assert sha(p)==record['sha256'];d=read(p)
    candidates=sorted(q for q in (N/'sequence_theory').glob('*.py') if sha(q)==d['entrypoint_sha256']);assert candidates
    result='{output}/'+name+'.json'
    jobs.append({'name':name,'script':rel(candidates[0]),'args':[rel(resolve(d['source']['path'])),result,'--parent-report',parent],'result':result})
prior='{output}/u_refinement_geometry.json'
for i,(_,layer) in enumerate(layers):
    keys=[k for k in ('monotone_rank_geometry','shared_moment_geometry','affine_geometry','affine_moment_geometry') if k in layer['bindings']];assert len(keys)==1
    name='final_geometry' if i==len(layers)-1 else f'layer_{i:02d}_geometry'
    geometry(layer['bindings'][keys[0]],name,prior);prior='{output}/'+name+'.json'
jobs.append({'name':'three_digit_rounding','script':rel(H/'replay_elementary_upper.py'),'args':['--output','{output}/three_digit_rounding.json','--lower-acceptance',rel(accept)],'result':'{output}/three_digit_rounding.json'})
gates=[];seen=set()
for rp in local_binding_reports:
    for item in read(rp)['complete_whole_q_bindings'].values():
        item=item.get('whole_q_premise',item)
        if 'accepted_component' not in item:continue
        cpath=resolve(item['accepted_component']['path']);c=read(cpath);w=read(resolve(c['whole_q_binding']['path']));entry=item['complete_component']
        key=(entry['profile']['sha256'],entry['reserve0'],entry['reserve6'])
        if key in seen:continue
        seen.add(key)
        tail=read(resolve(next(v for v in entry['whole_q_reports'] if v['part']=='tail')['path']));assert tail['status']=='PASS_INDEPENDENT_EXACT_R3_TAIL'
        for rec in entry['whole_q_reports']:
            d=read(resolve(rec['path']));label=d['profile'];branch=d['branch'];part=rec['part'];result='{output}/whole_q/independent_'+label+'_'+branch+'_'+part+'.json'
            gates.append({'name':label+'_'+part,'script':rel(resolve(w['adapter']['path'])),'args':['--profile',label,'--branch',branch,'--part',part,'--tail-R',tail['R'],'--data-dir',rel(cpath.parent/'input'),'--output-dir','{output}/whole_q'],'result':result,'expected_status':d['status'],'accepted':d.get('accepted')})
assert len(seen)==bank-29 and len(gates)==3*(bank-29)
cfg={'status':'FROZEN_PORTABLE_NEW_ROUND_REPLAY_CONFIG','original_root':str(R),'accepted_global':rel(accept),'operator_b':final['operator_b'],'norm_cap':final['norm_upper_bound'],'norm_cap_relation':final.get('norm_cap_relation','STRICT_BELOW_CLEAN_CAP'),
 'constructive_G_R_lower_bound':final['constructive_G_R_lower_bound'],'rounding_goal':'1.78','universal_strict_upper':'1760/987','baseline':rel(O),'jobs':jobs,'all_local_new_gates':gates,
 'mode_default':'Accepted immutable local gates; fresh complete numeric dependency graph; fresh exhaustive geometry, child partitions and protected baseline default replay.',
 'all_local_mode':f'Also reconstruct all {bank-29} current new whole-q criteria and the full protected baseline local/scalar mode.',
 'immutability':'Certificate bytes and original absolute paths are preserved; outputs are separate. Every criterion uses its own certified tail cutoff and literal budget.'}
cp.write_text(json.dumps(cfg,indent=2)+'\n')
files={};pending=[]
def add(p,h=None):
    p=resolve(p);assert p.is_file() and p.is_relative_to(R),p
    v=sha(p);assert h is None or h==v,('HASH',p);rr=rel(p)
    if rr not in files:files[rr]=v;pending.append(p)
def visit(v):
    if isinstance(v,dict):
        if isinstance(v.get('path'),str) and isinstance(v.get('sha256'),str):add(v['path'],v['sha256'])
        if isinstance(v.get('source'),str) and isinstance(v.get('alias'),str) and isinstance(v.get('sha256'),str):add(v['source'],v['sha256']);add(v['alias'],v['sha256'])
        for k in ('profile','acceptance'):
            if isinstance(v.get(k+'_path'),str) and isinstance(v.get(k+'_sha256'),str):add(v[k+'_path'],v[k+'_sha256'])
        for p,h in v.get('verified_premise_files',{}).items():add(p,h)
        for x in v.values():visit(x)
    elif isinstance(v,list):
        for x in v:visit(x)
for rr,h in read(O/'MANIFEST.json')['sha256'].items():add(O/rr,h)
for p in (O/'MANIFEST.json',accept,cp,Path(__file__),P/'verify_goal_latest_round.py',P/'Dockerfile',P/'README.md',P/'README_GOAL_REPLAY.md',H/'ELEMENTARY_KRIVINE_UPPER_INDEPENDENT_AUDIT.json',*visited):add(p)
for job in jobs+gates:
    add(R/job['script'])
    for arg in job['args']:
        if arg.startswith('new_round_2026-10-02/'):
            p=R/arg
            if p.is_file():add(p)
            elif p.is_dir():
                for q in p.iterdir():
                    if q.is_file():add(q)
while pending:
    p=pending.pop()
    if p.is_relative_to(O):continue
    if p.suffix=='.json' or p.name.endswith('.json.gz'):visit(read(p))
manifest={'status':'FROZEN_COMPLETE_NEW_ROUND_INPUT_WHITELIST','config':rel(cp),'extension_filter_used':False,'sha256':dict(sorted(files.items())),'file_count':len(files),'total_bytes':sum((R/rr).stat().st_size for rr in files),'excluded_mutable_outputs':'Fresh outputs, environments and unbound caches. Every explicit immutable input is included regardless of extension.'}
mp.write_text(json.dumps(manifest,indent=2)+'\n')
print(json.dumps({'config':rel(cp),'manifest':rel(mp),'sha256':sha(mp),'files':len(files),'bytes':manifest['total_bytes'],'numeric_and_geometry_jobs':len(jobs),'local_parts':len(gates)},indent=2))
