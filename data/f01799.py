#!/usr/bin/env python3
"""Freeze and freshly test the final coefficient transport after the713 replay."""
from pathlib import Path
import argparse,gzip,hashlib,json,os,re,shutil,subprocess,sys,tempfile,time
import build_preferred_portable_bundle as B
ROOT=B.ROOT;PREFIX=str(ROOT)+'/'
EXCLUDED={'replay_outputs','replay_outputs_targeted','replay_outputs_promoted','__pycache__'}
PATTERN=re.compile(r'(?<![0-9a-f])[0-9a-f]{64}(?![0-9a-f])')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(gzip.decompress(p.read_bytes())if p.suffix=='.gz'else p.read_text())
def write(p,v):p.write_text(json.dumps(v,indent=2)+'\n')

def main(stage):
    stage=stage.resolve();start=time.time();result=read(stage/'replay_outputs/result.json')
    assert result['status']=='PASS_PORTABLE_COMPLETE_FULLUPPER_PREFERRED_GLOBAL'
    assert result['operator_b']=='713/1000'and result['norm_upper_bound']=='35245449/62500000'
    old_manifest=read(stage/'MANIFEST.json');assert all(sha(stage/r)==h for r,h in old_manifest['sha256'].items())
    shutil.copy2(stage/'replay_outputs/result.json',stage/'PORTABLE_713_COMPLETE_REPLAY_RESULT.json')
    shutil.copy2(stage/'replay_outputs/preferred_complete_acceptance.json',stage/'PORTABLE_713_COMPLETE_REPLAY_ACCEPTANCE.json')
    provenance=read(stage/'PROVENANCE.json');pref=provenance['preferred'];config=read(stage/'REPLAY_CONFIG.json')
    audit='tonight/audit/global_b7171/';theory='tonight/theory/';parent=audit+'FINAL29_GLOBAL_ACCEPTANCE.json'
    scan='tonight/recovery/profile_cover/assigned_twenty-nine_coefficient_scan.json'
    lower=audit+'b711_complete_lower_target_transport.json';spectral=audit+'final29_b711_beta_spectral_audit.json'
    transport=theory+'assigned_b711_shifted_cells_audit.json';accepted=audit+'FINAL_B711_GLOBAL_ACCEPTANCE.json'
    jobs=[{'name':'final711_same_b_complete_lower_branch','args':[audit+'replay_lower_target_transport_exact_early.py','--b','711/1000','--output',lower],'bits':768},
          {'name':'final711_exact43_profile_beta_spectral_bindings','args':[audit+'replay_assigned_beta_bindings.py','--scan',scan,'--parent-global',parent,'--b','711/1000','--output',spectral],'bits':768},
          {'name':'final711_all272579_uniform_cell_transports','args':[theory+'audit_assigned_b711_shifted_cells.py'],'bits':768}]
    assembly=[audit+'assemble_b711_global_acceptance.py','--parent-global',parent,'--transport',transport,'--spectral',spectral,'--lower',lower,'--output',accepted]
    config['jobs'].append({'name':'complete_reference713_final29_global_assembly','args':config['assembly_args'],'bits':768})
    config['jobs'].extend(jobs);config['assembly_args']=assembly;config['reference_envelope_operator_b']='713/1000'
    pending=list(map(Path,[accepted,scan,lower,spectral,transport,*[j['args'][0]for j in jobs],assembly[0]]));seen=set();new=[]
    numeric=re.compile(r'^[+\-]?[0-9][0-9.eE+\-]*(?:/[0-9.eE+\-]+)?$')
    while pending:
        rel=pending.pop()
        if rel in seen:continue
        seen.add(rel)
        if(stage/rel).is_file():continue
        p=ROOT/rel;assert p.is_file(),rel;target=stage/rel;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,target);new.append(str(rel))
        if p.suffix=='.json'or p.name.endswith('.json.gz'):
            for value in B.strings(read(p)):
                if numeric.fullmatch(value)or value.startswith(('native:','R3:','native-contact:','R3-contact:','inherited:')):continue
                dep=B.candidate(value,p.parent)
                if dep is not None and dep not in seen:pending.append(dep)
    originals={};texts={};current={};ordered={}
    for rel in new:
        p=stage/rel
        if p.suffix=='.json':originals[rel]=p.read_bytes();texts[rel]=p.read_text();current[rel]=texts[rel]
        elif p.suffix in{'.py','.md','.tex','.txt'}:p.write_text(p.read_text().replace(PREFIX,''))
    oldmap={v['original_sha256']:v['portable_sha256']for v in pref['archived_operational_json']}
    rawhash={r:hashlib.sha256(v).hexdigest()for r,v in originals.items()}
    for iteration in range(24):
        mapping={**oldmap,**{rawhash[r]:sha(stage/r)for r in originals if sha(stage/r)!=rawhash[r]}};changed=0
        for rel in originals:
            pairs=[]
            def sub(m):
                before=m.group(0);after=mapping.get(before,before)
                if after!=before:pairs.append((before,after))
                return after
            text=PATTERN.sub(sub,texts[rel].replace(PREFIX,''));ordered[rel]=([(PREFIX,'')]if PREFIX in texts[rel]else[])+list(dict.fromkeys(pairs))
            if text!=current[rel]:(stage/rel).write_text(text);current[rel]=text;changed+=1
        if not changed:break
    else:raise RuntimeError('711 metadata hashes failed to stabilize')
    for rel,raw in originals.items():
        if not ordered[rel]:continue
        archive=stage/'original_data'/(rawhash[rel]+'.bin');archive.write_bytes(raw)
        transformed=texts[rel]
        for before,after in ordered[rel]:transformed=transformed.replace(before,after)
        assert json.loads(transformed)==read(stage/rel)
        pref['archived_operational_json'].append({'path':rel,'archive':str(archive.relative_to(stage)),'original_sha256':rawhash[rel],'portable_sha256':sha(stage/rel),'ordered_replacements':ordered[rel]})
    final=read(ROOT/accepted);assert final['status']=='PASS_INDEPENDENT_COMPLETE_B711_GLOBAL_SCALAR_CERTIFICATE'
    pref.update(operator_b='711/1000',acceptance=accepted,norm_cap=final['norm_upper_bound'],GR_lower_bound=final['constructive_G_R_lower_bound'],reference_envelope_operator_b='713/1000')
    provenance['proof_scope']='Complete frozen713 native-kappa68/29-profile proof followed by a uniform SAME coefficient711 transport of every selected upper cell and a separately fresh711 complete lower branch. The713 envelope is only a proof baseline.'
    write(stage/'PROVENANCE.json',provenance);write(stage/'REPLAY_CONFIG.json',config)
    shutil.copy2(ROOT/'tonight/audit/global_b7171/next_verify_bundle_template.py',stage/'verify_bundle.py')
    # Metadata describes the new target; every713 mathematical premise remains unchanged.
    checkpoint=read(stage/'BASE_REPLAY_CHECKPOINT.json')
    for name in('PROVENANCE.json','README.md','VALIDATION.json'):checkpoint['unchanged_base_inputs'].pop(name,None)
    checkpoint['baseline_provenance']={k:v for k,v in provenance.items()if k not in('preferred','proof_scope')}
    write(stage/'BASE_REPLAY_CHECKPOINT.json',checkpoint)
    text=(stage/'README.md').read_text().replace('N(713/1000)<35245449/62500000=0.563927184','N(711/1000)<112782853/200000000=0.563914265').replace('62500000/35245449=1.7732785869744488','200000000/112782853=1.7733192119195637').replace('1/N(713/1000)','1/N(711/1000)')
    text+='\nThe frozen envelope uses coefficient713/1000 as a proof baseline. The final four jobs apply one coefficient711/1000 uniformly to all272579 selected upper cells, verify43 exact profile/beta and target spectral bindings, reprice the complete same-coefficient lower branch, and assemble the final711 certificate. No coefficient choices are mixed.\n'
    (stage/'README.md').write_text(text)
    files=sorted(p for p in stage.rglob('*')if p.is_file()and not any(k in EXCLUDED for k in p.relative_to(stage).parts)and p.name!='MANIFEST.json')
    manifest={'status':'FROZEN_FINAL711_COMPLETE_PREFERRED_GLOBAL_INPUT_MANIFEST','extension_filter_used':False,'excluded_mutable_directories':sorted(EXCLUDED),'sha256':{str(p.relative_to(stage)):sha(p)for p in files}};write(stage/'MANIFEST.json',manifest)
    out=stage/'replay_outputs'
    with tempfile.TemporaryDirectory(prefix='preferred711_final_delta_')as temporary:
        work=Path(temporary)/'proof';shutil.copytree(stage,work,ignore=shutil.ignore_patterns(*EXCLUDED))
        env=dict(os.environ);env['GR_MPFR_BITS']='768';env['PYTHONINTMAXSTRDIGITS']='0';env['PYTHONPATH']=os.pathsep.join(str(work/p)for p in('tonight/audit','tonight/theory','tonight/operator/projection/certificate'))
        for job in [*jobs,{'name':'complete_final711_global_assembly','args':assembly}]:
            args=[work/v if v.startswith(('tonight/','reduced_search/'))else v for v in job['args']]
            with(out/(job['name']+'.log')).open('w')as log:r=subprocess.run([sys.executable,*map(str,args)],cwd=work,env=env,stdout=log,stderr=subprocess.STDOUT)
            assert r.returncode==0,(job['name'],'see replay_outputs log')
            print(json.dumps({'stage':'PASS','job':job['name']}),flush=True)
        fresh=read(work/accepted);assert fresh['norm_upper_bound']==pref['norm_cap']and fresh['operator_b']=='711/1000'
        for name,binding in fresh['bindings'].items():shutil.copy2(work/binding['path'],out/('preferred_'+name+'.json'))
        shutil.copy2(work/accepted,out/'preferred_complete_acceptance.json')
    assert all(sha(stage/r)==h for r,h in manifest['sha256'].items())
    summary={**result,'operator_b':'711/1000','norm_upper_bound':fresh['norm_upper_bound'],'constructive_G_R_lower_bound':fresh['constructive_G_R_lower_bound'],
        'reference_envelope_operator_b':'713/1000','final_uniform_target_coefficient_transport_freshly_replayed':True,
        'reference713_complete_portable_replay_result_sha256':sha(stage/'PORTABLE_713_COMPLETE_REPLAY_RESULT.json'),
        'uniform711_delta_seconds':time.time()-start,'seconds':result['seconds']+time.time()-start,
        'mode':'complete_fresh713_portable_global_replay_plus_fresh_uniform711_all_cell_and_lower_transport',
        'scope':'Complete single-coefficient711 global scalar norm and constructive lower bound. The713 envelope is only an independently accepted proof baseline; no three-digit enclosure of unrestricted G_R is asserted.'}
    write(out/'result.json',summary);print(json.dumps(summary,indent=2))

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--staging',type=Path,required=True);main(a.parse_args().staging)
