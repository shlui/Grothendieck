#!/usr/bin/env python3
"""Add the two frozen SAME-b monotone deltas to the completed portable base."""
from pathlib import Path
import argparse, gzip, hashlib, json, re, shutil
import build_preferred_portable_bundle as B

ROOT=B.ROOT; PREFIX=str(ROOT)+'/'
EXCLUDED={'replay_outputs','replay_outputs_targeted','replay_outputs_promoted','__pycache__'}
PATTERN=re.compile(r'(?<![0-9a-f])[0-9a-f]{64}(?![0-9a-f])')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(gzip.decompress(p.read_bytes()) if p.suffix=='.gz' else p.read_text())
def write(p,v):p.write_text(json.dumps(v,indent=2)+'\n')
def save(p,t):
    if p.suffix=='.gz':
        with p.open('wb') as f:
            with gzip.GzipFile(filename='',mode='wb',fileobj=f,mtime=0,compresslevel=3)as g:g.write(t.encode())
    else:p.write_text(t)

def main(out):
    out=out.resolve(); provenance=read(out/'PROVENANCE.json'); pref=provenance['preferred']
    assert pref['norm_cap']=='281976359/500000000'
    config=read(out/'REPLAY_CONFIG.json')
    config['jobs'].append({'name':'complete_native68_global_assembly','args':config['assembly_args'],'bits':768})
    base='tonight/recovery/profile_cover/'
    audit='tonight/audit/global_b7171/'
    theory='tonight/theory/'
    parent=base+'global_upper_b713_native_kappa68_refined_envelope.json.gz'
    parent_geometry=theory+'fullupper_b713_native_kappa68_refinement_audit.json'
    parent_global=audit+'FULLUPPER68_NATIVE_REFINED_GLOBAL_ACCEPTANCE.json'
    seeds=[Path(audit+'assemble_monotone_global_acceptance.py')]
    for tag,accepted in [('twentyeight','FULLUPPER68_TWENTYEIGHTBANK_GLOBAL_ACCEPTANCE.json'),('twentynine','FINAL29_GLOBAL_ACCEPTANCE.json')]:
        child=base+f'global_upper_b713_native_kappa68_{tag}bank_envelope.json.gz'
        geometry=theory+f'fullupper_b713_native_kappa68_{tag}bank_monotone_audit.json'
        assembly=[audit+'assemble_monotone_global_acceptance.py','--parent-global',parent_global,'--envelope',child,'--geometry',geometry,'--output',audit+accepted]
        config['jobs'].append({'name':f'final_{tag}_fresh_monotone_delta','args':[theory+'audit_monotone_R3_upgrade_fresh.py',parent,child,'--parent-audit',parent_geometry,'--output',geometry],'bits':768})
        if tag=='twentyeight':config['jobs'].append({'name':'complete_twentyeight_global_assembly','args':assembly,'bits':768})
        else:config['assembly_args']=assembly
        seeds.extend(map(Path,(child,geometry,audit+accepted)))
        parent=child;parent_geometry=geometry;parent_global=audit+accepted
    config['geometry_nodes']+=2
    pending=seeds[:];seen=set();new=[]
    numeric=re.compile(r'^[+\-]?[0-9][0-9.eE+\-]*(?:/[0-9.eE+\-]+)?$')
    while pending:
        rel=pending.pop()
        if rel in seen:continue
        seen.add(rel)
        if(out/rel).is_file():continue
        source=ROOT/rel;assert source.is_file(),rel
        target=out/rel;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(source,target);new.append(str(rel))
        if source.suffix=='.json'or source.name.endswith('.json.gz'):
            for value in B.strings(read(source)):
                if numeric.fullmatch(value)or value.startswith(('native:','R3:','native-contact:','R3-contact:','inherited:')):continue
                dependency=B.candidate(value,source.parent)
                if dependency is not None and dependency not in seen:pending.append(dependency)
    print(json.dumps({'stage':'final_delta_dependencies_copied','files':len(new)}),flush=True)
    originals={}; texts={};current={};ordered={}
    for rel in new:
        p=out/rel
        if p.suffix=='.json'or p.name.endswith('.json.gz'):
            originals[rel]=p.read_bytes();texts[rel]=gzip.decompress(originals[rel]).decode()if p.suffix=='.gz'else originals[rel].decode()
            current[rel]=texts[rel];ordered[rel]=[]
        elif p.suffix in{'.py','.md','.tex','.txt'}:
            p.write_text(p.read_text().replace(PREFIX,''))
    oldmap={v['original_sha256']:v['portable_sha256'] for v in pref['archived_operational_json']}
    rawhash={r:hashlib.sha256(v).hexdigest()for r,v in originals.items()}
    for iteration in range(32):
        mapping={**oldmap,**{rawhash[r]:sha(out/r)for r in originals if sha(out/r)!=rawhash[r]}}
        changed=0
        for rel in originals:
            pairs=[]
            def substitute(m):
                before=m.group(0);after=mapping.get(before,before)
                if after!=before:pairs.append((before,after))
                return after
            text=PATTERN.sub(substitute,texts[rel].replace(PREFIX,''))
            ordered[rel]=([(PREFIX,'')]if PREFIX in texts[rel]else[])+list(dict.fromkeys(pairs))
            if text!=current[rel]:save(out/rel,text);current[rel]=text;changed+=1
        print(json.dumps({'stage':'final_delta_metadata','iteration':iteration,'changed':changed}),flush=True)
        if not changed:break
    else:raise RuntimeError('Final delta metadata hashes did not stabilize')
    catalog=pref['archived_operational_json']
    # Retain the previous principal original source before replacing its name.
    previous=out/pref['original_source_archive']; dest=out/'original_data'/(pref['original_source_sha256']+'.bin')
    if not dest.exists():shutil.copy2(previous,dest)
    for item in catalog:
        if item['archive']==pref['original_source_archive']:item['archive']=str(dest.relative_to(out))
    for rel,raw in originals.items():
        if not ordered[rel]:continue
        digest=rawhash[rel]; archive=out/'original_data'/(digest+'.bin');archive.write_bytes(raw)
        transformed=texts[rel]
        for before,after in ordered[rel]:transformed=transformed.replace(before,after)
        assert json.loads(transformed)==read(out/rel),rel
        catalog.append({'path':rel,'archive':str(archive.relative_to(out)),'original_sha256':digest,'portable_sha256':sha(out/rel),'ordered_replacements':ordered[rel]})
    accepted=read(ROOT/parent_global);source=ROOT/parent;shutil.copy2(source,out/'original_preferred_envelope.bin')
    assert accepted['norm_upper_bound']=='35245449/62500000'
    if not any(v['path']==parent for v in catalog):catalog.append({'path':parent,'archive':'original_preferred_envelope.bin','original_sha256':sha(source),'portable_sha256':sha(out/parent),'ordered_replacements':[]})
    pref.update(acceptance=parent_global,source=parent,original_source_sha256=sha(source),portable_source_sha256=sha(out/parent),norm_cap=accepted['norm_upper_bound'],GR_lower_bound=accepted['constructive_G_R_lower_bound'])
    provenance['proof_scope']='Complete accepted SAME-b=.713 native-kappa68 global certificate followed by the frozen twenty-eight and twenty-nine profile monotone deltas; no later research premises.'
    write(out/'PROVENANCE.json',provenance);write(out/'REPLAY_CONFIG.json',config)
    files=sorted(p for p in out.rglob('*')if p.is_file()and not any(k in EXCLUDED for k in p.relative_to(out).parts)and p.name!='MANIFEST.json')
    write(out/'MANIFEST.json',{'status':'FROZEN_FINAL29_PREFERRED_GLOBAL_INPUT_MANIFEST','extension_filter_used':False,'excluded_mutable_directories':sorted(EXCLUDED),'sha256':{str(p.relative_to(out)):sha(p)for p in files}})
    print(json.dumps({'status':'FINAL29_PORTABLE_DELTA_APPENDED','new_files':len(new),'jobs':len(config['jobs']),'source_original_sha256':sha(source)},indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
