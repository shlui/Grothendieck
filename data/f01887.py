#!/usr/bin/env python3
"""Reconstruct every new metadata link from immutable originals in one graph."""
from pathlib import Path
import argparse,gzip,hashlib,json,re
import build_preferred_portable_bundle as B
ROOT=B.ROOT;PREFIX=str(ROOT)+'/'
EXCLUDED={'replay_outputs','replay_outputs_targeted','replay_outputs_promoted','__pycache__'}
PATTERN=re.compile(r'"([0-9a-f]{64})"')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(gzip.decompress(p.read_bytes())if p.suffix=='.gz'else p.read_text())
def write(p,v):p.write_text(json.dumps(v,indent=2)+'\n')
def save(p,t):
    if p.suffix=='.gz':
        with p.open('wb')as f:
            with gzip.GzipFile(filename='',mode='wb',fileobj=f,mtime=0,compresslevel=1)as g:g.write(t.encode())
    else:p.write_text(t)
def main(stage):
    stage=stage.resolve();base=read(B.BASE/'MANIFEST.json')['sha256'];prov=read(stage/'PROVENANCE.json');pref=prov['preferred']
    protected=list(map(Path,prov['primitive_packages']))+[Path('reduced_search/energy/633')]
    originals={};texts={};ordered={};current={};all_original_sha={}
    for p in stage.rglob('*'):
        if not p.is_file()or any(k in EXCLUDED for k in p.relative_to(stage).parts):continue
        rel=str(p.relative_to(stage));source=ROOT/rel
        if not source.is_file()or not(p.suffix=='.json'or p.name.endswith('.json.gz')):continue
        all_original_sha[rel]=sha(source)
        if rel in base or any(Path(rel).is_relative_to(q)for q in protected):continue
        originals[rel]=source.read_bytes();texts[rel]=gzip.decompress(originals[rel]).decode()if p.suffix=='.gz'else originals[rel].decode()
        current[rel]=gzip.decompress(p.read_bytes()).decode()if p.suffix=='.gz'else p.read_text()
    print(json.dumps({'stage':'entire_graph_originals_loaded','new_json_files':len(originals)}),flush=True)
    for iteration in range(32):
        mapping={h:sha(stage/r)for r,h in all_original_sha.items()if sha(stage/r)!=h};changed=0
        for rel in originals:
            pairs=[]
            def sub(m):
                before=m.group(1);after=mapping.get(before,before)
                if after!=before:pairs.append((before,after))
                return '"'+after+'"'
            text=PATTERN.sub(sub,texts[rel].replace(PREFIX,''))
            ordered[rel]=([(PREFIX,'')]if PREFIX in texts[rel]else[])+list(dict.fromkeys(pairs))
            if text!=current[rel]:save(stage/rel,text);current[rel]=text;changed+=1
        print(json.dumps({'stage':'entire_graph_hash_propagation','iteration':iteration,'changed':changed}),flush=True)
        if not changed:break
    else:raise RuntimeError('Entire ancestor graph failed to stabilize')
    catalog=[v for v in pref['archived_operational_json']if v['path']not in originals]
    for rel,raw in originals.items():
        if not ordered[rel]:continue
        digest=all_original_sha[rel];archive=stage/'original_data'/(digest+'.bin');archive.write_bytes(raw)
        transformed=texts[rel]
        for before,after in ordered[rel]:transformed=transformed.replace(before,after)
        assert json.loads(transformed)==read(stage/rel),('EXACT_STRUCTURED_REWRITE',rel)
        catalog.append({'path':rel,'archive':str(archive.relative_to(stage)),'original_sha256':digest,'portable_sha256':sha(stage/rel),'ordered_replacements':ordered[rel]})
    pref['archived_operational_json']=catalog;pref['portable_source_sha256']=sha(stage/pref['source'])
    write(stage/'PROVENANCE.json',prov)
    # The already freshly checked old base does not use the newly repaired
    # full-upper ancestor metadata. Preserve only its actual immutable inputs.
    checkpoint=read(stage/'BASE_REPLAY_CHECKPOINT.json')
    checkpoint['unchanged_base_inputs']={r:h for r,h in checkpoint['unchanged_base_inputs'].items()if(r in base or any(Path(r).is_relative_to(q)for q in protected))and r not in('PROVENANCE.json','README.md','VALIDATION.json','verify_bundle.py')}
    checkpoint['baseline_provenance']={k:v for k,v in prov.items()if k not in('preferred','proof_scope')}
    assert all(sha(stage/r)==h for r,h in checkpoint['unchanged_base_inputs'].items())
    write(stage/'BASE_REPLAY_CHECKPOINT.json',checkpoint)
    def records(v):
        if isinstance(v,dict):
            if isinstance(v.get('path'),str)and isinstance(v.get('sha256'),str):yield v
            for key,value in v.items():
                if key not in('cells','physically_empty_cells','cached_corner_constants','endpoint_prices'):yield from records(value)
        elif isinstance(v,list):
            for value in v:yield from records(value)
    links=0
    for rel in originals:
        data=read(stage/rel)
        for item in records(data):
            path=Path(item['path']);path=stage/path if not path.is_absolute()else path
            if path.is_file():assert sha(path)==item['sha256'],('FILE_BINDING',rel,item['path']);links+=1
        if isinstance(data,dict):
            for key in('monotone_upgrade','adaptive_amplitude_refinement','native_kappa_refinement'):
                if key in data:
                    binding=data[key];path=stage/binding['parent_path'];assert sha(path)==binding['parent_sha256'],('PARENT_SOURCE',rel);links+=1
    for item in catalog:
        raw=(stage/item['archive']).read_bytes();assert hashlib.sha256(raw).hexdigest()==item['original_sha256']
        text=gzip.decompress(raw).decode()if raw[:2]==b'\x1f\x8b'else raw.decode()
        for before,after in item['ordered_replacements']:text=text.replace(before,after)
        assert json.loads(text)==read(stage/item['path']),item['path']
    files=sorted(p for p in stage.rglob('*')if p.is_file()and not any(k in EXCLUDED for k in p.relative_to(stage).parts)and p.name!='MANIFEST.json')
    write(stage/'MANIFEST.json',{'status':'FROZEN_COMPLETE_REBOUND_FINAL29_INPUT_MANIFEST','extension_filter_used':False,'excluded_mutable_directories':sorted(EXCLUDED),'sha256':{str(p.relative_to(stage)):sha(p)for p in files}})
    print(json.dumps({'status':'PASS_ENTIRE_PORTABLE_ANCESTOR_GRAPH_REBOUND','checked_file_and_parent_links':links,'archived_exact_operational_rewrites':len(catalog),'byte_bound_base_checkpoint_inputs':len(checkpoint['unchanged_base_inputs'])},indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--staging',type=Path,required=True);main(p.parse_args().staging)
