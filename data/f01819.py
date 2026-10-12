#!/usr/bin/env python3
"""Build a separate full-upper/adaptive preferred package from frozen audits."""
from pathlib import Path
from fractions import Fraction as F
import argparse
import gzip
import hashlib
import json
import re
import shutil
import build_preferred_portable_bundle as B

ROOT=B.ROOT;BASE=B.BASE;PREFIX=str(ROOT)+'/'
EXCLUDED={'replay_outputs','replay_outputs_targeted','replay_outputs_promoted','__pycache__'}
SHA_PATTERN=re.compile(r'(?<![0-9a-f])[0-9a-f]{64}(?![0-9a-f])')

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):
    with(gzip.open(p,'rt')if p.suffix=='.gz'else p.open())as f:return json.load(f)
def write(p,v):p.write_text(json.dumps(v,indent=2)+'\n')
def raw_json(p):
    data=p.read_bytes();return gzip.decompress(data).decode()if p.suffix=='.gz'else data.decode()
def json_write(p,text):
    if p.suffix=='.gz':
        with p.open('wb')as f:
            with gzip.GzipFile(filename='',mode='wb',fileobj=f,mtime=0)as g:g.write(text.encode())
    else:p.write_text(text)

def main(x):
    out=x.output.resolve();assert not out.exists()or x.resume,'Use an empty separate staging path or --resume'
    accepted=read(x.acceptance);config=read(x.config)
    assert accepted['status'].startswith('PASS_INDEPENDENT_COMPLETE_')
    b=F(accepted['operator_b']);assert 0<b<=F(7171,10000)
    source_rel=Path(accepted['envelope']['path']);source=ROOT/source_rel
    assert sha(source)==accepted['envelope']['sha256']
    baseline_manifest=read(BASE/'MANIFEST.json')
    assert all(sha(BASE/rel)==h for rel,h in baseline_manifest['sha256'].items())
    shutil.copytree(BASE,out,dirs_exist_ok=x.resume,ignore=shutil.ignore_patterns(*EXCLUDED))
    shutil.copy2(out/'verify_bundle.py',out/'verify_strengthened_components.py')
    packages=set(read(BASE/'PROVENANCE.json')['primitive_packages'])
    for item in read(source)['native_endpoint_inputs']:packages.add(item['package'])
    protected=[Path(p)for p in packages]+[Path('reduced_search/energy/633')]
    def copy_file(rel):
        target=out/rel;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(ROOT/rel,target)
    for rel in protected:
        for p in(ROOT/rel).rglob('*'):
            if p.is_file()and'__pycache__'not in p.parts:copy_file(p.relative_to(ROOT))
    # Existing original/targeted whole-q replay helpers form the base; new
    # aliases, literal-reserve adapters and actual trees enter transitively.
    r3=ROOT/'tonight/operator/projection/replay_accepted_R3_with_targeted'
    for rel,h in read(r3/'MANIFEST.json')['sha256'].items():
        assert sha(r3/rel)==h
        if rel.startswith(('tonight/','reduced_search/')):
            if(ROOT/rel).is_file():copy_file(rel)
            else:
                target=out/rel;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(r3/rel,target)
    if(r3/'verify_b_transport_box.py').is_file():shutil.copy2(r3/'verify_b_transport_box.py',out/'verify_b_transport_box.py')
    templates={'verify_bundle.py':'next_verify_bundle_template.py',
        'verify_targeted_prerequisites.py':'preferred_targeted_verify_template.py',
        'verify_promoted_prerequisites.py':'next_verify_promoted_ranks_template.py'}
    for target,template in templates.items():shutil.copy2(ROOT/'tonight/audit/global_b7171'/template,out/target)
    seeds=[x.acceptance.resolve().relative_to(ROOT),source_rel,
        Path('tonight/audit/global_b7171/replay_generic_R3_cache_57.py'),
        Path('tonight/audit/global_b7171/replay_middlebank_R3_corners.py'),
        Path('tonight/audit/global_b7171/replay_robust_outer_boxes.py'),
        Path('tonight/theory/audit_general_r3_source_bands.py'),
        Path('tonight/theory/audit_general_r3_parameters_tail.py'),
        Path('tonight/theory/audit_global_native_R3_envelope.py'),
        Path('tonight/recovery/profile_cover/kappa_refinement_4_10/complete.json'),
        Path('tonight/recovery/profile_cover/kappa_refinement_10_40/complete.json')]
    for value in B.strings(config):
        dependency=B.candidate(value,x.config.resolve().parent)
        if dependency is not None:seeds.append(dependency)
    pending=seeds[:];seen=set();candidate_cache={}
    numeric=re.compile(r'^[+\-]?[0-9][0-9.eE+\-]*(?:/[0-9.eE+\-]+)?$')
    def candidate(value,parent):
        if not isinstance(value,str)or numeric.fullmatch(value):return None
        if value.startswith(('native:','R3:','native-contact:','R3-contact:','inherited:')):return None
        if '\n'in value or len(value)>2048:return None
        if not value.startswith((PREFIX,'tonight/','reduced_search/','../'))and'.'not in value:return None
        key=(str(parent),value)
        if key not in candidate_cache:candidate_cache[key]=B.candidate(value,parent)
        return candidate_cache[key]
    while pending:
        rel=pending.pop()
        if rel in seen:continue
        seen.add(rel);p=ROOT/rel;assert p.is_file(),rel;copy_file(rel)
        if p.suffix=='.json'or p.name.endswith('.json.gz'):
            for value in B.strings(read(p)):
                dependency=candidate(value,p.parent)
                if dependency is not None and dependency not in seen:pending.append(dependency)
    print(json.dumps(dict(stage='dependencies_copied',files=len(seen))),flush=True)
    # Record raw original bytes before any portable operational metadata edit.
    originals={};original_text={};ordered={}
    for p in out.rglob('*'):
        if not p.is_file()or any(k in EXCLUDED for k in p.relative_to(out).parts):continue
        rel=p.relative_to(out)
        if any(rel.is_relative_to(q)for q in protected):continue
        if len(rel.parts)==1 and p.name.startswith('original_'):continue
        if p.name in{'MANIFEST.json','PROVENANCE.json','VALIDATION.json'}:continue
        if p.suffix=='.json'or p.name.endswith('.json.gz'):
            originals[str(rel)]=p.read_bytes();original_text[str(rel)]=raw_json(p);ordered[str(rel)]=[]
            if PREFIX in original_text[str(rel)]:
                json_write(p,original_text[str(rel)].replace(PREFIX,''));ordered[str(rel)].append((PREFIX,''))
        elif p.suffix in{'.py','.md','.tex','.txt'}:
            content=p.read_text()
            if PREFIX in content:p.write_text(content.replace(PREFIX,''))
    # Path-only adapters for old accepted replay entrypoints.
    corner=out/'tonight/audit/global_b7171/replay_R3_corners.py';content=corner.read_text()
    if 'if not path.is_absolute():path=ROOT/path'not in content:
        content=content.replace("path=Path(item['path']);assert sha(path)==item['sha256'];checked_ledger.append(file_record(path))",
            "path=Path(item['path'])\n        if not path.is_absolute():path=ROOT/path\n        assert sha(path)==item['sha256'];checked_ledger.append(file_record(path))")
        corner.write_text(content)
    diffuse=out/'tonight/audit/diffuse_b7171_R3/independent_source_pricing.py';content=diffuse.read_text()
    if 'if not path.is_absolute():path=ROOT/path'not in content:
        content=content.replace('path=Path(filename);assert sha(path)==digest',
            'path=Path(filename)\n        if not path.is_absolute():path=ROOT/path\n        assert sha(path)==digest')
        diffuse.write_text(content)
    for filename in('check_polygon_contact_credits.py','check_adaptive_contact_credits.py'):
        contact=out/'tonight/audit/global_b7171'/filename
        if contact.exists():
            content=contact.read_text()
            if 'r3={}\nif not args.r3_report:'not in content:
                start=content.index("r3=json.loads(");end=content.index('if args.r3_report:',start)
                content=content[:start]+"r3={}\nif not args.r3_report:\n "+content[start:end].rstrip('\n').replace('\n','\n ')+'\n'+content[end:]
                contact.write_text(content)
    current_text={rel:raw_json(out/rel)for rel in originals}
    raw_hashes={rel:hashlib.sha256(raw).hexdigest()for rel,raw in originals.items()}
    for iteration in range(32):
        mapping={raw_hashes[rel]:sha(out/rel)for rel in originals if sha(out/rel)!=raw_hashes[rel]}
        changes=0
        for rel in originals:
            p=out/rel;old=current_text[rel];pairs=[]
            def replace(m):
                before=m.group(0);after=mapping.get(before,before)
                if after!=before:pairs.append((before,after))
                return after
            # Rebuild from original text every pass, so a dependency whose
            # portable SHA changes propagates through all ancestor references.
            new=SHA_PATTERN.sub(replace,original_text[rel].replace(PREFIX,''))
            ordered[rel]=([(PREFIX,'')]if PREFIX in original_text[rel]else[])+list(dict.fromkeys(pairs))
            if new!=old:json_write(p,new);current_text[rel]=new;changes+=1
        print(json.dumps(dict(stage='metadata_hashes',iteration=iteration,changed=changes)),flush=True)
        if not changes:break
    else:raise RuntimeError('Portable metadata hash references did not stabilize')
    archive_dir=out/'original_data';archive_dir.mkdir(exist_ok=True);catalog=[]
    for rel,raw in originals.items():
        if not ordered[rel]:continue
        digest=hashlib.sha256(raw).hexdigest();archive=archive_dir/(digest+'.bin')
        if not archive.exists():archive.write_bytes(raw)
        transformed=original_text[rel]
        for before,after in ordered[rel]:transformed=transformed.replace(before,after)
        assert json.loads(transformed)==read(out/rel),('EXACT_OPERATIONAL_REWRITE',rel)
        catalog.append(dict(path=rel,archive=str(archive.relative_to(out)),original_sha256=digest,
            portable_sha256=sha(out/rel),ordered_replacements=ordered[rel]))
    source_archive=out/'original_preferred_envelope.bin';shutil.copy2(source,source_archive)
    # Ensure the principal source is always represented, even when no absolute
    # operational metadata paths happened to require changing it.
    if not any(item['path']==str(source_rel)for item in catalog):
        catalog.append(dict(path=str(source_rel),archive='original_preferred_envelope.bin',
            original_sha256=sha(source),portable_sha256=sha(out/source_rel),ordered_replacements=[]))
    protected_count=0
    for rel in protected:
        for p in(ROOT/rel).rglob('*'):
            if p.is_file()and'__pycache__'not in p.parts:
                assert sha(p)==sha(out/p.relative_to(ROOT));protected_count+=1
    provenance=read(BASE/'PROVENANCE.json')
    provenance.update(kind='PREFERRED_FULLUPPER_ADAPTIVE_PORTABLE_PROVENANCE',
        primitive_packages=sorted(packages),proof_scope='Only the completely accepted chosen-coefficient full-upper/global proof; all later research is excluded as a global premise.',
        copied_code_adapter_scope='Relative path resolution and explicit fresh R3 input selection only; optional replay output files are isolated. Mathematical formulas unchanged.')
    provenance['preferred']=dict(operator_b=str(b),acceptance=str(x.acceptance.resolve().relative_to(ROOT)),
        source=str(source_rel),original_source_sha256=sha(source),portable_source_sha256=sha(out/source_rel),
        original_source_archive='original_preferred_envelope.bin',archived_operational_json=catalog,
        norm_cap=accepted['norm_upper_bound'],GR_lower_bound=accepted['constructive_G_R_lower_bound'],
        extension_filter_used=False,protected_native_package_bytes=True,protected_file_count=protected_count)
    write(out/'PROVENANCE.json',provenance);write(out/'REPLAY_CONFIG.json',config)
    files=sorted(p for p in out.rglob('*')if p.is_file()and not any(k in EXCLUDED for k in p.relative_to(out).parts)and p.name!='MANIFEST.json')
    write(out/'MANIFEST.json',dict(status='FROZEN_NEXT_PREFERRED_GLOBAL_INPUT_MANIFEST',extension_filter_used=False,excluded_mutable_directories=sorted(EXCLUDED),sha256={str(p.relative_to(out)):sha(p)for p in files}))
    # The copied default manifest will be regenerated after caller freezes the
    # README, chosen simple cap, validation metadata and all verifier commands.
    assert all(sha(BASE/rel)==h for rel,h in baseline_manifest['sha256'].items())
    print(json.dumps(dict(status='NEXT_PREFERRED_FROZEN_PAYLOAD_ASSEMBLED',output=str(out),
        operator_b=str(b),primitive_packages=len(packages),protected_byte_identical_files=protected_count,
        operational_json_archive_count=len(catalog),source_original_sha256=sha(source)),indent=2))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--acceptance',type=Path,required=True)
    ap.add_argument('--config',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--resume',action='store_true')
    main(ap.parse_args())
