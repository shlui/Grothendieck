#!/usr/bin/env python3
"""Freeze, verify and atomically promote the accepted b=.713 proof package."""
from pathlib import Path
import gzip
import hashlib
import json
import os
import zipfile

ROOT=Path(__file__).resolve().parents[3]
STAGE=ROOT/'tonight/operator/preferred_global_staging'
CANONICAL=ROOT/'tonight/operator/preferred_global'
ZIP=ROOT/'tonight/operator/preferred_global.zip'
PREFIX=str(ROOT)+'/'
EXCLUDED={'replay_outputs','replay_outputs_targeted','__pycache__'}

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):
    with(gzip.open(p,'rt')if p.suffix=='.gz'else p.open())as stream:return json.load(stream)
def write(p,value):p.write_text(json.dumps(value,indent=2)+'\n')
def records(v):
    if isinstance(v,dict):
        if isinstance(v.get('path'),str)and isinstance(v.get('sha256'),str):yield v
        for item in v.values():yield from records(item)
    elif isinstance(v,list):
        for item in v:yield from records(item)
def leaf_differences(a,b,path=''):
    assert type(a)==type(b),(path,'type')
    if isinstance(a,dict):
        assert a.keys()==b.keys(),(path,'keys')
        for k in a:yield from leaf_differences(a[k],b[k],path+'/'+k)
    elif isinstance(a,list):
        assert len(a)==len(b),(path,'length')
        for i,(x,y)in enumerate(zip(a,b)):yield from leaf_differences(x,y,path+'/'+str(i))
    elif a!=b:yield(path,a,b)

def main():
    result=read(STAGE/'replay_outputs/result.json')
    assert result['status']=='PASS_PORTABLE_COMPLETE_PREFERRED_GLOBAL_B7171'
    assert result['norm_upper_bound']=='564479481/1000000000'
    assert result['constructive_G_R_lower_bound']=='1000000000/564479481'
    assert result['middle_price_cells']==185761
    assert result['proved_physically_empty_middle_cells']==1574
    assert result['independent_R3_corner_prices']==2520
    immutable=[]
    for rel in('tonight/operator/native_global_b7171','tonight/operator/native_global_b7171_strengthened'):
        p=ROOT/rel;m=read(p/'MANIFEST.json')
        assert all(sha(p/name)==h for name,h in m['sha256'].items())
        immutable.append(dict(path=rel,manifest_sha256=sha(p/'MANIFEST.json'),unchanged=True))
    provenance=read(STAGE/'PROVENANCE.json');preferred=provenance['preferred']
    original=gzip.decompress((STAGE/preferred['original_source_archive']).read_bytes())
    archive=read(STAGE/preferred['source'])
    differences=list(leaf_differences(json.loads(original),archive))
    assert len(differences)==15
    for path,before,after in differences:
        assert path.startswith('/accepted_R3_binding/')and path.endswith('/path')
        assert before.startswith(PREFIX)and after==before[len(PREFIX):]
    protected=0
    for rel in provenance['primitive_packages']+['reduced_search/energy/633']:
        for p in(ROOT/rel).rglob('*'):
            if p.is_file()and'__pycache__'not in p.parts:
                assert sha(p)==sha(STAGE/p.relative_to(ROOT)),p
                protected+=1
    cache=read(STAGE/'tonight/audit/global_b7171/conditional_b713_fine46_R3_corner_replay.json')
    profiles={}
    for rec in records(cache['complete_whole_q_bindings']):
        rel=rec['path']
        if 'profile'in Path(rel).name:
            assert sha(STAGE/rel)==rec['sha256']==sha(ROOT/rel)
            profiles[rel]=rec['sha256']
    assert len(profiles)==14,len(profiles)
    provenance['proof_scope']='Accepted complete preferred b=713/1000 global bound with frozen fourteen-profile, forty-five-band inputs. The earlier original and strengthened proofs are included as prerequisites. Later fifteen-to-nineteen-profile research is excluded as global premises.'
    provenance['copied_code_path_adapters']=sorted(set(provenance.get('copied_code_path_adapters',[])+[
        'tonight/audit/global_b7171/check_polygon_contact_credits.py',
        'verify_targeted_prerequisites.py']))
    provenance['copied_code_adapter_scope']='Relative path resolution and selection of the explicit fresh complete R3 report. Optional prerequisite price checkers run in temporary copies so immutable accepted reports are not overwritten. Mathematical price and verification formulas are unchanged.'
    preferred.update(full_structured_source_equality_except_exact_prefix_rewrite=True,
        changed_source_leaf_paths=[x[0]for x in differences],
        deterministic_prefix_rewrite=dict(before=PREFIX,after='',scope='accepted_R3_binding/*/path only'),
        byte_identical_profiles=profiles,protected_file_count=protected)
    write(STAGE/'PROVENANCE.json',provenance)
    validation=dict(status='PASS_PREFERRED_GLOBAL_B713_FINAL46_FOURTEENBANK_COMPLETE_DEFAULT',
        default_result=result,combined_full_replayed_in_this_packaging_run=False,
        full_prerequisite_evidence='Previously fresh accepted complete scalar/local/R3 leaf proof reports are included with byte-bound original input profiles and trees. The optional full entrypoint and complete prerequisites are present.',
        post_default_changes='Provenance and validation metadata were appended. The optional targeted prerequisite wrapper now isolates historical default output files in a temporary copy. No default mathematical data or verification formulas changed.',
        preferred_default_command='python preferred_global/verify_bundle.py',
        complete_optional_command='python preferred_global/verify_bundle.py --full --workers 4',
        original_and_strengthened_packages=immutable,protected_byte_identical_files=protected,
        byte_identical_R3_profiles=len(profiles),full_original_source_only_prefix_rewrite=True)
    write(STAGE/'VALIDATION.json',validation)
    files=sorted(p for p in STAGE.rglob('*')if p.is_file()and not any(x in EXCLUDED for x in p.relative_to(STAGE).parts)and p.name!='MANIFEST.json')
    manifest=dict(status='FROZEN_FINAL_PREFERRED_GLOBAL_IMMUTABLE_INPUT_MANIFEST',
        extension_filter_used=False,excluded_mutable_directories=sorted(EXCLUDED),
        sha256={str(p.relative_to(STAGE)):sha(p)for p in files})
    write(STAGE/'MANIFEST.json',manifest)
    assert all(sha(STAGE/name)==h for name,h in manifest['sha256'].items())
    temporary_zip=ZIP.with_suffix('.zip.next')
    members=[STAGE/'MANIFEST.json',*files]
    with zipfile.ZipFile(temporary_zip,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=9)as z:
        for p in sorted(members):
            info=zipfile.ZipInfo('preferred_global/'+str(p.relative_to(STAGE)),date_time=(1980,1,1,0,0,0))
            info.compress_type=zipfile.ZIP_DEFLATED;info.external_attr=0o100644<<16
            z.writestr(info,p.read_bytes(),compresslevel=9)
    with zipfile.ZipFile(temporary_zip)as z:
        assert z.testzip()is None
        assert set(z.namelist())=={'preferred_global/'+str(p.relative_to(STAGE))for p in members}
        for rel,h in manifest['sha256'].items():assert hashlib.sha256(z.read('preferred_global/'+rel)).hexdigest()==h
        assert hashlib.sha256(z.read('preferred_global/MANIFEST.json')).hexdigest()==sha(STAGE/'MANIFEST.json')
    backup=ROOT/'tonight/operator/preferred_global_coarse11_tree'
    assert not backup.exists()
    if CANONICAL.exists():CANONICAL.rename(backup)
    STAGE.rename(CANONICAL);os.replace(temporary_zip,ZIP)
    summary=dict(status='PASS_ATOMIC_FINAL_PREFERRED_GLOBAL_PACKAGE',path=str(CANONICAL),zip=str(ZIP),
        zip_sha256=sha(ZIP),zip_bytes=ZIP.stat().st_size,immutable_manifest_files=len(files),
        manifest_sha256=sha(CANONICAL/'MANIFEST.json'),default_seconds=result['seconds'],
        norm_cap=result['norm_upper_bound'],constructive_G_R_lower_bound=result['constructive_G_R_lower_bound'],
        original_and_strengthened_unchanged=True,coarse_fallback_retained=True,
        byte_identical_R3_profiles=len(profiles),protected_byte_identical_files=protected,
        full_original_source_only_prefix_rewrite=True,crc_and_all_zip_member_hashes_pass=True)
    write(ROOT/'tonight/audit/global_b7171/PREFERRED_PACKAGE_CHECKPOINT.json',summary)
    print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
