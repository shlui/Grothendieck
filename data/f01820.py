#!/usr/bin/env python3
"""Assemble the accepted original global bound without optimizer dependencies."""
from pathlib import Path
import gzip
import hashlib
import json
import shutil

ROOT=Path(__file__).resolve().parents[3]
OUT=ROOT/'tonight/operator/native_global_b7171'
PREFIX=str(ROOT)+'/'


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def load(path):return json.loads(path.read_text())
def copy_file(rel):
    source=ROOT/rel;target=OUT/rel
    target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(source,target)
def copy_tree(rel):
    source=ROOT/rel
    for path in source.rglob('*'):
        if path.is_file()and '__pycache__'not in path.parts:copy_file(path.relative_to(ROOT))


def strings(value,key=''):
    if isinstance(value,dict):
        for k,v in value.items():
            if k not in {'origin','proposal_origin','optimizer_log','discovery_log'}:
                yield k
                yield from strings(v,k)
    elif isinstance(value,list):
        for v in value:yield from strings(v,key)
    elif isinstance(value,str):yield value


def candidate(value,parent):
    if value.startswith(PREFIX):path=Path(value)
    elif value.startswith(('tonight/','reduced_search/')):path=ROOT/value
    elif value.startswith('../'):path=(parent/value).resolve()
    else:return None
    try:rel=path.relative_to(ROOT)
    except ValueError:return None
    if 'source_package'in rel.parts or 'package'in rel.parts:return None
    if path.is_file()and path.suffix in {'.json','.gz','.py','.md','.tex','.txt'}:return rel
    return None


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    cover=ROOT/'tonight/recovery/profile_cover'
    package_names=set()
    ledgers=['native_union_b7171_target.json','native_upper_b7171_target_union.json',
             'native_lower_base_b7171_target_union.json','native_b7171_regional_transfer_assembly.json',
             'native_R3_b7171_u682_6905_union.json']
    for filename in ledgers:
        d=load(cover/filename)
        for binding in d.get('bindings',[]):
            package_names.add(binding['package'])
            mask=load(ROOT/binding['mask'])
            for alt in mask.get('alternative_scalar_bindings',[]):package_names.add(alt['package'])
        for row in d.get('rows',[]):package_names.add(row['package'])
    for rel in sorted(package_names):copy_tree(rel)
    # Existing standalone R3 prerequisite package supplies all local proof
    # code and four rank/diffuse trees. Copy only its proof payload topology.
    accepted=ROOT/'tonight/operator/projection/replay_accepted_R3'
    for rel in load(accepted/'MANIFEST.json')['sha256']:
        if rel.startswith(('tonight/','reduced_search/'))and(ROOT/rel).is_file():copy_file(rel)
    seeds=[
        'tonight/audit/global_b7171/GLOBAL_B7171_ACCEPTANCE.json',
        'tonight/audit/global_b7171/GLOBAL_B7171_REVIEW.md',
        'tonight/audit/global_b7171/GLOBAL_B7171_THEOREM.tex',
        'tonight/audit/global_b7171/replay_R3_corners.py',
        'tonight/audit/global_b7171/R3_CORNER_REPLAY_PROOF.md',
        'tonight/audit/global_b7171/assemble_independent_acceptance.py',
        'tonight/audit/replay_native_R3_unions.py','tonight/audit/check_b7171_first_transfers.py',
        'tonight/audit/mpfr_interval.py','tonight/audit/replay_R3_physical_q_trees.py',
        'tonight/theory/audit_same_b_native_R3_glue.py',
        'tonight/theory/same_b_native_R3_glue_audit.json',
        'tonight/theory/audit_global_native_R3_envelope.py',
        'tonight/theory/global_b7171_norm_and_sequence.tex',
        'tonight/theory/replay_r3_bernstein_center_cover.py',
        'tonight/theory/replay_r3_local_cover.py',
        'tonight/theory/replay_general_r3_cover.py',
        'tonight/theory/audit_general_r3_parameters_tail.py',
        'tonight/theory/GENERAL_R3_EXTENSION_AUDIT.md',
        'tonight/theory/general_r3_extension.tex',
        'tonight/theory/r3_physical_dual_criterion.tex',
        'tonight/audit/R3_PHYSICAL_DUAL_REVIEW.md',
        'tonight/audit/refit847/replay_complete_refit_cover.py',
        'tonight/audit/refit847/certify_refit_rows_independently.py',
        'tonight/audit/refit847/independent_row_68.json',
        'tonight/audit/refit847/independent_row_96.json',
        'reduced_search/energy/audit_stein_L1_independent.py',
        'reduced_search/energy/audit_stein_L1_independent.json',
        'reduced_search/energy/STEIN_L1_INDEPENDENT_AUDIT.md',
        'reduced_search/construct_reduced_sequence.py',
        'reduced_search/theory/universal_reduced_sequence.md',
        'tonight/paper_universal_sequence.tex',
        'tonight/recovery/profile_cover/global_middle_b7171_native_R3_envelope.json',
        *['tonight/recovery/profile_cover/'+x for x in ledgers],
    ]
    copy_tree('reduced_search/energy/633')
    if(ROOT/'reduced_search/audit/633_TRANSITION_REVIEW.md').is_file():seeds.append('reduced_search/audit/633_TRANSITION_REVIEW.md')
    # Frozen ten-row coefficients/table and previous eight-row binding ledger.
    copy_tree('tonight/operator/cubic_L1_refinement/refit847')
    for filename in ['kappa_shape_caps.json','refit_6_25_7_3_certificate.json','refit_8_85_10_3_certificate.json']:
        seeds.append('tonight/operator/cubic_L1_refinement/'+filename)
    for name in ('0','2','3','4','5','6','7','k8'):
        seeds.append('reduced_search/energy/stein_L1_row_'+name+'.json')
    seeds.append('reduced_search/energy/cubic_L1_certificate.json')
    pending=[Path(rel)for rel in seeds];seen=set()
    while pending:
        rel=pending.pop()
        if rel in seen:continue
        seen.add(rel);source=ROOT/rel
        if not source.is_file():raise FileNotFoundError(rel)
        copy_file(rel)
        if source.suffix=='.json':
            for value in strings(load(source)):
                dependency=candidate(value,source.parent)
                if dependency is not None and dependency not in seen:pending.append(dependency)
    # Explicit code dependencies named by imports rather than JSON records.
    for folder in ['tonight/operator/projection/certificate']:
        for path in(ROOT/folder).glob('*.py'):copy_file(path.relative_to(ROOT))
    for name in ['general_R3_cover.py','general_R3_price.py','price_finer_support_cells.py']:
        copy_file('tonight/operator/projection/lower_amp_certificate/'+name)
    originals={};rewrites=[]
    for target in OUT.rglob('*'):
        if not target.is_file()or '__pycache__'in target.parts:continue
        rel=target.relative_to(OUT);source=ROOT/rel
        if source.is_file():originals[str(rel)]=sha(source)
        if target.suffix in {'.json','.py','.md','.tex','.txt','.log'}:
            content=target.read_text()
            if PREFIX in content:
                transformed=content.replace(PREFIX,'')
                target.write_text(transformed)
                rewrites.append(dict(path=str(rel),original_sha256=originals[str(rel)],
                    portable_sha256=sha(target),deterministic_rewrite='Remove the original repository absolute prefix from metadata only.'))
    # Rebind only copied metadata references after deterministic path rewrites.
    # Mathematical profile/tree/scalar-row hashes remain untouched.
    rebound=[]
    for iteration in range(20):
        changed=False
        mapping={originals[str(path.relative_to(OUT))]:sha(path)
                 for path in OUT.rglob('*.json')if str(path.relative_to(OUT))in originals
                 and originals[str(path.relative_to(OUT))]!=sha(path)}
        for target in OUT.rglob('*.json'):
            rel=str(target.relative_to(OUT))
            if rel not in originals:continue
            content=target.read_text();updated=content
            for old,new in mapping.items():updated=updated.replace(old,new)
            if updated!=content:
                target.write_text(updated);changed=True
                if rel not in rebound:rebound.append(rel)
        if not changed:break
    else:raise RuntimeError('Copied metadata hash references did not stabilize')
    envelope_rel=Path('tonight/recovery/profile_cover/global_middle_b7171_native_R3_envelope.json')
    original=load(ROOT/envelope_rel);portable=load(OUT/envelope_rel)
    assert original['accepted_R3_binding']!=portable['accepted_R3_binding']
    mathematical_original={k:v for k,v in original.items()if k!='accepted_R3_binding'}
    mathematical_portable={k:v for k,v in portable.items()if k!='accepted_R3_binding'}
    assert mathematical_original==mathematical_portable
    for a,b in zip(original['accepted_R3_binding'],portable['accepted_R3_binding']):
        assert a['sha256']==b['sha256']and a['path'].replace(PREFIX,'')==b['path']
    archive=OUT/'original_global_envelope.json.gz'
    with archive.open('wb')as raw:
        with gzip.GzipFile(filename='',mode='wb',fileobj=raw,mtime=0)as stream:stream.write((ROOT/envelope_rel).read_bytes())
    with gzip.open(archive,'rb')as stream:assert hashlib.sha256(stream.read()).hexdigest()==sha(ROOT/envelope_rel)
    provenance=dict(kind='FROZEN_ORIGINAL_GLOBAL_B7171_PORTABLE_PROVENANCE',
        original_norm_cap='70756627/125000000',original_GR_lower_bound='125000000/70756627',
        original_envelope_sha256=sha(ROOT/envelope_rel),portable_envelope_sha256=sha(OUT/envelope_rel),
        archival_copy='original_global_envelope.json.gz',archive_decompressed_sha256=sha(ROOT/envelope_rel),
        mathematical_content_equal=True,only_envelope_field_changed='accepted_R3_binding path metadata',
        path_rewrites=rewrites,hash_rebound_metadata=rebound,original_source_sha256=originals,
        primitive_packages=sorted(package_names),
        proof_scope='Frozen original cap only; no newer source/contact/native refinements are substituted.')
    (OUT/'PROVENANCE.json').write_text(json.dumps(provenance,indent=2)+'\n')
    print(json.dumps({'copied_files':len(originals),'primitive_packages':len(package_names),
                      'metadata_path_rewrites':len(rewrites),'original_envelope_sha256':provenance['original_envelope_sha256']},indent=2))


if __name__=='__main__':main()
