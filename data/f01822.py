#!/usr/bin/env python3
"""Extend a copy of the original bundle with only accepted strengthened gates."""
from pathlib import Path
import gzip
import hashlib
import importlib.util
import json
import shutil

ROOT=Path(__file__).resolve().parents[3]
OLD=ROOT/'tonight/operator/native_global_b7171'
OUT=ROOT/'tonight/operator/native_global_b7171_strengthened'
spec=importlib.util.spec_from_file_location('portable_builder_definitions',Path(__file__).with_name('build_portable_bundle.py'))
B=importlib.util.module_from_spec(spec);spec.loader.exec_module(B);B.OUT=OUT


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def read(path):return json.loads(path.read_text())


def main():
    original_manifest=read(OLD/'MANIFEST.json')
    assert all(sha(OLD/rel)==h for rel,h in original_manifest['sha256'].items())
    shutil.copytree(OLD,OUT,dirs_exist_ok=True,ignore=shutil.ignore_patterns('replay_outputs','__pycache__'))
    shutil.copy2(OUT/'verify_bundle.py',OUT/'verify_original_components.py')
    original_provenance=read(OLD/'PROVENANCE.json');originals=dict(original_provenance['original_source_sha256'])
    B.copy_tree('tonight/native_frames/critical0750_0760')
    seeds=[
        'tonight/audit/global_b7171/STRENGTHENED_GLOBAL_B7171_ACCEPTANCE.json',
        'tonight/audit/global_b7171/assemble_strengthened_acceptance.py',
        'tonight/audit/global_b7171/replay_critical_native_endpoints.py',
        'tonight/audit/global_b7171/replay_ninebank_R3_corners.py',
        'tonight/audit/global_b7171/targeted_32_outer_vertices.json',
        'tonight/audit/global_b7171/replay_targeted_outer_boxes.py',
        'tonight/audit/global_b7171/check_physical_inverse_bounds.py',
        'tonight/audit/global_b7171/check_selected_contact_credits.py',
        'tonight/audit/verify_cubic_head_contact.py',
        'tonight/audit/verify_linear_contact_exclusion.py',
        'tonight/moments/cubic_head_contact/certify_direct_contact_majorant.py',
        'tonight/audit/replay_R3_generous_q_trees.py',
        'tonight/audit/targeted_lowr_R3_recheck/ACCEPTANCE.json',
        'tonight/theory/audit_clipped_contact_envelope.py',
        'tonight/theory/clipped_contact_envelope_audit.json',
        'tonight/theory/global_b7171_strengthened_certificate.tex',
        'tonight/theory/replay_targeted_highk.py',
        'tonight/theory/targeted_highk_whole_q_audit.json',
        'tonight/theory/cubic_head_contact_fragment.tex',
        'tonight/theory/linear_contact_source_fragment.tex',
        'tonight/recovery/profile_cover/global_middle_b7171_critical_ninebank_contact_envelope.json',
    ]
    seen=set();pending=list(map(Path,seeds))
    while pending:
        rel=pending.pop()
        if rel in seen:continue
        seen.add(rel);source=ROOT/rel
        if not source.is_file():raise FileNotFoundError(rel)
        B.copy_file(rel);originals[str(rel)]=sha(source)
        if source.suffix=='.json':
            for value in B.strings(read(source)):
                dependency=B.candidate(value,source.parent)
                if dependency is not None and dependency not in seen:pending.append(dependency)
    prefix=str(ROOT)+'/'
    rewrites=[]
    for path in OUT.rglob('*'):
        if not path.is_file()or path.suffix not in {'.json','.py','.md','.tex','.txt','.log'}:continue
        content=path.read_text()
        if prefix in content:
            path.write_text(content.replace(prefix,''));rewrites.append(str(path.relative_to(OUT)))
    # Map relative ledger references to the bundle root in copied audit code.
    own=OUT/'tonight/audit/global_b7171/replay_ninebank_R3_corners.py'
    # This source already handles both relative and absolute paths.
    assert "if not path.is_absolute():path=ROOT/path"in own.read_text()
    copied_R3=OUT/'tonight/audit/global_b7171/replay_R3_corners.py'
    content=copied_R3.read_text()
    if "if not path.is_absolute():path=ROOT/path"not in content:
        content=content.replace("path=Path(item['path']);assert sha(path)==item['sha256'];checked_ledger.append(file_record(path))",
            "path=Path(item['path'])\n        if not path.is_absolute():path=ROOT/path\n        assert sha(path)==item['sha256'];checked_ledger.append(file_record(path))")
        copied_R3.write_text(content)
    diffuse=OUT/'tonight/audit/diffuse_b7171_R3/independent_source_pricing.py';content=diffuse.read_text()
    if "if not path.is_absolute():path=ROOT/path"not in content:
        content=content.replace('path=Path(filename);assert sha(path)==digest',
            'path=Path(filename)\n        if not path.is_absolute():path=ROOT/path\n        assert sha(path)==digest')
        diffuse.write_text(content)
    rebound=[]
    for iteration in range(20):
        mapping={h:sha(OUT/rel)for rel,h in originals.items()if(OUT/rel).is_file()and(OUT/rel).suffix=='.json'and sha(OUT/rel)!=h}
        changed=False
        for path in OUT.rglob('*.json'):
            if path.name in {'PROVENANCE.json','MANIFEST.json','VALIDATION.json'}:continue
            content=path.read_text();updated=content
            for before,after in mapping.items():updated=updated.replace(before,after)
            if content!=updated:
                path.write_text(updated);changed=True
                rel=str(path.relative_to(OUT))
                if rel not in rebound:rebound.append(rel)
        if not changed:break
    else:raise RuntimeError('Strengthened metadata references failed to stabilize')
    rel=Path('tonight/recovery/profile_cover/global_middle_b7171_critical_ninebank_contact_envelope.json')
    original=read(ROOT/rel);portable=read(OUT/rel)
    assert {k:v for k,v in original.items()if k!='accepted_R3_binding'}=={
           k:v for k,v in portable.items()if k!='accepted_R3_binding'}
    assert len(original['accepted_R3_binding'])==len(portable['accepted_R3_binding'])
    for a,b in zip(original['accepted_R3_binding'],portable['accepted_R3_binding']):
        assert a['sha256']==b['sha256']and a['path'].replace(prefix,'')==b['path']
    archive=OUT/'original_strengthened_envelope.json.gz'
    with archive.open('wb')as raw:
        with gzip.GzipFile(filename='',mode='wb',fileobj=raw,mtime=0)as stream:stream.write((ROOT/rel).read_bytes())
    provenance={**original_provenance,'kind':'FROZEN_STRENGTHENED_GLOBAL_B7171_PORTABLE_PROVENANCE',
        'original_source_sha256':originals,
        'primitive_packages':[*original_provenance['primitive_packages'],'tonight/native_frames/critical0750_0760'],
        'strengthened':dict(original_envelope_sha256=sha(ROOT/rel),portable_envelope_sha256=sha(OUT/rel),
            archive='original_strengthened_envelope.json.gz',math_content_equal=True,
            path_rewrites=rewrites,hash_rebound_metadata=rebound,norm_cap='56513883/100000000',GR_lower_bound='100000000/56513883'),
        'proof_scope':'Accepted strengthened nine-profile/contact cap only. Newer joint-contact and middle-profile upgrades are excluded.'}
    (OUT/'PROVENANCE.json').write_text(json.dumps(provenance,indent=2)+'\n')
    # Install strengthened entrypoint and its standalone instructions after the
    # original payload copy, so a rebuild cannot silently restore the old wrapper.
    shutil.copy2(Path(__file__).with_name('strengthened_verify_bundle.py'),OUT/'verify_bundle.py')
    shutil.copy2(Path(__file__).with_name('strengthened_bundle_README.md'),OUT/'README.md')
    assert all(sha(OLD/p)==h for p,h in original_manifest['sha256'].items())
    print(json.dumps({'status':'ASSEMBLED_STRENGTHENED_PORTABLE_PAYLOAD','source_envelope_sha256':sha(ROOT/rel),
                      'additional_source_files':len(seen),'primitive_packages':13,'original_bundle_unchanged':True},indent=2))


if __name__=='__main__':main()
