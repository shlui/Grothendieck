#!/usr/bin/env python3
"""Assemble one complete preferred proof tree without extension exclusions."""
from pathlib import Path
import argparse
import gzip
import hashlib
import json
import shutil

ROOT=Path(__file__).resolve().parents[3]
BASE=ROOT/'tonight/operator/native_global_b7171_strengthened'
OUT=ROOT/'tonight/operator/preferred_global'
PREFIX=str(ROOT)+'/'


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def read(path):
    with(gzip.open(path,'rt')if path.suffix=='.gz'else path.open())as stream:return json.load(stream)
def strings(value):
    if isinstance(value,dict):
        for k,v in value.items():
            if k not in{'origin','proposal_origin','optimizer_log','discovery_log'}:
                yield k;yield from strings(v)
    elif isinstance(value,list):
        for child in value:yield from strings(child)
    elif isinstance(value,str):yield value
def candidate(value,parent):
    if '\n'in value or len(value)>2048:return None
    # Relative premise names contain an extension. Avoid treating every
    # large rational/decimal number in a proof tree as a filesystem path.
    if not value.startswith((PREFIX,'tonight/','reduced_search/','../'))and'.'not in value:return None
    if value and all(c in'0123456789/.-+eE'for c in value):return None
    if value.startswith(PREFIX):path=Path(value)
    elif value.startswith(('tonight/','reduced_search/')):path=ROOT/value
    else:path=(parent/value).resolve()
    try:rel=path.relative_to(ROOT)
    except ValueError:return None
    try:
        if path.is_file()and '__pycache__'not in rel.parts:return rel
    except OSError:return None
    return None
def copy_file(rel):
    target=OUT/rel;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(ROOT/rel,target)
def copy_tree(rel):
    for path in(ROOT/rel).rglob('*'):
        if path.is_file()and '__pycache__'not in path.parts:copy_file(path.relative_to(ROOT))


def main(acceptance_path,assembly_path):
    accepted=read(acceptance_path);assert accepted['status'].startswith('PASS_INDEPENDENT_COMPLETE_')
    base_manifest=read(BASE/'MANIFEST.json');assert all(sha(BASE/p)==h for p,h in base_manifest['sha256'].items())
    shutil.copytree(BASE,OUT,dirs_exist_ok=True,ignore=shutil.ignore_patterns('replay_outputs','__pycache__'))
    shutil.copy2(OUT/'verify_bundle.py',OUT/'verify_strengthened_components.py')
    originals=dict(read(BASE/'PROVENANCE.json')['original_source_sha256'])
    packages=set(read(BASE/'PROVENANCE.json')['primitive_packages'])
    source_rel=Path(accepted['envelope']['path']);original_source=ROOT/source_rel
    assert sha(original_source)==accepted['envelope']['sha256']
    for item in read(original_source)['native_endpoint_inputs']:packages.add(item['package'])
    # Every byte of each prerequisite package is included. This deliberately
    # includes bound.log, C sources, headers and any other hash-bound extension.
    for rel in sorted(packages):copy_tree(rel)
    copy_tree('reduced_search/energy/633')
    r3=ROOT/'tonight/operator/projection/replay_accepted_R3_with_targeted'
    manifest=read(r3/'MANIFEST.json')
    for rel,h in manifest['sha256'].items():
        assert sha(r3/rel)==h
        if rel.startswith(('tonight/','reduced_search/')):
            if(ROOT/rel).is_file():copy_file(rel)
            else:
                target=OUT/rel;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(r3/rel,target)
    shutil.copy2(r3/'verify_all.py',OUT/'verify_R3_prerequisites.py')
    shutil.copy2(r3/'verify_targeted.py',OUT/'verify_targeted_prerequisites.py')
    portable_targeted=ROOT/'tonight/audit/global_b7171/preferred_targeted_verify_template.py'
    if portable_targeted.is_file():shutil.copy2(portable_targeted,OUT/'verify_targeted_prerequisites.py')
    seeds=[acceptance_path.relative_to(ROOT),assembly_path.relative_to(ROOT),source_rel,
        Path('tonight/audit/global_b7171/replay_generic_R3_cache.py'),
        Path('tonight/audit/global_b7171/replay_middlebank_R3_corners.py'),
        Path('tonight/audit/global_b7171/replay_robust_outer_boxes.py'),
        Path('tonight/audit/global_b7171/replay_joint_contact_full.py'),
        Path('tonight/audit/global_b7171/replay_extended_critical_native_endpoints.py'),
        Path('tonight/audit/global_b7171/replay_final46_native_endpoints.py'),
        Path('tonight/audit/global_b7171/check_refined_native_transport.py'),
        Path('tonight/audit/global_b7171/transport_complete_outer_branches.py'),
        Path('tonight/recovery/profile_cover/kappa_refinement_2_4/complete.json'),
        Path('tonight/recovery/profile_cover/global_middle_b7171_extendedcritical_ninebank_planes_envelope.json'),
        Path('tonight/audit/global_b7171/check_polygon_plane_prices.py'),
        Path('tonight/audit/global_b7171/check_polygon_contact_credits.py'),
        Path('tonight/theory/audit_polygon_envelope.py'),
        Path('tonight/theory/audit_preferred_polygon_envelope.py'),
        Path('tonight/theory/source_polygon_prices.py'),
        Path('tonight/theory/final_preferred_polygon_certificate.tex'),
        Path('tonight/recovery/profile_cover/polygon_prices.py'),
        Path('tonight/audit/verify_linear_contact_bank.py')]
    seen=set();pending=seeds[:]
    while pending:
        rel=pending.pop()
        if rel in seen:continue
        seen.add(rel);source=ROOT/rel
        if not source.is_file():raise FileNotFoundError(rel)
        copy_file(rel)
        if source.suffix=='.json'or(source.suffix=='.gz'and source.name.endswith('.json.gz')):
            for value in strings(read(source)):
                dependency=candidate(value,source.parent)
                if dependency is not None and dependency not in seen:pending.append(dependency)
    # Preserve input bytes for every native package and every scalar/local
    # primitive. Only external operational metadata is rewritten below.
    protected=[Path(p)for p in packages]+[Path('reduced_search/energy/633')]
    for path in OUT.rglob('*'):
        if not path.is_file()or'__pycache__'in path.parts:continue
        rel=path.relative_to(OUT);source=ROOT/rel
        if source.is_file():originals[str(rel)]=sha(source)
        if any(rel.is_relative_to(p)for p in protected):continue
        if path.suffix in{'.json','.py','.md','.tex','.txt'}:
            content=path.read_text()
            if PREFIX in content:path.write_text(content.replace(PREFIX,''))
    source_target=OUT/source_rel
    corner=OUT/'tonight/audit/global_b7171/replay_R3_corners.py';content=corner.read_text()
    if 'if not path.is_absolute():path=ROOT/path'not in content:
        content=content.replace("path=Path(item['path']);assert sha(path)==item['sha256'];checked_ledger.append(file_record(path))",
            "path=Path(item['path'])\n        if not path.is_absolute():path=ROOT/path\n        assert sha(path)==item['sha256'];checked_ledger.append(file_record(path))")
        corner.write_text(content)
    diffuse=OUT/'tonight/audit/diffuse_b7171_R3/independent_source_pricing.py';content=diffuse.read_text()
    if 'if not path.is_absolute():path=ROOT/path'not in content:
        content=content.replace('path=Path(filename);assert sha(path)==digest',
            'path=Path(filename)\n        if not path.is_absolute():path=ROOT/path\n        assert sha(path)==digest')
        diffuse.write_text(content)
    if source_target.suffix=='.gz':
        operational=read(source_target);content=json.dumps(operational,indent=2).replace(PREFIX,'')+'\n'
        with source_target.open('wb')as raw:
            with gzip.GzipFile(filename='',mode='wb',fileobj=raw,mtime=0)as stream:stream.write(content.encode())
    # Hash references in copied metadata track only deterministic metadata
    # rewrites; raw scalar/local/native mathematical inputs retain their SHA.
    rebound=[]
    for _ in range(24):
        mapping={old:sha(OUT/rel)for rel,old in originals.items()if(OUT/rel).is_file()
            and((OUT/rel).suffix=='.json'or rel==str(source_rel))and sha(OUT/rel)!=old}
        changed=False
        for path in OUT.rglob('*.json'):
            rel=path.relative_to(OUT)
            if path.name in{'PROVENANCE.json','MANIFEST.json','VALIDATION.json'}or any(rel.is_relative_to(p)for p in protected):continue
            content=path.read_text();new=content
            for before,after in mapping.items():new=new.replace(before,after)
            if new!=content:path.write_text(new);changed=True;rebound.append(str(rel))
        if source_target.suffix=='.gz':
            content=json.dumps(read(source_target),indent=2)+'\n';new=content
            for before,after in mapping.items():new=new.replace(before,after)
            if new!=content:
                with source_target.open('wb')as raw:
                    with gzip.GzipFile(filename='',mode='wb',fileobj=raw,mtime=0)as stream:stream.write(new.encode())
                changed=True;rebound.append(str(source_rel))
        if not changed:break
    else:raise RuntimeError('Preferred metadata references did not stabilize')
    # The exact algorithm and all cell numbers must survive serialization.
    old=read(original_source);new=read(source_target)
    for field in('cells','physically_empty_cells','cached_corner_constants','profile_coefficients',
                  'maximum_cell_cap','strict_rational_cap','source_polygon_planes','physical_u_bounds'):
        assert old[field]==new[field],('MATHEMATICAL_SOURCE_CHANGED',field)
    archive=OUT/'original_preferred_envelope.bin';shutil.copy2(original_source,archive)
    provenance=read(BASE/'PROVENANCE.json')
    provenance.update(kind='PREFERRED_COMPLETE_GLOBAL_B7171_PORTABLE_PROVENANCE',
        original_source_sha256=originals,primitive_packages=sorted(packages),
        preferred=dict(acceptance=str(acceptance_path.relative_to(ROOT)),assembly=str(assembly_path.relative_to(ROOT)),
            source=str(source_rel),original_source_sha256=sha(original_source),portable_source_sha256=sha(source_target),
            original_source_archive='original_preferred_envelope.bin',mathematical_fields_exactly_equal=True,
            metadata_hash_rebound=sorted(set(rebound)),norm_cap=accepted['norm_upper_bound'],
            GR_lower_bound=accepted['constructive_G_R_lower_bound'],
            extension_filter_used=False,protected_native_package_bytes=True))
    (OUT/'PROVENANCE.json').write_text(json.dumps(provenance,indent=2)+'\n')
    assert all(sha(BASE/p)==h for p,h in base_manifest['sha256'].items())
    print(json.dumps(dict(status='PREFERRED_PAYLOAD_ASSEMBLED',files=len([p for p in OUT.rglob('*')if p.is_file()]),
        primitive_packages=len(packages),source=str(source_rel),immutable_strengthened_unchanged=True),indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--acceptance',type=Path,default=ROOT/'tonight/audit/global_b7171/POLYGON_GLOBAL_B7171_ACCEPTANCE.json')
    parser.add_argument('--assembly',type=Path,default=ROOT/'tonight/audit/global_b7171/assemble_polygon_acceptance.py')
    parser.add_argument('--output',type=Path,default=OUT);args=parser.parse_args();OUT=args.output.resolve()
    main(args.acceptance.resolve(),args.assembly.resolve())
