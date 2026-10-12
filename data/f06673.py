#!/usr/bin/env python3
"""Create a reproducible research archive from an explicit file whitelist."""
from pathlib import Path
import hashlib
import json
import os
import shutil
import zipfile

ROOT = Path(__file__).resolve().parent.parent
STAGE = ROOT / 'tonight/source_package'
OUT = ROOT / 'reduced_search/deliverables/Grothendieck_reduced_sequence_source_2026-10-01.zip'


def collect():
    files = {}
    old = ROOT / 'reduced_search/package'
    for src in old.rglob('*'):
        if src.is_file() and '__pycache__' not in src.parts and src.name != 'SHA256SUMS.json':
            rel = src.relative_to(old)
            if str(rel).startswith('tonight/recovery/profile_cover/continuous_'):
                continue
            live = ROOT / rel
            files[rel] = live if live.is_file() else src
    trees = ['reduced_search/operator/cubic7171', 'tonight/regional_certificates', 'tonight/moments',
             'tonight/native_frames', 'tonight/operator/preferred_global',
             'tonight/audit', 'tonight/theory', 'tonight/operator/projection',
             'tonight/operator/odd_tail_screen',
             'tonight/operator/cubic_L1_refinement', 'tonight/operator/adversarial3d']
    allowed = {'.py', '.json', '.md', '.tex', '.gz', '.npz', '.txt'}
    for tree in trees:
        for src in (ROOT / tree).rglob('*'):
            if src.is_file() and src.suffix in allowed and '__pycache__' not in src.parts:
                files[src.relative_to(ROOT)] = src
    # Frozen standalone proof manifests take precedence over extension filters.
    for bundle in ['tonight/operator/preferred_global']:
        base = ROOT / bundle
        binding = json.loads((base / 'MANIFEST.json').read_text())['sha256']
        for name, digest in binding.items():
            src = base / name
            assert hashlib.sha256(src.read_bytes()).hexdigest() == digest
            files[src.relative_to(ROOT)] = src
    for name in ['operator/literature_parallel_profiles.md', 'pdf_qa/final/BUILD_CHECK.json', 'main_preferred_constants.json', 'package_source.py', 'build_and_check_paper.py', 'adversarial_bivariate_signs.py', 'adversarial_bivariate_signs.json',
                 'adversarial_asymmetric_bivariate.py', 'adversarial_asymmetric_bivariate.json', 'optimize_cubic_witness_coefficients.py',
                 'optimize_cubic_witness_coefficients.json', 'certify_selected_witness_coefficients.py',
                 'selected_witness_coefficient_certificates.json']:
        src = ROOT / 'tonight' / name
        if src.is_file(): files[src.relative_to(ROOT)] = src
    # Published checks use the complete frozen primitive packages.  Their
    # source_profile inputs retain the accepted proposal data; omit large
    # duplicate discovery exports from the top-level research archive.
    return files


def main():
    if STAGE.exists(): shutil.rmtree(STAGE)
    STAGE.mkdir(parents=True)
    files = collect()
    for rel, src in files.items():
        dst = STAGE / rel
        dst.parent.mkdir(parents=True, exist_ok=True)
        # One read binds each member to a single byte snapshot.
        dst.write_bytes(src.read_bytes())
    readme = ROOT / 'tonight/SOURCE_README.md'
    (STAGE / 'README.md').write_bytes(readme.read_bytes())
    manifest = {str(src.relative_to(STAGE)): hashlib.sha256(src.read_bytes()).hexdigest()
                for src in sorted(STAGE.rglob('*')) if src.is_file() and src != STAGE / 'SHA256SUMS.json'}
    (STAGE / 'SHA256SUMS.json').write_text(json.dumps(manifest, indent=2) + '\n')
    temporary = OUT.with_suffix('.zip.pending')
    with zipfile.ZipFile(temporary, 'w', compression=zipfile.ZIP_DEFLATED, compresslevel=8) as archive:
        for src in sorted(STAGE.rglob('*')):
            if src.is_file(): archive.write(src, src.relative_to(STAGE))
    with zipfile.ZipFile(temporary) as archive:
        assert archive.testzip() is None
        for name, digest in manifest.items():
            assert hashlib.sha256(archive.read(name)).hexdigest() == digest
    with temporary.open('rb') as stream:
        os.fsync(stream.fileno())
    os.replace(temporary, OUT)
    print(json.dumps({'status': 'PASS_ARCHIVE_AND_ALL_MEMBER_HASHES',
                      'members': len(manifest) + 1, 'bytes': OUT.stat().st_size,
                      'output': str(OUT)}, indent=2))


if __name__ == '__main__': main()
