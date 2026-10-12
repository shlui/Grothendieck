#!/usr/bin/env python3
"""Build the logical verification inventory; reject changed archived bytes."""
from pathlib import Path
import argparse
import hashlib
import json

from package_paths import PackagePaths

root = Path(__file__).resolve().parents[1]
paths = PackagePaths(root)
manifest = json.loads(paths.frozen('new_round_2026-10-02/audit/portable_replay/FINAL_GOAL_INPUT_MANIFEST.json').read_text())

def sha(p):
    h = hashlib.sha256()
    with p.open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()

for name, expected in manifest['sha256'].items():
    assert sha(paths.frozen(name)) == expected, name
archive_path = paths.resolve('inventories/SOURCE_ARCHIVE_FILES.json')
archive = json.loads(archive_path.read_text())
assert archive['all_original_files_byte_identical'] is True
assert archive['file_count'] == 7647 and len(archive['files']) == 7647
assert len({item['path'] for item in archive['files']}) == 7647
for item in archive['files']:
    p = paths.frozen(item['path'])
    assert p.stat().st_size == item['size_bytes'] and sha(p) == item['sha256'], item['path']
names = sorted(name for name in paths.files if name.startswith('verification/'))
answer = dict(status='PACKAGED_VERIFICATION_FILE_INVENTORY',
              original_frozen_manifest_entries=len(manifest['sha256']),
              distinct_frozen_manifest_files=len({paths.frozen(n) for n in manifest['sha256']}),
              source_archive_inventory_sha256=sha(archive_path),
              frozen_files=sum(n.startswith('verification/frozen/') for n in names),
              strict_files=sum(n.startswith('verification/strict/') for n in names),
              total_bytes=sum(paths.resolve(n).stat().st_size for n in names),
              maximum_file_bytes=max(paths.resolve(n).stat().st_size for n in names),
              files_exceeding_github_100_mib=[n for n in names if paths.resolve(n).stat().st_size > 100 * 1024 * 1024],
              sha256={n: sha(paths.resolve(n)) for n in names})
ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument('--output', type=Path, help='Write a new inventory here; otherwise preserve the mapped inventory if already identical.')
args = ap.parse_args()
out = args.output or paths.resolve('inventories/VERIFICATION_FILES.json')
text = json.dumps(answer, indent=2) + '\n'
if args.output is None and out.read_text() != text:
    raise RuntimeError('Computed inventory differs from unchanged archived inventory; use --output to inspect it.')
if args.output is not None:
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(text)
print(json.dumps({k: v for k, v in answer.items() if k != 'sha256'}, indent=2))
