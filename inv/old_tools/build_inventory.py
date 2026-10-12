#!/usr/bin/env python3
"""Build the repository verification inventory; reject changed frozen bytes."""
from pathlib import Path
import hashlib
import json

root = Path(__file__).resolve().parents[1]
frozen = root / 'verification/frozen'
manifest = json.loads((frozen / 'new_round_2026-10-02/audit/portable_replay/FINAL_GOAL_INPUT_MANIFEST.json').read_text())
def sha(p):
    h = hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda: f.read(1024 * 1024), b''):
            h.update(b)
    return h.hexdigest()
for n, expected in manifest['sha256'].items():
    assert sha(frozen / n) == expected, n
archive_inventory = json.loads((root / 'inventories/SOURCE_ARCHIVE_FILES.json').read_text())
assert archive_inventory['all_original_files_byte_identical'] is True
assert archive_inventory['file_count'] == 7647
for item in archive_inventory['files']:
    p = frozen / item['path']
    assert p.stat().st_size == item['size_bytes'] and sha(p) == item['sha256'], item['path']
files = sorted(p for p in (root / 'verification').rglob('*') if p.is_file() and '__pycache__' not in p.parts)
answer = dict(status='PACKAGED_VERIFICATION_FILE_INVENTORY',
              original_frozen_manifest_entries=len(manifest['sha256']),
              distinct_frozen_manifest_files=len({(frozen / n).resolve() for n in manifest['sha256']}),
              source_archive_inventory_sha256=sha(root / 'inventories/SOURCE_ARCHIVE_FILES.json'),
              frozen_files=sum(p.is_relative_to(frozen) for p in files),
              strict_files=sum(p.is_relative_to(root / 'verification/strict') for p in files),
              total_bytes=sum(p.stat().st_size for p in files),
              maximum_file_bytes=max(p.stat().st_size for p in files),
              files_exceeding_github_100_mib=[str(p.relative_to(root)) for p in files if p.stat().st_size > 100 * 1024 * 1024],
              sha256={str(p.relative_to(root)): sha(p) for p in files})
out = root / 'inventories/VERIFICATION_FILES.json'
out.parent.mkdir(exist_ok=True)
out.write_text(json.dumps(answer, indent=2) + '\n')
print(json.dumps({k: v for k, v in answer.items() if k != 'sha256'}, indent=2))
