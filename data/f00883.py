#!/usr/bin/env python3
"""Verify the frozen proof whitelist and complete recovered historical sources.

This creates archive metadata only. It never changes a proof input and never
collects unrecorded runtime files, environments, or duplicate baseline copies.
"""
from pathlib import Path
import argparse
import hashlib
import json

ROOT = Path(__file__).resolve().parents[3]
HISTORY = ROOT / 'verified_2026-10-01'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


parser = argparse.ArgumentParser()
parser.add_argument('--manifest', type=Path, required=True)
parser.add_argument('--output', type=Path, required=True)
parser.add_argument('--extra', type=Path, action='append', default=[])
args = parser.parse_args()
manifest_path = args.manifest.resolve()
manifest = json.loads(manifest_path.read_text())
assert manifest['status'] == 'FROZEN_COMPLETE_NEW_ROUND_INPUT_WHITELIST'
assert manifest['extension_filter_used'] is False
records = {}


def add(path, expected=None):
    path = path.resolve()
    assert path.is_relative_to(ROOT) and path.is_file(), path
    relative = str(path.relative_to(ROOT))
    assert not relative.startswith(('new_round_2026-10-02/baseline_replay/',
                                    'new_round_2026-10-02/venv/')), relative
    digest = sha(path)
    assert expected is None or digest == expected, ('HASH_MISMATCH', relative)
    assert relative not in records or records[relative] == digest
    records[relative] = digest


for relative, digest in manifest['sha256'].items():
    add(ROOT / relative, digest)
add(manifest_path)
historical_map_path = HISTORY / 'SHA256SUMS.json'
historical_map = json.loads(historical_map_path.read_text())
for relative, digest in historical_map.items():
    add(HISTORY / relative, digest)
add(historical_map_path)
add(Path(__file__))
for extra in args.extra:
    add(extra)

output = args.output.resolve()
assert output.is_relative_to(ROOT) and not output.exists(), output
output.parent.mkdir(parents=True, exist_ok=True)
result = {
    'status': 'PASS_VERIFIED_FROZEN_PROOF_AND_COMPLETE_HISTORICAL_SOURCE_ALLOWLIST',
    'proof_input_manifest': {'path': str(manifest_path.relative_to(ROOT)),
                             'sha256': sha(manifest_path)},
    'historical_source_map': {'path': str(historical_map_path.relative_to(ROOT)),
                              'sha256': sha(historical_map_path)},
    'historical_recorded_members': len(historical_map),
    'historical_checksum_file_included_separately': True,
    'file_count': len(records),
    'total_bytes': sum((ROOT / relative).stat().st_size for relative in records),
    'extension_filter_used': False,
    'certificate_bytes_rewritten': False,
    'entrypoint_sha256': sha(Path(__file__)),
    'sha256': dict(sorted(records.items())),
    'archive_metadata_instruction': 'Include this allowlist JSON itself separately; '
        'its self hash cannot be a member of its own digest map. Archive every '
        'listed member under its exact relative path, verify hashes and ZIP CRCs '
        'after extraction, and preserve accepted certificate bytes.',
}
output.write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps({'path': str(output.relative_to(ROOT)), 'sha256': sha(output),
                  'files': len(records), 'bytes': result['total_bytes']}, indent=2))
