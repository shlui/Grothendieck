#!/usr/bin/env python3
"""Append a finished native replay's exact evidence to verified source inputs."""
from pathlib import Path
from fractions import Fraction
import argparse
import hashlib
import json

ROOT = Path(__file__).resolve().parents[2]


def sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b''):
            digest.update(chunk)
    return digest.hexdigest()


parser = argparse.ArgumentParser()
parser.add_argument('--original-allowlist', type=Path, required=True)
parser.add_argument('--replay-result', type=Path, required=True)
parser.add_argument('--replay-acceptance', type=Path, required=True)
parser.add_argument('--output', type=Path, required=True)
args = parser.parse_args()
original_path = args.original_allowlist.resolve()
result_path = args.replay_result.resolve()
acceptance_path = args.replay_acceptance.resolve()
output_path = args.output.resolve()
assert not output_path.exists() and output_path.is_relative_to(ROOT)
original = json.loads(original_path.read_text())
result = json.loads(result_path.read_text())
acceptance = json.loads(acceptance_path.read_text())
assert original['status'] == 'PASS_VERIFIED_FROZEN_PROOF_AND_COMPLETE_HISTORICAL_SOURCE_ALLOWLIST'
assert result['status'] == 'PASS_PORTABLE_COMPLETE_NEW_ROUND_GLOBAL_REPLAY'
assert result['all_local_replayed'] is True
assert acceptance['status'] == 'PASS_FRESH_END_TO_END_FINAL_GOAL_ALL_LOCAL_REPLAY'
assert acceptance['all_replay_stages'] == len(result['jobs']) == 219
assert acceptance['new_local_parts_replayed'] == 180
assert acceptance['fresh_default_numeric_geometry_upper_jobs'] == 38
assert Fraction(result['norm_upper_bound']) == Fraction(40, 71)
assert Fraction(result['constructive_G_R_lower_bound']) == Fraction(71, 40)
assert acceptance['bindings']['raw_all_local_result']['sha256'] == sha(result_path)
records = dict(original['sha256'])
for relative, digest in records.items():
    path = ROOT / relative
    assert path.resolve().is_relative_to(ROOT) and sha(path) == digest, relative


def add(path, expected=None):
    path = path.resolve()
    assert path.is_relative_to(ROOT) and path.is_file()
    relative = str(path.relative_to(ROOT))
    digest = sha(path)
    assert expected is None or digest == expected, relative
    assert relative not in records or records[relative] == digest
    records[relative] = digest


for binding in acceptance['bindings'].values():
    add(ROOT / binding['path'], binding['sha256'])
assert result['jobs'][0]['status'] == 'PASS_FRESH_PORTABLE_BASELINE'
assert result['jobs'][0]['all_local'] is True
for job in result['jobs'][1:]:
    path = Path(job['report'])
    assert job['status'].startswith('PASS')
    assert json.loads(path.read_text())['status'] == job['status']
    add(path, job['sha256'])
add(result_path.parent / 'baseline.log')
add(original_path)
add(result_path)
add(acceptance_path)
add(Path(__file__))
final = dict(original)
final['sha256'] = dict(sorted(records.items()))
final['file_count'] = len(records)
final['total_bytes'] = sum((ROOT / relative).stat().st_size for relative in records)
final['original_archive_allowlist'] = {
    'path': str(original_path.relative_to(ROOT)), 'sha256': sha(original_path)}
final['fresh_native_all_local_evidence'] = {
    'acceptance': str(acceptance_path.relative_to(ROOT)),
    'sha256': sha(acceptance_path), 'stages': 219,
    'exact_job_report_files': len(result['jobs']) - 1,
    'baseline_full_log_included': True,
    'runtime_or_temporary_directory_traversal_used': False,
    'execution_evidence_added_to_immutable_proof_manifest': False,
    'docker_execution_tested': False,
}
output_path.parent.mkdir(parents=True, exist_ok=True)
output_path.write_text(json.dumps(final, indent=2) + '\n')
print(json.dumps({'path': str(output_path.relative_to(ROOT)),
                  'sha256': sha(output_path), 'files': len(records),
                  'bytes': final['total_bytes']}, indent=2))
