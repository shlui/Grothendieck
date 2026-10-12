#!/usr/bin/env python3
"""Freeze a compact record only after the complete all-local run succeeds."""
from pathlib import Path
from fractions import Fraction as F
import argparse
import hashlib
import json

ROOT = Path(__file__).resolve().parents[3]
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
read = lambda p: json.loads(p.read_text())

parser = argparse.ArgumentParser()
parser.add_argument('--result', type=Path, required=True)
parser.add_argument('--manifest', type=Path, required=True)
parser.add_argument('--output', type=Path, required=True)
args = parser.parse_args()
result_path, manifest_path, output_path = (
    p.resolve() for p in (args.result, args.manifest, args.output))
assert not output_path.exists()
result, manifest = read(result_path), read(manifest_path)
config_path = ROOT / manifest['config']
config = read(config_path)
global_path = ROOT / config['accepted_global']
global_record = read(global_path)
assert result['status'] == 'PASS_PORTABLE_COMPLETE_NEW_ROUND_GLOBAL_REPLAY'
assert result['all_local_replayed'] is True
assert result['certificate_bytes_rewritten'] is False
assert result['operator_b'] == config['operator_b'] == global_record['operator_b'] == '711/1000'
assert F(result['norm_upper_bound']) == F(config['norm_cap']) == F(global_record['norm_upper_bound']) == F(40, 71)
assert F(result['constructive_G_R_lower_bound']) == F(71, 40)
assert result['three_significant_digit_rounding'] == '1.78'
assert F(result['universal_strict_upper']) == F(1760, 987) < F(357, 200)
assert result['immutable_input_files'] == manifest['file_count'] == len(manifest['sha256'])
assert all(sha(ROOT / path) == digest for path, digest in manifest['sha256'].items())
jobs = result['jobs']
assert len(jobs) == 1 + len(config['all_local_new_gates']) + len(config['jobs'])
assert jobs[0]['job'] == 'protected_baseline' and jobs[0]['all_local'] is True
assert jobs[0]['status'] == 'PASS_FRESH_PORTABLE_BASELINE'
expected = config['all_local_new_gates'] + config['jobs']
assert [job['job'] for job in jobs[1:]] == [job['name'] for job in expected]
for job in jobs[1:]:
    path = Path(job['report'])
    assert sha(path) == job['sha256']
    assert read(path)['status'] == job['status'] and job['status'].startswith('PASS')

def binding(path):
    assert path.is_relative_to(ROOT)
    return {'path': str(path.relative_to(ROOT)), 'sha256': sha(path)}

record = {
    'status': 'PASS_FRESH_END_TO_END_FINAL_GOAL_ALL_LOCAL_REPLAY',
    'operator_b': '711/1000', 'norm_upper_bound': '40/71',
    'norm_cap_relation': result['norm_cap_relation'],
    'constructive_G_R_lower_bound': '71/40',
    'constructive_G_R_lower_bound_decimal': '1.775',
    'universal_strict_upper': '1760/987',
    'three_significant_digit_rounding': '1.78',
    'immutable_input_files': len(manifest['sha256']),
    'protected_baseline_full_replayed': True,
    'new_local_criteria': len(config['all_local_new_gates']) // 3,
    'new_local_parts_replayed': len(config['all_local_new_gates']),
    'fresh_default_numeric_geometry_upper_jobs': len(config['jobs']),
    'all_replay_stages': len(jobs),
    'seconds': result['seconds'],
    'immutable_hashes_rechecked_after_run': True,
    'bindings': {'input_manifest': binding(manifest_path),
                 'replay_config': binding(config_path),
                 'complete_global_acceptance': binding(global_path),
                 'raw_all_local_result': binding(result_path),
                 'entrypoint': binding(Path(__file__).resolve())},
    'scope': 'The complete native Python all-local mode finished. It includes '
        'all default directed numerical, complete closed geometry and exact '
        'upper/rounding jobs, every current new whole-line criterion, and the '
        'protected baseline full local/scalar mode. No optimizer was called. '
        'The fixed Gaussian matrix sequence has actual ratio limit 1/N(.711), '
        'at least 1.775; equality to the cap or unrestricted constant is not asserted.',
    'detailed_reports': 'The raw summary retains SHA256 digests of every fresh '
        'report. Detailed output reports are reproducible from the frozen '
        'input manifest and are not immutable proof inputs.',
    'docker_execution_tested': False,
}
output_path.write_text(json.dumps(record, indent=2) + '\n')
print(json.dumps({'path': str(output_path.relative_to(ROOT)),
                  'sha256': sha(output_path), 'stages': len(jobs)}, indent=2))
