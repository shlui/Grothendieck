#!/usr/bin/env python3
"""Portable orchestration; frozen arithmetic and accepted inputs stay unchanged."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from fractions import Fraction
import argparse
import ctypes.util
import datetime
import hashlib
import importlib.util
import json
import os
import shlex
import shutil
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
FROZEN = ROOT / 'verification/frozen'
STRICT = ROOT / 'verification/strict'
MANIFEST = 'new_round_2026-10-02/audit/portable_replay/FINAL_GOAL_INPUT_MANIFEST.json'
CONFIG = 'new_round_2026-10-02/audit/portable_replay/FINAL_GOAL_REPLAY_CONFIG.json'
RECORDED = Path('/workspace/scratch/8242d848994f')
GEOMETRY = ('strong43_geometry', 'u_refinement_geometry', 'layer_00_geometry',
            'layer_01_geometry', 'layer_02_geometry', 'final_geometry')

def read(p):
    return json.loads(Path(p).read_text())

def sha(p):
    h = hashlib.sha256()
    with Path(p).open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()

def write(p, data):
    Path(p).parent.mkdir(parents=True, exist_ok=True)
    Path(p).write_text(json.dumps(data, indent=2) + '\n')

def package_inventory():
    m = read(FROZEN / MANIFEST)
    assert m['extension_filter_used'] is False
    failures = []
    for name, expected in m['sha256'].items():
        p = FROZEN / name
        if not p.is_file() or sha(p) != expected:
            failures.append(name)
    if failures:
        raise RuntimeError('Frozen manifest mismatch: ' + ', '.join(failures))
    inventory = read(ROOT / 'inventories/VERIFICATION_FILES.json')
    source_inventory = ROOT / 'inventories/SOURCE_ARCHIVE_FILES.json'
    assert sha(source_inventory) == inventory['source_archive_inventory_sha256']
    archive = read(source_inventory)
    assert archive['all_original_files_byte_identical'] is True and archive['file_count'] == 7647
    assert len(archive['files']) == archive['file_count']
    assert len({item['path'] for item in archive['files']}) == archive['file_count']
    for item in archive['files']:
        name = 'verification/frozen/' + item['path']
        assert inventory['sha256'].get(name) == item['sha256'], 'Archive identity differs: ' + name
        assert (ROOT / name).stat().st_size == item['size_bytes'], name
    for name, expected in inventory['sha256'].items():
        p = ROOT / name
        if not p.is_file() or sha(p) != expected:
            raise RuntimeError('Repository verification/support mismatch: ' + name)
    config = read(FROZEN / CONFIG)
    assert config['original_root'] == str(RECORDED)
    assert len(config['all_local_new_gates']) == 180 and len(config['jobs']) == 38
    for job in config['all_local_new_gates'] + config['jobs']:
        assert (FROZEN / job['script']).is_file(), job['script']
    result = dict(status='PASS_PACKAGED_VERIFICATION_INVENTORY',
                  original_archive_files_checked=archive['file_count'],
                  frozen_manifest_entries=len(m['sha256']),
                  distinct_frozen_manifest_files=len({(FROZEN / n).resolve() for n in m['sha256']}),
                  packaged_verification_files=len(inventory['sha256']),
                  frozen_manifest_sha256=sha(FROZEN / MANIFEST),
                  package_inventory_sha256=sha(ROOT / 'inventories/VERIFICATION_FILES.json'),
                  mathematical_replay_performed=False)
    return m, config, inventory, result

def dependencies(mode):
    required = ['flint', 'gmpy2', 'scipy'] if mode == 'full' else []
    missing = [n for n in required if importlib.util.find_spec(n) is None]
    if not ctypes.util.find_library('mpfr'):
        missing.append('system libmpfr')
    return dict(required_python_modules=required, missing=missing,
                mpfr_library=ctypes.util.find_library('mpfr'),
                python=sys.version.split()[0])

def plan(mode, config, workers):
    stages = []
    if mode == 'full':
        stages.append(dict(name='protected_baseline', all_local=True,
                           command=['python', '<disposable_baseline>/verify_bundle.py', '--full', '--workers', str(workers)]))
        stages += [dict(name=j['name'], script=j['script'], args=j['args'],
                        expected_status=j['expected_status'], accepted=j['accepted'])
                   for j in config['all_local_new_gates']]
    stages += [dict(name=j['name'], script=j['script'], args=j['args'])
               for j in config['jobs'] if '/audit/numeric/' in j['script']]
    stages += [dict(name=j['name'], script=j['script'], args=j['args'],
                   compatibility_copy=j['name'] in GEOMETRY[3:])
               for j in config['jobs'] if j['name'] in GEOMETRY]
    if mode == 'full':
        stages += [dict(name=j['name'], script=j['script'], args=j['args'])
                   for j in config['jobs'] if j['name'] == 'three_digit_rounding']
    assert len(stages) == (219 if mode == 'full' else 37)
    return dict(mode=mode, top_level_stages=len(stages), stages=stages,
                supplemental_checks=['exact uncapped 305039-leaf assembly',
                                     'finite modified-sequence arithmetic',
                                     'independent affine helper fixtures',
                                     'fresh report/hash acceptance'],
                materialization_root=str(RECORDED),
                source_root=str(FROZEN),
                full_local_replay=mode == 'full',
                mathematical_replay_performed=False)

def materialize(inventory):
    """Restore recorded paths without ever accepting different existing bytes."""
    owned = {n: h for n, h in inventory['sha256'].items()
             if n.startswith('verification/frozen/')}
    # Check all existing targets before copying any absent files.
    for name, expected in owned.items():
        relative = Path(name).relative_to('verification/frozen')
        target = RECORDED / relative
        if target.exists() and (not target.is_file() or sha(target) != expected):
            raise RuntimeError('Recorded-root file differs; use the provided clean Docker container: ' + str(target))
    for name in owned:
        relative = Path(name).relative_to('verification/frozen')
        target = RECORDED / relative
        if not target.exists():
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(ROOT / name, target)

def assert_materialized(m):
    for name, expected in m['sha256'].items():
        assert sha(RECORDED / name) == expected, 'Frozen input changed: ' + name

def execute(command, log, env, cwd):
    with Path(log).open('w') as stream:
        result = subprocess.run(command, cwd=cwd, env=env, stdout=stream, stderr=subprocess.STDOUT)
    if result.returncode:
        raise RuntimeError('Replay failed; inspect ' + str(log))

def full_local(config, out, env, workers):
    results = []
    # Baseline has mutable producers, so every nested baseline replay is disposable.
    with tempfile.TemporaryDirectory(prefix='grothendieck_baseline_') as temporary:
        work = Path(temporary) / 'proof'
        shutil.copytree(RECORDED / config['baseline'], work,
                        ignore=shutil.ignore_patterns('replay_outputs*', '__pycache__'))
        execute([sys.executable, str(work / 'verify_bundle.py'), '--full', '--workers', str(workers)],
                out / 'protected_baseline.log', env, work)
        shutil.copytree(work / 'replay_outputs', out / 'protected_baseline')
        baseline_report = read(out / 'protected_baseline/result.json')
        assert baseline_report['scalar_local_primitives_freshly_replayed_this_run'] is True
        results.append(dict(name='protected_baseline', status=baseline_report['status'],
                            report=str(out / 'protected_baseline/result.json'),
                            sha256=sha(out / 'protected_baseline/result.json')))
    def run(job):
        def arg(v):
            if '{output}' in v:
                return v.replace('{output}', str(out))
            return str(RECORDED / v) if v.startswith(('new_round_', 'verified_')) else v
        log = out / (job['name'] + '.log')
        execute([sys.executable, str(RECORDED / job['script']), *map(arg, job['args'])], log, env, RECORDED)
        source = Path(arg(job['result']))
        report = read(source)
        assert report['status'] == job['expected_status'] and report.get('accepted') == job['accepted']
        target = out / (job['name'] + '.json')
        if target != source:
            shutil.copy2(source, target)
        answer = dict(name=job['name'], status=report['status'], report=str(target), sha256=sha(target))
        print(json.dumps(answer), flush=True)
        return answer
    with ThreadPoolExecutor(workers) as pool:
        results += list(pool.map(run, config['all_local_new_gates']))
    write(out / 'FULL_LOCAL_REPLAYS.json', dict(status='PASS_FRESH_ALL_LOCAL_REPLAYS',
                                             top_level_stages=181, results=results))

def strict_replay(out, env):
    destination = out / 'strict'
    shutil.copytree(STRICT, destination, ignore=shutil.ignore_patterns('__pycache__'))
    # Each script writes new reports beside this copy. Archived repository reports remain untouched.
    for script in ('replay_all_numeric.py', 'replay_all_geometry.py', 'verify_strict_acceptance.py'):
        print('START', script, flush=True)
        execute([sys.executable, str(destination / script)], out / (script + '.log'), env, RECORDED)
        print('PASS', script, flush=True)
    acceptance = read(destination / 'STRICT_SEQUENCE_ACCEPTANCE.json')
    assert acceptance['status'] == 'PASS_STRICT_MODIFIED_SEQUENCE_CERTIFICATE'
    assert acceptance['fresh_directed_price_jobs'] == 31 and acceptance['fresh_closed_geometry_stages'] == 6
    return destination

def rounding(config, out, env):
    job = next(j for j in config['jobs'] if j['name'] == 'three_digit_rounding')
    def arg(v):
        if '{output}' in v:
            return v.replace('{output}', str(out))
        return str(RECORDED / v) if v.startswith(('new_round_', 'verified_')) else v
    execute([sys.executable, str(RECORDED / job['script']), *map(arg, job['args'])],
            out / 'three_digit_rounding.log', env, RECORDED)
    report = read(out / 'three_digit_rounding.json')
    assert report['three_significant_digit_rounding_proved'] is True
    assert report['rounding_result'] == '1.78'
    assert Fraction(report['universal_strict_upper']) == Fraction(config['universal_strict_upper'])

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--mode', choices=('inventory', 'strict', 'full'), default='inventory')
    ap.add_argument('--dry-run', action='store_true', help='Check packaged bytes and print complete stage plan; execute no mathematical replay.')
    ap.add_argument('--workers', type=int, default=4, help='Workers for full local and protected baseline replay; strict numeric script uses four.')
    ap.add_argument('--output', type=Path)
    ap.add_argument('--check-dependencies', action='store_true', help='Fail if mathematical replay dependencies are unavailable.')
    args = ap.parse_args()
    if args.workers < 1:
        ap.error('--workers must be positive')
    m, config, inventory, checked = package_inventory()
    if args.mode == 'inventory':
        if args.output:
            write(args.output / 'INVENTORY_RESULT.json', checked)
        print(json.dumps(checked, indent=2))
        return
    dep = dependencies(args.mode)
    if args.dry_run:
        p = plan(args.mode, config, args.workers)
        p.update(inventory=checked, dependencies=dep)
        if args.output:
            write(args.output / 'REPLAY_PLAN.json', p)
        print(json.dumps(p, indent=2))
        if args.check_dependencies and dep['missing']:
            raise RuntimeError('Missing dependencies: ' + ', '.join(dep['missing']))
        return
    if dep['missing']:
        raise RuntimeError('Missing dependencies: ' + ', '.join(dep['missing']) + '; install requirements.txt or use Dockerfile.')
    stamp = datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    out = (args.output or ROOT / 'replay_outputs' / stamp).resolve()
    if out.exists() and any(out.iterdir()):
        raise RuntimeError('Choose a fresh empty output directory: ' + str(out))
    out.mkdir(parents=True, exist_ok=True)
    write(out / 'REPLAY_PLAN.json', plan(args.mode, config, args.workers))
    materialize(inventory)
    assert_materialized(m)
    env = dict(os.environ, GR_FROZEN_BASELINE=str(RECORDED / config['baseline']),
               GR_MPFR_BITS='768', PYTHONDONTWRITEBYTECODE='1', PYTHONINTMAXSTRDIGITS='0')
    # Historical strict scripts call "python"; make it the same selected interpreter.
    shim = out / 'python_bin'
    shim.mkdir()
    launcher = shim / 'python'
    launcher.write_text('#!/bin/sh\nexec ' + shlex.quote(sys.executable) + ' "$@"\n')
    launcher.chmod(0o755)
    env['PATH'] = str(shim) + os.pathsep + env.get('PATH', '')
    if args.mode == 'full':
        full_local(config, out, env, args.workers)
    fresh = strict_replay(out, env)
    if args.mode == 'full':
        rounding(config, out, env)
        # The standalone strict checker describes its own targeted scope. In a
        # complete run, add the separately checked local-stage context to this
        # fresh report, without changing any arithmetic or archived report.
        acceptance = read(fresh / 'STRICT_SEQUENCE_ACCEPTANCE.json')
        acceptance['scope'] = ('All protected baseline and180 new local gate stages were freshly replayed by the portable orchestrator, followed by31 directed price jobs,6 exact geometry stages and the historical upper/rounding job. This strict acceptance additionally binds the uncapped complete cover, independent helper fixtures and finite modified-sequence arithmetic. No frozen input changed; no original byte identity is claimed for the reconstructed helper.')
        acceptance['full_replay_context_added_by_portable_orchestrator'] = {
            'local_stage_report': str(out / 'FULL_LOCAL_REPLAYS.json'),
            'local_stage_report_sha256': sha(out / 'FULL_LOCAL_REPLAYS.json'),
            'upper_rounding_report': str(out / 'three_digit_rounding.json'),
            'upper_rounding_report_sha256': sha(out / 'three_digit_rounding.json'),
            'arithmetic_changed': False,
        }
        write(fresh / 'STRICT_SEQUENCE_ACCEPTANCE.json', acceptance)
    assert_materialized(m)
    package_inventory()
    result = dict(status='PASS_PORTABLE_' + args.mode.upper() + '_REPLAY', mode=args.mode,
                  top_level_stages=219 if args.mode == 'full' else 37,
                  full_local_and_whole_line_replayed=args.mode == 'full',
                  strict_certificate=str(fresh / 'STRICT_SEQUENCE_ACCEPTANCE.json'),
                  strict_certificate_sha256=sha(fresh / 'STRICT_SEQUENCE_ACCEPTANCE.json'),
                  frozen_manifest_entries=len(m['sha256']), immutable_inputs_unchanged=True,
                  independent_reconstructed_affine_helper=True,
                  absent_original_helper_byte_identity_claimed=False,
                  portable_orchestrator_sha256=sha(Path(__file__)),
                  original_archive_byte_check_performed=False,
                  scope=('Complete protected baseline and all180 new local gate parts freshly replayed, then31 numeric,6 geometry and1 upper/rounding stages, followed by strict assembly and finite checks.'
                         if args.mode == 'full' else
                         '31 numeric and6 geometry stages freshly replayed, followed by strict assembly and finite checks. Frozen accepted whole-line/local inequalities remain premises.'))
    write(out / 'REPLAY_RESULT.json', result)
    print(json.dumps(result, indent=2))

if __name__ == '__main__':
    try:
        main()
    except (AssertionError, RuntimeError, OSError) as exc:
        print('ERROR:', exc, file=sys.stderr)
        sys.exit(1)
