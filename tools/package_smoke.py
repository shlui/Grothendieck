#!/usr/bin/env python3
"""Compact packaging checks only: no 31-job or 219-stage math replay."""
from pathlib import Path
import importlib.util
import json
import os
import shutil
import tempfile

from package_paths import PackagePaths

root = Path(__file__).resolve().parents[1]
paths = PackagePaths(root)
spec = importlib.util.spec_from_file_location('portable_replay', root / 'tools/replay.py')
driver = importlib.util.module_from_spec(spec)
spec.loader.exec_module(driver)
original_sha = driver.sha

def package_only_sha(p):
    p = Path(p).resolve()
    if p.is_relative_to(driver.RECORDED):
        raise RuntimeError('Packaging smoke attempted historical-root lookup: ' + str(p))
    return original_sha(p)

driver.sha = package_only_sha
manifest, config, inventory, checked = driver.package_inventory()
driver.sha = original_sha
full = driver.plan('full', config, 4)
strict = driver.plan('strict', config, 4)
assert full['top_level_stages'] == 219 and strict['top_level_stages'] == 37

# Hard links avoid a second 1.8GB copy. Unlink the selected member before
# tampering, so every byte of the real compact package stays unchanged.
with tempfile.TemporaryDirectory(prefix='grothendieck_compact_tamper_') as temporary:
    copy = Path(temporary) / 'G'
    shutil.copytree(root, copy, copy_function=os.link,
                    ignore=shutil.ignore_patterns('replay_outputs*', '__pycache__'))
    copied_paths = PackagePaths(copy)
    name = next(n for n in manifest['sha256']
                if sum(paths.frozen(v) == paths.frozen(n) for v in manifest['sha256']) == 1)
    target = copied_paths.frozen(name)
    target.unlink()
    target.write_bytes(b'intentionally invalid compact package fixture\n')
    saved = driver.ROOT, driver.FROZEN, driver.STRICT
    driver.ROOT, driver.FROZEN, driver.STRICT = copy, copy / 'data', copy / 's'
    try:
        driver.package_inventory()
    except RuntimeError as exc:
        assert str(exc) == 'Frozen manifest mismatch: ' + name
    else:
        raise AssertionError('Changed compact packaged input was accepted')
    finally:
        driver.ROOT, driver.FROZEN, driver.STRICT = saved
    assert original_sha(paths.frozen(name)) == manifest['sha256'][name]

guard_cases = {
    'parent_traversal': {'safe': '../outside'},
    'physical_target_collision': {'first': 'data/a', 'second': 'data/a'},
    'physical_case_collision': {'first': 'data/a', 'second': 'data/A'},
    'windows_reserved_name': {'safe': 'data/CON.json'},
    'physical_directory_case_collision': {'first': 'Data/a', 'second': 'data/b'},
}
rejected = []
for label, files in guard_cases.items():
    with tempfile.TemporaryDirectory(prefix='grothendieck_map_guard_') as temporary:
        test_root = Path(temporary)
        (test_root / 'paths.json').write_text(json.dumps({'format': 1, 'files': files}))
        try:
            PackagePaths(test_root)
        except RuntimeError:
            rejected.append(label)
        else:
            raise AssertionError('Invalid compact path map accepted: ' + label)
assert len(rejected) == len(guard_cases)
aliases = [name for name in manifest['sha256'] if '..' in Path(name).parts]
assert len(aliases) == 8
assert all(original_sha(paths.frozen(name)) == manifest['sha256'][name] for name in aliases)
try:
    paths.frozen('../../README.md')
except RuntimeError:
    pass
else:
    raise AssertionError('Frozen logical-root escape was accepted')

historical = driver.read(paths.frozen('new_round_2026-10-02/audit/final_goal_all_local_replay/REPLAY_RESULT.json'))
assert historical['all_local_replayed'] is True and len(historical['jobs']) == 219
assert all(item['status'].startswith('PASS') for item in historical['jobs'])
assert driver.dependencies('strict')['missing'] == []
files = ('tools/replay.py', 'tools/build_inventory.py', 'tools/package_smoke.py',
         'tools/package_paths.py', 'Dockerfile', 'requirements.txt',
         'inventories/VERIFICATION_FILES.json', 'inventories/SOURCE_ARCHIVE_FILES.json')
answer = {
    'status': 'PASS_COMPACT_PORTABLE_PACKAGE_BUILD_SMOKE',
    'inventory': checked,
    'package_only_inventory_historical_root_reads_blocked': True,
    'changed_packaged_input_rejected_while_historical_root_present': True,
    'path_map_guard_cases_rejected': rejected,
    'frozen_manifest_parent_aliases_checked': len(aliases),
    'frozen_logical_root_escape_rejected': True,
    'full_dry_run_top_level_stages': 219,
    'strict_dry_run_top_level_stages': 37,
    'saved_historical_successful_stage_records': len(historical['jobs']),
    'strict_dependencies': driver.dependencies('strict'),
    'full_dependencies': driver.dependencies('full'),
    'bindings': {name: original_sha(paths.resolve(name)) for name in files},
    'logical_path_map_sha256': original_sha(root / 'paths.json'),
    'new_complete_math_replay_performed': False,
    'new_31_job_strict_math_replay_performed': False,
    'scope': 'Compact map validity, original logical byte identities, complete stage plans, source lookup independence, and tamper rejection. Saved mathematical replay reports remain separate; this packaging smoke does not rerun them.'
}
# Preserve every archived inventory/report. Compact-specific new evidence has
# new names instead of overwriting the mapped original BUILD_VALIDATION/plans.
driver.write(root / 'inv/short_build_validation.json', answer)
driver.write(root / 'inv/short_full_plan.json', full)
driver.write(root / 'inv/short_strict_plan.json', strict)
print(json.dumps(answer, indent=2))
