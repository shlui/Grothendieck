#!/usr/bin/env python3
"""Packaging checks only: no 31-job or 219-stage mathematical replay."""
from pathlib import Path
import hashlib
import importlib.util
import json
import os
import shutil
import tempfile

root = Path(__file__).resolve().parents[1]
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

# Existing historical inputs cannot mask changed packaged bytes. Hard links
# avoid another 1.5GB copy; unlink the test file before tampering with its copy.
with tempfile.TemporaryDirectory(prefix='grothendieck_package_tamper_') as temporary:
    copy = Path(temporary) / 'frozen'
    shutil.copytree(driver.FROZEN, copy, copy_function=os.link)
    name = next(n for n in manifest['sha256'] if '..' not in Path(n).parts)
    target = copy / name
    target.unlink()
    target.write_bytes(b'intentionally invalid package fixture\n')
    saved = driver.FROZEN
    driver.FROZEN = copy
    try:
        driver.package_inventory()
    except RuntimeError as exc:
        assert str(exc) == 'Frozen manifest mismatch: ' + name
    else:
        raise AssertionError('Changed packaged input was accepted')
    finally:
        driver.FROZEN = saved
    assert original_sha(saved / name) == manifest['sha256'][name]

historical = driver.read(driver.FROZEN / 'new_round_2026-10-02/audit/final_goal_all_local_replay/REPLAY_RESULT.json')
assert historical['all_local_replayed'] is True and len(historical['jobs']) == 219
assert all(item['status'].startswith('PASS') for item in historical['jobs'])
assert driver.dependencies('strict')['missing'] == []
files = ('tools/replay.py', 'tools/build_inventory.py', 'tools/package_smoke.py',
         'Dockerfile', 'requirements.txt', 'inventories/VERIFICATION_FILES.json',
         'inventories/SOURCE_ARCHIVE_FILES.json')
answer = {
    'status': 'PASS_PORTABLE_PACKAGE_BUILD_SMOKE',
    'inventory': checked,
    'package_only_inventory_historical_root_reads_blocked': True,
    'changed_packaged_input_rejected_while_historical_root_present': True,
    'full_dry_run_top_level_stages': 219,
    'strict_dry_run_top_level_stages': 37,
    'saved_historical_successful_stage_records': len(historical['jobs']),
    'strict_dependencies': driver.dependencies('strict'),
    'full_dependencies': driver.dependencies('full'),
    'bindings': {name: original_sha(root / name) for name in files},
    'new_complete_math_replay_performed': False,
    'new_31_job_strict_math_replay_performed': False,
    'scope': 'Package byte identities, complete stage plans, source lookup independence, and tamper rejection. Saved mathematical replay reports are included and kept separate; this packaging smoke does not rerun them.'
}
driver.write(root / 'inventories/BUILD_VALIDATION.json', answer)
driver.write(root / 'inventories/FULL_REPLAY_PLAN.json', full)
driver.write(root / 'inventories/STRICT_REPLAY_PLAN.json', strict)
print(json.dumps(answer, indent=2))
