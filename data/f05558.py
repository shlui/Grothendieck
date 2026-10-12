"""Portable optimizer-free replay of one uniform Gaussian R3 source cell."""
from pathlib import Path
import hashlib,json,subprocess,sys

ROOT=Path(__file__).resolve().parent
def run(script):
    subprocess.run([sys.executable,str(ROOT/'code'/script)],cwd=ROOT,check=True)

manifest=ROOT/'MANIFEST.json'
if manifest.exists():
    d=json.loads(manifest.read_text())
    for name,digest in d['immutable_sha256'].items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest, name
run('replay_r3_local_cover.py')
run('replay_r3_bernstein_center_cover.py')
run('assemble_local_gate.py')
run('verify_moving_polygons.py')
run('audit_r3_moving_polygons.py')
print('PASS_PORTABLE_ALL_Q_R3_SOURCE_CELL; GLOBAL_OUTER_COVER_STILL_REQUIRED')
