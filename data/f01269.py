#!/usr/bin/env python3
"""Package a frozen source allowlist and verify extracted bytes and ZIP CRCs."""
from pathlib import Path
import argparse
import hashlib
import json
import zipfile

ROOT = Path(__file__).resolve().parents[2]


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


ap = argparse.ArgumentParser()
ap.add_argument('--allowlist', type=Path, required=True)
ap.add_argument('--archive', type=Path, required=True)
ap.add_argument('--extract', type=Path, required=True)
ap.add_argument('--report', type=Path, required=True)
a = ap.parse_args()
allow = a.allowlist.resolve()
d = json.loads(allow.read_text())
assert d['status'] == 'PASS_VERIFIED_FROZEN_PROOF_AND_COMPLETE_HISTORICAL_SOURCE_ALLOWLIST'
records = dict(d['sha256'])
records[str(allow.relative_to(ROOT))] = sha(allow)
assert not a.archive.exists() and not a.extract.exists() and not a.report.exists()
for rel, digest in records.items():
    path = ROOT / rel
    assert path.resolve().is_relative_to(ROOT) and sha(path) == digest, rel
with zipfile.ZipFile(a.archive, 'w', compression=zipfile.ZIP_DEFLATED,
                     compresslevel=6, allowZip64=True) as z:
    for rel in sorted(records):
        z.write(ROOT / rel, arcname=rel)
with zipfile.ZipFile(a.archive) as z:
    assert len(z.namelist()) == len(records)
    assert set(z.namelist()) == set(records)
    assert z.testzip() is None
    z.extractall(a.extract)
for rel, digest in records.items():
    assert sha(a.extract / rel) == digest, rel
out = {'status': 'PASS_ALL_ARCHIVE_CRCS_AND_EXTRACTED_SOURCE_HASHES',
       'archive': str(a.archive.resolve()), 'archive_sha256': sha(a.archive),
       'archive_bytes': a.archive.stat().st_size, 'member_count': len(records),
       'allowlist_sha256': sha(allow), 'extracted_to': str(a.extract.resolve()),
       'certificate_bytes_rewritten': False}
a.report.write_text(json.dumps(out, indent=2) + '\n')
print(json.dumps(out, indent=2))
