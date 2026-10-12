# Grothendieck verification repository

Repository: [shlui/Grothendieck](https://github.com/shlui/Grothendieck).

This directory supplies the fixed inputs, checking programs and saved
certificates for the real Grothendieck lower bound. It has the four folders
`data`, `inv`, `s`, `tools`, and nine root files. The current main paper and
the current Verification document are distributed separately. Neither has
a publication folder or a current source/PDF copy in this repository delivery.

## Files

| Folder | Contents |
|---|---|
| `data/` | Indexed frozen inputs and recovered historical verification material. |
| `inv/` | Byte inventories, audit records and repository layout checks. |
| `s/` | Strict supplement programs, support files and saved acceptance records. |
| `tools/` | Replay driver, original-path lookup and inventory helpers. |

| Root file | Purpose |
|---|---|
| `README.md` | File descriptions, dependencies and replay commands. |
| `PROVENANCE.md` | Origin and scope of the supplied mathematical evidence. |
| `CHECKSUMS.sha256` | SHA-256 of every repository file except this list itself. |
| `paths.json` | Original logical input labels mapped to stored short filenames. |
| `requirements.txt` | Recorded Python dependencies. |
| `Dockerfile` | Recorded replay environment and input location. |
| `.dockerignore` | Excludes generated outputs and local caches from the image. |
| `.gitattributes` | Preserves exact file bytes during Git checkout and commit. |
| `.gitignore` | Excludes generated replay outputs and Python caches. |

The stored names are short. Original input labels are resolved through
`paths.json`; mathematical inputs are not rewritten. The recovered frozen
material retains historical theory notes and older verification documents
as part of its original byte inventory. These are distinct from the current
paper and the separately supplied current Verification document.

## Run the checks

Run commands from the root of `Grothendieck`. The recorded environment is
Python 3.12 with `python-flint==0.9.0`, `gmpy2==2.3.1` and `scipy==1.17.0`,
using GMP and MPFR. The Docker specification includes the required system
libraries.

```bash
python -m pip install -r requirements.txt
python tools/replay.py --mode inventory
python tools/replay.py --mode strict
python tools/replay.py --mode full
```

| Mode | Scope |
|---|---|
| `inventory` | Checks 3300 manifest entries, representing 3292 distinct frozen inputs, and all recovered support files; no mathematical replay. |
| `strict` | Runs 31 directed price/plane jobs and six closed-cover geometry stages, followed by uncapped rational assembly, helper fixtures and finite-sequence arithmetic. Accepted local and whole-line comparisons remain premises. |
| `full` | Adds the protected baseline, 180 local component jobs and the historical upper/rounding job: 219 top-level stages, plus the strict supplemental checks. |

Use `--dry-run` to inspect a plan without executing its numerical stages,
`--output DIR` to choose a fresh empty output directory, and `--workers N`
to choose the worker count. Outputs include `REPLAY_RESULT.json`;
full mode also writes `FULL_LOCAL_REPLAYS.json`. Fresh strict acceptance
is under `strict/STRICT_SEQUENCE_ACCEPTANCE.json` in the output directory.
Without `--output`, the driver uses a timestamped `replay_outputs/` directory.

The inventory mode checks the fixed mathematical inputs and support files.
For all root files and current documentation as well, check the separate
repository-wide list with `sha256sum -c CHECKSUMS.sha256` on Linux.

## Docker

The accepted inputs record `/workspace/scratch/8242d848994f`. Replays restore
their original paths there without changing their bytes. The location must
be writable. The supplied Dockerfile creates it inside the container, so
the host does not need write access to `/workspace`.

```bash
docker build -t gr-proof-replay .
mkdir -p replay_outputs
docker run --rm -v "$PWD/replay_outputs:/output" gr-proof-replay --mode strict --output /output/strict
docker run --rm -v "$PWD/replay_outputs:/output" gr-proof-replay --mode full --output /output/full
```

Each run requires a fresh empty output directory. Container execution was
not repeated for this documentation revision.

## Saved evidence

- [Historical 219-stage acceptance](data/f00010.json)
- [Detailed historical replay result](data/f00147.json)
- [Strict sequence acceptance](s/STRICT_SEQUENCE_ACCEPTANCE.json)
- [Exact first-three arithmetic](s/finite_sequence_arithmetic.json)
- [Current repository layout checks](inv/REPOSITORY_LAYOUT.json)
- [Proof and notation audit](inv/TYPESETTING.json)
- [Verification provenance](PROVENANCE.md)

The modified sequence has certified limiting lower bound
`10^9/563380277 > 71/40`. Its first three certified lower bounds, rounded
downward, are `1.775000013258`, `1.775000014017`, and `1.775000014397`.
The matrices are specified by finite formulas; they are not materialized
entry by entry in this repository.

Reading saved reports does not execute a new verification. The full and
strict numerical workflows were not rerun for this documentation and
layout revision. Its fresh checks concern file identities, repository
structure, inventory mode and the two numerical-stage plans.
