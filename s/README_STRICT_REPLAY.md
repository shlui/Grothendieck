# Strict modified matrix sequence: reproduction

This supplement verifies a strict extension of the recovered proof. It needs
Python 3, MPFR with the dependencies used by the original verification package,
the original source ZIP, and this supplement.

The original source archive has SHA-256:

```
35f1119b00800862571051990e2e89eb97f6d32cbcf28d69842bccb4374b4d33
```

Extract its full contents under the recorded root
`/workspace/scratch/8242d848994f`. Full extraction is required: extracting only
files listed by its manifest omits some auxiliary modules. The recorded root
is required because the frozen inputs contain absolute paths. Do not modify
existing frozen proof inputs. Put this supplement's `strict_norm` directory
at any location accessible to Python; its own output paths are relative to
its scripts.

For example, after setting `GR_ORIGINAL_ZIP` to the original archive and
`GR_STRICT_DIR` to this supplement's `strict_norm` directory:

```bash
mkdir -p /workspace/scratch/8242d848994f
python -m zipfile -e "$GR_ORIGINAL_ZIP" /workspace/scratch/8242d848994f
python "$GR_STRICT_DIR/replay_all_numeric.py"
python "$GR_STRICT_DIR/replay_all_geometry.py"
python "$GR_STRICT_DIR/verify_strict_acceptance.py" --archive "$GR_ORIGINAL_ZIP"
```

`replay_all_numeric.py` freshly executes all 31 directed price/plane jobs
from the original current configuration at 768-bit precision. It preserves
frozen files and writes new reports. `replay_all_geometry.py` executes six
complete current geometry stages, with the transparent compatibility changes
below. `verify_strict_acceptance.py` checks report and input bindings, rechecks
all 3300 frozen manifest files, reruns the exact 305039-leaf uncapped assembly,
runs the independent helper fixtures and computes the first three finite
bounds by exact rational arithmetic. Its final report is
`STRICT_SEQUENCE_ACCEPTANCE.json`.

The complete 219-stage original whole-line/local replay is not repeated here.
Those established analytic certificates and the protected original complete
cover remain inherited premises. All current directed prices and current
closed-cover geometry are freshly checked; no numerical optimizer or sampled
maximum is used to prove an upper bound.

## Transparent archive compatibility changes

The original archive omitted `audit_affine_hybrid_geometry.py`, which two
frozen geometry checkers import. The supplement supplies a **new independent
checker**, `portable_geometry/audit_affine_hybrid_geometry.py`, reconstructed
from the paper's affine-envelope argument. It recomputes native chords,
source-piece coverage, quadratic secants or tangents, closed dominance
regions and every rational vertex value. It does not claim to recover the
absent original module's bytes or digest. `NEW_INDEPENDENT_AFFINE_HELPER_BINDING.json`
records its actual digest. The helper was separately reviewed and its manual
fixtures check ties, segments, points, knot boundaries and malformed inputs.

The archive also omitted an obsolete predecessor named
`audit_monotone_rank_children.py`. The current moment checker originally
hashed both that name and the existing `_v2.py` predecessor unconditionally.
A copy of that checker now hashes only predecessor files that actually exist.
The actual `_v2.py` predecessor digest still must match the fresh parent
report. This is the sole compatibility edit; it changes no mathematical
comparison or arithmetic. `GEOMETRY_COMPATIBILITY_PATCH.json` records the
exact original and added lines and both actual source hashes. The subsequent
copied checkers bind this patched checker and the new helper by their actual
hashes. Frozen source files are never overwritten.

## Mathematical result

Uncapping the accepted rational affine maxima gives the exact complete
upper-cover maximum

```
281690137769807/500000000000000 = 0.563380275539614 < 40/71.
```

The stronger inherited lower-amplitude bound is also below this cap. The
archived acceptance record used `40/71` as an intentionally coarser cap;
its equality did not identify the operator norm.

Round the cap upward to `B0=563380276/10^9` and change the coefficient to
`b*=711000001/10^9`. The elementary contraction
`||P3||_(L∞→L1) ≤ 1` gives the cap `B*=563380277/10^9` for the modified
operator. Consequently its rational reduced matrix sequence has actual
ratio limit at least

```
10^9/563380277 = 1.775000014776875... > 71/40.
```

Use `d_n=10^9*2^(n-1)` for `n>=1` in the paper's logarithmic-cutoff
construction. The first three Gaussian dimensions are `10^9`, `2*10^9`
and `4*10^9`; their grid radii are all 15. These are Gaussian dimensions,
not matrix orders. The matrices are fully specified implicitly, rather than
materialized entry by entry. The finite certificates, rounded downward to
12 decimal places, are respectively

```
1.775000013258
1.775000014017
1.775000014397
```

The exact fractions, radius certificates and arithmetic source appear in
`finite_sequence_arithmetic.json` and `finite_sequence_arithmetic.py`.
