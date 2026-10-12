# Verification provenance

## Current repository layout

The repository root `Grothendieck` has the four folders `data`, `inv`, `s`,
`tools`, and the nine root files listed in `README.md`. The current main
paper and current Verification source/PDF are distributed separately.
The previous `paper/` and `docs/` publication folders are excluded from
this delivery, and their unused entries are removed from `paths.json`.
All fixed numerical inputs, original checker programs and saved
mathematical certificates retain their previous bytes. Git attributes
preserve the exact bytes of all stored files.

## Frozen inputs

`data/` retains the 7647 recovered files from
`Grothendieck_reduced_sequence_source_2026-10-01.zip`, including historical
reports and theory notes. The original archive has 439387711 bytes and
SHA-256:

```text
35f1119b00800862571051990e2e89eb97f6d32cbcf28d69842bccb4374b4d33
```

The proof-input manifest is [data/f00878.json](data/f00878.json), with
original label
`new_round_2026-10-02/audit/portable_replay/FINAL_GOAL_INPUT_MANIFEST.json`.
It has 3300 entries covering 3292 distinct frozen files and SHA-256:

```text
2e49c298fc886e3fd1f4cc058f5a5194bb25ad616b4fadf9e0a194503ac2758a
```

[SOURCE_ARCHIVE_FILES.json](inv/SOURCE_ARCHIVE_FILES.json) records every
recovered original file. [VERIFICATION_FILES.json](inv/VERIFICATION_FILES.json)
records the fixed verification and strict support inventory. Both remain
unchanged. Inventory mode validates them and every corresponding file.
The outer replay and lookup programs are separate from those frozen inputs.

## Historical full evidence

The saved full run has 219 successful top-level stages: a protected
baseline, 180 local components and 38 default jobs. Its evidence is
[the complete acceptance](data/f00010.json) and
[the detailed result](data/f00147.json).
These retain the historical coefficient `711/1000`, coarse norm cap
`40/71` and historical upper bound `1760/987`. The historical upper-bound
stage does not verify the literature upper bound quoted in the current
paper. The saved full result is inherited evidence, not a fresh run of
this delivery.

## Strict modified-sequence evidence

The earlier targeted strict replay ran all 31 directed price/plane jobs
at 768-bit precision and six closed-cover geometry stages. It reconstructed
the exact uncapped 305039-leaf upper-cover assembly, checked independent
affine fixtures and recomputed the first-three certificates with rational
arithmetic. Previously accepted local and whole-line comparisons remain
premises in that strict workflow.

| Saved record in `s/` | Scope |
|---|---|
| `ALL_DEFAULT_NUMERIC_REPLAYS.json` | Directed price/plane jobs and report bindings. |
| `strong43_geometry.json`, `u_refinement_geometry.json`, `layer_00_geometry.json`, `layer_01_geometry.json`, `layer_02_geometry.json`, `final_geometry.json` | Closed geometry checks. |
| `STRICT_UNCAPPED_ASSEMBLY.json` | Exact retained upper-cover maximum. |
| `AFFINE_HELPER_MANUAL_FIXTURES.json` | Independent affine checks, ties and degenerate polygons. |
| `finite_sequence_arithmetic.json` | Exact first-three bounds and grid-radius arithmetic. |
| `STRICT_SEQUENCE_ACCEPTANCE.json` | Final report bindings and strict acceptance scope. |

The strict modified limit is `10^9/563380277 > 71/40`. Its strict margin
follows from the retained exact upper-cover value
`281690137769807/500000000000`, not from a rounded decimal.

## Fresh checks for this delivery

Neither full nor strict mathematical replay was repeated for the current
documentation revision. Inventory mode and both stage plans are checked
in the repository's current four-folder layout. The complete repository
checksum list is reconstructed, every current file is checked, and every
unchanged file is compared by SHA-256 with the previous delivered archive.
The current Verification source compiles independently; its opening names
the main-paper results, and Section E describes this repository directly.
All of its mathematical formulas and result statements are unchanged.

[REPOSITORY_LAYOUT.json](inv/REPOSITORY_LAYOUT.json) records these current
checks. [TYPESETTING.json](inv/TYPESETTING.json) retains the preceding proof
and notation audit and records the current, separately delivered
Verification source/PDF. Earlier inventory, smoke and mathematical records
retain their historical scopes and are not relabelled as fresh numerical
proofs. Publication to GitHub is not performed by this file preparation.
