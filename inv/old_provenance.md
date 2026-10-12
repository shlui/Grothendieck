# Verification provenance

The current 19-page paper is set in A4 at 11 points with 25 mm margins. It uses weighted Gaussian integrals and matrix notation throughout, defines tensors, contractions, `Sym`, Frobenius and operator norms, and positive-semidefinite comparisons, and cites the established Gaussian characterization instead of reproving it. The source upper-bound derivative argument uses complete monotonicity, the conditional residual estimate follows from two nonnegative squares, and the stronger cubic-shape restriction has a complete commuting-projection proof. Exact `PAQ` invariance and the universal rational balanced search-class corollary are included in the main text. Appendices A–D are subsections; the bibliography is alphabetical, and journal references omit arXiv identifiers when a published version is cited. The latest editorial pass also shortens headings and describes the Gaussian integration-by-parts and product-expansion calculations directly. Equation (1) uses automatically sized absolute-value delimiters. Mathematical displays, numerical constants, and reference-entry contents are otherwise unchanged.

These changes preserve the numerical constants, checker code, certificate premises, and frozen verification inputs. Neither a full nor a strict mathematical replay was performed again for this final editorial revision. The historical complete record and the earlier targeted strict record below retain their respective scopes.

## Correspondence with frozen notation

The manuscript's shape coordinate `σ = 9‖U‖₂²/κ − 1` is the literal field or variable `z` in the frozen certificates. Its alternative shape coordinate satisfies `σ = 2/5 + 3s/5`. The manuscript's scalar `χ` in the analytic head calculation, defined by `s = χ(A+Bχ)²`, is denoted `w` in those inputs and programs. These are notation changes, not replacements or transformations of the stored data.

For the longitudinal coordinate, the raw basis is `H = (1, z, z²−1, z³−3z)ᵀ` and the orthonormal basis is `Ĥ = (1, z, (z²−1)/√2, (z³−3z)/√6)ᵀ`. With `D_H = diag(1, 1, √2, √6)`, the correspondence is `H = D_H Ĥ` and `m̂ = D_H⁻¹m`. If `N₀` has the raw even block `E` and odd block `R`, its orthonormal quadratic matrix is `M = D_H N₀ D_H`. The orthonormal parameters satisfy `β̂ = √6β` and `ω̂ = √2ω`, while `ℓ` and `δ` are unchanged. In particular, the penalty `6β²/b` in raw coordinates is `β̂²/b` in orthonormal coordinates. The raw residual ellipsoid and its rank-one image remain in raw coordinates. This explicit correspondence preserves the accepted inequalities, their numerical normalization, and every frozen byte identity.

## Frozen proof inputs

`verification/frozen/` contains all 7647 files recovered from the original source archive, including its historical output reports and older drafts. The current manuscript is under `paper/`.

The 3300 proof-input manifest entries, covering 3292 distinct files, are listed in
`new_round_2026-10-02/audit/portable_replay/FINAL_GOAL_INPUT_MANIFEST.json`.
Its SHA-256 is:

```text
2e49c298fc886e3fd1f4cc058f5a5194bb25ad616b4fadf9e0a194503ac2758a
```

The authoritative original archive has 439387711 bytes and SHA-256:

```text
35f1119b00800862571051990e2e89eb97f6d32cbcf28d69842bccb4374b4d33
```

The manifest is checked by `python tools/replay.py --mode inventory`. It identifies proof inputs by their exact bytes, including compressed files. [SOURCE_ARCHIVE_FILES.json](inventories/SOURCE_ARCHIVE_FILES.json) records the restored archive-file inventory, and [VERIFICATION_FILES.json](inventories/VERIFICATION_FILES.json) records the packaged verification and support inventory. The new driver and compatibility support are separate from the frozen inputs. Recorded absolute paths are preserved by materialization at `/workspace/scratch/8242d848994f`.

## Two different replay records

The original complete `--all-local` run has a saved record of 219 successful top-level stages. Those stages comprise one protected full baseline, 180 new local components, and 38 default jobs. Its recovered evidence is:

- [Complete all-local acceptance](verification/frozen/new_round_2026-10-02/audit/FINAL_GOAL_ALL_LOCAL_REPLAY_ACCEPTANCE.json)
- [Detailed replay result and report bindings](verification/frozen/new_round_2026-10-02/audit/final_goal_all_local_replay/REPLAY_RESULT.json)

These records are historical evidence; they do not claim that all 219 stages were freshly repeated during this revision. They retain the former coefficient `711/1000`, coarse norm cap `40/71`, and historical universal upper bound `1760/987`. The revised paper uses the strict modified lower bound and Peng's external upper bound `1.7813`. The upper-bound job in a full historical replay verifies the former certificate; it does not verify Peng's literature theorem.

The earlier strict revision freshly replayed all 31 current directed price/plane jobs at 768-bit precision and six current closed-cover geometry stages. It also checked every manifest input, reconstructed the exact uncapped 305039-leaf upper-cover assembly, ran the independent helper fixtures, and recomputed the first three matrix certificates with rational arithmetic. Local and whole-line comparisons from the accepted full run remain inherited premises in that targeted replay.

`--mode full` requests the protected full baseline and 180 local components (181 stages), the 31 directed price/plane and six geometry jobs, and one historical upper/rounding job, totaling 219 stages. The strict uncapped assembly, helper fixtures, finite arithmetic, and acceptance checks are additional. A successful fresh run is established by that run's generated reports. Full mode was not rerun during packaging, and neither full nor strict mode was rerun during this final editorial revision. Reading the packaged saved reports alone does not execute either mode.

The supplied strict evidence is under `verification/strict/`:

| Record | What it establishes |
|---|---|
| `ALL_DEFAULT_NUMERIC_REPLAYS.json` | The 31 directed price/plane jobs and their report bindings. |
| `strong43_geometry.json`, `u_refinement_geometry.json`, `layer_00_geometry.json`, `layer_01_geometry.json`, `layer_02_geometry.json`, `final_geometry.json` | Six complete geometry stages. |
| `STRICT_UNCAPPED_ASSEMBLY.json` | Exact maximum across the retained upper-cover leaves. |
| `AFFINE_HELPER_MANUAL_FIXTURES.json` | Independent helper checks, including ties and degenerate polygons. |
| `finite_sequence_arithmetic.json` | Exact first-three bounds and grid-radius certificates. |
| `STRICT_SEQUENCE_ACCEPTANCE.json` | Final bindings and scope of the strict sequence certificate. |

The saved reports retain the paths of the machines on which they were produced. A new replay writes new paths and bindings into its own output directory: `REPLAY_RESULT.json`, `strict/STRICT_SEQUENCE_ACCEPTANCE.json`, and, for full mode, `FULL_LOCAL_REPLAYS.json`.

## Packaging validation

[BUILD_VALIDATION.json](inventories/BUILD_VALIDATION.json) records `PASS_PORTABLE_PACKAGE_BUILD_SMOKE`. The reproducible checker is `tools/package_smoke.py`. It verifies all 7647 original archive-file identities, all 3300 active manifest entries, and the packaged support inventory. It also checks that packaged source lookup succeeds while reads from the historical root are blocked, and that a deliberately changed packaged input is rejected even when the historical root exists.

The validation checks the complete [219-stage full plan](inventories/FULL_REPLAY_PLAN.json), the [37-stage strict plan](inventories/STRICT_REPLAY_PLAN.json), the presence of all 219 saved historical stage records, and dependency diagnostics. The final support inventory contains 7647 frozen archive files and 108 strict files. Its largest file is 98282817 bytes, below 100 MiB.

These are packaging and planning checks. They execute neither a new complete 219-stage mathematical replay nor a new 31-job strict mathematical replay. Docker execution is also outside this validation's scope. The supplied targeted strict mathematical reports described above were produced before packaging.

## Independent affine helper and compatibility patch

The recovered archive omitted `audit_affine_hybrid_geometry.py`, imported by later geometry checkers. The strict supplement supplies a new independent implementation in `verification/strict/portable_geometry/`, reconstructed from the paper's affine-envelope proof. It recomputes source-knot partitions, native chords, quadratic secants or tangents, closed dominance polygons, and rational vertex maxima. Its actual SHA-256 is:

```text
91225406ef522137d4813e1cced709192734e306e0c96f376ffe571131eebc31
```

This is the new helper's identity. It is not the absent original helper's digest and does not assert original byte identity. `NEW_INDEPENDENT_AFFINE_HELPER_BINDING.json` documents the reconstruction; a separate review and manual fixtures check its mathematical behavior.

An obsolete predecessor, `audit_monotone_rank_children.py`, was also absent. A copied moment checker now hashes predecessor files only when present; the active `audit_monotone_rank_children_v2.py` digest still must match its fresh parent report. No numerical or geometric comparison is changed. `GEOMETRY_COMPATIBILITY_PATCH.json` records the exact removed and added lines and both checker hashes. Later copied checkers bind the patched checker and new helper by their actual hashes. Frozen proof inputs are unchanged.

## Strict arithmetic result

The uncapped retained upper-cover maximum is

```text
281690137769807/500000000000000 = 0.563380275539614 < 40/71.
```

The lower-amplitude cap is smaller. Rounding the norm cap upward to `563380276/10^9` and changing the cubic coefficient from `711/1000` to `711000001/10^9` increases the cap by at most `10^−9`, since `||P3||_(L∞→L1) ≤ 1`. The modified cap is therefore `563380277/10^9`, with reciprocal strictly above `71/40`.

The certified sequence limit is a lower bound for the real Grothendieck constant. The certificate does not identify the true operator norm with its cap or identify this sequence's limit with the unrestricted constant. The literature upper bound used in the paper is cited separately.
