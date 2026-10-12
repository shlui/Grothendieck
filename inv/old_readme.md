# Real Grothendieck constant: paper and verification

This package contains the revised 19-page paper, set in A4 at 11 points with 25 mm margins, and the inputs and programs for its computer-assisted lower bound. The modified rational matrix sequence satisfies

$$
\lim_n G(H_n)\ge\frac{10^9}{563380277}>\frac{71}{40}=1.775.
$$

The paper combines this with Peng's literature upper bound, `1.7813`, to determine the real constant to three significant digits as `1.78`, which appears to be new in the literature.

## Read the paper

- [PDF](paper/Grothendieck_streamlined.pdf)
- [LaTeX source](paper/Grothendieck_streamlined.tex)
- [Verification provenance](PROVENANCE.md)

History, motivation, the established Naor–Regev relation between `α` and `G_R`, matrix construction, and analytical proofs are in the main text. The exposition uses weighted Gaussian integrals and matrix notation throughout, with explicit definitions of tensors, contractions, `Sym`, norms, and matrix comparisons. It includes exact `PAQ` invariance, the corollary that rational balanced matrices form a universal search class, and the complete projection proof of the stronger cubic-shape restriction. Numerical estimates and finite certificate checks are in appendices A–D, formatted as subsections. Appendix D writes out the first three matrices through explicit rational blocks and scalar entry formulas. The bibliography is alphabetical. Published journal references omit their arXiv identifiers; arXiv-only preprints retain them. Headings are concise, and the Gaussian integration-by-parts and product-expansion arguments are described directly.

## Paper notation and frozen certificates

The current notation corresponds to the unchanged frozen data as follows:

| Current paper | Literal frozen notation | Meaning |
|---|---|---|
| `σ = 9‖U‖₂²/κ − 1` | `z` | Cubic moment-shape coordinate; the alternative coordinate `s` satisfies `σ = 2/5 + 3s/5`. |
| `χ` in `s = χ(A+Bχ)²` | `w` | Scalar inverse variable in the analytic head calculation. |

The raw longitudinal Hermites are `H = (1, z, z²−1, z³−3z)ᵀ`. The orthonormal basis is `Ĥ = D_H⁻¹H`, where `D_H = diag(1, 1, √2, √6)`. Accordingly, raw moments `m` become `m̂ = D_H⁻¹m`, and a raw quadratic matrix `N₀` becomes `M = D_H N₀ D_H`. The normalized parameters are `β̂ = √6β` and `ω̂ = √2ω`; `ℓ` and `δ` are unchanged. Thus the raw penalty `6β²/b` is the normalized penalty `β̂²/b`. These correspondences describe the same inequalities and do not alter certificate fields, coefficients, or hashes.

To rebuild the PDF with a standard LaTeX installation:

```bash
cd paper
pdflatex -interaction=nonstopmode -halt-on-error Grothendieck_streamlined.tex
pdflatex -interaction=nonstopmode -halt-on-error Grothendieck_streamlined.tex
cd ..
```

## Run verification

Run from this directory. The recorded dependency environment is Python 3.12 with `python-flint==0.9.0`, `gmpy2==2.3.1`, and `scipy==1.17.0`; the numerical libraries use GMP and MPFR. The frozen Docker specification includes `libgmp10` and `libmpfr6`.

```bash
python -m pip install python-flint==0.9.0 gmpy2==2.3.1 scipy==1.17.0
python tools/replay.py --mode inventory
python tools/replay.py --mode strict
python tools/replay.py --mode full
```

| Mode | Scope |
|---|---|
| `inventory` | Verify all 3300 manifest entries, covering 3292 distinct frozen files, and the packaged support inventory; no mathematical replay. |
| `strict` | Replay 31 directed price/plane jobs and six closed-cover geometry stages, then check uncapped rational maxima, helper fixtures, and finite sequence arithmetic. Previously accepted local and whole-line inequalities remain premises. |
| `full` | Run the protected full baseline and 180 local-component jobs (181 stages), the same 31 price/plane and six geometry jobs, and one historical upper/rounding job: 219 top-level stages. Uncapped assembly, fixtures, finite arithmetic, and acceptance checks are additional. |

Use `--dry-run` to inspect the plan and `--output DIR` to choose the report directory. `--workers N` controls the full local and protected baseline replay (default `4`); the strict numerical script uses four workers. Otherwise reports go into a timestamped directory under `replay_outputs/`.

Each mathematical replay writes `REPLAY_RESULT.json`; its fresh strict acceptance is under `strict/STRICT_SEQUENCE_ACCEPTANCE.json` in that output directory. Full mode also writes `FULL_LOCAL_REPLAYS.json`. Full mode was not rerun during packaging; the saved complete record and the earlier targeted strict record are supplied separately. Neither a full nor a strict mathematical replay was rerun for this final editorial revision; numerical constants, checker code, and frozen proof inputs are unchanged.

The [package build validation](inventories/BUILD_VALIDATION.json) passes byte-inventory, source-lookup, tamper-rejection, and complete-stage planning checks. It does not execute a new mathematical replay. [PROVENANCE.md](PROVENANCE.md) describes these checks and the supplied mathematical evidence separately.

Strict and full replays materialize the frozen inputs, without rewriting their contents, at `/workspace/scratch/8242d848994f`, because accepted inputs contain that recorded absolute path. That location must be writable. Fresh strict support and outputs are copied into the selected output directory; the packaged evidence is preserved.

The provided Dockerfile creates the recorded path inside a container, so the host does not need write access to `/workspace`. Build from this directory and mount a report directory:

```bash
docker build -t gr-proof-replay .
mkdir -p replay_outputs
docker run --rm -v "$PWD/replay_outputs:/output" gr-proof-replay --mode strict --output /output/strict
docker run --rm -v "$PWD/replay_outputs:/output" gr-proof-replay --mode full --output /output/full
```

Each run requires a fresh empty output directory. Container execution is outside the scope of the packaged build validation.

## First three matrices

The sequence is `H_n = H*_(10^9·2^(n−1))(711000001/10^9)`. Its first three Gaussian dimensions are `10^9`, `2·10^9`, and `4·10^9`, with grid radius `15`. These dimensions are distinct from the much larger matrix orders. The paper specifies every matrix through finite entry and rounding formulas; the matrices have not been materialized entry by entry.

| Instance | Certified ratio lower bound, rounded downward |
|---|---:|
| `H_1` | `1.775000013258` |
| `H_2` | `1.775000014017` |
| `H_3` | `1.775000014397` |

Exact fractions and radius certificates are in [finite_sequence_arithmetic.json](verification/strict/finite_sequence_arithmetic.json). The complete sequence acceptance is [STRICT_SEQUENCE_ACCEPTANCE.json](verification/strict/STRICT_SEQUENCE_ACCEPTANCE.json). The [saved 219-stage acceptance](verification/frozen/new_round_2026-10-02/audit/FINAL_GOAL_ALL_LOCAL_REPLAY_ACCEPTANCE.json) and [detailed historical result](verification/frozen/new_round_2026-10-02/audit/final_goal_all_local_replay/REPLAY_RESULT.json) are also included. See [PROVENANCE.md](PROVENANCE.md) for their scope and the distinction between the historical upper-bound certificate and Peng's theorem used in the current paper.
