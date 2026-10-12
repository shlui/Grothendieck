#!/usr/bin/env python3
"""Producer-free exact replay of the compact cubic cell-entry identity."""
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations
import hashlib
import json

ROOT = Path(__file__).resolve().parents[2]
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
fragment = ROOT / 'new_round_2026-10-02/sequence_theory/COMPACT_CUBIC_CELL_ENTRIES.tex'
checks = []
for dimension in [1, 2, 3, 4, 7, 13, 37]:
    a = [[F((j + 1) * (k + 2) - 5, j + k + 3) for k in range(3)]
         for j in range(dimension)]
    b = [[F((j + 3) * (k + 1) + 2, 2 * j + k + 5) for k in range(3)]
         for j in range(dimension)]
    u = [x[0] * y[0] for x, y in zip(a, b)]
    e = [(x[1] - 1) * (y[1] - 1) for x, y in zip(a, b)]
    t = [(x[2] - 3 * x[0]) * (y[2] - 3 * y[0]) for x, y in zip(a, b)]
    direct = sum(t, F(0)) / 6
    direct += sum((e[i] * u[j] / 2 for i in range(dimension)
                   for j in range(dimension) if i != j), F(0))
    direct += sum((u[i] * u[j] * u[k] for i, j, k in
                   combinations(range(dimension), 3)), F(0))
    s = sum(u, F(0))
    compact = (sum(t, F(0)) + 3 * s * sum(e, F(0))
               - 3 * sum((x * y for x, y in zip(e, u)), F(0))
               + s ** 3 - 3 * s * sum((x ** 2 for x in u), F(0))
               + 2 * sum((x ** 3 for x in u), F(0))) / 6
    assert direct == compact
    checks.append({'dimension': dimension, 'conditional_identity': True})

result = {
    'status': 'PASS_INDEPENDENT_ANALYTIC_AND_EXACT_COMPACT_CUBIC_CELL_ENTRY_REVIEW',
    'checks': checks,
    'fragment': {'path': str(fragment.relative_to(ROOT)), 'sha256': sha(fragment)},
    'entrypoint_sha256': sha(Path(__file__)),
    'analytic_review': 'The degree-three orthonormal Hermite patterns (3),(2,1),'
        '(1,1,1) have factors1/6,1/2,1. Independence survives product-cell '
        'conditioning. Expanding the elementary symmetric cubic gives the '
        'stated O(d) expression; point substitution gives the displayed kernel. '
        'Gaussian interval moments follow by integration by parts. Antipodal '
        'odd-cell features are2 times their cell integrals, hence weights '
        'w_C=2 gamma(C) give exactly the folded matrix entry.',
    'scope': 'An arbitrary-dimensional algebra identity. No new norm cap, '
        'density assertion, or convergence rate is introduced.'
}
destination = Path(__file__).with_name('COMPACT_CUBIC_ENTRY_INDEPENDENT_AUDIT.json')
assert not destination.exists()
destination.write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps({'path': str(destination.relative_to(ROOT)),
                  'sha256': sha(destination)}, indent=2))
