"""Replay exact cubic-cell arithmetic against independent Hermite sums."""
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import hashlib
import json
import random

from compact_cubic_cells import conditional_cubic_kernel, projected_cubic_kernel


def main():
    rng = random.Random(20261002)
    cases = 0
    for dimension in range(1, 13):
        for _ in range(25):
            a = [tuple(F(rng.randrange(-20, 21), rng.randrange(1, 10))
                       for _ in range(3)) for _ in range(dimension)]
            b = [tuple(F(rng.randrange(-20, 21), rng.randrange(1, 10))
                       for _ in range(3)) for _ in range(dimension)]
            u = [x[0] * y[0] for x, y in zip(a, b)]
            e = [(x[1] - 1) * (y[1] - 1) for x, y in zip(a, b)]
            t = [(x[2] - 3 * x[0]) * (y[2] - 3 * y[0])
                 for x, y in zip(a, b)]
            # Three possible multi-index patterns: 3, 2+1, and 1+1+1.
            direct = sum(t) / 6 + sum(e[i] * u[j] / 2
                for i in range(dimension) for j in range(dimension) if i != j)
            direct += sum(u[i] * u[j] * u[k]
                          for i, j, k in combinations(range(dimension), 3))
            assert direct == conditional_cubic_kernel(a, b)
            x = [F(rng.randrange(-20, 21), rng.randrange(1, 10))
                 for _ in range(dimension)]
            y = [F(rng.randrange(-20, 21), rng.randrange(1, 10))
                 for _ in range(dimension)]
            assert projected_cubic_kernel(x, y) == conditional_cubic_kernel(
                [(v, v*v, v**3) for v in x],
                [(v, v*v, v**3) for v in y])
            cases += 1
    source = Path(__file__).with_name('compact_cubic_cells.py')
    result = dict(status='PASS_EXACT_FORMAL_MOMENT_AND_POINT_IDENTITIES',
                  formal_rational_cases=cases, point_cases=cases,
                  dimensions=[1, 12], seed=20261002,
                  source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                  scope='Arithmetic identity only; Gaussian cell moments and '
                  'the global scalar cap are separate certificate inputs.')
    recorded = Path(__file__).with_name('compact_cubic_cells_identity_check.json')
    if recorded.exists():
        assert result == json.loads(recorded.read_text())
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
