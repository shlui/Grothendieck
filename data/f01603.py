#!/usr/bin/env python3
"""Exact lossless matrix-class reduction, with an implicit clone representation.

The input is a rational finite matrix.  Clearing its denominator changes both
objectives by the same factor; cloning then preserves both objectives exactly.
The accompanying report supplies the analytic infinite sequence and its proof.
No output of this program is a new Grothendieck-constant lower bound.
"""
from __future__ import annotations

import argparse
from bisect import bisect_right
from dataclasses import dataclass
from fractions import Fraction as F
from itertools import product
import json
from math import lcm
from pathlib import Path


def rational_matrix(rows):
    a = tuple(tuple(F(x) for x in row) for row in rows)
    if not a or not a[0] or any(len(row) != len(a[0]) for row in a):
        raise ValueError("A nonempty rectangular matrix is required.")
    return a


def clear_denominators(a):
    a = rational_matrix(a)
    scale = lcm(*(x.denominator for row in a for x in row))
    z = tuple(tuple(int(scale*x) for x in row) for row in a)
    keep_r = [i for i, row in enumerate(z) if any(row)]
    keep_c = [j for j in range(len(z[0])) if any(row[j] for row in z)]
    if not keep_r:
        raise ValueError("The zero matrix has no Grothendieck ratio.")
    return tuple(tuple(z[i][j] for j in keep_c) for i in keep_r), scale


def offsets(counts):
    result = [0]
    for n in counts:
        result.append(result[-1]+n)
    return tuple(result)


@dataclass(frozen=True)
class CloneSpec:
    """B has order L, but its storage has the size of the integer seed."""
    z: tuple[tuple[int, ...], ...]
    denominator: int
    r: tuple[int, ...]
    c: tuple[int, ...]
    row_offsets: tuple[int, ...]
    column_offsets: tuple[int, ...]

    @classmethod
    def from_rational(cls, a):
        z, scale = clear_denominators(a)
        r = tuple(sum(abs(x) for x in row) for row in z)
        c = tuple(sum(abs(row[j]) for row in z) for j in range(len(z[0])))
        assert min(r+c) > 0 and sum(r) == sum(c)
        return cls(z, scale, r, c, offsets(r), offsets(c))

    @property
    def order(self):
        return sum(self.r)

    def entry(self, a, b):
        if not (0 <= a < self.order and 0 <= b < self.order):
            raise IndexError("Clone index outside matrix.")
        i = bisect_right(self.row_offsets, a)-1
        j = bisect_right(self.column_offsets, b)-1
        return F(self.z[i][j], self.r[i]*self.c[j])

    def materialize(self, maximum_order=32):
        if self.order > maximum_order:
            raise ValueError("Keep this matrix implicit; clone order exceeds limit.")
        return tuple(tuple(self.entry(a, b) for b in range(self.order))
                     for a in range(self.order))

    def absolute_row_sum(self, original_row):
        i = original_row
        return sum((self.c[j]*F(abs(self.z[i][j]), self.r[i]*self.c[j])
                    for j in range(len(self.c))), F(0))

    def absolute_column_sum(self, original_column):
        j = original_column
        return sum((self.r[i]*F(abs(self.z[i][j]), self.r[i]*self.c[j])
                    for i in range(len(self.r))), F(0))

    def symmetric_balanced_entry(self, a, b):
        """H=J0 tensor [[0,B],[B^T,0]], of order4L, without materializing."""
        size = 2*self.order
        if not (0 <= a < 2*size and 0 <= b < 2*size):
            raise IndexError("Balanced index outside matrix.")
        block_a, a = divmod(a, size)
        block_b, b = divmod(b, size)
        sign = 1 if block_a == block_b else -1
        if a < self.order <= b:
            value = self.entry(a, b-self.order)
        elif b < self.order <= a:
            value = self.entry(b, a-self.order)
        else:
            value = F(0)
        return F(sign, 2)*value

    def to_json(self):
        return {
            "integer_seed": self.z,
            "cleared_denominator": self.denominator,
            "row_clone_counts": self.r,
            "column_clone_counts": self.c,
            "absolute_doubly_stochastic_order": self.order,
            "symmetric_zero_diagonal_zero_sum_order": 4*self.order,
            "formula": "B[(i,p),(j,q)]=Z[i,j]/(r[i]*c[j])",
            "exact_factors": {
                "C(B)/C(input)": self.denominator,
                "Q(B)/Q(input)": self.denominator,
                "C(H)/C(B)": 4,
                "Q(H)/Q(B)": 4,
                "G(B)/G(input)": 1,
                "G(H)/G(input)": 1,
            },
        }


def scalar_norm(a):
    """Exhaustive exact scalar norm; enumeration fixes a global sign."""
    a = rational_matrix(a)
    if len(a) > len(a[0]):
        a = tuple(zip(*a))
    if len(a) > 20:
        raise ValueError("Exhaustive scalar norm limited to the smaller side<=20.")
    maximum = F(-1)
    best = None
    for tail in product((-1, 1), repeat=len(a)-1):
        s = (1,)+tail
        columns = [sum((s[i]*a[i][j] for i in range(len(a))), F(0))
                   for j in range(len(a[0]))]
        value = sum(map(abs, columns), F(0))
        if value > maximum:
            maximum, best = value, s
    return maximum, best


def dot(x, y):
    return sum((F(a)*F(b) for a, b in zip(x, y)), F(0))


def vector_objective(a, u, v):
    a = rational_matrix(a)
    if len(a) != len(u) or len(a[0]) != len(v):
        raise ValueError("Vector-list lengths must match the matrix.")
    return sum((a[i][j]*dot(u[i], v[j]) for i in range(len(a))
                for j in range(len(a[0]))), F(0))


def exact_audit():
    """Check actual small objectives, not only implementation identities."""
    inputs = [
        [[1, 1], [1, -1]],
        [[2, -1], [1, -3]],
        [[F(1, 2), F(-1, 3)], [F(1, 6), F(1, 2)]],
        [[0, 0, 0], [1, -1, 1], [-1, 0, 1]],
    ]
    results = []
    for data in inputs:
        spec = CloneSpec.from_rational(data)
        b = spec.materialize(maximum_order=20)
        c_z, _ = scalar_norm(spec.z)
        c_b, _ = scalar_norm(b)
        assert c_z == c_b
        assert all(spec.absolute_row_sum(i) == 1 for i in range(len(spec.r)))
        assert all(spec.absolute_column_sum(j) == 1 for j in range(len(spec.c)))
        # Rational unit vectors give independent cross-objective checks.
        circle = [(F(1), F(0)), (F(3, 5), F(4, 5)), (F(-5, 13), F(12, 13))]
        u = [circle[i % 3] for i in range(len(spec.r))]
        v = [circle[(j+1) % 3] for j in range(len(spec.c))]
        assert all(dot(w, w) == 1 for w in u+v)
        uu = [u[i] for i in range(len(u)) for _ in range(spec.r[i])]
        vv = [v[j] for j in range(len(v)) for _ in range(spec.c[j])]
        q_z = vector_objective(spec.z, u, v)
        assert vector_objective(b, uu, vv) == q_z
        n = 4*spec.order
        h = tuple(tuple(spec.symmetric_balanced_entry(i, j) for j in range(n))
                  for i in range(n))
        assert all(h[i][j] == h[j][i] for i in range(n) for j in range(n))
        assert all(h[i][i] == 0 and sum(h[i]) == 0
                   and sum(map(abs, h[i])) == 1 for i in range(n))
        # K's independent two terms can attain C(B) simultaneously;
        # grouping gives a finite upper bound4*C(B), as proved in the report.
        results.append({"input": [[str(x) for x in row] for row in data],
                        "seed_scalar_norm": str(c_z),
                        "clone_scalar_norm": str(c_b),
                        "vector_objective_both": str(q_z),
                        "clone_order": spec.order,
                        "balanced_order": n,
                        "absolute_row_and_column_sums": "1 exactly"})
    return {"status": "PASS_EXACT_RATIONAL_REDUCTION_CHECKS",
            "scope": "Exact finite examples and vector assignments; the general proof is analytic.",
            "new_constant_bound": None, "cases": results}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path,
                        help="JSON matrix entries, with rational strings permitted")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--audit", action="store_true")
    args = parser.parse_args()
    if args.audit:
        result = exact_audit()
    elif args.input:
        result = CloneSpec.from_rational(json.loads(args.input.read_text())).to_json()
    else:
        parser.error("Supply --input or --audit.")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps({"status": result.get("status", "EXACT_IMPLICIT_REDUCTION"),
                      "output": str(args.output)}, indent=2))


if __name__ == "__main__":
    main()
