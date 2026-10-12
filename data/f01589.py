#!/usr/bin/env python3
"""Independent exact finite checks of cloning and balancing reductions.

Classical norms and structural identities are exact. Halmos numerical checks
are diagnostics only; the functional-calculus proof is in REDUCED_CLASS_REVIEW.
"""
from fractions import Fraction as F
from math import lcm
from pathlib import Path
import json
import numpy as np


def transpose(a):
    return [list(row) for row in zip(*a)]


def classical(a):
    q = lcm(*(v.denominator for row in a for v in row))
    z = [[int(v*q) for v in row] for row in a]
    n = len(z[0])
    row_values = [sum(row) for row in z]
    best = sum(abs(x) for x in row_values)
    signs = [1]*n
    previous = 0
    for k in range(1, 1 << (n-1)):
        code = k ^ (k >> 1)
        flip = (code ^ previous).bit_length()
        old = signs[flip]
        signs[flip] = -old
        for i in range(len(z)):
            row_values[i] -= 2*old*z[i][flip]
        best = max(best, sum(abs(x) for x in row_values))
        previous = code
    return F(best, q)


def clone(z):
    r = [sum(abs(v) for v in row) for row in z]
    c = [sum(abs(v) for v in col) for col in zip(*z)]
    assert all(r) and all(c)
    rows = [i for i, mass in enumerate(r) for _ in range(mass)]
    cols = [j for j, mass in enumerate(c) for _ in range(mass)]
    assert len(rows) == len(cols)
    return [[F(z[i][j], r[i]*c[j]) for j in cols] for i in rows]


def dilation(a):
    m, n = len(a), len(a[0])
    return [[F(0)]*m+row for row in a] + [row+[F(0)]*n for row in transpose(a)]


def balance(k):
    j = [[F(1, 2), F(-1, 2)], [F(-1, 2), F(1, 2)]]
    return [[j[a][b]*k[i][t] for b in range(2) for t in range(len(k))]
            for a in range(2) for i in range(len(k))]


def abs_ds(a):
    assert all(sum(abs(v) for v in row) == 1 for row in a)
    assert all(sum(abs(v) for v in row) == 1 for row in transpose(a))


def check(z):
    original = [[F(v) for v in row] for row in z]
    b = clone(z)
    k = dilation(b)
    h = balance(k)
    abs_ds(b); abs_ds(k); abs_ds(h)
    assert h == transpose(h)
    assert all(h[i][i] == 0 for i in range(len(h)))
    assert all(sum(row) == 0 for row in h)
    norms = [classical(a) for a in (original, b, k, h)]
    assert norms[1] == norms[0]
    assert norms[2] == 2*norms[0]
    assert norms[3] == 4*norms[0]
    a = np.array([[float(v) for v in row] for row in h])
    vals, vecs = np.linalg.eigh(a)
    assert np.max(np.abs(vals)) < 1+1e-13
    root = (vecs*np.sqrt(np.maximum(0., 1-vals*vals)))@vecs.T
    u = np.block([[a, root], [root, -a]])
    p = (np.eye(len(u))-u)/2
    diag = np.r_[np.full(len(h), 2/len(h)), np.zeros(len(h))]
    w = np.sqrt(diag)[:, None]*u*np.sqrt(diag)[None, :]
    target = np.zeros_like(w); target[:len(h), :len(h)] = 2*a/len(h)
    errors = {
        "U_squared": float(np.max(np.abs(u@u-np.eye(len(u))))),
        "projection_squared": float(np.max(np.abs(p@p-p))),
        "projection_half_diagonal": float(np.max(np.abs(np.diag(p)-.5))),
        "projection_zero_row_sums": float(np.max(np.abs(p.sum(axis=1)))),
        "weighted_identity": float(np.max(np.abs(w-target))),
        "trace_D": float(abs(diag.sum()-2)),
    }
    assert max(errors.values()) < 5e-8
    return {"integer_seed": z, "orders": [len(original), len(b), len(k), len(h)],
            "exact_classical_norms": [str(c) for c in norms],
            "exact_structural_checks": "PASSED",
            "floating_Halmos_diagnostic_errors": errors}


def main():
    result = {"status": "PASSED",
              "cases": [check([[1, 1], [1, -1]]),
                        check([[1, 1], [2, -1]])],
              "scope": "Exact finite classical and structural checks, floating Halmos diagnostics; no new endpoint."}
    out = Path(__file__).with_name("reduced_maps_checks.json")
    out.write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
