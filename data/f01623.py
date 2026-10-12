#!/usr/bin/env python3
"""Exact counterexample to the scalar degree-six moment relaxation.

This does not claim the discrete law is third Gaussian chaos.
"""
from fractions import Fraction as F
import json
from pathlib import Path


def multiply(a, b):
    c = [F(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            c[i + j] += x * y
    return c


def symmetric_two_atom(u, v):
    q = (1 - u * u) / (v * v - u * u)
    p = 1 - q
    moments = {j: p * u**j + q * v**j for j in (1, 2, 4, 6)}
    return p, q, moments


def verify_dual_factorization(u, v):
    # D[P(x^2)-x] = (x-u)^2(x-v)^2 [x^2+2(u+v)x+(u+v)^2+uv].
    S, P = u + v, u * v
    D = 2 * P * S**3
    lhs = [P**2 * (S**2 + P), -D,
           S**4 - P * S**2 + 3 * P**2, F(0),
           -2 * S**2 + 3 * P, F(0), F(1)]
    rhs = multiply(multiply([u*u, -2*u, F(1)],
                            [v*v, -2*v, F(1)]),
                   [S*S + P, 2*S, F(1)])
    assert lhs == rhs
    assert D > 0 and all(x > 0 for x in (S*S + P, 2*S, F(1)))
    return {"denominator": str(D), "numerator_coefficients": list(map(str, lhs))}


u, v = F(9, 10), F(5)
p, q, moments = symmetric_two_atom(u, v)
assert p == F(2400, 2419) and q == F(19, 2419)
assert moments[1] == F(55, 59)
assert moments[2] == 1
assert moments[4] == F(139, 25)
assert moments[6] == F(154067, 1250)
kappa = moments[4] - 3
kappa6 = moments[6] - 15 * moments[4] + 30
assert kappa == F(64, 25) and kappa6 == F(87317, 1250)
assert 0 <= kappa <= 90

# For every permitted Stein-variance number V, the deficit budget is minimized
# at V=2*kappa/9. The resulting sixth-cumulant allowance is
# (194/sqrt(10)-211/5)*kappa^(3/2).
# sqrt(10)<16/5 is certified by this rational square comparison.
assert F(10) < F(16, 5)**2
K_lower = F(194) / F(16, 5)
strong_coefficient_lower = K_lower - F(211, 5)
assert strong_coefficient_lower == F(737, 40)
sqrt_kappa = F(8, 5)
assert sqrt_kappa**2 == kappa
strong_cumulant_lower = strong_coefficient_lower * kappa * sqrt_kappa
slack = strong_cumulant_lower - kappa6
assert strong_cumulant_lower == F(47168, 625)
assert slack == F(7019, 1250) and slack > 0

target = F(1769, 2000)  # 0.8845 exactly
assert moments[1] > target
assert moments[1]**2 > F(500, 639)

# A second exact example sits near the source's small-cumulant range.
p2, q2, m2 = symmetric_two_atom(F(9, 10), F(4))
k2 = m2[4] - 3
c2 = m2[6] - 15*m2[4] + 30
assert m2[1] == F(46, 49) and k2 == F(17, 20)
assert c2 >= 0 and c2*c2 < F(194**2, 10) * k2**3

# At kappa=0 the exact two-node Gaussian quadrature has relaxed L1 squared
# (3+sqrt(3))/6. sqrt(3)>17/10 already implies >47/60>500/639.
assert F(17, 10)**2 < 3
assert F(47, 60) > F(500, 639)
assert F(47, 60) > target*target

record = {
    "status": "PASS",
    "interpretation": "Scalar moment relaxation barrier; no third-chaos realization claimed",
    "main_example": {
        "absolute_atoms": [str(u), str(v)],
        "absolute_atom_probabilities": [str(p), str(q)],
        "symmetric_atom_probabilities": [str(p/2), str(q/2)],
        "moments": {str(j): str(x) for j, x in moments.items()},
        "kappa4": str(kappa), "kappa6": str(kappa6),
        "deficit_valid_for_every_variance_in": [str(7*kappa/45), str(2*kappa/9)],
        "strict_rational_slack_at_worst_variance": str(slack),
        "absolute_first_moment_decimal": float(moments[1]),
        "squared_absolute_first_moment": str(moments[1]**2),
    },
    "dual_factorization": verify_dual_factorization(u, v),
    "small_cumulant_example": {
        "absolute_atoms": ["9/10", "4"],
        "absolute_atom_probabilities": [str(p2), str(q2)],
        "absolute_first_moment": str(m2[1]),
        "kappa4": str(k2), "kappa6": str(c2),
        "satisfies_sharp_bound_at_minimum_variance": True,
    },
    "kappa_zero_relaxation": {
        "squared_atoms": ["3-sqrt(6)", "3+sqrt(6)"],
        "absolute_first_moment_squared": "(3+sqrt(3))/6",
        "strict_rational_lower_bound_on_squared_value": "47/60",
    },
}
out = Path(__file__).with_suffix('.json')
out.write_text(json.dumps(record, indent=2) + '\n')
print(json.dumps(record, indent=2))
