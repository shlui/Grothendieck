#!/usr/bin/env python3
"""Directed dyadic seed and implicit reduced matrix for P1-bP3.

This builds a deliberately small illustration of the analytically defined
sequence.  It certifies the entry approximation, exact scalar norm of the
rational seed, and one vector witness.  It does not estimate the global
all-dimension operator norm and does not certify1.775.
"""
from fractions import Fraction as F
from itertools import product
from math import factorial
from pathlib import Path
import argparse
import json
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE/'operator'/'cubic720'))
from mpfr_interval import I, selftest
from construct_reduced_sequence import CloneSpec, scalar_norm


def multiindices(k, d):
    if d == 1:
        yield (k,)
    else:
        for x in range(k+1):
            for rest in multiindices(k-x, d-1):
                yield (x,)+rest


def hermite(x, k):
    if k == 0:
        return I(1)
    a, b = I(1), x
    for j in range(1, k):
        a, b = b, x*b-j*a
    return b


def gaussian_bin_moments(left, right):
    def boundary(x, k):
        if x is None:
            return I(0)
        x = I(x)
        return (-x.square()/2).exp()/((2*I.pi()).sqrt())*hermite(x, k-1)/I(factorial(k)).sqrt()
    def cdf(x, is_left):
        if x is None:
            return I(0 if is_left else 1)
        return (1+(I(x)/I(2).sqrt()).erf())/2
    return [cdf(right, False)-cdf(left, True)]+[
        boundary(left, k)-boundary(right, k) for k in range(1, 4)]


def multiply(items):
    result = I(1)
    for x in items:
        result *= x
    return result


def certified_dyadic(x, bits):
    endpoints = x.record()
    lower, upper = F(endpoints['lower']), F(endpoints['upper'])
    scale = 1 << bits
    midpoint = (lower+upper)/2
    q = F((midpoint*scale+F(1, 2)).numerator // (midpoint*scale+F(1, 2)).denominator, scale)
    error = max(abs(q-lower), abs(q-upper))
    return q, error


def build(d=2, radius=F(1), mesh=F(1), coefficient=F(18, 25), bits=20):
    count = 2*radius/mesh
    if d < 1 or radius <= 0 or mesh <= 0 or count.denominator != 1 or count.numerator % 2:
        raise ValueError('Require d>=1, radius>0, mesh>0, and even integral2radius/mesh.')
    edges = [None]+[-radius+j*mesh for j in range(count.numerator+1)]+[None]
    bins = list(zip(edges[:-1], edges[1:]))
    moments = [gaussian_bin_moments(a, b) for a, b in bins]
    representative = [-radius if a is None else radius if b is None else (a+b)/2
                      for a, b in bins]
    base = len(bins)
    order = base**d//2
    if order > 64:
        raise ValueError('Keep the infinite sequence symbolic; demo limited to64cells.')
    cells = list(product(range(base//2), *[range(base) for _ in range(d-1)]))
    features = {}
    for k in (1, 3):
        indices = list(multiindices(k, d))
        features[k] = [[2*multiply(moments[cell[t]][exponent[t]] for t in range(d))
                        for exponent in indices] for cell in cells]
    def inner(x, y):
        return sum((a*b for a, b in zip(x, y)), I(0))
    interval_a = [[inner(features[1][i], features[1][j])-
                   I(coefficient)*inner(features[3][i], features[3][j])
                   for j in range(order)] for i in range(order)]
    rational_a = [[F(0)]*order for _ in range(order)]
    error = F(0)
    for i in range(order):
        for j in range(i+1):
            q, e = certified_dyadic(interval_a[i][j], bits)
            rational_a[i][j] = rational_a[j][i] = q
            error += e*(1 if i == j else 2)
    c, _ = scalar_norm(rational_a)
    points = [tuple(representative[k] for k in cell) for cell in cells]
    lengths = [I(sum(x*x for x in p)).sqrt() for p in points]
    q = sum((I(rational_a[i][j])*I(sum(x*y for x, y in zip(points[i], points[j])))/
             (lengths[i]*lengths[j]) for i in range(order) for j in range(order)), I(0))
    spec = CloneSpec.from_rational(rational_a)
    assert all(spec.absolute_row_sum(i) == 1 for i in range(len(spec.r)))
    assert all(spec.absolute_column_sum(j) == 1 for j in range(len(spec.c)))
    return {
        'status': 'PASS_DIRECTED_FINITE_SEED_AND_EXACT_CLASS_CONSTRUCTION',
        'scope': 'Small finite illustration; not an all-dimension norm upper certificate.',
        'engine': selftest(), 'dimension': d, 'radius': str(radius), 'mesh': str(mesh),
        'coefficient': str(coefficient), 'dyadic_bits': bits, 'base_order': order,
        'rational_seed': [[str(x) for x in row] for row in rational_a],
        'certified_entrywise_l1_error_upper': str(error),
        'exact_rational_scalar_norm': str(c),
        'certified_vector_assignment_value': q.record(),
        'certified_finite_ratio_lower': (q/I(c)).record(),
        'reduced_class': spec.to_json(),
        'sequence_global_norm_cap': None,
        'new_constant_endpoint': None,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--coefficient', default='18/25')
    parser.add_argument('--bits', type=int, default=20)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = build(coefficient=F(args.coefficient), bits=args.bits)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k: result[k] for k in ['status', 'coefficient', 'base_order',
                      'exact_rational_scalar_norm', 'certified_finite_ratio_lower']}, indent=2))


if __name__ == '__main__':
    main()
