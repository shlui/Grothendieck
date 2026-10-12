#!/usr/bin/env python3
"""Independent exact tail-polynomial and directed moment-price reconstruction.

The rational row parameters are read without importing the archived verifier.
Compact interval certificates are replayed separately in stein_halfline_replay.
The analytic sharp633 theorem is separately reconstructed in energy/633;
it is not certified by this script. Run with the standard library and
adjacent mpfr_interval.py.
"""
from fractions import Fraction as F
from math import comb
from pathlib import Path
import hashlib
import json
from mpfr_interval import I, selftest


def trim(a):
    while len(a) > 1 and not a[-1]:
        a.pop()
    return a


def add(a, b):
    return trim([(a[i] if i < len(a) else F(0)) +
                 (b[i] if i < len(b) else F(0))
                 for i in range(max(len(a), len(b)))])


def scale(a, q):
    return trim([x*q for x in a])


def mul(a, b):
    c = [F(0)]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            c[i+j] += x*y
    return trim(c)


def power(a, n):
    c = [F(1)]
    for _ in range(n):
        c = mul(c, a)
    return c


def derivative(a):
    return trim([i*a[i] for i in range(1, len(a))] or [F(0)])


def translated(a, r):
    return trim([sum((a[j]*comb(j, k)*r**(j-k)
                      for j in range(k, len(a))), F(0))
                 for k in range(len(a))])


def tail(row):
    assert list(map(F, row['odd_weights'])) == [F(1), F(0), F(0)]
    assert F(row['sixth_offset']) == 0
    n, length = row['degree'], F(row['scale'])
    v = list(map(F, row['phi_coefficients']))
    eta = F(row['eta'])
    assert eta > 0
    x, d, z = [F(0), F(1)], [length, F(1)], [-length, F(1)]
    cheb = [[F(1)], z]
    for _ in range(2, n+1):
        cheb.append(add(scale(mul(z, cheb[-1]), F(2)),
                        scale(mul(power(d, 2), cheb[-2]), F(-1))))
    primitive = [F(0)] + [v[j]/(j+1) for j in range(4)]
    numerator = mul(primitive, power(d, n+1))
    for j, coefficient in enumerate(v[4:]):
        numerator = add(numerator, scale(mul(mul(x, cheb[j]),
                                            power(d, n-j)), coefficient))
    denominator = power(d, n+2)
    deriv_num = add(mul(derivative(numerator), d),
                    scale(numerator, F(-(n+1))))
    q = [F(0)]*7
    for j, coefficient in enumerate(map(F, row['majorant'])):
        q[2*j] = coefficient
    gap_num = add(add(mul(q, denominator), scale(deriv_num, F(-1))),
                  mul(mul(x, numerator), d))
    base = scale(mul(mul(power(x, 4), denominator), gap_num), 2*eta)
    cu = F(row['tail']['c_upper'])
    assert (I(2)/I.pi()).sqrt().below(cu)
    results = []
    for endpoint in (0, 1):
        residual = add(scale(mul(power(x, 2), deriv_num), F(-1)),
                       scale(denominator, -endpoint*cu))
        polynomial = add(base, scale(power(residual, 2), F(-1)))
        shift = F(row['tail']['domain'][0])
        coefficients = translated(polynomial, shift)
        assert coefficients[0] > 0 and all(c >= 0 for c in coefficients)
        recorded = row['tail']['checks'][endpoint]
        assert F(recorded['minimum_translated_coefficient']) == min(coefficients)
        assert recorded['translated_coefficient_count'] == len(coefficients)
        digest = hashlib.sha256('|'.join(map(str, polynomial)).encode()).hexdigest()
        assert digest == recorded['polynomial_sha256']
        results.append({'endpoint': endpoint, 'coefficient_count': len(coefficients),
                        'minimum_coefficient': str(min(coefficients)),
                        'polynomial_sha256': digest})
    return results


def prices(row, lo, hi, cap):
    q0, q2, q4, q6 = map(F, row['majorant'])
    eta = F(row['eta'])
    assert q6 > 0 and eta > 0
    c = (I(2)/I.pi()).sqrt()
    big = I(194)/I(10).sqrt()
    assert (big-I(F(633, 15))).above(0)
    values = []
    for kappa in (lo, hi):
        for offset, halfvariance in ((F(0), F(7, 90)),
                                    (F(633, 15), F(1, 9))):
            k = I(kappa)
            price = (c+I(q0+q2+q4*(3+kappa))+
                     I(q6)*(15+15*k+(big-I(offset))*k*k.sqrt())+
                     I(eta*kappa*halfvariance))
            assert price.below(cap)
            values.append({'kappa': str(kappa), 'sixth_coefficient_offset': str(offset),
                           'price_interval': price.record()})
    return values


def main():
    source = Path(__file__).resolve().parents[2]/'reduced_search'/'energy'
    certificate = json.loads((source/'cubic_L1_certificate.json').read_text())
    cap = F(certificate['L1_upper_rational'])
    assert cap == F(849, 1000)
    intervals = [(item['row_path'], *map(F, item['kappa_interval']))
                 for item in certificate['rows']]
    assert len(intervals) == 8
    assert intervals[0][1] == 0 and intervals[-1][2] == 90
    assert all(lo < hi for _, lo, hi in intervals)
    assert all(intervals[i][2] == intervals[i+1][1]
               for i in range(len(intervals)-1))
    results = []
    for filename, lo, hi in intervals:
        path = source/filename
        row = json.loads(path.read_text())
        results.append({'row_path': filename, 'row_sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
                        'kappa_interval': [str(lo), str(hi)],
                        'exact_independent_tail': tail(row),
                        'independent_directed_prices': prices(row, lo, hi, cap)})
    result = {'status': 'PASS_INDEPENDENT_TAIL_AND_PRICE_RECONSTRUCTION',
              'engine': selftest(), 'L1_bound': str(cap), 'rows': results,
              'scope': 'Exact tail and directed endpoint prices; compact certificate replay and fresh sharp633 theorem proof are separate.'}
    out = Path(__file__).with_name('stein_prices_and_tails_replay.json')
    out.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'status': result['status'], 'rows': len(results),
                      'L1_bound': result['L1_bound']}, indent=2))


if __name__ == '__main__':
    main()
