#!/usr/bin/env python3
"""Independent Fraction algebra and MPFR replay of the cubic-head gate.

Imports no producer polynomial or Arb code. The analytic bridge is separately
reviewed in CUBIC_HEAD_POLYNOMIAL_REVIEW.md. Run beside mpfr_interval.py.
"""
from collections import defaultdict
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
from mpfr_interval import I, selftest


def add(*polys):
    result = defaultdict(F)
    for polynomial in polys:
        for exponent, coefficient in polynomial.items():
            result[exponent] += coefficient
    return {e: c for e, c in result.items() if c}


def mul(left, right):
    result = defaultdict(F)
    for le, lc in left.items():
        for re, rc in right.items():
            result[tuple(a+b for a, b in zip(le, re))] += lc*rc
    return {e: c for e, c in result.items() if c}


def scale(polynomial, coefficient):
    return {e: c*coefficient for e, c in polynomial.items() if c*coefficient}


def derivative_w(polynomial):
    return {(a, b, w-1): c*w for (a, b, w), c in polynomial.items() if w}


def main():
    unit = {(0, 0, 0): F(1)}
    a_var, b_var, w_var = ({(1, 0, 0): F(1)},
                          {(0, 1, 0): F(1)},
                          {(0, 0, 1): F(1)})
    d1 = add(a_var, mul(b_var, w_var))
    d3 = add(a_var, scale(mul(b_var, w_var), 3))
    p, r, s = unit, 0, 1
    records = []
    for order, count in enumerate((5, 12, 22), 1):
        p = add(mul(add(scale(p, F(1, 2)), scale(derivative_w(p), -1)),
                    mul(d1, d3)),
                scale(mul(mul(b_var, p), d3), r),
                scale(mul(mul(b_var, p), d1), 3*s))
        r, s = r+2, s+2
        assert len(p) == count and all(c > 0 for c in p.values())
        if order == 1:
            assert {e: c for e, c in p.items() if e[2] == 0} == {
                (2, 0, 0): F(1, 2), (1, 1, 0): F(3)}
        if order == 2:
            assert {e: c for e, c in p.items() if e[2] == 0} == {
                (4, 0, 0): F(1, 4), (3, 1, 0): F(5),
                (2, 2, 0): F(30)}
        records.append({'order': order, 'denominator_powers': [r, s],
                        'coefficients': [{'powers': list(e), 'coefficient': str(c)}
                                         for e, c in sorted(p.items())]})

    a0, b0 = F(97, 40), F(1, 6)
    longitudinal_a = a0+3*b0
    q = [a0+2*b0, 1/(2*a0), -(a0+6*b0)/(24*a0**4),
         (a0*a0+20*a0*b0+120*b0*b0)/(240*a0**7)]
    assert q == [F(331, 120), F(20, 97), F(-1096000, 265587843),
                 F(106274560000, 727184560303017)]
    assert longitudinal_a == F(117, 40) and q[3] > 0
    producer = Path(__file__).resolve().parents[1]/'moments'/'cubic_head_polynomial_certificate.json'
    certificate = json.loads(producer.read_text())
    assert list(map(F, certificate['polynomial_majorant_coefficients_without_c'])) == q
    kappa, u, v = F(53, 100), F(151, 200), F(-1, 4)
    c = (I(2)/I.pi()).sqrt()
    big = I(194)/I(10).sqrt()
    longitudinal_d = I(6).sqrt()*b0
    prices = []
    for shape, cap in ((0, F(252844293, 10**9)),
                       (1, F(250945617, 10**9))):
        variance = kappa*(F(7, 45)+F(shape, 15))
        sixth = 15+15*I(kappa)+(big-I(F(633*shape, 15)))*I(kappa)*I(kappa).sqrt()
        price = c*(I(q[0]+q[1]+q[2]*(3+kappa))+I(q[3])*sixth)
        source_cap = price-I(longitudinal_a*u)-longitudinal_d*v
        assert source_cap.below(cap)
        recorded = certificate['example']['endpoint_prices'][shape]
        assert F(recorded['V']) == variance
        assert F(recorded['source_r_simple_upper']) == cap
        assert price.overlaps(I(F(recorded['A_upper'])))
        assert source_cap.overlaps(I(F(recorded['source_r_upper_at_u755_vminus025'])))
        prices.append({'shape': shape, 'variance': str(variance),
                       'price_interval': price.record(),
                       'source_cap_interval': source_cap.record(),
                       'verified_simple_cap': str(cap)})
    result = {'status': 'PASS_INDEPENDENT_FRACTION_ALGEBRA_AND_MPFR_PRICES',
              'engine': selftest(), 'derivative_numerators': records,
              'majorant_coefficients_without_c': list(map(str, q)),
              'prices': prices,
              'producer_certificate_sha256': hashlib.sha256(producer.read_bytes()).hexdigest(),
              'scope': 'All-halfline polynomial source gate; no global operator endpoint.'}
    Path(__file__).with_name('cubic_head_polynomial_independent.json').write_text(
        json.dumps(result, indent=2)+'\n')
    print(json.dumps({'status': result['status'], 'counts': [5, 12, 22],
                      'simple_caps': [p['verified_simple_cap'] for p in prices]}, indent=2))


if __name__ == '__main__':
    main()
