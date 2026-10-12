#!/usr/bin/env python3
"""Finite matrix ratios under the explicitly retained .699 norm premise.

The numerator certificate is unconditional. The denominator 1000/1773
is an inherited project premise, not freshly proved by this calculation.
"""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
from mpfr_interval import I, selftest

HERE = Path(__file__).resolve().parent


def main():
    source = HERE / 'numerator_rate_certificate.json'
    raw = source.read_bytes()
    data = json.loads(raw)
    assert data['status'] == 'PASS_DIRECTED_UNCONDITIONAL_NUMERATOR_CERTIFICATE'
    rows = []
    for row in data['profiles']['retained_cubic_699']['rows']:
        n = row['n']
        # A lower endpoint remains a lower bound after subtracting eta.
        numerator = F(row['finite_vector_numerator_lower_bound']['lower'])
        eta = F(1, n**4)
        value = I(numerator - eta) / (I(F(1000, 1773)) + I(eta))
        lower = F(value.record()['lower'])
        floor = lower * 10**12
        integer = floor.numerator // floor.denominator
        rows.append({'n': n, 'dimension': n*n, 'entry_error_budget': str(eta),
                     'conditional_finite_ratio': value.record(),
                     'downfloor12': f'{integer//10**12}.{integer%10**12:012d}'})
    result = {'status': 'PASS_CONDITIONAL_ON_RETAINED_NORM_PREMISE',
              'premise': 'N(699/1000) <= 1000/1773',
              'source_sha256': hashlib.sha256(raw).hexdigest(),
              'engine': selftest(), 'rows': rows}
    (HERE / 'retained_finite_ratio_certificate.json').write_text(
        json.dumps(result, indent=2) + '\n')
    print(json.dumps({'status': result['status'], 'premise': result['premise'],
                      'rows': [{'n': r['n'], 'ratio_lower': r['downfloor12']}
                               for r in rows]}, indent=2))


if __name__ == '__main__':
    main()
