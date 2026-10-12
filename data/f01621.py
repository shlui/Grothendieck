#!/usr/bin/env python3
"""Complete-halfline Stein certificate for a normalized cubic-chaos L1 cap.

Uses the standard fixed-chaos Stein identity and the included fresh coupled
sixth-cumulant theorem with coefficient 633 (see 633/PROOF.md and its checks).
This script certifies the scalar part; the tensor part has separate checkers.
Discovery coefficients are frozen as exact rationals; no grid sampling enters
acceptance. Directed Taylor bounds cover a compact interval and exact translated
polynomials cover the unbounded tail. See the imported archived verifier.
"""
from pathlib import Path
from fractions import Fraction as F
import json
import importlib.util
import hashlib
from flint import arb, ctx

HERE = Path(__file__).resolve().parent
VERIFIER = HERE/'halfline_stein_verifier.py'
spec = importlib.util.spec_from_file_location('archived_mean_verifier', VERIFIER)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
ctx.prec = 512


def a(x):
    x = F(x)
    return arb(x.numerator)/x.denominator


def certify_interval(row, lo, hi):
    q = list(map(F, row['majorant']))
    eta = F(row['eta'])
    assert q[3] > 0 and eta > 0
    c = (2/arb.pi()).sqrt()
    K = 194/arb(10).sqrt()
    assert K-a(F(633, 15)) > 0
    values = []
    for k in map(F, (lo, hi)):
        for c6, hv in [(K, F(7, 90)), (K-a(F(633, 15)), F(1, 9))]:
            value = (c+a(q[0]+q[1]+q[2]*(3+k))+
                     a(q[3])*(15+15*a(k)+c6*a(k)*a(k).sqrt())+
                     a(eta*k*hv))
            values.append(F(str(value.upper().fmpq())))
    high = max(values)
    assert high < F(849, 1000), ('PRICE_FAIL', lo, hi, float(high))
    return {'kappa_interval': [str(lo), str(hi)], 'price_upper': str(high),
            'price_display': float(high),
            'proof': 'Affine in coupled shape, so use its two endpoints; '
                     'each resulting price is convex in kappa, so use '
                     'the two interval endpoints.'}


if __name__ == '__main__':
    raw = json.loads((HERE/'stein_mean_discovery.json').read_text())['rows']
    specs = [(0, '0', '39/100'), (2, '39/100', '39/25'),
             (3, '39/25', '4'), (4, '4', '34/5'),
             ('k8', '34/5', '48/5'), (5, '48/5', '76/5'),
             (6, '76/5', '28'), (7, '28', '90')]
    rows = []
    for index, lo, hi in specs:
        output = HERE/f'stein_L1_row_{index}.json'
        new_row = not output.exists()
        if output.exists():
            row = json.loads(output.read_text())
        else:
            discovery = (json.loads((HERE/'stein_k8_discovery.json').read_text())
                         ['rows'][0] if index == 'k8' else raw[index])
            row = module.freeze(discovery, F(1, 10**7))
            # A generous constant reserve avoids relying on cancellation in
            # the floating discovery regression, while its price is only
            # 0.0001. The sixth-degree reserve remains much smaller.
            row['majorant'][0] = str(F(row['majorant'][0])+F(1, 10000))
            row['additional_constant_reserve'] = '1/10000'
        # Replay the complete compact and exterior tests even for frozen rows.
        row = module.certify(row)
        if new_row:
            output.write_text(json.dumps(row, indent=2)+'\n')
        interval = certify_interval(row, lo, hi)
        rows.append({'row_path': output.name,
                     'row_sha256': hashlib.sha256(output.read_bytes()).hexdigest(),
                     'halfline_status': row['status'], **interval})
        print(json.dumps({'row': index, 'compact_leaves': row['compact_leaves'],
                          **interval}), flush=True)
    assert [r['kappa_interval'] for r in rows] == [
        [lo, hi] for _, lo, hi in specs]
    assert F(rows[0]['kappa_interval'][0]) == 0
    assert F(rows[-1]['kappa_interval'][1]) == 90
    assert all(F(left['kappa_interval'][1]) == F(right['kappa_interval'][0])
               for left, right in zip(rows, rows[1:]))
    maximum = max(F(r['price_upper']) for r in rows)
    out = {'status': 'PASS_DIRECTED_COMPLETE_CUBIC_L1_BOUND',
           'scope': 'All normalized third Gaussian chaoses, using the '
                    'standard Stein identity and included fresh sharp633 chain. '
                    'Not a full Grothendieck operator norm certificate.',
           'precision_bits': ctx.prec, 'rows': rows,
           'L1_upper_rational': '849/1000', 'L1_upper_display': .849,
           'stronger_maximum_price_display': float(maximum),
           'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
           'archived_verifier_sha256': hashlib.sha256(VERIFIER.read_bytes()).hexdigest()}
    (HERE/'cubic_L1_certificate.json').write_text(json.dumps(out, indent=2)+'\n')
    print(json.dumps({k: out[k] for k in ('status', 'L1_upper_rational',
                                       'stronger_maximum_price_display')}), flush=True)
