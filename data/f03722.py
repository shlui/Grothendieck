#!/usr/bin/env python3
"""Independent MPFR replay of the cubic L1 refinement ledger.

No producer module is imported. Scalar half-line inequalities and the sharp633
moment theorem are accepted premises, not replayed by this price-only checker.
All coefficients and partitions are read as exact rational numbers. The price
is evaluated in a regrouping distinct from the producer, using directed MPFR.
"""
from __future__ import annotations

from fractions import Fraction as F
from pathlib import Path
from collections import Counter
from functools import lru_cache
import hashlib
import json
import os
import sys

os.environ.setdefault('GR_MPFR_BITS', '384')
HERE = Path(__file__).resolve().parent
AUDIT = HERE.parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(AUDIT))
from mpfr_interval import I, selftest

SOURCE = ROOT / 'tonight/operator/cubic_L1_refinement/kappa_shape_caps.json'
ENERGY = ROOT / 'reduced_search/energy'
NAMES = ('0', '2', '3', '4', 'k8', '5', '6', '7')
BIN_SPECS = (
    ('coarse0_1', '0', '1'), ('coarse1_10', '1', '10'),
    ('coarse10_90', '10', '90'), ('fine0_03', '0', '3/10'),
    ('fine03_11', '3/10', '11/10'), ('fine11_4', '11/10', '4'),
    ('fine4_10', '4', '10'), ('high10_20', '10', '20'),
    ('high20_40', '20', '40'), ('high40_60', '40', '60'),
    ('high60_90', '60', '90'))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def decimal_ceiling(bound, places=12):
    # Real.text(UP) is outward rounded; convert its decimal expansion exactly.
    x = F(bound.hi.text(2))
    scale = 10**places
    return F(-((-x.numerator * scale) // x.denominator), scale)


def run():
    engine = selftest()
    recorded = json.loads(SOURCE.read_text())
    assert recorded['status'] == 'PASS_CERTIFIED_KAPPA_SHAPE_CUBIC_L1_CAPS'
    accepted_path = ENERGY / 'cubic_L1_certificate.json'
    accepted = json.loads(accepted_path.read_text())
    assert accepted['status'] == 'PASS_DIRECTED_COMPLETE_CUBIC_L1_BOUND'
    assert sha(accepted_path) == recorded['accepted_energy_certificate_sha256']
    assert len(accepted['rows']) == 8

    rows = {}
    row_audit = {}
    for entry in accepted['rows']:
        path = ENERGY / entry['row_path']
        assert sha(path) == entry['row_sha256']
        raw = json.loads(path.read_text())
        assert raw['status'] == 'PASS_DIRECTED_SHARP_GENERIC_ODD_CHEBYSHEV_RATIONAL_HALF_LINE'
        assert tuple(map(F, raw['odd_weights'])) == (F(1), F(0), F(0))
        assert F(raw['sixth_offset']) == 0
        q = tuple(map(F, raw['majorant']))
        eta = F(raw['eta'])
        assert len(q) == 4 and q[3] > 0 and eta > 0
        name = path.name.removeprefix('stein_L1_row_').removesuffix('.json')
        saved = recorded['rows'][name]
        assert saved['sha256'] == sha(path)
        assert tuple(map(F, saved['majorant'])) == q
        assert F(saved['eta']) == eta
        rows[name] = q, eta
        row_audit[name] = {
            'path': str(path.relative_to(ROOT)), 'sha256': sha(path),
            'majorant': list(map(str, q)), 'eta': str(eta),
            'q6_positive': True, 'eta_positive': True,
            'old_pricing_interval_is_not_a_domain_restriction': True}
    assert set(rows) == set(NAMES) == set(recorded['rows'])

    normal_mean = (I(2) / I.pi()).sqrt()
    sixth_leading = I(194) / I(10).sqrt()
    convex_factor = sixth_leading - I(F(633, 15))
    assert convex_factor.above(0)

    @lru_cache(None)
    def k_threehalves(k):
        return I(k) * I(k).sqrt()

    @lru_cache(None)
    def price(name, k, shape):
        q0, q2, q4, q6 = rows[name][0]
        eta = rows[name][1]
        # Exact rational terms first; only c and K*k^(3/2) need MPFR.
        rational = (q0 + q2 + 3*q4 + 15*q6
                    + k*(q4 + 15*q6 + F(7, 90)*eta)
                    + shape*k*eta/F(30))
        irrational_factor = sixth_leading - I(F(633, 15)*shape)
        return normal_mean + I(rational) + I(q6)*irrational_factor*k_threehalves(k)

    @lru_cache(None)
    def rounded(name, k, shape):
        return decimal_ceiling(price(name, k, shape))

    def four_corners(name, lo, hi):
        return tuple(rounded(name, k, shape)
                     for k in (lo, hi) for shape in (F(0), F(1)))

    leaves = recorded['leaves']
    assert len(leaves) == 2780
    previous = F(0)
    audited_leaves = []
    for idx, leaf in enumerate(leaves):
        lo, hi = map(F, leaf['kappa_interval'])
        assert previous == lo and lo < hi <= 90
        previous = hi
        name = leaf['row']
        assert name in rows
        corners = four_corners(name, lo, hi)
        assert tuple(map(F, leaf['corners'])) == corners
        assert F(leaf['cap']) == max(corners)
        selected_prices = []
        for corner_index, (k, shape) in enumerate(
                (k, s) for k in (lo, hi) for s in (F(0), F(1))):
            interval = price(name, k, shape)
            assert interval.below(corners[corner_index])
            selected_prices.append({'kappa': str(k), 'shape': str(shape),
                                    'directed_interval': interval.record(),
                                    'ceiling_12_places': str(corners[corner_index])})
        candidate_caps = {n: max(four_corners(n, lo, hi)) for n in NAMES}
        best_cap, best_name = min((v, n) for n, v in candidate_caps.items())
        assert (F(leaf['cap']), name) == (best_cap, best_name)
        audited_leaves.append({
            'index': idx, 'kappa_interval': [str(lo), str(hi)], 'row': name,
            'cap': str(best_cap), 'selected_corner_prices': selected_prices,
            'all_eight_candidate_caps': {n: str(v) for n, v in candidate_caps.items()}})
    assert previous == 90

    tables = []
    assert len(recorded['table']) == len(BIN_SPECS)
    for table, (bin_name, a, b) in zip(recorded['table'], BIN_SPECS):
        lo, hi = F(a), F(b)
        assert table['name'] == bin_name
        assert tuple(map(F, table['kappa_interval'])) == (lo, hi)
        pieces = [(i, p) for i, p in enumerate(leaves)
                  if F(p['kappa_interval'][0]) >= lo
                  and F(p['kappa_interval'][1]) <= hi]
        assert pieces
        assert F(pieces[0][1]['kappa_interval'][0]) == lo
        assert F(pieces[-1][1]['kappa_interval'][1]) == hi
        assert all(F(pieces[j][1]['kappa_interval'][1]) ==
                   F(pieces[j+1][1]['kappa_interval'][0])
                   for j in range(len(pieces)-1))
        maximum, active = max((F(p['cap']), i) for i, p in pieces)
        assert (F(table['rho']), table['active_leaf'], table['pieces']) == (
            maximum, active, len(pieces))
        tables.append({'name': bin_name, 'interval': [str(lo), str(hi)],
                       'maximum': str(maximum), 'active_leaf': active,
                       'pieces': len(pieces)})

    points = []
    for item in recorded['points']:
        k = F(item['kappa'])
        cap, name = min((max(four_corners(n, k, k)), n) for n in NAMES)
        assert (F(item['rho']), item['row']) == (cap, name)
        assert tuple(map(F, item['shape_corners'])) == four_corners(name, k, k)
        points.append({'kappa': str(k), 'rho': str(cap), 'row': name,
                       'shape0_interval': price(name, k, F(0)).record(),
                       'shape1_interval': price(name, k, F(1)).record()})
    assert any(F(p['kappa']) == F(53, 100) for p in points)

    maximum, active = max((F(p['cap']), i) for i, p in enumerate(leaves))
    assert maximum == F(recorded['global_rho'])
    assert maximum < F(848521, 10**6) < F(849, 1000)
    out = {
        'status': 'PASS_INDEPENDENT_DIRECTED_REFINED_CUBIC_L1_PRICE_REPLAY',
        'engine': engine, 'source_path': str(SOURCE.relative_to(ROOT)),
        'source_sha256': sha(SOURCE),
        'accepted_energy_certificate_sha256': sha(accepted_path),
        'rows': row_audit,
        'leaf_count': len(leaves), 'selected_corner_count': 4*len(leaves),
        'all_eight_candidate_corner_checks': 8*4*len(leaves),
        'unique_mpfr_price_evaluations': price.cache_info().currsize,
        'coverage': {'interval': ['0', '90'], 'exact_rational': True,
                     'no_gaps': True, 'strictly_positive_widths': True},
        'analytic_validity': {
            'convex_factor_K_minus_633_over_15': convex_factor.record(),
            'formula': 'c+q0+q2+3q4+15q6+k(q4+15q6+7eta/90+s eta/30)+q6(K-633s/15)k^(3/2)',
            'convexity': 'For k>0, d2L/dk2=(3/4)q6(K-633s/15)/sqrt(k)>0. Continuity at k=0 gives convexity on [0,90].',
            'shape': 'L is affine in s; both s=0 and s=1 are checked.',
            'whole_line_row': 'The fixed scalar inequality contains no kappa domain, so old discovery/pricing intervals do not restrict repricing.',
            'premises': 'Accepted eight global scalar inequalities and separately proved sharp633 coupled moments.'},
        'table': tables, 'points': points, 'leaves': audited_leaves,
        'global_rho': str(maximum), 'active_leaf': active,
        'simplified_safe_rho': '848521/1000000',
        'scope': 'Fresh directed moment repricing and rational interval cover only; no scalar half-line replay and no global Grothendieck operator bound.'}
    output = HERE / 'independent_refined_l1_price_ledger.json'
    output.write_text(json.dumps(out, indent=2)+'\n')
    print(json.dumps({k: out[k] for k in (
        'status', 'engine', 'leaf_count', 'selected_corner_count',
        'all_eight_candidate_corner_checks', 'unique_mpfr_price_evaluations',
        'global_rho', 'active_leaf', 'table', 'points')}, indent=2))


if __name__ == '__main__':
    run()
