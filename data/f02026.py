#!/usr/bin/env python3
"""Fresh direct-Chebyshev whole-line certification of two refit L1 rows.

Imports the independently written replay definitions, never a producer module.
Uses new adaptive compact partitions on [0,80] at 768 bits, Taylor order six,
and independent Fraction tail polynomials on [80,infinity), with c<4/5.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse
import hashlib
import importlib.util
import json
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
ENGINE = ROOT/'reduced_search/energy/audit_stein_L1_independent.py'
BASE = ROOT/'tonight/operator/cubic_L1_refinement'
ROWS = {
    '68': BASE/'refit_6_25_7_3_certificate.json',
    '96': BASE/'refit_8_85_10_3_certificate.json'}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(name):
    start = time.time()
    spec = importlib.util.spec_from_file_location('independent_chebyshev_engine', ENGINE)
    engine = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(engine)
    assert engine.ctx.prec == 768 and engine.ctx.cap == 8
    engine_hash = sha(ENGINE)
    path = ROWS[name]
    content = path.read_bytes()
    digest = hashlib.sha256(content).hexdigest()
    row = json.loads(content)
    copy_path = BASE/'refit847'/f'row_{name}.json'
    assert sha(copy_path) == digest
    assert row['degree'] == 12 and F(row['scale']) == 2
    assert tuple(map(F, row['odd_weights'])) == (F(1), F(0), F(0))
    assert F(row['sixth_offset']) == 0
    assert len(row['phi_coefficients']) == row['degree']+5
    assert len(row['majorant']) == 4
    assert F(row['eta']) > 0 and F(row['majorant'][3]) > 0
    # Alternating arctangent bounds in Machin's identity are exact Fractions.
    pi_lower = (16*engine.arctan(F(1,5),41)
                -4*engine.arctan(F(1,239),40))
    assert pi_lower > F(25,8)
    compact = engine.compact(row, F(80))
    tail = engine.exact_tail(row, F(80))
    prices = engine.prices(row, *map(F, row['proposed_kappa_interval']))
    assert sha(path) == digest and sha(copy_path) == digest
    assert sha(ENGINE) == engine_hash
    output = {
        'status': 'PASS_FRESH_INDEPENDENT_WHOLE_LINE_CUBIC_L1_REFIT',
        'row': name, 'source_path': str(path.relative_to(ROOT)),
        'source_sha256': digest,
        'frozen_copy_path': str(copy_path.relative_to(ROOT)),
        'frozen_copy_matches': True,
        'independent_engine_path': str(ENGINE.relative_to(ROOT)),
        'independent_engine_sha256': engine_hash,
        'producer_modules_imported': False,
        'precision_bits': engine.ctx.prec,
        'degree': row['degree'], 'scale': row['scale'],
        'majorant': row['majorant'], 'eta': row['eta'],
        'phi_coefficients': row['phi_coefficients'],
        'compact': compact, 'tail': tail,
        'tail_Mills_multiplier': {
            'c': 'sqrt(2/pi)', 'upper_rational': '4/5',
            'exact_Machin_pi_lower': str(pi_lower),
            'pi_lower_exceeds_25_over_8': True},
        'proposed_interval_prices': prices,
        'seconds': time.time()-start,
        'scope': 'Fresh full-line scalar comparison only. Moment pricing accepts the separate Stein and sharp633 premises; this row does not assert a Grothendieck operator bound.'}
    output_path = HERE/f'independent_row_{name}.json'
    output_path.write_text(json.dumps(output, indent=2)+'\n')
    print(json.dumps({k:output[k] for k in (
        'status', 'row', 'source_sha256', 'compact', 'tail',
        'proposed_interval_prices', 'seconds')}, indent=2), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--row', choices=tuple(ROWS), required=True)
    args = parser.parse_args()
    run(args.row)
