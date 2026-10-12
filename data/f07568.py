#!/usr/bin/env python3
"""Directed finite-ratio formula evaluation for an explicitly supplied norm cap.

This program verifies arithmetic, not the all-dimension operator norm cap.
An accepted norm proof must be supplied separately before its displayed
finite ratios become unconditional Grothendieck lower bounds.
"""
from argparse import ArgumentParser
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
from mpfr_interval import I, selftest
from certify_numerator_rate import row, floor12


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--b', required=True, type=F)
    parser.add_argument('--cap', required=True, type=F)
    parser.add_argument('--norm-proof', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    assert 0 <= args.b < 1 and args.cap > 0
    rows = []
    for n in (32, 128, 1024):
        numerator = row(n, F(1), args.b, F(0))
        eta = F(1, n**4)
        lower = F(numerator['finite_vector_numerator_lower_bound']['lower'])
        ratio = I(lower - eta) / (I(args.cap) + I(eta))
        rows.append({'n': n, 'dimension': n*n,
                     'rationalization_tolerance': str(eta),
                     'numerator': numerator,
                     'ratio_lower_enclosure': ratio.record(),
                     'ratio_lower_downfloor12': floor12(ratio)})
    proof = None
    if args.norm_proof:
        raw = args.norm_proof.read_bytes()
        proof = {'file': str(args.norm_proof),
                 'sha256': hashlib.sha256(raw).hexdigest()}
    result = {'status': 'PASS_DIRECTED_FINITE_FORMULA_ARITHMETIC',
              'scope': 'Arithmetic only; requires the separately proved global norm cap.',
              'operator_b': str(args.b), 'required_norm_cap': str(args.cap),
              'limit_of_displayed_lower_bounds': str(1 / args.cap),
              'norm_proof_binding': proof, 'engine': selftest(), 'rows': rows}
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({'status': result['status'], 'operator_b': str(args.b),
                      'required_norm_cap': str(args.cap),
                      'limit_of_lower_bounds': str(1 / args.cap),
                      'rows': [{'n': r['n'], 'lower': r['ratio_lower_downfloor12']}
                               for r in rows]}, indent=2))


if __name__ == '__main__':
    main()
