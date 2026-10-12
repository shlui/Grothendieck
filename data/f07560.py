#!/usr/bin/env python3
"""Check arithmetic for the cited upper bound; its theorem is external input."""
from fractions import Fraction as F
from pathlib import Path
import json
from mpfr_interval import I, selftest

def main():
    upper = I.pi()/(2*(1+I(2).sqrt()).log()) - I(F(347,10**6))
    assert upper.below(F(1781867,10**6))
    assert upper.below(F(1785,1000))
    result = dict(status='PASS_DIRECTED_PUBLISHED_UPPER_FORMULA_ARITHMETIC',
                  formula='pi/(2*log(1+sqrt(2)))-347/1000000',
                  published_primary_source='https://arxiv.org/html/2608.11158v2',
                  theorem_scope='Cited external upper-bound theorem; its numerical proof is not replayed here.',
                  directed_enclosure=upper.record(),
                  outward_decimal_upper='1781867/1000000',
                  three_significant_digit_rounding_interval=['1775/1000','1785/1000'],
                  lower_bound_still_required='1775/1000',
                  engine=selftest())
    (Path(__file__).parent/'published_upper_rounding_check.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
