#!/usr/bin/env python3
"""Isolated exact-slack adapter for the no-optimizer whole-line verifier."""
import argparse
from fractions import Fraction
import sys
import certify_continuous_profile_scalars as verifier
p=argparse.ArgumentParser(add_help=False)
p.add_argument('--q0',type=Fraction,default=Fraction(1,5000))
p.add_argument('--q6',type=Fraction,default=Fraction(1,10**9))
a,remaining=p.parse_known_args()
assert a.q0>=0 and a.q6>=0
verifier.Q0_SLACK=a.q0;verifier.Q6_SLACK=a.q6
sys.argv=[sys.argv[0]]+remaining
verifier.main()
