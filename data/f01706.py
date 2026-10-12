#!/usr/bin/env python3
"""Isolated small-slack complete-halfline replay for the frozen high-u profile."""
from fractions import Fraction
from pathlib import Path
import sys
import certify_continuous_profile_scalars as verifier

HERE=Path(__file__).resolve().parent
verifier.Q0_SLACK=Fraction(1,5000)
verifier.Q6_SLACK=Fraction(1,10**9)
if len(sys.argv)==1:
    sys.argv += ['--source',str(HERE.parents[1]/'tonight/operator/high_u_literature/source_profile.json'),
                 '--output-dir',str(HERE/'high_u_scalars')]
verifier.main()
