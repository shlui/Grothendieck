#!/usr/bin/env python3
"""Exact rational finite bounds for the established logarithmic-cutoff sequence.

The scalar cap is an input: this arithmetic does not certify that cap.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path


def floor_decimal(value: F, places: int = 12) -> str:
    multiplier = 10**places
    scaled = value.numerator*multiplier//value.denominator
    return f'{scaled//multiplier}.{scaled%multiplier:0{places}d}'


def evaluate(d: int, b: F, cap: F) -> dict:
    assert isinstance(d,int) and d>=4 and 0<b<1 and cap>0
    q=(1-F(1,2*d))*(1-b/F(2*(d+2)))
    bulk=F(1,d)*((1-b/2)/F(d-1)+b/F(2*(d+2))+b/F(d-2))
    loss=bulk+F(4,d*d)*(1+b/F(2*(d+2)))+4*b/F(d**4)
    eta=F(1,d*d)
    lower=max(F(0),q-loss-eta)/(cap+eta)
    return {'dimension':d,'q_lower':str(q),'loss_upper':str(loss),
            'rationalization_error':str(eta),'ratio_lower':str(lower),
            'ratio_lower_floor12':floor_decimal(lower)}


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--b',required=True,type=F)
    ap.add_argument('--cap',required=True,type=F)
    ap.add_argument('--dimension',nargs='+',type=int,default=[1000,10000,1000000,100000000])
    ap.add_argument('--norm-proof',type=Path)
    ap.add_argument('--output',required=True,type=Path)
    a=ap.parse_args()
    assert not a.output.exists(), 'Preserve the previous arithmetic record.'
    proof=None
    if a.norm_proof:
        data=json.loads(a.norm_proof.read_text())
        assert data['status'].startswith('PASS_INDEPENDENT_COMPLETE_')
        assert F(data['operator_b'])==a.b and F(data['norm_upper_bound'])==a.cap
        proof={'path':str(a.norm_proof),'sha256':hashlib.sha256(a.norm_proof.read_bytes()).hexdigest()}
    result={'status':'PASS_EXACT_RATIONAL_FINITE_ARITHMETIC',
            'operator_b':str(a.b),'scalar_cap_input':str(a.cap),'norm_proof':proof,
            'endpoint_lower':str(1/a.cap),'endpoint_lower_floor12':floor_decimal(1/a.cap),
            'finite_bounds':[evaluate(d,a.b,a.cap) for d in a.dimension],
            'arithmetic':'fractions.Fraction','entrypoint_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'scope':'Finite lower bounds conditional on the stated scalar cap; actual family limit1/N is distinct from1/cap.'}
    a.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
