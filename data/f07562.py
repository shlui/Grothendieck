#!/usr/bin/env python3
"""Directed MPFR256 finite bounds for the exactly soluble b=1/2 seed.

The scalar theorem is a separately attributed primary-source input. This
checker evaluates only the explicit finite transfer formulas, with exact
rational input and total entrywise rationalization tolerance eta_n=n^-4.
"""
from fractions import Fraction as F
from pathlib import Path
import json
from mpfr_interval import I, selftest

def floor12(x):
    k=(F(x.record()['lower'])*10**12).__floor__()
    return f'{k//10**12}.{k%10**12:012d}'

def row(n):
    d=n*n
    phi=(-I(d)/2).exp()/(2*I.pi()).sqrt()
    p=2*n*phi
    loss=I(F(1,d))*(I(F(3,4))/(d-1)+I(F(1,4))/(d+2)+I(F(1,2))/(d-2))
    loss+=4*(I(F(3,4))+I(F(1,4))*I(F(d+4,d+2)).sqrt())*p.sqrt()+2*p
    radial=(I(1)-I(F(1,2*d)))*(I(1)-I(F(1,4*(d+2))))
    q=radial-loss
    scalar=11/(6*I.pi())
    endpoint=6*I.pi()/11
    eta=F(1,d*d)
    analytic=q/scalar
    rational=(q-I(eta))/(scalar+I(eta))
    rational_upper=(I(1)+I(eta))/(scalar-I(eta))
    assert q.above(eta)
    return {'n':n,'dimension':d,'entrywise_rationalization_tolerance':str(eta),
            'radial_numerator_lower_bound':radial.record(),
            'quadratic_midpoint_objective_loss_bound':loss.record(),
            'finite_vector_numerator_lower_bound':q.record(),
            'exact_analytic_scalar_norm':scalar.record(),
            'analytic_seed_ratio_lower_bound':analytic.record(),
            'analytic_seed_ratio_lower_downfloor12':floor12(analytic),
            'analytic_seed_ratio_upper_bound':endpoint.record(),
            'reduced_rational_ratio_lower_bound':rational.record(),
            'reduced_rational_ratio_lower_downfloor12':floor12(rational),
            'reduced_rational_ratio_upper_bound':rational_upper.record()}

def main():
    result={'status':'PASS_DIRECTED_FINITE_BENCHMARK_FORMULAS',
            'engine':selftest(),
            'source_input':'Saha et al., arXiv:2608.11158v2, Sections 8--9, affine all-odd-sign inequality',
            'scope':'Finite formula evaluation; the external sign theorem is invoked, not independently replayed.',
            'operator':'P1-(1/2)P3',
            'actual_analytic_and_reduced_ratio_limit':(6*I.pi()/11).record(),
            'rows':[row(n) for n in (32,128,1024)]}
    Path(__file__).with_name('saha_benchmark_certificate.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':result['status'],'endpoint':result['actual_analytic_and_reduced_ratio_limit'],
        'rows':[{k:r[k] for k in ('n','analytic_seed_ratio_lower_downfloor12',
                                  'reduced_rational_ratio_lower_downfloor12')}
                for r in result['rows']]},indent=2))

if __name__=='__main__':main()
