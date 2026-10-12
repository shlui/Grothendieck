#!/usr/bin/env python3
"""Exact rational arithmetic for the self-contained Krivine upper proof."""
from fractions import Fraction as F
from math import factorial
from pathlib import Path
import json,hashlib


def multiply(a,b):
    c=[F(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):c[i+j]+=x*y
    return c


def main():
    x=F(141,160)
    terms=[x**k/F(factorial(k))for k in range(13)]
    ratio=x/14;assert 0<ratio<1
    tail=x**13/F(factorial(13))/(1-ratio)
    U=sum(terms,F(0))+tail
    gap=F(239,99)-U;assert gap>0
    intermediary_gap=F(1207,500)-U;assert intermediary_gap>0
    sqrt_gap=F(2)-F(140,99)**2;assert sqrt_gap==F(2,9801)>0
    numerator=list(map(F,[0,0,0,0,1,-4,6,-4,1]))
    denominator=list(map(F,[1,0,1]))
    quotient=list(map(F,[4,0,-4,0,5,-4,1]));remainder=F(-4)
    recovered=multiply(denominator,quotient);recovered[0]+=remainder
    assert recovered==numerator
    integral_polynomial=sum((c/F(i+1)for i,c in enumerate(quotient)),F(0))
    assert integral_polynomial==F(22,7)
    upper=F(22,7)/(2*x);assert upper==F(1760,987)
    rounding_gap=F(357,200)-upper;assert rounding_gap==F(359,197400)>0
    report=dict(status='PASS_EXACT_RATIONAL_ELEMENTARY_KRIVINE_UPPER_COMPONENTS',
        arithmetic='fractions.Fraction; no floating point or external producer imports',
        exponent=str(x),taylor_terms=list(map(str,terms)),geometric_tail_upper=str(tail),
        exponential_upper=str(U),rational_exp_comparator='239/99',positive_exp_comparison_gap=str(gap),
        optional_exp_comparator='1207/500',optional_positive_gap=str(intermediary_gap),
        sqrt_comparator='140/99',positive_sqrt_square_gap=str(sqrt_gap),
        pi_integrand_numerator=list(map(str,numerator)),pi_integrand_denominator=list(map(str,denominator)),
        polynomial_quotient=list(map(str,quotient)),remainder=str(remainder),
        polynomial_integral=str(integral_polynomial),analytic_integral_identity='Integral = 22/7-pi > 0',
        c_lower=str(x),pi_upper='22/7',universal_constant_strict_upper=str(upper),
        three_significant_digit_lower='71/40',rounding_interval_strict_upper='357/200',
        rounding_upper_positive_gap=str(rounding_gap),
        scope='Arithmetic components only; the signed tensor/Gaussian rounding proof is analytic. Rounding to1.78 additionally requires a separately certified G_R>=71/40.',
        entrypoint_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    output=Path(__file__).with_name('ELEMENTARY_KRIVINE_UPPER_EXACT_CERTIFICATE.json')
    output.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:report[k]for k in ('status','exponential_upper','positive_exp_comparison_gap','universal_constant_strict_upper','rounding_upper_positive_gap')},indent=2))

if __name__=='__main__':main()
