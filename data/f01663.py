#!/usr/bin/env python3
"""Exact replay of a feasible Gram relaxation point, not a sign norm bound."""
from fractions import Fraction as F
from decimal import Decimal, localcontext, ROUND_FLOOR, ROUND_CEILING
from math import factorial, isqrt
from pathlib import Path
import json

M = F(1000, 1773)
L = 20001

def atan_interval(q, n):
    s = sum((F((-1)**j, (2*j+1)*q**(2*j+1)) for j in range(n)), F(0))
    t = s + F((-1)**n, (2*n+1)*q**(2*n+1))
    return min(s,t), max(s,t)

def cos_sum(x,n):
    return sum((F((-1)**j,factorial(2*j))*x**(2*j) for j in range(n+1)),F(0))

def sin_sum(x,n):
    return sum((F((-1)**j,factorial(2*j+1))*x**(2*j+1) for j in range(n+1)),F(0))

def sqrt_lower(x, places=60):
    scale=10**places
    a=isqrt(x.numerator*scale*scale//x.denominator)
    value=F(a,scale)
    assert value*value<=x<(value+F(1,scale))**2
    return value

def dec(x, upper=False):
    with localcontext() as ctx:
        ctx.prec=55
        ctx.rounding=ROUND_CEILING if upper else ROUND_FLOOR
        return str(Decimal(x.numerator)/Decimal(x.denominator))

def floor_gap(x):
    value=F((x*10**30).__floor__(),10**30)
    assert value>0
    return dec(value)

def main():
    a,A=atan_interval(5,102);b,B=atan_interval(239,34)
    pi_lo,pi_hi=16*a-4*B,16*A-4*b
    c_lo,c_hi=2/pi_hi,2/pi_lo
    nu_lo=sqrt_lower(c_lo)
    assert M<c_lo
    u_hi=cos_sum(pi_lo/7,10)
    v_lo=sin_sum(pi_lo/14,9)
    assert u_hi<F(901,1000) and v_lo>F(222,1000)
    phase_upper=7*M*F(679,1000)+7*(1-M)*F(901,1000)**9
    assert phase_upper<4
    noise_records=[]
    for r in (F(1,10),F(1,2),F(9,10),F(99,100),F(999,1000),F(9999,10000)):
        gap=c_lo-M-(1-M)*r**(L-1)
        assert gap>0
        noise_records.append({'r':str(r),'normalized_positive_gap_floor':floor_gap(gap)})
    mh,mk=(1+M)/2,(1-M)/2
    # P(|Z|<=.3) >= nu*.3*exp(-.045) >= nu*.3*.955.
    assert nu_lo>F(797,1000)
    center_mass_lower=F(797,1000)*F(3,10)*F(191,200)
    assert center_mass_lower>1-mh
    h_first_cap_lower=c_lo*F(91,100)
    assert h_first_cap_lower>M
    # Earlier independent norm caps.
    cap697=F(20000,35457)
    refined_scalar_cap_lower=c_lo*F(2290180883,2500000000)
    assert cap697>M and refined_scalar_cap_lower>M
    # Song--Yue theorem, with positive odd-identity coefficient nu*kappa.
    al=F('1.433927284243');ka=F('.0079949058')
    C=al-F(1,6)+F('.00000000862054')
    sy_margin_lower=c_lo*C-al*M+nu_lo*ka*(1-M)
    assert sy_margin_lower>0
    result={
        'status':'PASS_EXACT_RELAXATION_LIMITATION_NOT_NEW_LOWER_BOUND',
        'baseline_operator':'P1-(699/1000)P3',
        'baseline_uniform_cap':str(M),
        'baseline_endpoint':str(1/M),
        'tail_degree':L,
        'formal_gram':{'e1_f':str(M),'e1_g':str(M),'c1':str(M),
                       'eL_f':str(1-M),'eL_g':str(1-M),'cL':'0',
                       'odd_identity_correlation':str(M)},
        'first_energy_cap_gap_floor':floor_gap(c_lo-M),
        'full_seven_phase_lhs_rational_upper':str(phase_upper),
        'full_seven_phase_gap_floor':floor_gap(4-phase_upper),
        'finite_borell_noise_parameters':noise_records,
        'support_masses':{'h':str(mh),'k':str(mk)},
        'center_mass_threshold_gap_floor':floor_gap(center_mass_lower-(1-mh)),
        'h_first_moment_cap_gap_floor':floor_gap(h_first_cap_lower-M),
        'earlier_cap_697_gap_floor':floor_gap(cap697-M),
        'earlier_refined_cap_gap_floor':floor_gap(refined_scalar_cap_lower-M),
        'song_yue_cap_margin_floor':floor_gap(sy_margin_lower),
        'pi_interval':[dec(pi_lo),dec(pi_hi,True)],
        'limitation':'For a finite head plus constant odd tail, relaxed pairing equals first_eigenvalue*(1000/1773). No strict improvement follows from this relaxation.'}
    path=Path(__file__).with_name('joint_cap_limitation_certificate.json')
    path.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ['status','baseline_endpoint','full_seven_phase_gap_floor','song_yue_cap_margin_floor']},indent=2))

if __name__=='__main__':main()
