#!/usr/bin/env python3
"""Exact rational checks for the strict Grothendieck sequence.

No external mathematical certificate is replayed here. The norm caps are
inputs; this file checks their arithmetic and the finite-sequence formulas.
Run from any directory. The JSON is written beside this script.
"""
from fractions import Fraction as F
from pathlib import Path
import json


def frac(x):
    return {"numerator": str(x.numerator), "denominator": str(x.denominator)}


def down_decimal(x, digits=12):
    scaled = x.numerator * 10**digits // x.denominator
    return str(scaled // 10**digits) + "." + str(scaled % 10**digits).zfill(digits)


def log_bounds(t, terms=40):
    # log((1+t)/(1-t)) = 2 sum_{j>=0} t^(2j+1)/(2j+1).
    lo = 2 * sum((t**(2*j+1)/F(2*j+1) for j in range(terms)), F(0))
    tail = 2*t**(2*terms+1)/F(2*terms+1)/(1-t*t)
    return lo, lo+tail


ln2_lo, ln2_hi = log_bounds(F(1, 3))
ln125_lo, ln125_hi = log_bounds(F(1, 9))
ln10_lo = 3*ln2_lo+ln125_lo
ln10_hi = 3*ln2_hi+ln125_hi

b = F(711000001, 10**9)
B0 = F(563380276, 10**9)
B = F(563380277, 10**9)
Braw = F(281690137769807, 500000000000000)
assert 0 < b < 1
assert Braw < B0
assert B == B0 + F(1, 10**9)
assert B < F(40, 71)


def concise_bound(d):
    return (1-(1+b)/F(2*d)-F(10, d*d))/(B+F(1, d*d))


# With a=(1+b)/2, the derivative with respect to real d has numerator
# a*B*d^2+(20*B+2)*d-a. It is positive for d>=4 because 0<a<1
# and (20*B+2)*d-a > 2*d-1 > 0. Thus all subsequent terms exceed
# the already strictly verified first term.
a = (1+b)/2
assert 0 < a < 1


def original_rational_bound(d, cap):
    q = (1-F(1, 2*d))*(1-b/F(2*(d+2)))
    loss = (F(1, d)*((1-b/2)/F(d-1)+b/F(2*(d+2))+b/F(d-2))
            +F(4, d*d)*(1+b/F(2*(d+2)))+4*b/F(d**4))
    assert loss <= F(9, d*d)
    assert q >= 1-(1+b)/F(2*d)
    return (q-loss-F(1, d*d))/(cap+F(1, d*d))


records = []
for n in (1, 2, 3):
    d = 10**9*2**(n-1)
    # ln(d)=9 ln(10)+(n-1)ln(2); prove 14<sqrt(10 ln(d))<15.
    log_lo = 9*ln10_lo+(n-1)*ln2_lo
    log_hi = 9*ln10_hi+(n-1)*ln2_hi
    assert F(196) < 10*log_lo < 10*log_hi < F(225)
    radius = 15
    intervals = 2*d*radius+2
    value = concise_bound(d)
    sharper = original_rational_bound(d, B)
    assert value > F(71, 40)
    assert sharper >= value
    records.append({
        "index": n,
        "gaussian_dimension": str(d),
        "grid_radius": radius,
        "grid_mesh": frac(F(1, d)),
        "coordinate_interval_count": str(intervals),
        "seed_order_formula": f"({intervals})^{d}/2",
        "final_matrix_order_formula": "4 L_d, where L_d=sum_{i,j}|(Z_d)_{ij}|",
        "lower_bound": frac(value),
        "lower_bound_downrounded_12_decimals": down_decimal(value),
        "comparison_40_numerator_minus_71_denominator":
            str(40*value.numerator-71*value.denominator),
        "radius_certificate": {
            "10_log_lower": frac(10*log_lo),
            "10_log_upper": frac(10*log_hi),
            "between_integer_squares": [196, 225],
        },
        "sharper_original_formula_lower_bound": frac(sharper),
    })

# Optional sharper cap and power-of-ten sequence from the initial request.
Bsharp = F(563380276257023, 10**15)
assert Bsharp == Braw + F(1, 10**9)*F(847, 1000)**2
optional = []
for d, radius in ((10**8, 14), (10**9, 15), (10**10, 16)):
    value = original_rational_bound(d, Bsharp)
    assert value > F(71, 40)
    exponent = len(str(d))-1
    assert F((radius-1)**2) < 10*exponent*ln10_lo
    assert 10*exponent*ln10_hi < F(radius**2)
    optional.append({"gaussian_dimension": str(d), "grid_radius": radius,
                     "lower_bound": frac(value),
                     "lower_bound_downrounded_12_decimals": down_decimal(value)})

report = {
    "description": "Exact rational arithmetic for the strict modified matrix sequence",
    "sequence": "H_n = H^*_{d_n}(b), d_n=10^9 2^(n-1), n>=1",
    "b": frac(b),
    "unperturbed_raw_cap": frac(Braw),
    "unperturbed_rounded_cap": frac(B0),
    "perturbed_norm_cap": frac(B),
    "lower_formula": "[1-(1+b)/(2d)-10/d^2]/[B+1/d^2]",
    "finite_lower_bound_limit": frac(1/B),
    "finite_lower_bound_limit_downrounded_12_decimals": down_decimal(1/B),
    "monotonicity": "For real d>=4: F'(d)=[a B d^2+(20 B+2)d-a]/(B d^2+1)^2>0, a=(1+b)/2",
    "norm_cap_gap_from_40_over_71": frac(F(40, 71)-B),
    "conditions": {"b_between_zero_and_one": True, "d_at_least_four": True,
                   "loss_bound_9_over_d_squared": "For d>=4, 0<b<1: Lbar_d <= (101/12)/d^2 < 9/d^2"},
    "first_three": records,
    "optional_sharper_cap": frac(Bsharp),
    "optional_first_three_power_ten": optional,
    "scope": "Arithmetic verification only; the operator norm cap is an independently supplied premise.",
}
destination = Path(__file__).with_suffix('.json')
destination.write_text(json.dumps(report, indent=2)+"\n")
for record in records:
    print(record['gaussian_dimension'], record['lower_bound_downrounded_12_decimals'],
          'R='+str(record['grid_radius']),
          'strict comparison='+record['comparison_40_numerator_minus_71_denominator'])
print(destination)
