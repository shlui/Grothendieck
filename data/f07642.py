"""Directed inverse bounds for the accepted support--source envelope.

For rho in (0,1), c=sqrt(2/pi), physical r<=R_rho(u) implies
u<=U_rho(r)=c*(1-rho^2+rho^2*exp(-erfinv(r/rho)^2)).
No inverse special function is assumed: a directed erf bisection provides
an exact rational bracket.  Inputs and outputs are Fractions; no binary
floating-point discovery input is used.  Uses the accepted LP64 MPFR
adapter, available at the same relative path in the portable proof bundle.
"""
from fractions import Fraction as F
from pathlib import Path
from functools import lru_cache
import json
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tonight/operator/projection/certificate"))
from mpfr_interval import I, UP, DOWN


@lru_cache(None)
def inverse_source(r, rho=F(847, 1000), steps=80):
    if not isinstance(r, (F, int, str)) or not isinstance(rho, (F, int, str)):
        raise TypeError("Use exact Fraction, integer or decimal-string inputs")
    r, rho = F(r), F(rho)
    assert 0 < rho < 1 and 0 <= r <= rho and steps > 0
    c = (I(2)/I.pi()).sqrt()
    if r == 0:
        return {"r": str(r), "rho": str(rho), "s_bracket": None,
                "u_lower": str(F(c.lo.text(DOWN, 80))),
                "u_upper": str(F(c.hi.text(UP, 80)))}
    if r == rho:
        value = c*I(1-rho*rho)
        return {"r": str(r), "rho": str(rho), "s_bracket": "infinity",
                "u_lower": str(F(value.lo.text(DOWN, 80))),
                "u_upper": str(F(value.hi.text(UP, 80)))}
    target = I(r/rho)
    lo, hi = F(0), F(1)
    while not I(hi).erf().above(target):
        hi *= 2
    for _ in range(steps):
        mid = (lo+hi)/2
        image = I(mid).erf()
        if image.below(target):
            lo = mid
        elif image.above(target):
            hi = mid
        else:
            # Stopping preserves the previous certified bracket.  At
            # MPFR256 this cannot affect an 80-step dyadic bisection.
            break
    assert I(lo).erf().below(target) and I(hi).erf().above(target)
    upper = c*(I(1-rho*rho)+I(rho*rho)*I(-lo*lo).exp())
    lower = c*(I(1-rho*rho)+I(rho*rho)*I(-hi*hi).exp())
    return {"r": str(r), "rho": str(rho), "s_bracket": [str(lo), str(hi)],
            "u_lower": str(F(lower.lo.text(DOWN, 80))),
            "u_upper": str(F(upper.hi.text(UP, 80)))}


def physical_u_clip(ulo, uhi, rlo, rho=F(847, 1000)):
    """Outward clipped amplitude range; None means provably empty."""
    ulo, uhi, rlo = map(F, (ulo, uhi, rlo))
    assert ulo <= uhi
    record = inverse_source(rlo, rho)
    high = min(uhi, F(record["u_upper"]))
    return None if high < ulo else (ulo, high)


if __name__ == "__main__":
    records = [inverse_source(r) for r in (F(0), F(129, 500), F(1, 4), F(847, 1000))]
    target = Path(__file__).with_name("source_inverse_bounds_examples.json")
    target.write_text(json.dumps({"status": "PASS_DIRECTED_ERF_BISECTION",
                                  "records": records}, indent=2)+"\n")
    print(json.dumps({"r": records[1]["r"], "u_upper_decimal": float(F(records[1]["u_upper"])),
                      "worst_cell_clip": list(map(str, physical_u_clip(F(151, 200), F(19, 25), F(129, 500))))}))
