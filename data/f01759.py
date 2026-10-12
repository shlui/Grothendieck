#!/usr/bin/env python3
"""Exact amplitude extension of the already accepted broad scalar profile."""
from fractions import Fraction as F
from pathlib import Path
import json

HERE=Path(__file__).resolve().parent


def main():
    ledger=json.loads((HERE/'regional_u2100_3000_endpoint_ledger.json').read_text())
    assert ledger['status']=='PASS_EXACT_REGIONAL_ASSEMBLY'
    grid=list(map(F,ledger['r_grid']))
    rho=F(849,1000)
    assert grid[-1]==rho
    ell=F(ledger['profile']['ell'])
    constant=sum((F(ledger[k]) for k in ('maximum_endpoint_with_gram',
        'quadratic_chord_reserve','uniform_height_inflation','tail_slope_charge')),F(0))
    low_endpoint=F(7,50)
    energy=low_endpoint**2+F(18,25)*rho*rho
    assert energy==F(53857672,10**8)
    extensions=[]
    for end,cap in ((F(21,50),F(561523,10**6)),
                    (F(10661,25000),F(5633801,10**7))):
        values=[u*u-2*ell*u+constant for u in (low_endpoint,end)]
        assert max(values)<cap<F(40,71)
        assert energy<cap
        extensions.append({'physical_profile_interval':[str(low_endpoint),str(end)],
                           'combined_energy_profile_interval':['0',str(end)],
                           'endpoint_values':list(map(str,values)),
                           'strict_rational_cap':str(cap),
                           'cap_display':float(cap),
                           'maximum_profile_price_display':float(max(values))})
    result={'status':'PASS_EXACT_AMPLITUDE_EXTENSION',
            'source_bound_reason':'||k0||infinity<=1 and the fresh projection dual bound give r<=849/1000 for every amplitude.',
            'fixed_outer_constant':str(constant),'ell':str(ell),
            'low_energy_interval':['0',str(low_endpoint)],'low_energy_cap':str(energy),
            'extensions':extensions,
            'scope':'Continuous amplitude range only; remaining larger amplitudes are unproved.'}
    (HERE/'extended_broad_amplitude_ledger.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':result['status'],'extensions':[
        {'range':x['combined_energy_profile_interval'],'cap':x['strict_rational_cap']}
        for x in extensions]},indent=2))


if __name__=='__main__':
    main()
