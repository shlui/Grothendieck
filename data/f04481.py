#!/usr/bin/env python3
"""Independent expanded MPFR repricing of the ten-row .847 cover."""
from pathlib import Path
from fractions import Fraction as F
from functools import lru_cache
import hashlib,json,os,sys
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
os.environ['GR_MPFR_BITS']='768'
from mpfr_interval import I,selftest


def run():
    data=json.loads((HERE/'complete_price_cover.json').read_text())
    assert data['status']=='PASS_COMPLETE_CUBIC_L1_REFIT_PRICE_COVER'
    rows={}
    for name,record in data['rows'].items():
        path=(HERE/record['path']).resolve()
        assert hashlib.sha256(path.read_bytes()).hexdigest()==record['sha256']
        row=json.loads(path.read_text())
        assert row['status']=='PASS_DIRECTED_SHARP_GENERIC_ODD_CHEBYSHEV_RATIONAL_HALF_LINE'
        q=list(map(F,row['majorant']));eta=F(row['eta'])
        assert q==list(map(F,record['majorant'])) and eta==F(record['eta'])
        assert q[3]>0 and eta>0
        rows[name]=(q,eta)
    c=(I(2)/I.pi()).sqrt();K=I(194)/I(10).sqrt()
    assert (K-F(633,15)).above(0)
    @lru_cache(None)
    def price(name,k,s):
        q,eta=rows[name]
        constant=c+q[0]+q[1]+3*q[2]+15*q[3]
        linear=q[2]+15*q[3]+eta*(F(7,90)+s/F(30))
        return constant+I(linear)*k+(K-F(633,15)*s)*q[3]*k*I(k).sqrt()
    last=F(0);maximum=F(0);checks=0
    for leaf in data['leaves']:
        lo,hi=map(F,leaf['kappa_interval']);cap=F(leaf['cap'])
        assert lo==last and lo<hi<=90;last=hi
        listed=list(map(F,leaf['corners']));j=0
        for k in (lo,hi):
            for s in (F(0),F(1)):
                assert price(leaf['row'],k,s).below(listed[j]) and listed[j]<=cap
                j+=1;checks+=1
        assert cap==max(listed)
        maximum=max(maximum,cap)
    assert last==90 and maximum==F(data['maximum_price']) and maximum<F(847,1000)
    for record in data['table']:
        lo,hi=map(F,record['kappa_interval'])
        leaves=[p for p in data['leaves']if F(p['kappa_interval'][0])>=lo and F(p['kappa_interval'][1])<=hi]
        assert leaves and F(leaves[0]['kappa_interval'][0])==lo and F(leaves[-1]['kappa_interval'][1])==hi
        assert len(leaves)==record['pieces'] and max(F(p['cap'])for p in leaves)==F(record['cap'])
    out=dict(status='PASS_INDEPENDENT_REFIT847_PRICE_COVER',engine=selftest(),
        rows=len(rows),leaves=len(data['leaves']),corner_prices=checks,target='847/1000',
        maximum_price=str(maximum),maximum_display=float(maximum),
        source_cover_sha256=hashlib.sha256((HERE/'complete_price_cover.json').read_bytes()).hexdigest(),
        scope='Independent expanded prices and exhaustive rational kappa cover; complete-halfline acceptance of all ten rows and fresh633/Stein premises are separate proof inputs.')
    (HERE/'independent_price_verification.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))


if __name__=='__main__':run()
