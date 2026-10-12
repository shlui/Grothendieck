#!/usr/bin/env python3
"""Independent MPFR verification of the complete cubic-L1 price cover.

Imports no producer, optimization or scalar verifier implementation.
The eight accepted whole-line inequalities are explicit proof inputs.
"""
from pathlib import Path
from fractions import Fraction as F
from functools import lru_cache
import hashlib,json,os
os.environ['GR_MPFR_BITS']='768'
from mpfr_interval import I,UP,selftest

HERE=Path(__file__).resolve().parent
ENERGY=HERE.parents[2]/'reduced_search/energy'


def run():
    data=json.loads((HERE/'kappa_shape_caps.json').read_text())
    assert data['status']=='PASS_CERTIFIED_KAPPA_SHAPE_CUBIC_L1_CAPS'
    accepted=json.loads((ENERGY/'cubic_L1_certificate.json').read_text())
    assert accepted['status']=='PASS_DIRECTED_COMPLETE_CUBIC_L1_BOUND'
    assert hashlib.sha256((ENERGY/'cubic_L1_certificate.json').read_bytes()).hexdigest()==data['accepted_energy_certificate_sha256']
    rows={}
    for record in accepted['rows']:
        path=ENERGY/record['row_path'];sha=hashlib.sha256(path.read_bytes()).hexdigest()
        assert sha==record['row_sha256']
        frozen=json.loads(path.read_text())
        assert frozen['status']=='PASS_DIRECTED_SHARP_GENERIC_ODD_CHEBYSHEV_RATIONAL_HALF_LINE'
        name=path.name[len('stein_L1_row_'):-len('.json')]
        listed=data['rows'][name]
        assert listed['sha256']==sha
        q=list(map(F,frozen['majorant']));eta=F(frozen['eta'])
        assert q==list(map(F,listed['majorant'])) and eta==F(listed['eta'])
        assert q[3]>0 and eta>0
        rows[name]=(q,eta)
    c=(I(2)/I.pi()).sqrt();K=I(194)/I(10).sqrt()
    assert (K-F(633,15)).above(0)

    @lru_cache(None)
    def price(name,k,shape):
        q,eta=rows[name]
        constant=c+q[0]+q[1]+3*q[2]+15*q[3]
        linear=q[2]+15*q[3]+eta*(F(7,90)+shape/F(30))
        root=(K-F(633,15)*shape)*q[3]
        # Independent expanded form, rather than assembled moment prices.
        return constant+I(linear)*k+root*k*I(k).sqrt()

    leaves=data['leaves'];last=F(0);maximum=F(0);corners=0
    for leaf in leaves:
        lo,hi=map(F,leaf['kappa_interval']);cap=F(leaf['cap']);name=leaf['row']
        assert lo==last and lo<hi<=90;last=hi
        recorded=list(map(F,leaf['corners']))
        index=0
        for k in (lo,hi):
            for shape in (F(0),F(1)):
                assert price(name,k,shape).below(recorded[index])
                assert recorded[index]<=cap
                index+=1;corners+=1
        assert cap==max(recorded)
        maximum=max(maximum,cap)
    assert last==90 and maximum==F(data['global_rho'])
    assert maximum<F(848515,10**6)<F(849,1000)
    table=[]
    for record in data['table']:
        lo,hi=map(F,record['kappa_interval'])
        pieces=[p for p in leaves if F(p['kappa_interval'][0])>=lo and F(p['kappa_interval'][1])<=hi]
        assert pieces and F(pieces[0]['kappa_interval'][0])==lo and F(pieces[-1]['kappa_interval'][1])==hi
        assert len(pieces)==record['pieces']
        assert max(F(p['cap'])for p in pieces)==F(record['rho'])
        table.append(dict(kappa_interval=record['kappa_interval'],rho=record['rho']))
    for record in data['points']:
        k=F(record['kappa']);name=record['row'];cap=F(record['rho'])
        assert all(price(name,k,s).below(cap)for s in (F(0),F(1)))
    out=dict(status='PASS_INDEPENDENT_COMPLETE_KAPPA_SHAPE_L1_PRICES',
        engine=selftest(),leaves=len(leaves),corner_prices_checked=corners,
        global_cap='848515/1000000',stronger_cover_maximum=str(maximum),
        interval_caps=table,
        source_price_sha256=hashlib.sha256((HERE/'kappa_shape_caps.json').read_bytes()).hexdigest(),
        verifier_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        scope='Independent expanded MPFR formulas and exact exhaustive kappa cover. Accepted whole-line row inequalities and fresh633 moment/Stein identities are proof inputs; no new Grothendieck operator norm endpoint.')
    (HERE/'independent_price_verification.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))


if __name__=='__main__':
    run()
