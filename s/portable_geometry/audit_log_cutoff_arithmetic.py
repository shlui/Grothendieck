#!/usr/bin/env python3
"""Independent MPFR arithmetic check of the root Arb log-cutoff table.

The scalar norm proof and chi-mean inequality remain analytic inputs.
Imports no root formula helper and uses the accepted portable MPFR adapter.
"""
import argparse, hashlib, json, os
from pathlib import Path
from fractions import Fraction as F
os.environ['GR_MPFR_BITS']='768'
from mpfr_interval import I, selftest


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main(path,output):
    data=json.loads(Path(path).read_text());assert data['status']=='PASS_DIRECTED_FINITE_FORMULA_ARITHMETIC'
    b,B=F(data['operator_b']),F(data['required_global_norm_cap'])
    assert 0<b<1 and B>0
    rows=[]
    for record in data['rows']:
        d,R=record['dimension'],record['integer_cutoff'];assert d>=4 and R>=1
        h,eta=F(1,d),F(1,d*d)
        assert F(record['mesh'])==h and F(record['rationalization_error'])==eta
        assert record['intervals_per_coordinate']==2*R*d+2
        accepted_difference=I(R*R)-10*I(d).log()
        previous_difference=I((R-1)**2)-10*I(d).log()
        assert accepted_difference.above(0) and previous_difference.below(0)
        phi=I(F(-R*R,2)).exp()/(2*I.pi()).sqrt()
        p=(I(2)*d*phi)/R;assert p.below(1)
        # Deliberately group the bulk differently from the root Arb helper.
        bulk=I(h*h)*(I(F(d*(2-b),2*(d-1)))+
                       I(F(d*b,2*(d+2)))+I(F(d*b,d-2)))
        tail=4*p.sqrt()*(I(1-b/2)+I(b/2)*I(F(d+4,d+2)).sqrt())+4*I(b)*p
        loss=bulk+tail
        radial=I(F(2*d-1,2*d))*I(1-F(b,2*(d+2)))
        numerator=radial-loss
        ratio=(numerator-I(eta))/(I(B)+I(eta))
        assert I(F(record['tail_probability_upper'])).above(p)
        assert I(F(record['quadratic_loss_upper'])).above(loss)
        assert I(F(record['radial_lower'])).below(radial)
        assert I(F(record['numerator_lower'])).below(numerator)
        assert I(F(record['finite_ratio_lower'])).below(ratio)
        assert F(record['finite_ratio_lower'])==F(record['finite_ratio_decimal'])
        rows.append({'n':record['n'],'d':d,'minimal_integer_cutoff':R,
                     'cutoff_difference':accepted_difference.record(),
                     'previous_cutoff_difference':previous_difference.record(),
                     'loss_interval':loss.record(),'radial_lower_formula_interval':radial.record(),
                     'finite_ratio_interval':ratio.record(),
                     'root_floor':record['finite_ratio_decimal']})
    result={'status':'PASS_INDEPENDENT_MPFR768_LOG_CUTOFF_ARITHMETIC',
            'source':{'path':str(Path(path).resolve()),'sha256':sha(path)},
            'mpfr_selftest':selftest(),'operator_b':str(b),'norm_cap':str(B),
            'endpoint_interval':(1/I(B)).record(),'rows':rows,
            'scope':'Directed arithmetic and minimal cutoff only; global norm proof and chi-mean lower inequality are separate analytic inputs.',
            'producer_imports':False,'entrypoint_sha256':sha(__file__),
            'adapter_sha256':sha(Path(__file__).with_name('mpfr_interval.py'))}
    Path(output).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':result['status'],'checked_rows':len(rows),
                      'floors':[r['root_floor'] for r in rows]},indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('source');p.add_argument('output')
    a=p.parse_args();main(a.source,a.output)
