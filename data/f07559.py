#!/usr/bin/env python3
"""Directed MPFR256 verification of the unconditional finite numerator rate.

No scalar-norm cap or Grothendieck endpoint is assumed or output.
"""
from fractions import Fraction as F
from pathlib import Path
import json
from mpfr_interval import I, selftest

def floor12(x):
    lower=F(x.record()['lower']);scaled=lower*10**12
    k=scaled.numerator//scaled.denominator
    assert k>=0
    return f'{k//10**12}.{k%10**12:012d}'

def row(n,alpha,epsilon,tau):
    d=n*n
    phi=(-I(d)/2).exp()/(2*I.pi()).sqrt()
    p=2*n*phi
    assert F(p.record()['upper'])<1
    direction=(I(F(1,d*(d-2)))+8*n*phi).sqrt()
    loss=I(F(1,d))*(I(alpha-epsilon/2)/(d-1)+I(epsilon)/(2*(d+2))+I(epsilon)/(d-2))
    loss+=4*(I(alpha-epsilon/2)+I(epsilon)/2*I(F(d+4,d+2)).sqrt())*p.sqrt()+4*I(epsilon)*p
    radial=(I(1)-I(F(1,2*d)))*(I(alpha)-I(epsilon)/(2*(d+2)))+I(tau)
    lower=radial-loss
    return {'n':n,'dimension':d,'direction_error_upper_bound':direction.record(),
            'quadratic_objective_loss_upper_bound':loss.record(),
            'radial_value_lower_bound':radial.record(),
            'finite_vector_numerator_lower_bound':lower.record(),
            'finite_vector_numerator_downfloor12':floor12(lower)}

def main():
    profiles={'pure_cubic_720':(F(1),F(18,25),F(0)),
              'summed_tail':(F(126,125),F(699,1000),F(-1,125)),
              'retained_cubic_699':(F(1),F(699,1000),F(0))}
    result={'status':'PASS_DIRECTED_UNCONDITIONAL_NUMERATOR_CERTIFICATE',
      'engine':selftest(),'scope':'Vector numerator only; the global scalar cap remains a separate proof obligation.',
      'profiles':{name:{'alpha':str(a),'epsilon':str(e),'tau':str(t),
                       'rows':[row(n,a,e,t) for n in [8,16,32,64,128,256,512,1024]]}
                  for name,(a,e,t) in profiles.items()}}
    Path(__file__).with_name('numerator_rate_certificate.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':result['status'],'display':
      {name:[{'n':r['n'],'numerator_floor12':r['finite_vector_numerator_downfloor12']}
             for r in profile['rows']]for name,profile in result['profiles'].items()}},indent=2))

if __name__=='__main__':main()
