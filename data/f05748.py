#!/usr/bin/env python3
"""Independent exact source pricing of tight and finer lower-amplitude cells.

Imports no producer code and executes no local q gate. The q-gate artifacts
are bound by hashes and frozen-profile metadata only; their proof replay is a
separate audit. Outer formulas are expanded over exact rational coefficients
and evaluated with independently implemented directed MPFR intervals.
"""
from fractions import Fraction as F
from pathlib import Path
from itertools import product
import gzip
import hashlib
import json
import os
import sys

os.environ.setdefault('GR_MPFR_BITS','384')
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
BASE = ROOT/'tonight/operator/projection/lower_amp_certificate'
sys.path.insert(0,str(HERE.parent))
from mpfr_interval import I,selftest


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def exact_profile(path):
    raw = json.loads(path.read_text())
    assert raw['kind'] == 'EXACT_RATIONAL_R3_PROFILE'
    return {k:[F(x) for x in v] if isinstance(v,list) else F(v)
            for k,v in raw['profile'].items()}


def upward_decimal(interval,places):
    x = F(interval.hi.text(2))
    scale = 10**places
    return F(-((-x.numerator*scale)//x.denominator),scale)


def coefficients(p,reserve0,reserve6):
    assert p['eps']>0 and p['p2']!=0 and p['eta']>0 and p['G']>0
    assert p['kappa']==F(53,100)
    theta = F(53,759)
    assert theta==p['kappa']/(3*(p['kappa']+2))
    q0,q2,q4,q6 = p['major']
    q0+=reserve0
    q6+=reserve6
    zeta2 = p['zeta']**2
    A = q6+zeta2
    B = p['G']-p['d0']-zeta2
    assert A>=0 and B>=0
    vy = p['y']/p['p2']-p['p0']*p['r']
    constant = (6*p['beta']**2/p['eps']+vy**2
                +q0+q2+3*q4+15*q6+p['eta']*theta**2
                -2*p['nuisance2']*theta+6*zeta2)
    return {'theta':theta,'q0':q0,'q2':q2,'q4':q4,'q6':q6,
            'A':A,'B':B,'vy':vy,'constant':constant,
            'r_linear':-2*(p['gamma']-p['p0']*vy)}


def price(p,c,u,r,k,z):
    # Expanded outer price. All terms except K and the 3/2 powers are rational.
    theta,q4,q6 = c['theta'],c['q4'],c['q6']
    linear_k = (q4+15*q6+p['eta']/2*((1+z)/9-2*theta/3+theta**2)
                +p['nuisance2']*(F(1,3)-theta)
                +9*p['zeta']**2-p['c'])
    rational = (c['constant']+u*u-2*p['ell']*u
                +p['eps']*r*r+c['r_linear']*r
                +linear_k*k+c['B']*k*k)
    nonlinear = (I(c['A'])*(I(194)/I(10).sqrt()-I(F(633,9)*(z-F(2,5))))
                 +I(F(9,4)*c['B'])*I(z)*I(z).sqrt())
    assert nonlinear.above(0)
    return I(rational)+I(k)*I(k).sqrt()*nonlinear


def verify_hash_ledger(source):
    checked=[]
    for filename,digest in source['sha256'].items():
        path=ROOT/filename
        assert sha(path)==digest,('HASH_MISMATCH',filename)
        checked.append({'path':str(path.relative_to(ROOT)),'sha256':digest})
    return checked


def bind_local_metadata(label,p,suffix,reserve0,reserve6):
    # Deliberately no q-gate calculation or producer import.
    expected={k:[str(x) for x in v] if isinstance(v,list) else str(v)
              for k,v in p.items()}
    expected['major'][0]=str(p['major'][0]+reserve0)
    expected['major'][3]=str(p['major'][3]+reserve6)
    result=[]
    for prefix,a,b in [('body','0','8'),('exterior','8','64')]:
        path=BASE/f'{prefix}_u{label}_{suffix}.json.gz'
        with gzip.open(path,'rt') as stream:
            record=json.load(stream)
        assert record['status']=='PASS_PHYSICAL_BERNSTEIN_LOCAL_SUPREMUM'
        assert record['pending']==0
        assert (record['qa'],record['qb'])==(a,b)
        assert F(record['reserve0'])==reserve0
        assert F(record['reserve'])==reserve6
        assert record['profile']==expected
        result.append({'path':str(path.relative_to(ROOT)),'sha256':sha(path),
                       'q_interval':[a,b],'recorded_leaves':record['accepted'],
                       'profile_binding':'exactly matches inflated frozen profile',
                       'q_gate_replayed_here':False})
    return result


def source_bridge(source):
    rho=F(source['rho'])
    assert rho==F(830644,10**6)
    a=F(source['support_gate']['a'])
    upper=F(source['support_gate']['upper'])
    assert a==F(9,5) and upper==F(1649921,10**6)
    l1path=ROOT/'tonight/operator/cubic_L1_refinement/kappa_shape_caps.json'
    l1=json.loads(l1path.read_text())
    row=next(v for v in l1['table'] if v['name']=='coarse0_1')
    assert tuple(map(F,row['kappa_interval']))==(F(0),F(1))
    assert F(row['rho'])==F(830643632499,10**12)<rho
    auditpath=HERE.parent/'refined_l1/independent_refined_l1_price_ledger.json'
    accepted=json.loads(auditpath.read_text())
    assert accepted['status']=='PASS_INDEPENDENT_DIRECTED_REFINED_CUBIC_L1_PRICE_REPLAY'
    assert accepted['source_sha256']==sha(l1path)
    assert next(v for v in accepted['table'] if v['name']=='coarse0_1')['maximum']==row['rho']
    # Quadratic Hermite interpolant at 0 and a double contact t=1/rho.
    # J'''(x)=-(sqrt(2/pi)/a^3)*x*exp(-x^2/(2a^2))<=0 on x>=0.
    c=(I(2)/I.pi()).sqrt()
    t=F(1)/rho
    exponential=(-(I(t)**2)/(2*I(a)**2)).exp()
    derivative=(I(t)/(I(2).sqrt()*I(a))).erf()
    J=I(a)*c*exponential+I(t)*derivative
    alpha=(2*(J-I(a)*c)-I(t)*derivative)/I(t)
    beta=(I(t)*derivative-J+I(a)*c)/(I(t)**2)
    assert alpha.above(0) and beta.above(0)
    bridge=I(a)*c+I(rho)*alpha+beta
    alternative=(1-I(rho)**2)*I(a)*c+I(rho)**2*J
    assert bridge.overlaps(alternative)
    assert bridge.below(upper)
    recorded=F(source['support_gate']['directed_price_upper'])
    assert bridge.below(recorded)
    assert upward_decimal(bridge,49)==recorded
    return {'rho':str(rho),'accepted_l1_cap':row['rho'],
            'accepted_l1_source_sha256':sha(l1path),
            'independent_l1_audit_sha256':sha(auditpath),
            'a':str(a),'source_bridge_upper':str(upper),
            'contact_t':str(t),'J_contact':J.record(),
            'linear_coefficient':alpha.record(),'quadratic_coefficient':beta.record(),
            'bridge_price':bridge.record(),'recorded_upper':str(recorded),
            'majorant_justification':'Hermite remainder J-P=J[0,t,t,x]*x*(x-t)^2<=0 because J third derivative is nonpositive.',
            'mass_justification':'a*u+r=E[a*Z*h+Q*k]<=E max(a|Z|,|Q|), using the four physical point values and independence of Z,Q.',
            'l1_use':'Positive linear interpolant coefficient permits E|Q|<=rho; E Q^2=1 supplies its quadratic expectation.'}


def audit_branch(name,profiles):
    finer=name=='finer'
    source_path=BASE/('lower_source_domain_finer.json' if finer else 'lower_source_cells_tight.json')
    initial_sha=sha(source_path)
    source=json.loads(source_path.read_text())
    assert source['status']==('PASS_ALL_Q_R3_SMALL_CUMULANT_HIGH_SOURCE_DOMAIN'
                              if finer else 'PASS_LOWER_AMPLITUDE_SOURCE_CELLS')
    reserve0,reserve6=map(F,(source['reserve0'],source['reserve6']))
    assert (reserve0,reserve6)==((F(1,50000),F(1,5000000)) if finer
                                else (F(1,20000),F(1,10**6)))
    hashes=verify_hash_ledger(source)
    bridge=source_bridge(source) if finer else None
    results=[]
    for cell in source['cells']:
        label=cell['profile'];p=profiles[label]
        profile_path=BASE/f'profile_u.{label}.json'
        assert source['sha256'][str(profile_path.relative_to(ROOT))]==sha(profile_path)
        c=coefficients(p,reserve0,reserve6)
        uends=(F(341,500),F(69,100))
        assert tuple(map(F,cell['u']))==uends
        zends=(F(2,5),F(1))
        assert tuple(map(F,cell['shape']))==zends
        if finer:
            kends=(F(0),F(1))
            if label=='70':
                r_ends={u:(F(1,5),F(199,500)) for u in uends}
                assert cell['r']=='[1/5,199/500]'
                cap=F(5631,10000)
            else:
                assert label=='68'
                r_ends={u:(F(199,500),F(1649921,10**6)-F(9,5)*u) for u in uends}
                assert cell['r']=='[199/500,1649921/1000000-(9/5)u]'
                cap=F(14083,25000)
            assert F(cell['rounded_exact_cap'])==cap
        else:
            kends=(F(0),F(53,100)) if label=='68' else (F(0),F(1))
            rend=(F(37,100),F(43,100)) if label=='68' else (F(1,5),F(199,500))
            assert tuple(map(F,cell['r']))==rend
            r_ends={u:rend for u in uends}
            cap=F(563,1000) if label=='68' else F(40,71)
            assert F(cell['rounded_cap'])==cap
        assert tuple(map(F,cell['kappa']))==kends
        assert all(lo<hi for lo,hi in r_ends.values())
        expected_vertices=[(u,r,k,z) for u in uends
                           for r,k,z in product(r_ends[u],kends,zends)]
        assert len(cell['vertices'])==len(expected_vertices)==16
        vertices=[]
        for old,(u,r,k,z) in zip(cell['vertices'],expected_vertices):
            assert tuple(F(old[x]) for x in ('u','r','kappa','shape'))==(u,r,k,z)
            actual=price(p,c,u,r,k,z)
            saved=F(old['upper'])
            assert actual.below(saved) and actual.below(cap)
            assert upward_decimal(actual,50)==saved
            # The plus sign is mandatory, although covariance itself may be negative.
            covariance=k/3-c['theta']*(k+2)
            assert covariance==2*(k-F(53,100))/(3*(F(53,100)+2))
            y=p['y']+p['p0']*p['p2']*(r-p['r'])
            assert y-p['p0']*p['p2']*r==p['y']-p['p0']*p['p2']*p['r']
            original_r=(p['eps']-p['p0']**2)*r*r-2*p['gamma']*r+y*y/p['p2']**2
            transported_r=p['eps']*r*r+c['r_linear']*r+c['vy']**2
            assert original_r==transported_r
            vertices.append({'u':str(u),'r':str(r),'kappa':str(k),'shape':str(z),
                              'price_interval':actual.record(),'producer_upper':str(saved),
                              'positive_sign_nuisance_term':str(p['nuisance2']*covariance),
                              'frozen_theta':str(c['theta']),'y_transport':str(y)})
        maximum=max(F(v['producer_upper']) for v in vertices)
        assert maximum==F(cell['maximum_upper'])
        assert F(cell['margin_to_40_71'])==F(40,71)-maximum
        metadata=bind_local_metadata(label,p,name,reserve0,reserve6)
        results.append({'profile':label,'profile_sha256':sha(profile_path),
                        'theta':str(c['theta']),'coefficients':{k:str(v) for k,v in c.items()},
                        'u_interval':list(map(str,uends)),
                        'ur_polygon_vertices':[[str(u),str(r)] for u in uends for r in r_ends[u]],
                        'kappa_interval':list(map(str,kends)),'shape_interval':list(map(str,zends)),
                        'cap':str(cap),'maximum_upper':str(maximum),
                        'margin_to_cap':str(cap-maximum),'margin_to_40_71':str(F(40,71)-maximum),
                        'local_gate_metadata':metadata,'vertices':vertices})
    assert sha(source_path)==initial_sha
    producer=BASE/('price_finer_support_cells.py' if finer else 'price_source_cells.py')
    return {'branch':name,'source_path':str(source_path.relative_to(ROOT)),
            'source_sha256':initial_sha,'producer_path':str(producer.relative_to(ROOT)),
            'producer_sha256':sha(producer),'checked_hash_ledger':hashes,
            'reserve0':str(reserve0),'reserve6':str(reserve6),
            'source_bridge':bridge,'cells':results,
            'scope':('Conditional on separately replayed q gates: all physical sources with u in [.682,.69], kappa in [0,1], and r>=.2, via the audited bridge. Does not cover r<.2 or kappa>1.'
                     if finer else 'Conditional on separately replayed q gates: only the two specified source rectangles; no unrestricted amplitude-strip coverage.')}


def audit_extension(profiles,primary_output):
    results=[]
    for label in ('70','68'):
        p=profiles[label]
        c=coefficients(p,F(1,50000),F(1,5000000))
        vertices=[]
        for u in (F(341,500),F(1381,2000)):
            rs=((F(1,5),F(199,500)) if label=='70' else
                (F(199,500),F(1649921,10**6)-F(9,5)*u))
            assert rs[0]<rs[1]
            for r,k,z in product(rs,(F(0),F(1)),(F(2,5),F(1))):
                actual=price(p,c,u,r,k,z)
                assert actual.below(F(563371,10**6))
                vertices.append({'u':str(u),'r':str(r),'kappa':str(k),
                                  'shape':str(z),'price':actual.record()})
        maximum=max(F(v['price']['upper']) for v in vertices)
        results.append({'profile':label,'profile_sha256':sha(BASE/f'profile_u.{label}.json'),
                        'maximum_upper':str(maximum),'maximum_display':float(maximum),
                        'vertices':vertices})
    assert F(563371,10**6)<F(40,71)
    out={'status':'PASS_INDEPENDENT_OUTER_ONLY_AMPLITUDE_EXTENSION',
         'u':['341/500','1381/2000'],'kappa':['0','1'],'r_lower':'1/5',
         'bridge':'(9/5)u+r<=1649921/1000000',
         'primary_independent_source_audit_sha256':sha(primary_output),
         'script_sha256':sha(Path(__file__)),
         'uniform_cap':'563371/1000000','profiles':results,
         'scope':'Conditional on same separately replayed finer all-q gates and accepted source bridge; source cells only, no kappa>1 or r<.2 coverage.'}
    (HERE/'independent_finer_amplitude_extension.json').write_text(json.dumps(out,indent=2)+'\n')
    return out


def main():
    profiles={label:exact_profile(BASE/f'profile_u.{label}.json') for label in ('68','70')}
    # All shape coefficients multiplying k^(3/2) are nonnegative.
    Kmin=I(194)/I(10).sqrt()-I(F(211,5))
    assert Kmin.above(0)
    result={'status':'PASS_INDEPENDENT_LOWER_AMPLITUDE_R3_SOURCE_PRICING',
            'engine':selftest(),'producer_modules_imported':False,
            'q_gates_executed':False,
            'expanded_price':'C+u^2-2ell*u+eps*r^2-2(gamma-p0*vy)r+L(z)k+B*k^2+[A(K-633(z-2/5)/9)+(9/4)B*z^(3/2)]k^(3/2)',
            'convexity':{
                'minimum_K_shape':Kmin.record(),
                'u_r':'Hessian diag(2,2eps) is positive. Polygon vertex bounds are valid jointly in (u,r).',
                'kappa':'Second derivative =2B+(3/4)[A*K(z)+(9/4)B*z^(3/2)]/sqrt(k)>=0; continuity extends to k=0.',
                'shape':'Second derivative =(27/16)B*k^(3/2)/sqrt(z)>=0; all remaining shape dependence is affine.',
                'scope':'Separate convexity in kappa and shape suffices by repeated endpoint interpolation; joint convexity in these two variables is not required.'},
            'frozen_theta':'53/759 for both profiles, at fitting kappa0=53/100.',
            'covariance_sign':'PLUS nuisance2*[k/3-theta*(k+2)]; covariance may be negative for k<53/100.',
            'nuisance0':'Its expectation is exactly zero because E Utheta=0.',
            'transport_identity':'y(r)-p0*p2*r=y0-p0*p2*r0 and r cost =eps*r^2-2(gamma-p0*vy)r+vy^2.',
            'branches':[audit_branch(name,profiles) for name in ('tight','finer')],
            'scope':'Independent outer-source pricing, domain algebra, source bridge, and file binding. All-q local proof replay is separate; no global operator norm or G_R endpoint is asserted.'}
    output=HERE/'independent_source_price.json'
    output.write_text(json.dumps(result,indent=2)+'\n')
    extension=audit_extension(profiles,output)
    print(json.dumps({'status':result['status'],'engine':result['engine'],
                      'outer_only_extension_cap':extension['uniform_cap'],
                      'branches':[{'branch':b['branch'],
                                   'hash_files_checked':len(b['checked_hash_ledger']),
                                   'cells':[{'profile':c['profile'],'maximum_upper':c['maximum_upper'],'cap':c['cap']} for c in b['cells']]} for b in result['branches']]},indent=2))


if __name__=='__main__':
    main()
