#!/usr/bin/env python3
"""Independent diffuse R3 source-price algebra and directed evaluations.

No producer imports or q-gate execution. The profile is kept frozen at its
fitting cumulant; b may be lowered from its validated reference value without
altering the local potential or scalar comparison.
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
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
BASE=ROOT/'tonight/operator/projection/diffuse_certificate'
PROFILE=BASE/'profile_interval_u7478_k0_01.json'
sys.path.insert(0,str(HERE.parent))
from mpfr_interval import I,selftest


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def upward_decimal(value,places):
    q=F(value.hi.text(2));scale=10**places
    return F(-((-q.numerator*scale)//q.denominator),scale)


def profile():
    raw=json.loads(PROFILE.read_text())
    assert raw['kind']=='EXACT_RATIONAL_R3_PROFILE'
    p={k:[F(t)for t in v]if isinstance(v,list)else F(v)
       for k,v in raw['profile'].items()}
    assert p['kappa']==F(1,20) and p['eps']==F(18,25)
    assert p['eta']>0 and p['G']>0 and p['p2']!=0
    return p


def coefficients(p,b,reserve0=F(1,50000),reserve6=F(1,5000000)):
    assert 0<b<=p['eps']
    theta=p['kappa']/(3*(p['kappa']+2))
    assert theta==F(1,123)
    q0,q2,q4,q6=p['major'];q0+=reserve0;q6+=reserve6
    A=q6+p['zeta']**2
    B=p['G']-p['d0']-p['zeta']**2
    assert A>0 and B>0
    vy=p['y']/p['p2']-p['p0']*p['r']
    C=(6*p['beta']**2/b+vy**2+q0+q2+3*q4+15*q6
       +p['eta']*theta**2-2*p['nuisance2']*theta+6*p['zeta']**2)
    return {'theta':theta,'q0':q0,'q2':q2,'q4':q4,'q6':q6,
            'A':A,'B':B,'vy':vy,'constant':C,
            'r_linear':-2*(p['gamma']-p['p0']*vy),'b':b}


def price(p,c,u,r,k,z):
    assert k>=0 and F(2,5)<=z<=1
    theta=c['theta']
    L=(c['q4']+15*c['q6']+p['eta']/2*((1+z)/9-2*theta/3+theta**2)
       +p['nuisance2']*(F(1,3)-theta)+9*p['zeta']**2-p['c'])
    rational=(c['constant']+u*u-2*p['ell']*u+c['b']*r*r
              +c['r_linear']*r+L*k+c['B']*k*k)
    nonlinear=(I(c['A'])*(I(194)/I(10).sqrt()-I(F(633,9)*(z-F(2,5))))
               +I(F(9,4)*c['B'])*I(z)*I(z).sqrt())
    assert nonlinear.above(0)
    return I(rational)+I(k)*I(k).sqrt()*nonlinear


def spectral_gates(p,b):
    # Finite rational prerequisites only, not local q proofs.
    even=p['p2']**2-b*(1+p['p0']**2)/2
    assert even>=0
    matrices=[]
    for a,d,e in [(p['R11'],p['R13'],p['R33']),
                  (p['R11']+1,p['R13'],p['R33']-b/6),
                  (p['R11']-b,p['R13'],p['R33'])]:
        determinant=a*e-d*d
        assert a>=0 and e>=0 and determinant>=0
        matrices.append({'diagonal':[str(a),str(e)],'off_diagonal':str(d),
                         'determinant':str(determinant)})
    return {'b':str(b),'even_gate_margin':str(even),'odd_gate_matrices':matrices}


def metadata(p):
    expected={k:[str(t)for t in v]if isinstance(v,list)else str(v)
              for k,v in p.items()}
    expected['major'][0]=str(p['major'][0]+F(1,50000))
    expected['major'][3]=str(p['major'][3]+F(1,5000000))
    rows=[]
    for prefix,a,b in [('body','0','8'),('exterior','8','64')]:
        path=BASE/f'{prefix}_interval_u7478_k0_01.json.gz'
        with gzip.open(path,'rt')as stream:record=json.load(stream)
        assert record['status']=='PASS_PHYSICAL_BERNSTEIN_LOCAL_SUPREMUM'
        assert record['pending']==0
        assert (record['qa'],record['qb'])==(a,b)
        assert F(record['reserve0'])==F(1,50000)
        assert F(record['reserve'])==F(1,5000000)
        assert record['profile']==expected
        rows.append({'path':str(path.relative_to(ROOT)),'sha256':sha(path),
                     'q_interval':[a,b],'recorded_leaves':record['accepted'],
                     'q_gate_replayed_here':False})
    return rows


def anchor_prices(p):
    rows=[]
    for b in (F(18,25),F(7171,10000)):
        c=coefficients(p,b)
        for label,u,r in [('fitted',F(3739,5000),F(1391,5000)),
                          ('witness',F(747698091,10**9),F(278354057,10**9))]:
            vertices=[]
            for k,z in product((F(0),F(1,10)),(F(2,5),F(1))):
                value=price(p,c,u,r,k,z)
                assert value.below(F(40,71))
                vertices.append({'kappa':str(k),'shape':str(z),'price':value.record(),
                                  'plus_sign_covariance_correction':str(p['nuisance2']*(k/3-c['theta']*(k+2)))})
            maximum=max(F(v['price']['upper'])for v in vertices)
            rows.append({'label':label,'b':str(b),'u':str(u),'r':str(r),
                         'kappa_interval':['0','1/10'],'shape_interval':['2/5','1'],
                         'maximum_upper':str(maximum),'maximum_display':float(maximum),
                         'vertices':vertices})
    return rows


def support_bridge(source):
    l1path=ROOT/'tonight/operator/cubic_L1_refinement/kappa_shape_caps.json'
    l1=json.loads(l1path.read_text())
    row=next(x for x in l1['table']if x['name']=='coarse0_1')
    assert tuple(map(F,row['kappa_interval']))==(F(0),F(1))
    rho=F(830644,10**6)
    assert F(row['rho'])<rho
    auditpath=HERE.parent/'refined_l1/independent_refined_l1_price_ledger.json'
    accepted=json.loads(auditpath.read_text())
    assert accepted['status']=='PASS_INDEPENDENT_DIRECTED_REFINED_CUBIC_L1_PRICE_REPLAY'
    assert accepted['source_sha256']==sha(l1path)
    assert next(x for x in accepted['table']if x['name']=='coarse0_1')['maximum']==row['rho']
    a=F(3);t=F(1)/rho;c=(I(2)/I.pi()).sqrt()
    derivative=(I(t)/(I(2).sqrt()*I(a))).erf()
    J=I(a)*c*(-I(t)**2/(2*I(a)**2)).exp()+I(t)*derivative
    alpha=(2*(J-I(a)*c)-I(t)*derivative)/I(t)
    beta=(I(t)*derivative-J+I(a)*c)/(I(t)**2)
    assert alpha.above(0) and beta.above(0)
    bridge=I(a)*c+I(rho)*alpha+beta
    saved=F(source['support_upper'])
    assert bridge.below(saved) and upward_decimal(bridge,49)==saved
    residual=bridge-3*I(F(743,1000))
    assert residual.below(F(3,10))
    return {'a':'3','rho':str(rho),'accepted_kappa_interval':['0','1'],
            'accepted_l1_cap':row['rho'],'accepted_l1_source_sha256':sha(l1path),
            'independent_l1_audit_sha256':sha(auditpath),
            'contact':str(t),'linear_coefficient':alpha.record(),
            'quadratic_coefficient':beta.record(),'bridge':bridge.record(),
            'recorded_upper':str(saved),'r_upper_at_u743':residual.record(),
            'proof':'J_a third derivative is nonpositive, so its quadratic Hermite interpolant at 0,t,t is a global upper bound; positive linear coefficient permits E|Q|<=rho. Physical diamond and independence give 3u+r<=E J_a(|Q|).'}


def audit_masks(p):
    path=BASE/'critical_source_upper_masks.json'
    initial_sha=sha(path);source=json.loads(path.read_text())
    assert source['status']=='PASS_ALL_Q_CONTINUOUS_CRITICAL_HIGH_SOURCE_MASKS'
    assert F(source['reserve0'])==F(1,50000)
    assert F(source['reserve6'])==F(1,5000000)
    hashes=[]
    for filename,digest in source['sha256'].items():
        file=ROOT/filename;assert sha(file)==digest,('HASH_MISMATCH',filename)
        hashes.append({'path':str(file.relative_to(ROOT)),'sha256':digest})
    local=metadata(p)
    assert len(source['regions'])==len(local)==2
    for record,ours in zip(source['regions'],local):
        assert record['q']==ours['q_interval']
        assert record['leaves']==ours['recorded_leaves']
    tail_path=BASE/'tail_interval_u7478_k0_01.json'
    assert json.loads(tail_path.read_text())==source['tail']
    assert source['tail']['status']=='PASS' and F(source['tail']['R'])==64
    bridge=support_bridge(source)
    masks=[]
    assert len(source['masks'])==2
    for old in source['masks']:
        b,cap=F(old['operator_b']),F(old['cap'])
        assert (b,cap) in ((F(18,25),F(56337,10**5)),(F(7171,10**4),F(3521,6250)))
        assert tuple(map(F,old['kappa']))==(F(0),F(1,10))
        assert tuple(map(F,old['shape']))==(F(2,5),F(1))
        spectral=spectral_gates(p,b)
        c=coefficients(p,b)
        source_corners=[]
        for k,z in product((F(0),F(1,10)),(F(2,5),F(1))):
            value=price(p,c,F(0),F(0),k,z)
            source_corners.append({'kappa':str(k),'shape':str(z),'price':value.record(),
                                   'upward_60_significant_digits':str(upward_decimal(value,59))})
        constant=max(F(x['upward_60_significant_digits'])for x in source_corners)
        assert constant==F(old['source_constant_upper'])
        gamma=p['gamma']-p['p0']*c['vy']
        assert gamma==F(old['effective_gamma']) and gamma>b*F(3,10)
        def Fbound(u,r):
            return u*u-2*p['ell']*u+b*r*r-2*gamma*r+constant
        cells=[];previous=F(743,1000);previous_r=None
        assert len(old['cells'])==12
        for index,cell in enumerate(old['cells']):
            ua,ub=map(F,cell['u']);ra,rb=map(F,cell['lower_r_vertices'])
            assert ua==previous and ub-ua==F(1,1000)
            if previous_r is not None:assert ra==previous_r
            previous,previous_r=ub,rb
            assert F(cell['upper_r'])==F(3,10) and F(cell['cap'])==cap
            assert 0<ra<F(3,10) and 0<rb<F(3,10)
            # Verify exact integer threshold proposals independently. They are
            # safe boundary endpoints, not numerical proof by root sampling.
            threshold=[]
            for u,r in ((ua,ra),(ub,rb)):
                assert (r*10**12).denominator==1
                assert Fbound(u,r)<cap
                assert Fbound(u,r-F(1,10**12))<cap
                assert Fbound(u,r-F(2,10**12))>=cap
                threshold.append({'u':str(u),'r':str(r),'bound':str(Fbound(u,r)),
                                  'preceding_grid_point_passes':True,
                                  'second_preceding_grid_point_fails':True})
            ur=((ua,ra),(ub,rb),(ua,F(3,10)),(ub,F(3,10)))
            maximum=max(Fbound(u,r)for u,r in ur)
            assert maximum==F(cell['maximum_upper']) and maximum<cap
            vertex_prices=[]
            for u,r in ur:
                rational_upper=Fbound(u,r)
                for k,z in product((F(0),F(1,10)),(F(2,5),F(1))):
                    actual=price(p,c,u,r,k,z)
                    assert actual.below(rational_upper) and actual.below(cap)
                    covariance=k/3-c['theta']*(k+2)
                    assert covariance==2*(k-F(1,20))/(3*(F(1,20)+2))
                    y=p['y']+p['p0']*p['p2']*(r-p['r'])
                    assert y-p['p0']*p['p2']*r==p['y']-p['p0']*p['p2']*p['r']
                    oldr=(b-p['p0']**2)*r*r-2*p['gamma']*r+y*y/p['p2']**2
                    assert oldr==b*r*r-2*gamma*r+c['vy']**2
                    b0=p['eps']
                    change=(b-b0)*(r*r-6*p['beta']**2/(b*b0))
                    assert change==(b-b0)*r*r+6*p['beta']**2*(1/b-1/b0)
                    vertex_prices.append({'u':str(u),'r':str(r),'kappa':str(k),'shape':str(z),
                                          'price':actual.record(),'rational_bound':str(rational_upper),
                                          'plus_sign_covariance_term':str(p['nuisance2']*covariance),
                                          'exact_b_price_change':str(change)})
            cells.append({'index':index,'u':[str(ua),str(ub)],
                          'lower_r_vertices':[str(ra),str(rb)],'upper_r':'3/10',
                          'threshold_checks':threshold,'maximum_rational_bound':str(maximum),
                          'margin_to_cap':str(cap-maximum),'vertices':vertex_prices})
        assert previous==F(151,200)
        example=old['example_rational_point']
        eu,er=F(example['u']),F(example['r'])
        actual=max(upward_decimal(price(p,c,eu,er,k,z),50)
                   for k,z in product((F(0),F(1,10)),(F(2,5),F(1))))
        assert actual==F(example['upper']) and actual<cap
        containing=next(x for x in cells if F(x['u'][0])<=eu<=F(x['u'][1]))
        ua,ub=map(F,containing['u']);ra,rb=map(F,containing['lower_r_vertices'])
        lower=ra+(rb-ra)*(eu-ua)/(ub-ua)
        assert lower<=er<=F(3,10)
        # A later proposed source witness is deliberately outside the masks.
        # This checks coverage only and asserts no physical counterexample.
        gu,gr=F(750515,10**6),F(270685,10**6)
        gapcell=next(x for x in cells if F(x['u'][0])<=gu<=F(x['u'][1]))
        ga,gb=map(F,gapcell['u']);gra,grb=map(F,gapcell['lower_r_vertices'])
        gaplower=gra+(grb-gra)*(gu-ga)/(gb-ga)
        assert gr<gaplower
        gapprices=[price(p,c,gu,gr,k,z)for k,z in product((F(0),F(1,10)),(F(2,5),F(1)))]
        gapmaximum=max(gapprices,key=lambda x:F(x.hi.text(2)))
        assert gapmaximum.above(F(40,71))
        masks.append({'operator_b':str(b),'cap':str(cap),'spectral_gates':spectral,
                      'source_constant':str(constant),'source_constant_corner_checks':source_corners,
                      'effective_gamma':str(gamma),'cells':cells,
                      'example':{'u':str(eu),'r':str(er),'mask_lower_r':str(lower),
                                 'inside_mask':True,'maximum_upper':str(actual)},
                      'later_proposed_source_point':{'u':str(gu),'r':str(gr),
                          'mask_lower_r':str(gaplower),'below_mask_by':str(gaplower-gr),
                          'inside_mask':False,'uniform_profile_price':gapmaximum.record(),
                          'scope':'Outside the certified source domain. The profile price exceeds the target, which does not prove an actual operator-norm counterexample or physical realizability of this point.'}})
    assert sha(path)==initial_sha
    out={'status':'PASS_INDEPENDENT_DIFFUSE_R3_CONTINUOUS_SOURCE_MASK_PRICING',
         'engine':selftest(),'source_path':str(path.relative_to(ROOT)),
         'source_sha256':initial_sha,'profile_sha256':sha(PROFILE),
         'checked_hash_ledger':hashes,'local_q_metadata':local,
         'tail_record_sha256':sha(tail_path),'tail_proof_replayed_here':False,
         'frozen_theta':'1/123','reserve0':'1/50000','reserve6':'1/5000000',
         'support_bridge':bridge,'mask_count':2,'cell_count':24,
         'polygon_kappa_shape_vertex_checks':384,'source_constant_corner_checks':8,
         'convexity':{
             'ur':'Hessian diag(2,2b)>0 gives polygon vertex bounds.',
             'kappa':'2B+(3/4)[A*K(z)+(9/4)B*z^(3/2)]/sqrt(k)>=0, continuously at k=0.',
             'shape':'(27/16)B*k^(3/2)/sqrt(z)>=0; remaining shape terms affine.',
             'r_monotonicity':'effective_gamma>b*.3 ensures the rational bound decreases in r throughout the masks.'},
         'masks':masks,
         'scope':'Only these high-r polygons for u[.743,.755], kappa[0,.1], shape[.4,1], conditional on separately replayed all-q gates. Source bridge bounds every physical r above by .3, but does not prove it lies above the mask lower boundary. No full-amplitude coverage or global norm claim.'}
    (HERE/'independent_continuous_source_masks.json').write_text(json.dumps(out,indent=2)+'\n')
    return out


def main():
    p=profile()
    gates=[spectral_gates(p,b)for b in (F(18,25),F(7171,10000))]
    assert (I(194)/I(10).sqrt()-I(F(211,5))).above(0)
    out={'status':'PASS_INDEPENDENT_DIFFUSE_R3_ANCHOR_PRICING_ONLY',
         'engine':selftest(),'profile_path':str(PROFILE.relative_to(ROOT)),
         'profile_sha256':sha(PROFILE),'frozen_theta':'1/123',
         'reserve0':'1/50000','reserve6':'1/5000000',
         'finite_spectral_gates':gates,'local_metadata':metadata(p),
         'coefficients_b7171':{k:str(v)for k,v in coefficients(p,F(7171,10000)).items()},
         'anchors':anchor_prices(p),
         'b_transfer':'For 0<b<=.72 all three b-dependent spectral gates weaken, while the local potential and scalar majorant are exactly unchanged.',
         'b_price_difference':'P_b-P_b0=(b-b0)*[r^2-6*beta^2/(b*b0)].',
         'scope':'Anchor and witness prices for kappa[0,.1] only, conditional on separately replayed whole-q gates. No source-domain or amplitude-strip coverage is asserted by this anchor artifact.'}
    (HERE/'independent_anchor_prices.json').write_text(json.dumps(out,indent=2)+'\n')
    masks=audit_masks(p)
    print(json.dumps({'status':out['status'],'engine':out['engine'],
                      'source_mask_status':masks['status'],
                      'mask_polygon_vertices':masks['polygon_kappa_shape_vertex_checks'],
                      'anchors':[{k:v for k,v in r.items()if k!='vertices'}for r in out['anchors']]},indent=2))


if __name__=='__main__':main()
