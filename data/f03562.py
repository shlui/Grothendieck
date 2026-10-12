#!/usr/bin/env python3
"""Independent pricing of the new b7171 diffuse mask and remainder box.

Reuses only independently written expanded-price definitions. No producer
module is imported and no local q gate or native-domain union is executed.
"""
from fractions import Fraction as F
from pathlib import Path
from itertools import product
import gzip
import hashlib
import importlib.util
import json

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
BASE=ROOT/'tonight/operator/projection/diffuse_certificate'
LOWER=ROOT/'tonight/operator/projection/lower_amp_certificate'
PROFILE=BASE/'profile_interval_b7171_newworst.json'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def independent_module(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


ENGINE=HERE.parent/'diffuse_R3/independent_diffuse_source_price.py'
OLD=independent_module('independent_expanded_diffuse_price',ENGINE)
LOW_ENGINE=HERE.parent/'lower_amp_R3/independent_source_price.py'
LP=independent_module('independent_expanded_lower_price',LOW_ENGINE)


def exact_profile(path):
    row=json.loads(path.read_text())
    assert row['kind']=='EXACT_RATIONAL_R3_PROFILE'
    return {k:[F(t)for t in v]if isinstance(v,list)else F(v)
            for k,v in row['profile'].items()}


def hash_ledger(source):
    results=[]
    for filename,digest in source['sha256'].items():
        path=Path(filename)
        if not path.is_absolute():path=ROOT/path
        assert sha(path)==digest
        results.append({'path':str(path.relative_to(ROOT)),'sha256':digest})
    return results


def diffuse_metadata(p):
    expected={k:[str(t)for t in v]if isinstance(v,list)else str(v)for k,v in p.items()}
    expected['major'][0]=str(p['major'][0]+F(1,50000))
    expected['major'][3]=str(p['major'][3]+F(1,5000000))
    results=[]
    for prefix,a,b in [('body','0','8'),('exterior','8','64')]:
        path=BASE/f'{prefix}_interval_b7171_newworst.json.gz'
        with gzip.open(path,'rt')as stream:record=json.load(stream)
        assert record['status']=='PASS_PHYSICAL_BERNSTEIN_LOCAL_SUPREMUM'
        assert record['pending']==0 and (record['qa'],record['qb'])==(a,b)
        assert F(record['reserve0'])==F(1,50000)
        assert F(record['reserve'])==F(1,5000000)
        assert record['profile']==expected
        results.append({'path':str(path.relative_to(ROOT)),'sha256':sha(path),
                        'q':[a,b],'recorded_leaves':record['accepted'],
                        'q_gate_replayed_here':False})
    return results


def audit_diffuse():
    path=BASE/'critical_source_upper_masks_b7171_newworst.json'
    initial_sha=sha(path);source=json.loads(path.read_text())
    assert source['status']=='PASS_ALL_Q_CONTINUOUS_CRITICAL_HIGH_SOURCE_MASKS'
    assert F(source['reserve0'])==F(1,50000) and F(source['reserve6'])==F(1,5000000)
    hashes=hash_ledger(source)
    p=exact_profile(PROFILE)
    assert p['eps']==F(7171,10000) and p['kappa']==F(1,20)
    b=p['eps'];cap=F(40,71)-F(1,10**6)
    c=OLD.coefficients(p,b)
    spectral=OLD.spectral_gates(p,b)
    metadata=diffuse_metadata(p)
    for old,new in zip(source['regions'],metadata):
        assert old['q']==new['q'] and old['leaves']==new['recorded_leaves']
    assert len(source['regions'])==2
    tail_path=BASE/'tail_interval_b7171_newworst.json'
    assert json.loads(tail_path.read_text())==source['tail']
    assert source['tail']['status']=='PASS' and F(source['tail']['R'])==64
    bridge=OLD.support_bridge(source)
    assert len(source['masks'])==1
    mask=source['masks'][0]
    assert F(mask['operator_b'])==b and F(mask['cap'])==cap
    assert tuple(map(F,mask['kappa']))==(F(0),F(1,10))
    assert tuple(map(F,mask['shape']))==(F(2,5),F(1))
    corners=[]
    for k,z in product((F(0),F(1,10)),(F(2,5),F(1))):
        actual=OLD.price(p,c,F(0),F(0),k,z)
        corners.append({'kappa':str(k),'shape':str(z),'price':actual.record(),
                        'ceiling_60_significant_digits':str(OLD.upward_decimal(actual,59))})
    constant=max(F(x['ceiling_60_significant_digits'])for x in corners)
    assert constant==F(mask['source_constant_upper'])
    gamma=p['gamma']-p['p0']*c['vy']
    assert gamma==F(mask['effective_gamma']) and gamma>b*F(3,10)
    def bound(u,r):return u*u-2*p['ell']*u+b*r*r-2*gamma*r+constant
    previous=F(743,1000);previous_r=None;cells=[]
    assert len(mask['cells'])==12
    for index,cell in enumerate(mask['cells']):
        ua,ub=map(F,cell['u']);ra,rb=map(F,cell['lower_r_vertices'])
        assert ua==previous and ub-ua==F(1,1000)
        if previous_r is not None:assert previous_r==ra
        previous,previous_r=ub,rb
        assert F(cell['upper_r'])==F(3,10) and F(cell['cap'])==cap
        assert 0<ra<F(3,10) and 0<rb<F(3,10)
        for u,r in ((ua,ra),(ub,rb)):
            assert (r*10**12).denominator==1
            assert bound(u,r)<cap and bound(u,r-F(1,10**12))<cap
            assert bound(u,r-F(2,10**12))>=cap
        ur=((ua,ra),(ub,rb),(ua,F(3,10)),(ub,F(3,10)))
        maximum=max(bound(u,r)for u,r in ur)
        assert maximum==F(cell['maximum_upper']) and maximum<cap
        vertices=[]
        for u,r in ur:
            for k,z in product((F(0),F(1,10)),(F(2,5),F(1))):
                actual=OLD.price(p,c,u,r,k,z)
                assert actual.below(bound(u,r)) and actual.below(cap)
                covariance=k/3-c['theta']*(k+2)
                assert covariance==(40*k-2)/123
                y=p['y']+p['p0']*p['p2']*(r-p['r'])
                assert y-p['p0']*p['p2']*r==p['y']-p['p0']*p['p2']*p['r']
                oldr=(b-p['p0']**2)*r*r-2*p['gamma']*r+y*y/p['p2']**2
                assert oldr==b*r*r-2*gamma*r+c['vy']**2
                vertices.append({'u':str(u),'r':str(r),'kappa':str(k),'shape':str(z),
                                 'price':actual.record(),'rational_bound':str(bound(u,r)),
                                 'plus_covariance_term':str(p['nuisance2']*covariance)})
        cells.append({'index':index,'u':[str(ua),str(ub)],
                      'lower_r_vertices':[str(ra),str(rb)],'upper_r':'3/10',
                      'maximum_rational_bound':str(maximum),'margin_to_cap':str(cap-maximum),
                      'vertices':vertices})
    assert previous==F(151,200)
    example=mask['example_rational_point']
    eu,er=F(example['u']),F(example['r'])
    assert eu==F('0.750515528144') and er==F('0.270685360072')
    actual=max(OLD.upward_decimal(OLD.price(p,c,eu,er,k,z),50)
               for k,z in product((F(0),F(1,10)),(F(2,5),F(1))))
    assert actual==F(example['upper']) and actual<cap
    containing=next(x for x in cells if F(x['u'][0])<=eu<=F(x['u'][1]))
    ua,ub=map(F,containing['u']);ra,rb=map(F,containing['lower_r_vertices'])
    lower=ra+(rb-ra)*(eu-ua)/(ub-ua)
    assert lower<=er<F(3,10)
    assert sha(path)==initial_sha
    out={'status':'PASS_INDEPENDENT_NEW_B7171_DIFFUSE_SOURCE_MASK_PRICING',
         'engine':OLD.selftest(),'source_path':str(path.relative_to(ROOT)),
         'source_sha256':initial_sha,'profile_sha256':sha(PROFILE),
         'independent_expanded_price_engine_sha256':sha(ENGINE),
         'checked_hash_ledger':hashes,'local_q_metadata':metadata,
         'tail_record_sha256':sha(tail_path),'tail_proof_replayed_here':False,
         'operator_b':str(b),'cap':str(cap),'frozen_theta':'1/123',
         'coefficients':{k:str(v)for k,v in c.items()},'finite_spectral_gates':spectral,
         'source_bridge':bridge,'source_constant':str(constant),
         'source_constant_corners':corners,'effective_gamma':str(gamma),
         'cell_count':12,'polygon_kappa_shape_vertex_checks':192,'cells':cells,
         'newer_witness':{'u':str(eu),'r':str(er),'mask_lower_r':str(lower),
                         'inside_mask':True,'margin_above_lower_boundary':str(er-lower),
                         'maximum_uniform_price':str(actual)},
         'convexity':'Expanded price is jointly convex in (u,r), separately convex in kappa and shape because A>0,B>0,K(z)>0. Polygon vertices and moment endpoints therefore cover the exact domains.',
         'scope':'Only the twelve high-r polygons for u[.743,.755], kappa[0,.1], all admissible shape; conditional on separate all-q replay. No sources below the lower boundary, whole-amplitude, or global norm claim.'}
    (HERE/'independent_diffuse_source_masks.json').write_text(json.dumps(out,indent=2)+'\n')
    return out


def audit_remainder_box():
    path=LOWER/'b7171_native_remainder_extended_box.json'
    initial_sha=sha(path);source=json.loads(path.read_text())
    assert source['status']=='PASS_ALL_Q_B7171_NATIVE_REMAINDER_BOX'
    hashes=hash_ledger(source)
    p=exact_profile(LOWER/'profile_u.70.json')
    assert p['eps']==F(18,25) and F(source['baseline_b'])==p['eps']
    b=F(7171,10**4);assert F(source['operator_b'])==b
    pb={**p,'eps':b}
    c=LP.coefficients(pb,F(1,50000),F(1,5000000))
    assert c['theta']==F(53,759)
    spectral=OLD.spectral_gates(p,b)
    metadata=LP.bind_local_metadata('70',p,'finer',F(1,50000),F(1,5000000))
    assert len(source['regions'])==2
    for old,new in zip(source['regions'],metadata):
        assert old['q']==new['q_interval'] and old['leaves']==new['recorded_leaves']
    finer=json.loads((LOWER/'lower_source_domain_finer.json').read_text())
    fine70=next(x for x in finer['cells']if x['profile']=='70')
    assert fine70['tail']==source['tail']
    bounds=((F(137,200),F(1381,2000)),(F(47,125),F(39,100)),
            (F(0),F(1)),(F(2,5),F(1)))
    assert tuple(tuple(map(F,source[k]))for k in ('u','r','kappa','shape'))==bounds
    assert F(source['reserve0'])==F(1,50000) and F(source['reserve6'])==F(1,5000000)
    cap=F(5633,10000);assert F(source['rounded_exact_cap'])==cap
    vertices=[]
    assert len(source['vertices'])==16
    for old,(u,r,k,z)in zip(source['vertices'],product(*bounds)):
        assert tuple(F(old[x])for x in ('u','r','kappa','shape'))==(u,r,k,z)
        actual=LP.price(pb,c,u,r,k,z)
        saved=F(old['upper'])
        assert actual.below(saved) and actual.below(cap)
        assert OLD.upward_decimal(actual,60)==saved
        baseline=LP.coefficients(p,F(1,50000),F(1,5000000))
        change=(b-p['eps'])*(r*r-6*p['beta']**2/(b*p['eps']))
        assert change==(b-p['eps'])*r*r+6*p['beta']**2*(1/b-1/p['eps'])
        assert c['constant']-baseline['constant']+(b-p['eps'])*r*r==change
        vertices.append({'u':str(u),'r':str(r),'kappa':str(k),'shape':str(z),
                         'price':actual.record(),'saved_upper':str(saved),
                         'exact_b_price_change':str(change),
                         'plus_covariance_term':str(p['nuisance2']*(k/3-c['theta']*(k+2)))})
    maximum=max(F(v['saved_upper'])for v in vertices)
    assert maximum==F(source['maximum_upper']) and maximum<cap<F(40,71)
    assert F(source['margin_to_target'])==F(40,71)-maximum
    assert sha(path)==initial_sha
    out={'status':'PASS_INDEPENDENT_B7171_EXTENDED_REMAINDER_BOX_PRICING',
         'engine':OLD.selftest(),'source_path':str(path.relative_to(ROOT)),
         'source_sha256':initial_sha,'independent_lower_price_engine_sha256':sha(LOW_ENGINE),
         'profile_sha256':sha(LOWER/'profile_u.70.json'),'checked_hash_ledger':hashes,
         'local_q_metadata':metadata,'operator_b':str(b),'baseline_b':str(p['eps']),
         'frozen_theta':'53/759','finite_spectral_gates':spectral,
         'u':list(map(str,bounds[0])),'r':list(map(str,bounds[1])),
         'kappa':list(map(str,bounds[2])),'shape':list(map(str,bounds[3])),
         'cap':str(cap),'maximum_upper':str(maximum),'vertices':vertices,
         'scope':'Exact restricted source box only, conditional on previously accepted finer profile70 all-q gates. Native complementary-domain union is not replayed or asserted here.'}
    (HERE/'independent_extended_remainder_box.json').write_text(json.dumps(out,indent=2)+'\n')
    return out


def main():
    diffuse=audit_diffuse()
    box=audit_remainder_box()
    print(json.dumps({'diffuse_status':diffuse['status'],
                      'diffuse_polygon_vertices':diffuse['polygon_kappa_shape_vertex_checks'],
                      'diffuse_cap':diffuse['cap'],'newer_witness':diffuse['newer_witness'],
                      'remainder_status':box['status'],'remainder_cap':box['cap'],
                      'remainder_maximum':box['maximum_upper']},indent=2))


if __name__=='__main__':main()
