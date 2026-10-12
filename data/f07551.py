"""Independent directed audit of three constrained-r physical polygons.

The polygons are affine images of (u,t) rectangles.  This program proves
their bounds, not all-r amplitude coverage or a global operator norm.
"""
from fractions import Fraction as F
from pathlib import Path
from itertools import product
import gzip,hashlib,json
from mpfr_interval import I,UP,DOWN

ROOT=Path(__file__).resolve().parents[1]
PROFILE=ROOT/'operator/projection/certificate/frozen_R3_profile_exact.json'
RESULTS=ROOT/'operator/projection/results'


def main():
    data=json.load(open(PROFILE));assert data['kind']=='EXACT_RATIONAL_R3_PROFILE'
    p={k:[F(x) for x in v] if isinstance(v,list) else F(v)
       for k,v in data['profile'].items()}
    q0,q2,q4,q6=p['major'];q0+=F(1,50000);q6+=F(1,1000000)
    local=json.load(open(Path(__file__).parent/'r3_bernstein_center_cover_audit.json'))
    assert local['status']=='PASS_FRESH_PHYSICAL_CENTER_PARTITION_AND_WHOLE_Q_GLUE'
    coverpath=RESULTS/'directed_bernstein_region_0_8.json.gz'
    assert local['certificate_sha256']==hashlib.sha256(coverpath.read_bytes()).hexdigest()
    cp=json.load(gzip.open(coverpath,'rt'))['profile']
    assert all(cp[k]==data['profile'][k] for k in cp if k!='major')
    assert [F(x) for x in cp['major']]==[q0,q2,q4,q6]
    theta=p['kappa']/(3*(p['kappa']+2));vy=p['y']/p['p2']-p['p0']*p['r']
    A=q6+p['zeta']**2;B=p['G']-p['d0']-p['zeta']**2
    assert A>=0 and B>=0
    assert (I(194)/I(10).sqrt()-I(F(633,15))).above(0)
    assert p['r']==F(5151,20000)
    names=['target_polygon','fallback_polygon','fallback_wide_kappa_polygon']
    cases=[]
    for name in names:
        sourcepath=RESULTS/(name+'.json');source=json.load(open(sourcepath))
        assert source['reserve0']=='1/50000' and source['reserve6']=='1/1000000'
        assert source['shape']==['2/5','1']
        ue=[F(x) for x in source['u']];te=[F(x) for x in source['t']]
        ke=[F(x) for x in source['kappa']];cap=F(source['cap'])
        saved={(F(row['u']),F(row['t']),F(row['kappa']),F(row['shape'])):row
               for row in source['rows']}
        assert len(saved)==16
        vertices=[]
        for u,t,k,z in product(ue,te,ke,[F(2,5),F(1)]):
            r=F(5151,20000)-3*(u-F(151,200))+t
            assert r>0
            row=saved[u,t,k,z];assert F(row['r'])==r
            K=I(194)/I(10).sqrt()-I(F(633,9)*(z-F(2,5)))
            c6=K*I(k)*I(k).sqrt()
            sourceD=I(k*k)+I(F(9,4))*I(k*z)*I(k*z).sqrt()
            variance=k*(1+z)/9-2*theta*k/3+theta**2*(k+2)
            covariance=k/3-theta*(k+2)
            rational=(u*u-2*p['ell']*u+6*p['beta']**2/p['eps']+
                      p['eps']*r*r-2*(p['gamma']-p['p0']*vy)*r+vy*vy+
                      q0+q2+(3+k)*q4+q6*(15+15*k)+
                      p['eta']*variance/2+p['nuisance2']*covariance+
                      p['zeta']**2*(6+9*k)-p['c']*k)
            value=I(rational)+I(A)*c6+I(B)*sourceD
            upper=F(value.hi.text(UP,70))
            assert upper<=F(row['upper'])
            vertices.append({'u':str(u),'t':str(t),'r':str(r),'kappa':str(k),
                             'shape':str(z),'upper':str(upper)})
        maximum=max(F(v['upper']) for v in vertices)
        assert maximum<F(source['maximum_upper'])<cap
        rvalues=[F(5151,20000)-3*(u-F(151,200))+t for u,t in product(ue,te)]
        cases.append({'name':name,'status':'PASS_CONSTRAINED_R_POLYGON',
                      'u':[str(x) for x in ue],'t':[str(x) for x in te],
                      'kappa':[str(x) for x in ke],'shape':['2/5','1'],
                      'r_global_extrema':[str(min(rvalues)),str(max(rvalues))],
                      'cap':str(cap),'maximum_upper':str(maximum),
                      'maximum_upper_decimal':str(float(maximum)),
                      'margin_to_cap':str(cap-maximum),'vertices':vertices,
                      'producer_record_sha256':hashlib.sha256(sourcepath.read_bytes()).hexdigest()})
    # Joint convexity of the transported u,r quadratic is retained after
    # the affine map: Hessian/2=[[1+9b,-3b],[-3b,b]], determinant=b>0.
    b=p['eps'];hessian=[[1+9*b,-3*b],[-3*b,b]]
    assert hessian[0][0]>0 and hessian[0][0]*hessian[1][1]-hessian[0][1]**2==b>0
    out={'status':'PASS_INDEPENDENT_MOVING_POLYGON_PAIRING_BOUNDS',
         'scope':'Constrained-r polygons only, not an all-r amplitude cover or global operator norm.',
         'parameterization':'r=5151/20000-3*(u-151/200)+t',
         'affine_quadratic_hessian_divided_by_2':[[str(x) for x in row] for row in hessian],
         'affine_quadratic_determinant':str(b),
         'shape_endpoints':'Shared dimensionless cubic shape z=1296*s/kappa in [2/5,1].',
         'source_inputs':'Fresh633 sixth cumulant bound; exact cubic D upper bound9/4; Gaussian product identities.',
         'profile_sha256':hashlib.sha256(PROFILE.read_bytes()).hexdigest(),
         'whole_q_local_certificate_sha256':local['certificate_sha256'],
         'cases':cases}
    output=Path(__file__).parent/'r3_moving_polygons_audit.json'
    output.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({**{k:v for k,v in out.items() if k!='cases'},
                      'cases':[{k:v for k,v in c.items() if k!='vertices'} for c in cases]},indent=2))


if __name__=='__main__':main()
