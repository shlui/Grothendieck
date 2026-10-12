"""Independent directed price and separate-convexity checks for one region.

The regional result becomes available only after whole-q local coverage for
the same exact profile and the recorded reserves has passed.
"""
from fractions import Fraction as F
from pathlib import Path
from itertools import product
import hashlib,json
import gzip
from mpfr_interval import I,UP,DOWN

ROOT=Path(__file__).resolve().parents[1]
PROFILE=ROOT/'operator/projection/certificate/frozen_R3_profile_exact.json'


def main():
    d=json.load(open(PROFILE));assert d['kind']=='EXACT_RATIONAL_R3_PROFILE'
    p={k:[F(x) for x in v] if isinstance(v,list) else F(v)
       for k,v in d['profile'].items()}
    q0,q2,q4,q6=p['major'];q0+=F(1,50000);q6+=F(1,1000000)
    theta=p['kappa']/(3*(p['kappa']+2))
    A=q6+p['zeta']**2;B=p['G']-p['d0']-p['zeta']**2
    assert A>=0 and B>=0
    Kmin=I(194)/I(10).sqrt()-I(F(633,15))
    assert Kmin.above(0)
    vy=p['y']/p['p2']-p['p0']*p['r']
    vertices=[]
    for u,r,k,z in product([F(7549,10000),F(7551,10000)],
                           [F(2574,10000),F(2577,10000)],
                           [F(51,100),F(55,100)],[F(2,5),F(1)]):
        K=I(194)/I(10).sqrt()-I(F(633,9)*(z-F(2,5)))
        power=I(k)*I(k).sqrt()
        c6=K*power
        sourceD=I(k*k)+I(F(9,4))*I(k*z)*I(k*z).sqrt()
        variance=k*(1+z)/9-2*theta*k/3+theta**2*(k+2)
        covariance=k/3-theta*(k+2)
        # The transported-r quadratic is simplified independently.
        rational=(u*u-2*p['ell']*u+6*p['beta']**2/p['eps']+
                  p['eps']*r*r-2*(p['gamma']-p['p0']*vy)*r+vy*vy+
                  q0+q2+(3+k)*q4+q6*(15+15*k)+
                  p['eta']*variance/2+p['nuisance2']*covariance+
                  p['zeta']**2*(6+9*k)-p['c']*k)
        value=I(rational)+I(A)*c6+I(B)*sourceD
        vertices.append({'u':str(u),'r':str(r),'kappa':str(k),'shape':str(z),
                         'lower':value.lo.text(DOWN,70),
                         'upper':value.hi.text(UP,70)})
    maximum=max(F(v['upper']) for v in vertices)
    cap=F(56334,100000)
    assert maximum<cap<F(40,71)
    localpath=Path(__file__).parent/'r3_bernstein_center_cover_audit.json'
    local=json.load(open(localpath))
    assert local['status']=='PASS_FRESH_PHYSICAL_CENTER_PARTITION_AND_WHOLE_Q_GLUE'
    coverpath=ROOT/'operator/projection/results/directed_bernstein_region_0_8.json.gz'
    assert local['certificate_sha256']==hashlib.sha256(coverpath.read_bytes()).hexdigest()
    cd=json.load(gzip.open(coverpath,'rt'))
    cp={k:[F(x) for x in v] if isinstance(v,list) else F(v)
        for k,v in cd['profile'].items()}
    assert all(cp[k]==p[k] for k in p if k!='major')
    assert cp['major']==[q0,q2,q4,q6]
    out={'status':'PASS_INDEPENDENT_REGIONAL_PAIRING_BOUND',
         'scope':'One explicit rectangle and all shapes; not a global norm bound.',
         'rectangle':{'u':['.7549','.7551'],'r':['.2574','.2577'],
                      'kappa':['.51','.55'],'shape':['2/5','1']},
         'constant_reserve':'1/50000','sixth_degree_reserve':'1/1000000',
         'vertex_count':len(vertices),'maximum_upper':str(maximum),
         'maximum_upper_decimal':str(float(maximum)),
         'rational_cap':str(cap),'margin_to_cap':str(cap-maximum),
         'separate_convexity':{
             'u_quadratic_coefficient':'1','r_quadratic_coefficient':str(p['eps']),
             'kappa_squared_coefficient':str(B),
             'kappa_three_halves_coefficient':'A*K(shape)+(9/4)*B*shape^(3/2)>0',
             'shape_three_halves_coefficient':'(9/4)*B*kappa^(3/2)>=0'},
         'profile_sha256':hashlib.sha256(PROFILE.read_bytes()).hexdigest(),
         'whole_q_local_certificate_sha256':local['certificate_sha256'],
         'source_moment_inputs':'Fresh633 sixth-cumulant comparison; Gaussian cubic product identities; exact D upper bound 9/4.',
         'vertices':vertices}
    output=Path(__file__).parent/'r3_outer_rectangle_audit.json'
    output.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k!='vertices'},indent=2))


if __name__=='__main__':main()
