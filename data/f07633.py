"""Independently traverse and recheck the recorded local partition.

The MPFR endpoint integrals and curvature enclosure are the separately audited
arithmetic primitives. Partition geometry, interpolation, and tail coefficient
algebra are reconstructed here rather than invoking the producer's cover.
"""
from fractions import Fraction as F
from pathlib import Path
from itertools import product
from math import comb
import gzip, hashlib, json, sys, time

ROOT=Path(__file__).resolve().parents[1]
CERT=ROOT/'operator/projection/certificate'
RESULTS=ROOT/'operator/projection/results'
sys.path.insert(0,str(CERT))
from certify_local_region import second_lower
from certify_local_point import context,value,majorant
from mpfr_interval import I,UP


def decode(d):
    return {k:[F(x) for x in v] if isinstance(v,list) else F(v)
            for k,v in d.items()}


def tail_coefficients(p,R):
    Aphi=sum(abs(x) for x in p['phi'])
    Ader=sum(j*j*abs(x) for j,x in enumerate(p['phi']))
    theta=p['kappa']/(3*(p['kappa']+2))
    derivative=Aphi/R**2+4*Ader/R**4
    stein=Aphi+theta*Aphi+Aphi/R**2+4*theta*Ader/R**2+4*Ader/R**4
    aa=abs(p['nuisance0']-p['nuisance2'])+derivative
    bb=abs(p['nuisance2'])
    rho=1+p['zeta']**2/p['G']
    c=p['c']-3*p['d0'];d=p['d0']
    l0=p['gamma']+p['zeta']*(3+c/(2*p['G']))*p['p0']
    l2=p['y']-p['p0']*p['p2']*p['r']+p['zeta']*(3+c/(2*p['G']))*p['p2']
    l3=abs(p['zeta']*(-1+d/(2*p['G'])))*(abs(p['p0'])+abs(p['p2']))
    trace=rho*(p['p0']**2+2*p['p2']**2)+p['R11']+6*p['R33']
    constant=2*abs(p['ell'])+5*abs(p['beta'])+trace
    pol=[F(0)]*7
    for j,cj in enumerate(p['major']):pol[2*j]=cj
    pol[0]-=constant+stein+aa**2/(2*p['eta'])
    pol[1]-=2*(abs(l0)+abs(l2))
    pol[2]-=aa*bb/p['eta']+c*c/(4*p['G'])
    pol[3]-=2*l3
    pol[4]-=bb*bb/(2*p['eta'])+abs(c*d)/(2*p['G'])
    pol[6]-=d*d/(4*p['G'])
    shifted=[sum(comb(j,k)*pol[j]*R**(j-k) for j in range(k,7))
             for k in range(7)]
    return pol,shifted


def main():
    path=Path(sys.argv[1]) if len(sys.argv)>1 else RESULTS/'directed_region_8_64.json.gz'
    d=json.load(gzip.open(path,'rt'));p=decode(d['profile'])
    assert d['status']=='PASS_COMPACT_LOCAL_SUPREMUM'
    assert p['eta']>0 and p['G']>0
    assert p['p0']**2+2*p['p2']**2<1
    for a,b,c in [(p['R11'],p['R13'],p['R33']),
                  (p['R11']+1,p['R13'],p['R33']-p['eps']/6),
                  (p['R11']-p['eps'],p['R13'],p['R33'])]:
        assert a>=0 and c>=0 and a*c>=b*b
    assert p['p2']**2>=p['eps']*(1+p['p0']**2)/2
    qa,qb=F(d['qa']),F(d['qb']);assert 0<=qa<qb
    stack=[[(qa,qb),(F(-1),F(1)),(F(-1),F(1)),(F(-5,2),F(5,2))]]
    rho=1+p['zeta']**2/p['G']
    C3=abs(p['zeta']*(-1+p['d0']/(2*p['G'])))*(abs(p['p0'])+abs(p['p2']))
    diag=[rho,p['R11'],p['R33']]
    corner={};C={};maj={};curvature={};accepted=0;largest=F(-10**20)
    start=time.time()
    for event in d['events']:
        assert stack,'Event sequence continued after covering root.'
        box=stack.pop()
        if event!='C':
            axis=int(event);assert 0<=axis<4
            a,b=box[axis];mid=(a+b)/2
            left=list(box);right=list(box)
            left[axis]=(a,mid);right[axis]=(mid,b)
            stack.append(right);stack.append(left)
            continue
        a,b=box[0]
        if (a,b) not in curvature:curvature[a,b]=max(F(0),-second_lower(p,a,b))
        endpoints=[]
        for q,m0,m1,m3 in product(*box):
            key=(q,m0,m1,m3)
            if key not in corner:
                if q not in C:C[q]=context(p,q);maj[q]=majorant(p,q)
                corner[key]=value(p,C[q],(m0,m1,m3))-I(maj[q])
            endpoints.append(corner[key])
        allowance=(F(3,2)*b*C3+curvature[a,b]/8)*(b-a)**2
        allowance+=sum(v*(r-l)**2/4 for v,(l,r) in zip(diag,box[1:]))
        upper=I(max(x.hi for x in endpoints),max(x.hi for x in endpoints))+I(allowance)
        assert upper.below(0),('Rejected recorded leaf',box,upper)
        largest=max(largest,F(upper.hi.text(UP,50)));accepted+=1
    assert not stack,'Partition has unvisited cells.'
    assert accepted==d['accepted'] and len(d['events'])==d['nodes']
    tailpath=RESULTS/'directed_local_tail.json';td=json.load(open(tailpath))
    pol,shift=tail_coefficients(p,F(td['R']))
    assert pol==[F(x) for x in td['difference_polynomial']]
    assert shift==[F(x) for x in td['shifted_coefficients']]
    assert all(x>=0 for x in shift)
    out={'status':'PASS_RECORDED_COMPACT_PARTITION_AND_RATIONAL_TAIL',
         'scope':'Frozen profile only; no outer source-parameter cover.',
         'compact_interval':[str(qa),str(qb)],'leaves_rechecked':accepted,
         'partition_events':len(d['events']),'directed_corners':len(corner),
         'tail_from':td['R'],'tail_nonnegative_shifted_coefficients':len(shift),
         'largest_leaf_residual':str(largest),'seconds':time.time()-start,
         'certificate_sha256':hashlib.sha256(path.read_bytes()).hexdigest()}
    outpath=Path(__file__).parent/'r3_local_partition_audit.json'
    outpath.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))


if __name__=='__main__':main()
