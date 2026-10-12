#!/usr/bin/env python3
"""Independent full-four-moment partition replay, using audited primitives.

No general_R3 producer is imported. The tree walker, raw-Hermite matrix,
endpoint objective, radius, exclusions and glue are reconstructed here.
Audited shared primitives are MPFR256 polynomial Gaussian integration,
rational Stein majorant and its certified second derivative lower bound.
"""
from fractions import Fraction as F
from pathlib import Path
from multiprocessing import Pool
import gzip,hashlib,json,sys,time,argparse

ROOT=Path(__file__).resolve().parents[2]
BASE=ROOT/'tonight/operator/projection/lower_amp_certificate'
CERT=ROOT/'tonight/operator/projection/certificate'
sys.path.insert(0,str(CERT))
from mpfr_interval import I,UP
from gaussian_max_certificate import max_integral
from certify_local_point import majorant
from certify_local_bernstein_region import second_lower
from audit_general_r3_parameters_tail import gates,tail,load

def decode(d):
    return {k:[F(x) for x in v] if isinstance(v,list) else F(v) for k,v in d.items()}

def build_matrix(p):
    h=p['zeta']**2/p['G'];N=[[F(0) for _ in range(4)] for _ in range(4)]
    N[0][0]=p['E00']+h*p['p0']**2
    N[0][2]=N[2][0]=p['E02']+h*p['p0']*p['p2']
    N[2][2]=p['E22']+h*p['p2']**2
    N[1][1]=p['R11'];N[1][3]=N[3][1]=p['R13'];N[3][3]=p['R33']
    return N

def init_worker(p,N,major):
    global WP,WN,WM
    WP,WN,WM=p,N,major

def directed_endpoint(key):
    q,*center=key;p=WP;N=WN
    gg=p['c']*q+p['d0']*(q*q*q-3*q)
    aa=p['zeta']*(-(q*q*q-3*q)+gg/(2*p['G']))
    raw=[q*p['gamma']+aa*p['p0'],F(0),
         q*(p['y']-p['E02']*p['r'])+aa*p['p2'],F(0)]
    raw=[raw[j]+sum(N[j][k]*center[k] for k in range(4)) for j in range(4)]
    # Raw Hermites -> ordinary cubic coefficients.
    poly=[raw[0]-raw[2],raw[1]-3*raw[3],raw[2],raw[3]]
    head=[F(0),p['ell']-3*p['beta'],F(0),p['beta']]
    quad=sum(center[j]*N[j][k]*center[k] for j in range(4) for k in range(4))
    value=2*max_integral([head,[-x for x in head],poly,[-x for x in poly]])
    value+=I(gg*gg/(4*p['G'])-quad-WM[q])
    return key,F(value.hi.text(UP,65))

def reconstruct(path):
    d=json.load(gzip.open(path,'rt'))
    assert d['status']=='PASS_GENERAL_PHYSICAL_R3_COMPACT_COVER'
    assert d['pending']==0 and not d['pending_boxes']
    p=decode(d['profile']);gates(p)
    a,b=F(d['qa']),F(d['qb']);assert 0<=a<b
    stack=[[(a,b),(F(-1),F(1)),(F(0),F(1)),(F(-3,2),F(3,2)),(F(-5,2),F(5,2))]]
    leaves=[];excluded=0
    for event in d['events']:
        assert stack,'Events beyond exhausted root'
        box=stack.pop()
        if event not in ('C','N'):
            axis=int(event);assert 0<=axis<5
            lo,hi=box[axis];mid=(lo+hi)/2
            left=list(box);right=list(box);left[axis]=(lo,mid);right[axis]=(mid,hi)
            stack.extend([right,left]);continue
        if event=='N':
            v=[lo if lo>0 else hi if hi<0 else F(0) for lo,hi in box[1:]]
            assert v[0]**2+v[1]**2+v[2]**2/2+v[3]**2/6>1
            excluded+=1
        else:leaves.append(box)
    assert not stack,'Unvisited partition cells'
    assert len(leaves)==d['accepted'] and excluded==d['excluded']
    assert len(d['events'])==d['nodes']
    return d,p,leaves,excluded

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--workers',type=int,default=4)
    ap.add_argument('--label',default='70');args=ap.parse_args()
    start=time.time()
    sources={str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest()
      for path in [BASE/'general_R3_cover.py',BASE/'general_R3_price.py']+
       [CERT/name for name in ('mpfr_interval.py','gaussian_max_certificate.py',
         'certify_local_point.py','certify_local_bernstein_region.py','rational_phi_bounds.py','bernstein_integer_bounds.py')]}
    paths=[BASE/f'general_u{args.label}_high_{part}.json.gz' for part in ('body','exterior')]
    cert_hash={str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest() for path in paths}
    groups=[reconstruct(path) for path in paths]
    d,p,_,_=groups[0]
    assert F(d['qa'])==0 and F(d['qb'])==8
    assert F(groups[1][0]['qa'])==8 and F(groups[1][0]['qb'])==64
    assert groups[1][1]==p
    rawpath=BASE/f'general_profile_u{args.label}_high_exact.json';raw=load(rawpath)
    assert all(p[k]==raw[k] for k in p if k!='major')
    assert p['major']==[raw['major'][j]+(F(d['reserve0']) if j==0 else
        F(d['reserve6']) if j==3 else 0) for j in range(4)]
    N=build_matrix(p);leaves=[bx for _,_,lx,_ in groups for bx in lx]
    keys=set();intervals=set()
    for bx in leaves:
        center=tuple((a+b)/2 for a,b in bx[1:]);intervals.add(bx[0])
        keys.update((q,)+center for q in bx[0])
    qnodes={key[0] for key in keys}
    majors={q:majorant(p,q) for q in qnodes}
    curvature={ab:max(F(0),-second_lower(p,*ab)) for ab in intervals}
    print(json.dumps({'stage':'PARTITION_AND_CURVATURE_PASS','events':sum(x[0]['nodes'] for x in groups),
      'leaves':len(leaves),'excluded':sum(x[3] for x in groups),'centers':len(keys),'qnodes':len(qnodes),
      'seconds':time.time()-start}),flush=True)
    endpoints={}
    with Pool(args.workers,initializer=init_worker,initargs=(p,N,majors)) as pool:
        for count,(key,value) in enumerate(pool.imap_unordered(directed_endpoint,keys,chunksize=100),1):
            endpoints[key]=value
            if count%4000==0:print(json.dumps({'stage':'DIRECTED_ENDPOINT_REPLAY','centers':count,
                'seconds':time.time()-start}),flush=True)
    C3=abs(p['zeta']*(-1+p['d0']/(2*p['G'])))*(abs(p['p0'])+abs(p['p2']))
    maximum=F(-10**20)
    for bx in leaves:
        a,b=bx[0];center=tuple((x+y)/2 for x,y in bx[1:]);width=[(y-x)/2 for x,y in bx[1:]]
        radius=sum(N[j][j]*width[j]**2 for j in range(4))
        radius+=2*abs(N[0][2])*width[0]*width[2]+2*abs(N[1][3])*width[1]*width[3]
        qreserve=(F(3,2)*b*C3+curvature[a,b]/8)*(b-a)**2
        bound=max(endpoints[(a,)+center],endpoints[(b,)+center])+radius+qreserve
        assert bound<0,('Recorded leaf does not certify',bx,bound)
        maximum=max(maximum,bound)
    t=tail(raw,F(64),F(d['reserve0']),F(d['reserve6']))
    assert all(hashlib.sha256((ROOT/f).read_bytes()).hexdigest()==h for f,h in sources.items())
    assert all(hashlib.sha256((ROOT/f).read_bytes()).hexdigest()==h for f,h in cert_hash.items())
    out={'status':'PASS_INDEPENDENT_GENERAL_R3_WHOLE_Q_PROFILE',
      'scope':'Whole scalar q line for one frozen full-E profile; global outer physical coverage remains separate.',
      'label':args.label,'events_replayed':sum(x[0]['nodes'] for x in groups),
      'leaves_rechecked':len(leaves),'excluded_boxes_rechecked':sum(x[3] for x in groups),
      'directed_endpoints':len(keys),'q_nodes':len(qnodes),'largest_leaf_residual':str(maximum),
      'tail':t,'profile_sha256':hashlib.sha256(rawpath.read_bytes()).hexdigest(),
      'certificate_sha256':cert_hash,'shared_source_hashes':sources,'seconds':time.time()-start}
    Path(__file__).with_name(f'general_r3_u{args.label}_whole_q_audit.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k not in ('tail','shared_source_hashes')},indent=2))

if __name__=='__main__':main()
