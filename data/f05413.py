"""Fresh independent partition replay for the physical-center R3 cover.

Each accepted box bounds original physical assignments whose moments lie in
that box, by expanding the PSD quadratic about its center.  It does not assume
the dual potential is bounded by its value at a center.
"""
from fractions import Fraction as F
from pathlib import Path
import gzip,hashlib,json,sys,time

ROOT=Path(__file__).resolve().parents[1]
CERT=ROOT/'operator/projection/certificate'
RESULTS=ROOT/'operator/projection/results'
sys.path.insert(0,str(CERT))
from certify_local_bernstein_region import second_lower
from certify_local_point import context,value,majorant
from mpfr_interval import I,UP
from replay_r3_local_cover import decode,tail_coefficients


def main():
    path=RESULTS/'directed_bernstein_region_0_8.json.gz'
    d=json.load(gzip.open(path,'rt'));p=decode(d['profile'])
    assert d['status']=='PASS_PHYSICAL_BERNSTEIN_LOCAL_SUPREMUM'
    assert p['eta']>0 and p['G']>0
    pnorm=p['p0']**2+2*p['p2']**2
    assert 0<pnorm<1
    for a,b,c in [(p['R11'],p['R13'],p['R33']),
                  (p['R11']+1,p['R13'],p['R33']-p['eps']/6),
                  (p['R11']-p['eps'],p['R13'],p['R33'])]:
        assert a>=0 and c>=0 and a*c>=b*b
    assert p['p2']**2>=p['eps']*(1+p['p0']**2)/2
    qa,qb=F(d['qa']),F(d['qb']);assert qa==0 and qb==8
    files=['certify_local_bernstein_region.py','certify_local_point.py',
           'rational_phi_bounds.py','bernstein_integer_bounds.py',
           'gaussian_max_certificate.py','mpfr_interval.py']
    initial_hashes={f:hashlib.sha256((CERT/f).read_bytes()).hexdigest() for f in files}
    stack=[[(qa,qb),(F(-1),F(1)),(F(0),F(1)),(F(-5,2),F(5,2))]]
    rho=1+p['zeta']**2/p['G']
    C3=abs(p['zeta']*(-1+p['d0']/(2*p['G'])))*(abs(p['p0'])+abs(p['p2']))
    cache={};contexts={};majorants={};curvature={}
    accepted=excluded=0;largest=F(-10**20);start=time.time()
    for count,event in enumerate(d['events'],1):
        assert stack,'Events continued after root was exhausted.'
        box=stack.pop()
        if event not in ['C','N']:
            axis=int(event);assert 0<=axis<4
            a,b=box[axis];mid=(a+b)/2
            left=list(box);right=list(box)
            left[axis]=(a,mid);right[axis]=(mid,b)
            stack.append(right);stack.append(left);continue
        nearest=[a if a>0 else b if b<0 else F(0) for a,b in box[1:]]
        outside=nearest[0]**2/pnorm+nearest[1]**2+nearest[2]**2/6>1
        if event=='N':
            assert outside,('Unsound ellipsoid exclusion',box)
            excluded+=1;continue
        a,b=box[0]
        if (a,b) not in curvature:curvature[a,b]=max(F(0),-second_lower(p,a,b))
        center=tuple((x+y)/2 for x,y in box[1:])
        endpoints=[]
        for q in (a,b):
            key=(q,)+center
            if key not in cache:
                if q not in contexts:contexts[q]=context(p,q);majorants[q]=majorant(p,q)
                cache[key]=value(p,contexts[q],center)-I(majorants[q])
            endpoints.append(cache[key])
        # Exact maximum of the PSD center remainder over the rectangular
        # radius box: diagonal radii and the one odd cross term.
        rs,r1,r3=[(y-x)/2 for x,y in box[1:]]
        radius=rho*rs*rs+p['R11']*r1*r1+p['R33']*r3*r3+2*abs(p['R13'])*r1*r3
        qreserve=(F(3,2)*b*C3+curvature[a,b]/8)*(b-a)**2
        upper=I(max(x.hi for x in endpoints),max(x.hi for x in endpoints))+I(radius+qreserve)
        assert upper.below(0),('Rejected recorded physical-center leaf',box,upper)
        largest=max(largest,F(upper.hi.text(UP,50)));accepted+=1
        if count%1000==0:
            print(json.dumps({'events':count,'leaves':accepted,'centers':len(cache),
                              'seconds':time.time()-start}),flush=True)
    assert not stack,'Partition has unvisited cells.'
    assert accepted==d['accepted'] and excluded==d['excluded']
    assert len(d['events'])==d['nodes']
    # The previously replayed [8,64] certificate and tail use a smaller
    # constant reserve.  Their majorant is therefore dominated by this one.
    oldpath=RESULTS/'directed_region_8_64.json.gz'
    old=json.load(gzip.open(oldpath,'rt'));op=decode(old['profile'])
    assert old['status']=='PASS_COMPACT_LOCAL_SUPREMUM'
    assert all(p[k]==op[k] for k in p if k!='major')
    assert p['major'][0]>=op['major'][0] and p['major'][3]>=op['major'][3]
    assert p['major'][1:3]==op['major'][1:3]
    previous=json.load(open(Path(__file__).parent/'r3_local_partition_audit.json'))
    assert previous['status']=='PASS_RECORDED_COMPACT_PARTITION_AND_RATIONAL_TAIL'
    assert previous['certificate_sha256']==hashlib.sha256(oldpath.read_bytes()).hexdigest()
    td=json.load(open(RESULTS/'directed_local_tail.json'))
    tail_pol,tail_shift=tail_coefficients(op,F(td['R']))
    assert [F(x) for x in td['difference_polynomial']]==tail_pol
    assert [F(x) for x in td['shifted_coefficients']]==tail_shift
    assert all(x>=0 for x in tail_shift)
    final_hashes={f:hashlib.sha256((CERT/f).read_bytes()).hexdigest() for f in files}
    assert initial_hashes==final_hashes,'Producer source changed during replay.'
    out={'status':'PASS_FRESH_PHYSICAL_CENTER_PARTITION_AND_WHOLE_Q_GLUE',
         'scope':'Whole scalar line for one frozen profile; outer source coverage is separate.',
         'compact_interval':[str(qa),str(qb)],'leaves_rechecked':accepted,
         'excluded_boxes_rechecked':excluded,'partition_events':len(d['events']),
         'directed_centers':len(cache),'q_nodes':len(contexts),
         'largest_leaf_residual':str(largest),'seconds':time.time()-start,
         'constant_reserve':d['reserve0'],'sixth_degree_reserve':d['reserve'],
         'larger_q_glue':'Previously replayed [8,64] and >=64 tail, dominated constant majorant.',
         'certificate_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
         'source_hashes':initial_hashes}
    output=Path(__file__).parent/'r3_bernstein_center_cover_audit.json'
    output.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))


if __name__=='__main__':main()
