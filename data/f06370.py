"""Directed physical four-moment R3 cover for a general even PSD matrix.

This isolated extension leaves the independently audited rank-one producer
unchanged.  It accepts finite compact scalar regions only; global source
parameter coverage and independent audit remain separate.
"""
from fractions import Fraction as Q
from pathlib import Path
from itertools import product
import argparse,gzip,json,sys,time,math
CORE=Path(__file__).resolve().parents[1]/'certificate'
sys.path.insert(0,str(CORE))
from mpfr_interval import I,UP
from gaussian_max_certificate import max_integral
from certify_local_point import rational,majorant
from certify_local_bernstein_region import second_lower

def psd(a,b,c):return a>=0 and c>=0 and a*c>=b*b

def validate(p):
    assert p['eps']>0 and p['kappa']>=0 and p['eta']>0 and p['G']>0
    p0,p2=p['p0'],p['p2']
    assert p0*p0+2*p2*p2<1
    for a,b,c in [(p['E00'],p['E02'],p['E22']),(p['E00']+1,p['E02'],p['E22']-p['eps']/2),
                  (p['E00']-p0*p0,p['E02']-p0*p2,p['E22']-p2*p2),
                  (p['R11'],p['R13'],p['R33']),(p['R11']+1,p['R13'],p['R33']-p['eps']/6),
                  (p['R11']-p['eps'],p['R13'],p['R33'])]:assert psd(a,b,c),(a,b,c)
    return p

def freeze(path):
    d=json.loads(Path(path).read_text())
    if d.get('kind')=='EXACT_RATIONAL_GENERAL_R3_PROFILE':
        return validate({k:[Q(t) for t in v] if isinstance(v,list) else Q(v) for k,v in d['profile'].items()})
    x=d['history'][-1]['x'];n=d['order'];ia=14+n
    p0=rational(d['p'][0]);p2=rational(d['p'][1]/math.sqrt(2))
    p=dict(ell=rational(x[0]),beta=rational(x[1]/math.sqrt(6)),
           p0=p0,p2=p2,E00=rational(x[2])+Q(1,10**7),E02=rational(x[3]/math.sqrt(2)),E22=rational(x[4]/2)+Q(1,10**7),
           R11=rational(x[5])+Q(1,10**7),R13=rational(x[6]/math.sqrt(6)),R33=rational(x[7]/6)+Q(1,10**7),
           gamma=rational(x[8]),y=rational(x[9]/math.sqrt(2)),r=Q(str(d['r'])),kappa=Q(str(d['kappa'])),eps=Q(18,25),
           major=[rational(t)for t in x[10:14]],phi=[rational(t)for t in x[14:14+n]],
           nuisance0=rational(x[ia]),nuisance2=rational(x[ia+1]),eta=rational(x[ia+2]),zeta=rational(x[ia+3]),
           c=rational(x[ia+4]),d0=rational(x[ia+5]),G=rational(x[ia+6]))
    return validate(p)

def matrix(p):
    h=p['zeta']**2/p['G']
    return [[p['E00']+h*p['p0']**2,Q(0),p['E02']+h*p['p0']*p['p2'],Q(0)],
            [Q(0),p['R11'],Q(0),p['R13']],
            [p['E02']+h*p['p0']*p['p2'],Q(0),p['E22']+h*p['p2']**2,Q(0)],
            [Q(0),p['R13'],Q(0),p['R33']]]

def value(p,N,q,m):
    g=p['c']*q+p['d0']*(q**3-3*q);coef=p['zeta']*(-(q**3-3*q)+g/(2*p['G']))
    B0=[q*p['gamma']+coef*p['p0'],Q(0),q*(p['y']-p['E02']*p['r'])+coef*p['p2'],Q(0)]
    h=[B0[j]+sum(N[j][k]*m[k]for k in range(4))for j in range(4)]
    B=[h[0]-h[2],h[1]-3*h[3],h[2],h[3]];L=[Q(0),p['ell']-3*p['beta'],Q(0),p['beta']]
    quad=sum(m[j]*N[j][k]*m[k]for j in range(4)for k in range(4))
    return 2*max_integral([L,[-t for t in L],B,[-t for t in B]])-I(quad)+I(g*g/(4*p['G']))

def cover(p,qa,qb,reserve0,reserve6,maxnodes,output):
    p['major'][0]+=reserve0;p['major'][3]+=reserve6;N=matrix(p)
    C3=abs(p['zeta']*(-1+p['d0']/(2*p['G'])))*(abs(p['p0'])+abs(p['p2']))
    stack=[[(qa,qb),(Q(-1),Q(1)),(Q(0),Q(1)),(Q(-3,2),Q(3,2)),(Q(-5,2),Q(5,2))]]
    events=[];cache={};maj={};Ks={};nodes=accepted=excluded=0;start=time.time();largest=Q(-10**9)
    while stack and nodes<maxnodes:
        box=stack.pop();a,b=box[0];nodes+=1
        near=[l if l>0 else r if r<0 else Q(0)for l,r in box[1:]]
        if near[0]**2+near[1]**2+near[2]**2/2+near[3]**2/6>1:
            events.append('N');excluded+=1;continue
        if (a,b)not in Ks:Ks[a,b]=max(Q(0),-second_lower(p,a,b))
        center=tuple((l+r)/2 for l,r in box[1:]);radius=[(r-l)/2 for l,r in box[1:]]
        vals=[]
        for q in [a,b]:
            key=(q,*center)
            if key not in cache:
                if q not in maj:maj[q]=majorant(p,q)
                cache[key]=value(p,N,q,center)-I(maj[q])
            vals.append(cache[key])
        extras=[(Q(3,2)*b*C3+Ks[a,b]/8)*(b-a)**2]+[N[j][j]*radius[j]**2 for j in range(4)]
        cross=2*abs(N[0][2])*radius[0]*radius[2]+2*abs(N[1][3])*radius[1]*radius[3]
        U=I(max(v.hi for v in vals),max(v.hi for v in vals))+I(sum(extras)+cross)
        if U.below(0):
            events.append('C');accepted+=1;largest=max(largest,Q(U.hi.text(UP,40)))
        else:
            split=list(extras);split[1]+=abs(N[0][2])*radius[0]*radius[2];split[3]+=abs(N[0][2])*radius[0]*radius[2]
            split[2]+=abs(N[1][3])*radius[1]*radius[3];split[4]+=abs(N[1][3])*radius[1]*radius[3]
            axis=max(range(5),key=lambda j:split[j]);l,r=box[axis];mid=(l+r)/2;left=list(box);right=list(box);left[axis]=(l,mid);right[axis]=(mid,r);stack.extend([right,left]);events.append(str(axis))
        if nodes%2000==0:print(json.dumps(dict(nodes=nodes,accepted=accepted,excluded=excluded,pending=len(stack),centers=len(cache),qnodes=len(maj),seconds=time.time()-start)),flush=True)
    out=dict(status='PASS_GENERAL_PHYSICAL_R3_COMPACT_COVER'if not stack else'INCOMPLETE_NO_CERTIFICATE',qa=str(qa),qb=str(qb),reserve0=str(reserve0),reserve6=str(reserve6),nodes=nodes,accepted=accepted,excluded=excluded,pending=len(stack),largest_residual=str(largest),seconds=time.time()-start,
             profile={k:[str(t)for t in v]if isinstance(v,list)else str(v)for k,v in p.items()},events=events,pending_boxes=[[[str(a),str(b)]for a,b in bx]for bx in stack])
    with gzip.open(output,'wt')as f:json.dump(out,f)
    print(json.dumps({k:v for k,v in out.items()if k not in ['profile','events','pending_boxes']},indent=2),flush=True)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('profile');ap.add_argument('--qa',default='0');ap.add_argument('--qb',default='8');ap.add_argument('--reserve0',default='.001');ap.add_argument('--reserve6',default='.000001');ap.add_argument('--max-nodes',type=int,default=100000);ap.add_argument('--output',required=True);args=ap.parse_args();cover(freeze(args.profile),Q(args.qa),Q(args.qb),Q(args.reserve0),Q(args.reserve6),args.max_nodes,Path(args.output))
