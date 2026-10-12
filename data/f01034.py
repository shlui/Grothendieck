"""Directed q-and-moment cover for one frozen rank-one R3 profile.

Acceptance uses exact rational algebra and outward MPFR enclosures.  A finite
compact cover is not an operator certificate until the q tail and all outer
source parameters are also covered.
"""
from fractions import Fraction as Q
from pathlib import Path
from itertools import product
import argparse,json,time,gzip
from mpfr_interval import I,UP,DOWN,BITS
from certify_local_point import freeze,context,value,majorant,ROOT
from rational_phi_bounds import phi_bounds

class Jet:
    def __init__(self,v):self.v=[I(t) for t in v]+[I(0)]*(4-len(v))
    @staticmethod
    def cast(x):return x if isinstance(x,Jet) else Jet([x])
    def __add__(self,x):x=self.cast(x);return Jet([a+b for a,b in zip(self.v,x.v)])
    __radd__=__add__
    def __neg__(self):return Jet([-t for t in self.v])
    def __sub__(self,x):return self+-self.cast(x)
    def __rsub__(self,x):return self.cast(x)+-self
    def __mul__(self,x):
        x=self.cast(x);return Jet([sum((self.v[j]*x.v[k-j] for j in range(k+1)),I(0)) for k in range(4)])
    __rmul__=__mul__
    def inv(self):
        out=[self.v[0].inv()]
        for k in range(1,4):out.append(-out[0]*sum((self.v[j]*out[k-j] for j in range(1,k+1)),I(0)))
        return Jet(out)
    def __truediv__(self,x):return self*self.cast(x).inv()
    def __rtruediv__(self,x):return self.cast(x)*self.inv()

def phi_jet(p,a,b):
    q=Jet([I(I(a).lo,I(b).hi),1]);q2=q*q;t=(q2-1)/(q2+1);base=q/(q2+1)
    Ts=[Jet([1])]
    if len(p['phi'])>1:Ts.append(t)
    for j in range(2,len(p['phi'])):Ts.append(2*t*Ts[-1]-Ts[-2])
    ph=base*sum((c*t for c,t in zip(p['phi'],Ts)),Jet([0]))
    return ph.v[0],ph.v[1],2*ph.v[2],6*ph.v[3]

def second_lower(p,a,b):
    q=I(I(a).lo,I(b).hi);q2=q.square();q4=q2.square();_,fp,fpp,fppp=phi_bounds(p,a,b)
    theta=p['kappa']/(3*(p['kappa']+2));hh=q2-1;lo=-1-theta*hh
    res=fp-p['nuisance0']-p['nuisance2']*hh
    rp=fpp-2*p['nuisance2']*q;rpp=fppp-2*p['nuisance2'];eta=p['eta']
    inside=(rp.square()+res*rpp)/eta
    outside=-2*theta*res-4*theta*q*rp+lo*rpp-eta*(4*theta*theta*q2-2*theta*lo)
    gate=res-eta*lo
    if I(0).lo<=gate.lo:pen=inside
    elif gate.hi<=I(0).hi:pen=outside
    else:pen=I(min(inside.lo,outside.lo),max(inside.hi,outside.hi))
    cc=p['c']-3*p['d0'];g=cc*q+p['d0']*q*q2;gp=cc+3*p['d0']*q2;gpp=6*p['d0']*q
    poly2=2*p['major'][1]+12*p['major'][2]*q2+30*p['major'][3]*q4
    sec=(gp.square()+g*gpp)/(2*p['G'])-poly2-2*fp-q*fpp+2*theta*fp+4*theta*q*fpp+(1+theta*hh)*fppp+pen
    return Q(sec.lo.text(DOWN,40))

def tail(p,R):
    assert R>=1
    Ap=sum(abs(c) for c in p['phi']);Ad=sum(j*j*abs(c) for j,c in enumerate(p['phi']))
    theta=p['kappa']/(3*(p['kappa']+2));dp=Ap/R**2+4*Ad/R**4
    Dtheta=theta*Ap+Ap/R**2+4*theta*Ad/R**2+4*Ad/R**4
    A=abs(p['nuisance0']-p['nuisance2'])+dp;B=abs(p['nuisance2']);eta=p['eta']
    cc=p['c']-3*p['d0'];d=p['d0'];rho=1+p['zeta']**2/p['G']
    h=p['y']-p['p0']*p['p2']*p['r'];l1=p['gamma']+p['zeta']*(3+cc/(2*p['G']))*p['p0'];l2=h+p['zeta']*(3+cc/(2*p['G']))*p['p2']
    a3=abs(p['zeta']*(-1+d/(2*p['G'])))*(abs(p['p0'])+abs(p['p2']))
    Ctrace=rho*(p['p0']**2+2*p['p2']**2)+p['R11']+6*p['R33']
    # 2 E|L| <= 2(|ell|+sqrt(6)|beta|), sqrt(6)<5/2.
    constant=2*(abs(p['ell'])+Q(5,2)*abs(p['beta']))+Ctrace
    pol=[Q(0)]*7
    for j,c in enumerate(p['major']):pol[2*j]=c
    pol[0]-=constant+Ap+Dtheta+A*A/(2*eta)
    pol[1]-=2*(abs(l1)+abs(l2));pol[2]-=A*B/eta+cc*cc/(4*p['G'])
    pol[3]-=2*a3;pol[4]-=B*B/(2*eta)+abs(cc*d)/(2*p['G']);pol[6]-=d*d/(4*p['G'])
    # Exact nonnegative coefficients after q=R+t prove every q>=R.
    import math
    shifted=[sum(pol[j]*math.comb(j,k)*R**(j-k) for j in range(k,7)) for k in range(7)]
    return dict(status='PASS' if all(t>=0 for t in shifted) else 'FAIL',R=str(R),difference_polynomial=[str(t) for t in pol],shifted_coefficients=[str(t) for t in shifted])

def cover(p,qa,qb,capreserve,maxnodes,output,resume=None,reserve0=None):
    reserve0=capreserve if reserve0 is None else reserve0
    p['major'][0]+=reserve0;p['major'][3]+=capreserve
    rho=1+p['zeta']**2/p['G'];C3=abs(p['zeta']*(-1+p['d0']/(2*p['G'])))*(abs(p['p0'])+abs(p['p2']))
    diag=[rho,p['R11'],p['R33']]
    cache={};contexts={};majorcache={};Kcache={};events=[];count=accepted=0;start=time.time();largest=Q(-10**9)
    stack=[[(qa,qb),(Q(-1),Q(1)),(Q(0),Q(1)),(Q(-5,2),Q(5,2))]]
    excluded=0
    pnorm=p['p0']**2+2*p['p2']**2
    oldseconds=0
    if resume:
        with gzip.open(resume,'rt') as f:old=json.load(f)
        assert old['qa']==str(qa) and old['qb']==str(qb) and old['reserve']==str(capreserve)
        expected={k:[str(t) for t in v] if isinstance(v,list) else str(v) for k,v in p.items()}
        assert expected==old['profile'], 'Changed mathematical profile cannot resume.'
        stack=[[(Q(a),Q(b)) for a,b in box] for box in old['pending_boxes']]
        events=old['events'];count=old['nodes'];accepted=old['accepted'];excluded=old.get('excluded',0);largest=Q(old['largest_residual']);oldseconds=old['seconds']
    stopcount=count+maxnodes
    while stack and count<stopcount:
        box=stack.pop();a,b=box[0];count+=1
        nearest=[l if l>0 else r if r<0 else Q(0) for l,r in box[1:]]
        if nearest[0]**2/pnorm+nearest[1]**2+nearest[2]**2/6>1:
            excluded+=1;events.append('N');continue
        if (a,b) not in Kcache:Kcache[a,b]=max(Q(0),-second_lower(p,a,b))
        K=Kcache[a,b];vals=[]
        center=[(l+r)/2 for l,r in box[1:]]
        for q in (a,b):
            key=(q,*center)
            if key not in cache:
                if q not in contexts:contexts[q]=context(p,q);majorcache[q]=majorant(p,q)
                cache[key]=value(p,contexts[q],tuple(center))-I(majorcache[q])
            vals.append(cache[key])
        Cq=Q(3,2)*b*C3+K/8
        extras=[Cq*(b-a)**2]+[d*(r-l)**2/4 for d,(l,r) in zip(diag,box[1:])]
        radius_cross=abs(p['R13'])*(box[2][1]-box[2][0])*(box[3][1]-box[3][0])/2
        upper=I(max(v.hi for v in vals),max(v.hi for v in vals))+I(sum(extras)+radius_cross)
        if upper.below(0):
            accepted+=1;U=Q(upper.hi.text(UP,35));largest=max(largest,U);events.append('C')
        else:
            axis=max(range(4),key=lambda j:extras[j]);l,r=box[axis];mid=(l+r)/2;left=list(box);right=list(box);left[axis]=(l,mid);right[axis]=(mid,r);stack.extend([right,left]);events.append(str(axis))
        if count%1000==0:print(json.dumps(dict(nodes=count,accepted=accepted,pending=len(stack),corners=len(cache),q_nodes=len(contexts),q_intervals=len(Kcache),seconds=time.time()-start)),flush=True)
    out=dict(status='PASS_PHYSICAL_BERNSTEIN_LOCAL_SUPREMUM' if not stack else 'INCOMPLETE_NO_CERTIFICATE',scope='One frozen source point, compact q interval; actual physical momentboxes boundedbyquadraticcentercompletion, L2 ball and odd-vector reflection. No operator or outer-source certificate.',qa=str(qa),qb=str(qb),reserve=str(capreserve),reserve0=str(reserve0),nodes=count,accepted=accepted,pending=len(stack),excluded=excluded,corners=len(cache),q_nodes=len(contexts),seconds=oldseconds+time.time()-start,largest_residual=str(largest),derivative_backend='EXACT_RATIONAL_BERNSTEIN',bits=BITS,profile={k:[str(t) for t in v] if isinstance(v,list) else str(v) for k,v in p.items()},events=events,pending_boxes=[[[str(a),str(b)] for a,b in box] for box in stack])
    with gzip.open(output,'wt') as f:json.dump(out,f)
    print(json.dumps({k:v for k,v in out.items() if k not in ['profile','events','pending_boxes']},indent=2),flush=True)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('profile');ap.add_argument('--qa',default='0');ap.add_argument('--qb',default='1');ap.add_argument('--reserve',default='0.000001');ap.add_argument('--reserve0');ap.add_argument('--max-nodes',type=int,default=100000);ap.add_argument('--tail',default='');ap.add_argument('--output');ap.add_argument('--resume');args=ap.parse_args();p=freeze(args.profile)
    if args.tail:
        p['major'][0]+=Q(args.reserve0 or args.reserve);p['major'][3]+=Q(args.reserve);out=tail(p,Q(args.tail));print(json.dumps(out,indent=2));Path(args.output or ROOT/'results'/'directed_local_tail.json').write_text(json.dumps(out,indent=2))
    else:cover(p,Q(args.qa),Q(args.qb),Q(args.reserve),args.max_nodes,Path(args.output or ROOT/'results'/'directed_local_region.json.gz'),args.resume,Q(args.reserve0) if args.reserve0 else None)
