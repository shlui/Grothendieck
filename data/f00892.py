#!/usr/bin/env python3
"""Independent reconstruction of physical R3 whole-q certificate trees.

Does not import the producer's local values, derivatives, Bernstein code,
tree code, or tail code. Gaussian absolute-cubic integrals use the separately
audited Arb root-isolation primitive. All frozen algebra is exact rational.
"""
from pathlib import Path
from fractions import Fraction as F
from functools import lru_cache
from math import comb
import argparse,gzip,hashlib,importlib.util,json,time,os,sys
from flint import arb,ctx
sys.dont_write_bytecode=True
BASELINE=Path(os.environ["GR_FROZEN_BASELINE"]).resolve()
sys.path.insert(0,str(BASELINE/"tonight/audit"))
from mpfr_interval import I,DOWN,UP

HERE=Path(__file__).resolve().parent
BASE=BASELINE
DATA=BASE/'tonight/operator/projection/lower_amp_certificate'
BACKEND=BASE/'tonight/recovery/profile_cover/directed_cubic_integral.py'
spec=importlib.util.spec_from_file_location('independent_audited_arb_cubic',BACKEND)
C=importlib.util.module_from_spec(spec);spec.loader.exec_module(C)
ctx.prec=384
TAIL_R=F(64)

def aa(x):
    x=F(x);return arb(x.numerator)/x.denominator

def plus(a,b):
    return [(a[i] if i<len(a) else F(0))+(b[i] if i<len(b) else F(0)) for i in range(max(len(a),len(b)))]
def times(a,b):
    z=[F(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):z[i+j]+=x*y
    return z
def scale(a,c):return [x*c for x in a]
def diff(a):return [i*a[i] for i in range(1,len(a))] or [F(0)]
def ev(a,q):
    z=F(0)
    for c in a[::-1]:z=z*q+c
    return z

@lru_cache(None)
def derivatives(coeff):
    # q/(1+q²) sum a_j T_j((q²-1)/(q²+1)).
    D=[F(1),F(0),F(1)];T=[F(-1),F(0),F(1)]
    powers=[[F(1)]]
    for i in range(len(coeff)):powers.append(times(powers[-1],D))
    cheb=[[F(1)],T]
    for i in range(2,len(coeff)):
        cheb.append(plus(scale(times(T,cheb[-1]),2),scale(times(powers[2],cheb[-2]),-1)))
    n=len(coeff);num=[F(0)]
    for j,a in enumerate(coeff):num=plus(num,scale(times(cheb[j],powers[n-j-1]),a))
    num=times([F(0),F(1)],num);out=[]
    for m in range(n,n+4):
        out.append((tuple(num),m))
        num=plus(times(diff(num),D),scale(times(num,[F(0),F(2)]),-m))
    return tuple(out)

@lru_cache(None)
def polybounds(num,a,b):
    # Expand p(a+(b-a)t), then convert monomials to Bernstein form.
    m=len(num)-1;w=b-a
    ap=[F(1)];wp=[F(1)]
    for i in range(m):ap.append(ap[-1]*a);wp.append(wp[-1]*w)
    power=[wp[i]*sum((num[j]*comb(j,i)*ap[j-i] for j in range(i,m+1)),F(0)) for i in range(m+1)]
    bs=[sum((power[i]*F(comb(k,i),comb(m,i)) for i in range(k+1)),F(0)) for k in range(m+1)]
    return min(bs),max(bs)

def phi_ranges(p,a,b):
    vals=[]
    for num,m in derivatives(tuple(p['phi'])):
        lo,hi=polybounds(num,a,b)
        vals.append(I(I(lo).lo,I(hi).hi)/I(I((1+a*a)**m).lo,I((1+b*b)**m).hi))
    return vals

def phi_point(p,q):
    return [ev(num,q)/(1+q*q)**m for num,m in derivatives(tuple(p['phi']))]

def lower_second(p,a,b):
    q=I(I(a).lo,I(b).hi);q2=q.square();q4=q2.square()
    _,fp,fpp,fppp=phi_ranges(p,a,b)
    theta=p['kappa']/(3*(p['kappa']+2));s=q2-1;low=-1-theta*s
    e=fp-p['nuisance0']-p['nuisance2']*s
    ep=fpp-2*p['nuisance2']*q;epp=fppp-2*p['nuisance2'];eta=p['eta']
    pin=(ep.square()+e*epp)/eta
    pout=-2*theta*e-4*theta*q*ep+low*epp-eta*(4*theta*theta*q2-2*theta*low)
    test=e-eta*low
    if I(0).lo<=test.lo:pen=pin
    elif test.hi<=I(0).hi:pen=pout
    else:pen=I(min(pin.lo,pout.lo),max(pin.hi,pout.hi))
    c=p['c']-3*p['d0'];g=c*q+p['d0']*q*q2
    gp=c+3*p['d0']*q2;gpp=6*p['d0']*q
    z=(gp.square()+g*gpp)/(2*p['G'])-2*p['major'][1]-12*p['major'][2]*q2-30*p['major'][3]*q4
    z+=-2*fp-q*fpp+2*theta*fp+4*theta*q*fpp+(1+theta*s)*fppp+pen
    return F(z.lo.text(DOWN,40))

def scalar_majorant(p,q):
    ph,fp,_,_=phi_point(p,q);theta=p['kappa']/(3*(p['kappa']+2));s=q*q-1
    low=-1-theta*s;e=fp-p['nuisance0']-p['nuisance2']*s
    pen=e*e/(2*p['eta']) if e>=p['eta']*low else low*e-p['eta']*low*low/2
    return sum((v*q**(2*i) for i,v in enumerate(p['major'])),F(0))+q*ph-(1+theta*s)*fp-pen

def physical_value(p,q,m):
    s,v1,v3=m;rho=1+p['zeta']**2/p['G'];g=p['c']*q+p['d0']*(q**3-3*q)
    coeff=p['zeta']*(-q**3+3*q+g/(2*p['G']))
    h0=q*p['gamma']+coeff*p['p0']+rho*p['p0']*s
    h2=q*(p['y']-p['p0']*p['p2']*p['r'])+coeff*p['p2']+rho*p['p2']*s
    h1=p['R11']*v1+p['R13']*v3;h3=p['R13']*v1+p['R33']*v3
    z=arb(0)
    for sign in (-1,1):
        z+=C.abs_cubic([aa(h0),aa(h1+sign*p['ell']),arb(2).sqrt()*aa(h2),arb(6).sqrt()*aa(h3+sign*p['beta'])])
    quad=rho*s*s+p['R11']*v1*v1+2*p['R13']*v1*v3+p['R33']*v3*v3
    return z+aa(g*g/(4*p['G'])-quad)

def freeze(which,branch):
    path=DATA/(f'profile_{which}.json'if which.startswith('interval_')else f'profile_u.{which}.json');d=json.loads(path.read_text())
    assert d['kind']=='EXACT_RATIONAL_R3_PROFILE'
    p={k:[F(t)for t in v] if isinstance(v,list) else F(v) for k,v in d['profile'].items()}
    assert p['eps']>0 and p['G']>0 and p['eta']>0 and p['kappa']>=0
    assert p['p2']**2>=p['eps']*(1+p['p0']**2)/2 and p['p0']**2+2*p['p2']**2<1
    for x,y,z in [(p['R11'],p['R13'],p['R33']),(p['R11']+1,p['R13'],p['R33']-p['eps']/6),(p['R11']-p['eps'],p['R13'],p['R33'])]:
        assert x>=0 and z>=0 and x*z>=y*y
    reserves={'tight':(F(1,20000),F(1,1000000)),'finer':(F(1,50000),F(1,5000000)),'generous':(F(1,1000),F(1,1000000)),'middle':(F(1,10000),F(1,1000000)),'robust':(F(1,200000),F(1,100000000)),'robustcurrent':(F(1,50000),F(1,10000000)),'highk':(F(1,1000),F(1,10000000)),'lowerface':(F(1,10000),F(1,100000000)),'nextface':(F(1,50000),F(1,100000000)),'q50':(F(1,20000),F(1,100000000)),'q25':(F(1,40000),F(1,100000000))}
    q0,q6=reserves[branch];p['major'][0]+=q0;p['major'][3]+=q6
    return p,path,q0,q6

def replay(which,branch,part):
    p,path,q0,q6=freeze(which,branch)
    tree=DATA/(f'{part}_{which}.json.gz'if which.startswith('interval_')else f'{part}_u{which}_{branch}.json.gz')
    with gzip.open(tree,'rt')as f:r=json.load(f)
    assert r['status']=='PASS_PHYSICAL_BERNSTEIN_LOCAL_SUPREMUM' and r['pending']==0
    assert F(r['reserve'])==q6 and F(r['reserve0'])==q0
    assert {k:[str(t)for t in v]if isinstance(v,list)else str(v)for k,v in p.items()}==r['profile']
    qa,qb=F(r['qa']),F(r['qb']);assert (qa,qb)==((F(0),F(8))if part=='body'else(F(8),TAIL_R))
    stack=[[(qa,qb),(F(-1),F(1)),(F(0),F(1)),(F(-5,2),F(5,2))]]
    rho=1+p['zeta']**2/p['G'];pnorm=p['p0']**2+2*p['p2']**2
    C3=abs(p['zeta']*(-1+p['d0']/(2*p['G'])))*(abs(p['p0'])+abs(p['p2']))
    derivatives_cache={};points={};accepted=excluded=splits=0;maxgap=None;start=time.time()
    for index,event in enumerate(r['events']):
        assert stack;box=stack.pop();a,b=box[0]
        if event=='N':
            near=[l if l>0 else h if h<0 else F(0)for l,h in box[1:]]
            assert near[0]**2/pnorm+near[1]**2+near[2]**2/6>1
            excluded+=1
        elif event=='C':
            center=tuple((l+h)/2 for l,h in box[1:]);vals=[]
            for q in(a,b):
                key=(q,*center)
                if key not in points:points[key]=physical_value(p,q,center)-aa(scalar_majorant(p,q))
                vals.append(points[key])
            if(a,b)not in derivatives_cache:derivatives_cache[a,b]=max(F(0),-lower_second(p,a,b))
            K=derivatives_cache[a,b]
            radius=(F(3,2)*b*C3+K/8)*(b-a)**2
            radius+=sum((z*(h-l)**2/4 for z,(l,h)in zip((rho,p['R11'],p['R33']),box[1:])),F(0))
            radius+=abs(p['R13'])*(box[2][1]-box[2][0])*(box[3][1]-box[3][0])/2
            upper=max(F(str(v.upper().fmpq()))for v in vals)+radius
            assert upper<0,('FAILED_PHYSICAL_Q_BOX',which,branch,part,index,box,str(upper))
            maxgap=upper if maxgap is None else max(maxgap,upper);accepted+=1
        else:
            axis=int(event);assert 0<=axis<4
            l,h=box[axis];mid=(l+h)/2;left,right=list(box),list(box)
            left[axis]=(l,mid);right[axis]=(mid,h);stack.extend((right,left));splits+=1
        if index and index%2000==0:print(json.dumps({'profile':which,'branch':branch,'part':part,'events':index,'accepted':accepted,'seconds':time.time()-start}),flush=True)
    assert not stack and accepted==r['accepted'] and excluded==r['excluded']
    assert accepted+excluded+splits==r['nodes']==len(r['events'])
    return {'status':'PASS_INDEPENDENT_PHYSICAL_Q_TREE','profile':which,'branch':branch,'part':part,'nodes':r['nodes'],'accepted':accepted,'excluded':excluded,'distinct_points':len(points),'q_intervals':len(derivatives_cache),'largest_directed_gap':str(maxgap),'seconds':time.time()-start,'tree_sha256':hashlib.sha256(tree.read_bytes()).hexdigest(),'profile_sha256':hashlib.sha256(path.read_bytes()).hexdigest()}

def replay_tail(which,branch):
    p,path,q0,q6=freeze(which,branch);R=TAIL_R
    Aphi=sum(abs(v)for v in p['phi']);Dphi=sum(j*j*abs(v)for j,v in enumerate(p['phi']))
    theta=p['kappa']/(3*(p['kappa']+2))
    dphi=Aphi/R**2+4*Dphi/R**4
    stein=theta*Aphi+Aphi/R**2+4*theta*Dphi/R**2+4*Dphi/R**4
    e0=abs(p['nuisance0']-p['nuisance2'])+dphi;e2=abs(p['nuisance2']);eta=p['eta']
    cc=p['c']-3*p['d0'];dd=p['d0'];rho=1+p['zeta']**2/p['G']
    h=p['y']-p['p0']*p['p2']*p['r']
    l1=p['gamma']+p['zeta']*(3+cc/(2*p['G']))*p['p0'];l2=h+p['zeta']*(3+cc/(2*p['G']))*p['p2']
    a3=abs(p['zeta']*(-1+dd/(2*p['G'])))*(abs(p['p0'])+abs(p['p2']))
    const=2*(abs(p['ell'])+F(5,2)*abs(p['beta']))+rho*(p['p0']**2+2*p['p2']**2)+p['R11']+6*p['R33']
    pol=[F(0)]*7
    for j,v in enumerate(p['major']):pol[2*j]=v
    pol[0]-=const+Aphi+stein+e0*e0/(2*eta)
    pol[1]-=2*(abs(l1)+abs(l2));pol[2]-=e0*e2/eta+cc*cc/(4*p['G'])
    pol[3]-=2*a3;pol[4]-=e2*e2/(2*eta)+abs(cc*dd)/(2*p['G']);pol[6]-=dd*dd/(4*p['G'])
    shifted=[sum((pol[j]*comb(j,k)*R**(j-k)for j in range(k,7)),F(0))for k in range(7)]
    assert all(v>=0 for v in shifted)
    if which.startswith('interval_'):
        optional=DATA/f'tail_{which}.json'
        record=json.loads(optional.read_text())if optional.exists()else None
    elif branch=='finer':
        ledger=json.loads((DATA/'lower_source_domain_finer.json').read_text())
        record=next(cell['tail'] for cell in ledger['cells'] if cell['profile']==which)
    else:record=json.loads((DATA/f'tail_u{which}_{branch}.json').read_text())
    if record is not None:
        assert record['status']=='PASS' and F(record['R'])==R
        assert list(map(F,record['difference_polynomial']))==pol and list(map(F,record['shifted_coefficients']))==shifted
    return {'status':'PASS_INDEPENDENT_EXACT_R3_TAIL','profile':which,'branch':branch,'R':str(R),'difference_polynomial':list(map(str,pol)),'shifted_coefficients':list(map(str,shifted))}

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--profile',required=True);ap.add_argument('--branch',choices=['tight','finer','generous','middle','robust','robustcurrent','highk','lowerface','nextface','q50','q25'],default='finer');ap.add_argument('--part',choices=['body','exterior','tail'],required=True);ap.add_argument('--data-dir',type=Path);ap.add_argument('--output-dir',type=Path);ap.add_argument('--tail-R',default='64');args=ap.parse_args();TAIL_R=F(args.tail_R);assert TAIL_R>=1
    if args.data_dir:DATA=args.data_dir.resolve()
    result=replay_tail(args.profile,args.branch)if args.part=='tail'else replay(args.profile,args.branch,args.part)
    result.update({'precision_bits':ctx.prec,'cubic_integral_backend_sha256':hashlib.sha256(BACKEND.read_bytes()).hexdigest(),'scope':'Frozen R3 physical moment and complete q interval only. Outer source parameters require a separate convex price proof.'})
    out=args.output_dir or HERE/'lower_amp_R3';out.mkdir(parents=True,exist_ok=True)
    (out/f'independent_{args.profile}_{args.branch}_{args.part}.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True)
