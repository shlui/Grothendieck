#!/usr/bin/env python3
"""Independent selected-endpoint and complete native/R3 domain replay.

Imports no native producer/pricer/union module. Uses directed MPFR for
Gaussian target means and exact Fraction quadratic/chord/partition algebra.
"""
from pathlib import Path
from fractions import Fraction as F
from functools import lru_cache
import argparse,bisect,hashlib,json,time
from mpfr_interval import I,UP,BITS

HERE=Path(__file__).resolve().parent;BASE=HERE.parents[1]
COVER=BASE/'tonight/recovery/profile_cover'
L1FILE=BASE/'tonight/operator/cubic_L1_refinement/refit847/complete_price_cover.json'
L1=json.loads(L1FILE.read_text());assert L1['status']=='PASS_COMPLETE_CUBIC_L1_REFIT_PRICE_COVER'
PI=I.pi();NU=(I(2)/PI).sqrt();K=I(194)/I(10).sqrt()
PACKAGES={};FRAMES={};POINTS=0

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def f(x):return F(str(x))
def high(x):return F(x.hi.text(UP,235))
def decimal_ceiling(x,digits):
    x=F(x);s=10**digits;return F(-((-x.numerator*s)//x.denominator),s)
def qmax(c,a,g,lo,hi):
    values=[c+a*lo+g*lo*lo,c+a*hi+g*hi*hi]
    if g<0:
        x=-a/(2*g)
        if lo<x<hi:values.append(c+a*x+g*x*x)
    return max(values)
def psd(a,o,c):assert a>=0 and c>=0 and a*c>=o*o

def package(name):
    if name in PACKAGES:return PACKAGES[name]
    p=BASE/name;manifest=json.loads((p/'SHA256.json').read_text())
    for rel,digest in manifest.items():assert sha(p/rel)==digest,('PACKAGE_HASH',name,rel)
    source=json.loads((p/'data/source_profile.json').read_text());local=json.loads((p/'data/local_summary.json').read_text());replay=json.loads((p/'replay.json').read_text())
    assert replay['status'] in ['PASS_COMPLETE_LOCAL_AND_SCALAR_PRIMITIVES','PASS_COMPLETE_REGIONAL_BOUND']
    assert local['status']=='PASS' and len(local['nodes'])==51 and all(n['status']=='PASS'and n['pending_boxes']==0 for n in local['nodes'])
    digest=sha(p/'data/source_profile.json');assert local['source_sha256']==digest
    if 'source_sha256'in replay:assert replay['source_sha256']==digest
    profile=local['profile'];M=[[F(v)for v in row]for row in profile['M']]
    assert all(M[i][j]==M[j][i]for i in range(4)for j in range(4))
    assert all(M[i][j]==0 for i in range(4)for j in range(4)if(i-j)%2)
    b=F(profile['epsilon']);assert b==F(18,25)
    for a,o,c in [(M[0][0],M[0][2],M[2][2]),(M[1][1],M[1][3],M[3][3]),(M[0][0]+1,M[0][2],M[2][2]-b),(M[1][1]+1,M[1][3],M[3][3]-b),(M[1][1]-b,M[1][3],M[3][3])]:psd(a,o,c)
    t=list(map(f,source['t']));h=list(map(f,source['profile']['nodal_majorant']));rg=list(map(f,source['r_grid']))
    assert t[0]==rg[0]==0 and len(t)==len(h)==51 and all(a<c for a,c in zip(t,t[1:]))and all(a<c for a,c in zip(rg,rg[1:]))
    slopes=[(c-a)/(y-x)for x,y,a,c in zip(t,t[1:],h,h[1:])]
    tau=f(source['profile']['tail_slope']);assert slopes[0]>=0 and all(a<=c for a,c in zip(slopes,slopes[1:]))and tau>=slopes[-1]
    inf=F(profile['nodal_majorant'][0])-h[0];tail=F(profile['tail_slope'])-tau;assert inf>=0 and tail>=0
    assert all(F(a)-c==inf for a,c in zip(profile['nodal_majorant'],h))
    result={'path':p,'source':source,'local':local,'profile':profile,'M':M,'rg':rg,'t':t,'h':h,'slopes':slopes,'tau':tau,'inflation':inf,'tail':tail,'source_sha256':digest,'manifest_sha256':sha(p/'SHA256.json'),'replay_sha256':sha(p/'replay.json')}
    PACKAGES[name]=result;return result

def sameframe(a,b):
    assert a['M']==b['M'] and a['t']==b['t']
    for k in['ell','beta','delta','nu','epsilon','tail_slope']:assert F(a['profile'][k])==F(b['profile'][k])
    assert list(map(F,a['profile']['nodal_majorant']))==list(map(F,b['profile']['nodal_majorant']))

@lru_cache(None)
def target_mean(name,r):
    d=package(name);s=d['slopes'];weights=[r*s[0]]+[r*(b-a)for a,b in zip(s,s[1:])]+[r*(d['tau']-s[-1])]
    assert len(weights)==len(d['t'])and all(w>=0 for w in weights)
    mean=I(d['h'][0])
    for knot,w in zip(d['t'],weights):
        x=I(knot/r);excess=NU*(-x.square()/2).exp()-x*(x/I(2).sqrt()).erfc()
        mean+=w*excess
    return mean

@lru_cache(None)
def scalar(name,index):
    d=package(name);path=d['path']/f'data/scalars/row{index:02d}.json';cert=json.loads(path.read_text());rr=d['source']['rows'][index];fr=cert['frozen_row']
    assert cert['status']=='PASS_FRESH_COMPLETE_SCALAR_HALF_LINE'and cert['source_sha256']==d['source_sha256']and cert['row_index']==index
    q0,q6=F(fr['added_q0']),F(fr['added_q6']);assert q0>=0 and q6>=0
    q=list(map(F,fr['majorant']));assert q==[f(x)+(q0 if j==0 else q6 if j==3 else 0)for j,x in enumerate(rr['majorant'])]
    for key,raw in[('r','r'),('lambda_L1','lambda_L1'),('eta','eta')]:assert F(fr[key])==f(rr[raw])
    assert list(map(F,fr['phi_coefficients']))==list(map(f,rr['phi_coefficients']))
    assert F(cert['tail']['domain'][0])>=max(d['t'][-1]/F(fr['r']),F(12))
    assert F(cert['tail']['domain'][0])==F(cert['compact']['domain'][1])
    assert q[3]>=0 and F(fr['lambda_L1'])>=0 and F(fr['eta'])>0
    return {'r':F(fr['r']),'mu':F(fr['lambda_L1']),'eta':F(fr['eta']),'q':q,'certificate_sha256':sha(path)}

@lru_cache(None)
def pure(name):
    record=L1['rows'][name];path=L1FILE.parent/record['path'];d=json.loads(path.read_text())
    assert sha(path)==record['sha256']and d['status'].startswith('PASS_DIRECTED_')
    assert list(map(F,d['majorant']))==list(map(F,record['majorant']))and F(d['eta'])==F(record['eta'])
    assert list(map(F,d.get('odd_weights',['1','0','0'])))==[F(1),F(0),F(0)]and F(d.get('sixth_offset','0'))==0
    q=list(map(F,record['majorant']));eta=F(record['eta']);assert q[3]>=0 and eta>0
    return q,eta

def price(name,index,k,s,rho,l1):
    global POINTS
    rr=scalar(name,index);q=rr['q'][:];eta=rr['eta'];mu=rr['mu']
    if l1!='constant':
        qq,ee=pure(l1);q=[a+mu*c for a,c in zip(q,qq)];eta+=mu*ee
    value=target_mean(name,rr['r'])+q[0]+q[1]+q[2]*(3+k)
    value+=I(q[3])*(15+15*k+(K-F(633,15)*s)*I(k)*I(k).sqrt())
    value+=I(eta*k*(F(7,90)+s/30))
    if l1=='constant':value+=I(mu)*(I(rho)-NU)
    POINTS+=1;return value

def frame(binding,b):
    name=binding['mask']
    if name in FRAMES:return FRAMES[name]
    path=BASE/name;d=json.loads(path.read_text());assert sha(path)==binding['mask_sha256']
    assert d['status']=='PASS_DIRECTED_NATIVE_CELL_PRICING_BASE_LEAF_CERTIFICATE_REQUIRED'and F(d['operator_b'])==b
    primary=package(binding['package']);assert d['source_sha256']==primary['source_sha256']==binding['source_sha256']
    origins={binding['package']}
    for alt in d['alternative_scalar_bindings']:
        secondary=package(alt['package']);assert secondary['source_sha256']==alt['source_sha256'];sameframe(primary,secondary);origins.add(alt['package'])
    M=primary['M'];p=primary['profile'];old=F(p['epsilon']);assert 0<b<=old
    g=(F(p['nu'])-M[0][2])**2/M[2][2]-M[0][0]+b-2*F(p['delta']);shift=F(p['beta'])**2*(1/b-1/old)
    assert F(d['exact_gram_cost'])==g and F(d['local_constant_transport'])==shift and F(d['ell'])==F(p['ell'])
    if d.get('l1_rows_file_sha256'):assert d['l1_rows_file_sha256']==sha(L1FILE)
    records={};keylist=[]
    for row in d['endpoint_prices']:
        lo,hi=map(F,row['kappa_interval']);assert lo<hi;keylist.append((lo,hi));rho=F(row['rho']);assert 0<rho<=F(847,1000)
        endpoints=list(map(lambda rr:list(map(F,rr)),row['prices']));selected=row['selected_rows'];rg=primary['rg']
        assert len(endpoints)==len(selected)==len(rg)
        assert selected[0]=='r=0'and endpoints[0]==[F(p['nodal_majorant'][0])+shift]*4
        for j,choice in enumerate(selected[1:],1):
            origin=choice['source_package'];assert origin in origins
            rr=scalar(origin,choice['native_row']);assert rr['r']==rg[j]
            local=package(origin)
            for index,(k,s)in enumerate(( (k,s)for k in(lo,hi)for s in(F(0),F(1)) )):
                value=price(origin,choice['native_row'],k,s,rho,choice['l1_row'])+local['inflation']+local['tail']*F(849,1000)*rg[j]+shift
                assert value.below(endpoints[j][index]),('ENDPOINT_REPRICE',name,j,k,s,value.record())
        records[lo,hi]=endpoints
    keylist.sort();assert keylist[0][0]==0 and keylist[-1][1]==90 and all(a[1]==c[0]for a,c in zip(keylist,keylist[1:]))
    result={'name':name,'rg':primary['rg'],'ell':F(p['ell']),'g':g,'endpoints':records,'keys':keylist,'mask_sha256':sha(path),'b':b}
    FRAMES[name]=result;return result

def frame_bound(d,ulo,uhi,key,rlo,rhi):
    rg=d['rg'];assert rhi<=rg[-1];j=bisect.bisect_right(rg,rlo)-1
    if j==len(rg)-1:j-=1
    assert rg[j]<=rlo<rhi<=rg[j+1]
    ep=d['endpoints'][key];outer=max(ulo*ulo-2*d['ell']*ulo,uhi*uhi-2*d['ell']*uhi);values=[]
    for corner in range(4):
        a=(ep[j+1][corner]-ep[j][corner])/(rg[j+1]-rg[j]);c=ep[j][corner]-a*rg[j]
        values.append(qmax(c,a,d['g'],rlo,rhi)+outer)
    return decimal_ceiling(max(values),12)

def source_upper(u,rho):
    z=(I(u)/NU-1+rho*rho)/(rho*rho)
    if z.hi<=I(0).lo:return rho
    assert z.above(0)and z.below(1)
    return decimal_ceiling(high(I(rho)*(-z.log()).sqrt().erf()),15)

def replay(path,coverage,r3=None,envelope=False):
    n=json.loads(path.read_text());assert n['status']=='PASS_EXACT_CONTINUOUS_NATIVE_FRAME_UNION_PRICING'
    b=F(n['operator_b']);assert b==F(7171,10000);threshold=F(n['threshold']);assert threshold==F(40,71)or threshold==F(451,800)
    rho=F(n['source_rho']);assert rho==F(847,1000)
    ds={d['name']:d for d in(frame(binding,b)for binding in n['bindings'])};keys=next(iter(ds.values()))['keys'];assert all(d['keys']==keys for d in ds.values())
    up=list(map(F,n['u_grid']));assert all(a<c for a,c in zip(up,up[1:]));lo,hi=coverage;assert lo in up and hi in up
    summaries={(tuple(map(F,s['u_interval'])),tuple(map(F,s['kappa_interval']))):s for s in n['summary']};assert len(summaries)==len(n['summary'])
    cells={};native=other=0;maximum=F(0);groups=0;missing=0
    for c in n['cells']:
        key=(tuple(map(F,c['u_interval'])),tuple(map(F,c['kappa_interval'])),tuple(map(F,c['r_interval'])))
        assert key not in cells;cells[key]=c
    for ulo,uhi in zip(up,up[1:]):
        if not(lo<=ulo<uhi<=hi):continue
        upper=source_upper(ulo,rho)
        for band in keys:
            summary=summaries[((ulo,uhi),band)];assert F(summary['source_upper'])==upper;groups+=1
            group=sorted([(r,c)for (u,k,r),c in cells.items()if u==(ulo,uhi)and k==band],key=lambda x:x[0][0])
            assert group and group[0][0][0]==0 and group[-1][0][1]==upper
            assert all(a[0][1]==c[0][0]for a,c in zip(group,group[1:]));closed=True
            for(rlo,rhi),cell in group:
                assert 0<=rlo<rhi<=upper
                if cell['frame']is not None:
                    assert cell['frame']in ds;cap=frame_bound(ds[cell['frame']],ulo,uhi,band,rlo,rhi);assert F(cell['cap'])==cap
                    assert bool(cell['accepted'])==(cap<threshold)
                else:assert cell['cap']is None and not cell['accepted']
                if cell['accepted']:native+=1;maximum=max(maximum,cap)
                elif envelope:
                    closed=False
                    if cell['frame']is not None:native+=1;maximum=max(maximum,cap)
                    else:missing+=1
                else:
                    closed=False;assert r3 is not None,('UNCLOSED_NATIVE_CELL',path,cell)
                    assert F(r3['operator_b'])==b and F(r3['maximum_upper'])<F(r3['rounded_exact_cap'])<threshold
                    for arr,bounds in[((ulo,uhi),r3['u']),((rlo,rhi),r3['r']),(band,r3['kappa'])]:
                        x,y=map(F,bounds);assert x<=arr[0]<arr[1]<=y
                    assert list(map(F,r3['shape']))==[F(2,5),F(1)]
                    other+=1;maximum=max(maximum,F(r3['rounded_exact_cap']))
            assert summary['whole_source_closed']==closed
    return {'status':'PASS_INDEPENDENT_NATIVE_ENVELOPE_INPUT_REPLAY'if envelope else'PASS_INDEPENDENT_COMPLETE_SAME_B_AMPLITUDE_UNION','path':str(path.relative_to(BASE)),'input_sha256':sha(path),'b':str(b),'interval':list(map(str,coverage)),'threshold':str(threshold),'native_cells':native,'R3_cells':other,'missing_native_price_cells':missing,'source_band_groups':groups,'maximum_rounded_cap':str(maximum),'source_rho':str(rho)}

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--ledger',type=Path,required=True);ap.add_argument('--u',type=F,nargs=2,required=True);ap.add_argument('--r3',type=Path);ap.add_argument('--envelope',action='store_true');ap.add_argument('--output',type=Path,required=True);args=ap.parse_args();start=time.time()
    r3=json.loads(args.r3.read_text())if args.r3 else None
    result=replay(args.ledger.resolve(),args.u,r3,args.envelope)
    result.update({'precision_bits':BITS,'selected_endpoint_prices_recomputed':POINTS,'frames_repriced':len(FRAMES),'accepted_package_bindings':[{'package':name,'source_sha256':d['source_sha256'],'manifest_sha256':d['manifest_sha256'],'replay_sha256':d['replay_sha256']}for name,d in PACKAGES.items()],'producer_modules_imported':False,'seconds':time.time()-start,'scope':'Complete regional all-source amplitude theorem for ONE b. No unrestricted norm or Grothendieck endpoint.'})
    args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
