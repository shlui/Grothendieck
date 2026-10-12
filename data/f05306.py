"""Replay a packaged regional theorem with no optimizer or discovery imports.

Dependencies: Python 3.12+ and python-flint. All acceptance uses exact
Fractions or directed Arb arithmetic. Stored partitions merely nominate
finite leaves; their inequalities and exhaustive coverage are rechecked.
"""
from pathlib import Path
from fractions import Fraction as F
import json,gzip,hashlib,argparse,time,sys
from flint import arb,ctx
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
import scalar_halfline as S
from directed_cubic_integral import local_point,exact

def qa(x):
    x=F(str(x));return arb(x.numerator)/x.denominator

def excluded(box):
    return sum(a*a if a>0 else b*b if b<0 else F(0) for a,b in box)>1

def local_replay(root,data,summary):
    assert summary['status']=='PASS' and len(summary['nodes'])==len(data['t'])
    p=summary['profile'];M=[[F(z) for z in row]for row in p['M']]
    raw=data['profile'];oldM=[[F(str(z))for z in row]for row in raw['M']]
    repair=M[0][0]-oldM[0][0]
    assert repair>=0
    assert all(M[i][j]-oldM[i][j]==(repair if i==j else 0)for i in range(4)for j in range(4))
    for key in ['ell','beta','delta','nu','epsilon']:assert F(p[key])==F(str(raw[key]))
    inflation=F(p['nodal_majorant'][0])-F(str(raw['nodal_majorant'][0]))
    assert inflation>=0 and all(F(a)-F(str(b))==inflation for a,b in zip(p['nodal_majorant'],raw['nodal_majorant']))
    tail_change=F(p['tail_slope'])-F(str(raw['tail_slope']));assert tail_change>=0
    eps=F(p['epsilon']);delta=F(p['delta']);nu=F(p['nu'])
    for a,b,c in [(M[0][0],M[0][2],M[2][2]),(M[1][1],M[1][3],M[3][3]),
                  (M[0][0]+1,M[0][2],M[2][2]-eps),
                  (M[1][1]+1,M[1][3],M[3][3]-eps),
                  (M[1][1]-eps,M[1][3],M[3][3])]:
        assert a>=0 and c>=0 and a*c>=b*b
    assert F(p['tail_slope'])**2>=4*(delta*delta+nu*nu)
    initial=[tuple(F(z)for z in pair)for pair in summary['initial_box']]
    C=(2/arb.pi()).sqrt()
    norms=[arb(1),C,2*C*(-arb(1)/2).exp()/arb(2).sqrt(),
           C*(1+4*(-arb(3)/2).exp())/arb(6).sqrt()]
    for j,(a,b) in enumerate(initial):assert qa(b)>=norms[j] and qa(-a)>=norms[j]
    total=0
    for node in summary['nodes']:
        index=node['node_index'];cert=json.load(gzip.open(root/'local'/f'node{index:03d}.json.gz','rt'))
        assert cert['status']=='PASS' and F(cert['q'])==F(str(data['t'][index]))
        assert F(cert['target'])==F(p['nodal_majorant'][index])
        stack=[initial];accepted=outside=0
        for event in cert['events']:
            assert stack
            box=stack.pop();total+=1
            if event=='N':assert excluded(box);outside+=1
            elif event=='C':
                center=[(a+b)/2 for a,b in box];w=[(b-a)/2 for a,b in box]
                radius=sum(abs(M[i][j])*w[i]*w[j]for i in range(4)for j in range(4))
                bound=local_point(p,cert['q'],list(map(str,center)))+qa(radius)
                assert bound<=qa(cert['target']),(index,center,bound,cert['target'])
                accepted+=1
            else:
                assert event in ['S0','S1','S2','S3'];axis=int(event[1])
                a,b=box[axis];mid=(a+b)/2
                left=list(box);right=list(box);left[axis]=(a,mid);right[axis]=(mid,b)
                stack.extend([right,left])
        assert not stack and accepted==cert['accepted'] and outside==cert['excluded']
        if index%10==0:print(json.dumps(dict(event='local_replay',node=index,total_boxes=total)),flush=True)
    return p,M,inflation,tail_change,total

def compact_leaf_replay(row,proof):
    cutoff=F(proof['domain'][1]);leaves=proof['leaves']
    assert leaves and F(leaves[0][0])==0 and F(leaves[-1][1])==cutoff
    eta=S.aa(row['eta']);prior=F(0)
    all_breaks={b for b,_ in row['target']}|{b for b,_ in row['regression']}
    for leaf in leaves:
        lo,hi=F(leaf[0]),F(leaf[1]);assert lo==prior and hi>lo;prior=hi
        assert not any(lo<b<hi for b in all_breaks)
        constants=S.segment_constants(row,lo);mid,rad=(lo+hi)/2,(hi-lo)/2
        whole=S.jets(row,S.aa(lo).union(S.aa(hi)),S.mills(lo).union(S.mills(hi)),constants)
        center=S.jets(row,S.aa(mid),S.mills(mid),constants)
        linear=leaf[3]=='linear'
        assert leaf[3] in ['linear','quadratic']
        if linear:assert whole[2][0]>=eta
        def gap(parts):
            q,reg,res=parts
            return q-reg-(res-eta/2 if linear else res*res/(2*eta))
        cc,ww=gap(center),gap(whole)
        lower=cc[0]-sum((abs(cc[j])*S.aa(rad**j)for j in range(1,4)),arb(0))-abs(ww[4])*S.aa(rad**4)
        assert lower>0,(row['index'],lo,hi,lower)
    return len(leaves)

def scalar_replay(root,data):
    prices=[];total_leaves=0
    source_hash=hashlib.sha256((root/'source_profile.json').read_bytes()).hexdigest()
    for index in range(len(data['rows'])):
        cert=json.load(open(root/'scalars'/f'row{index:02d}.json'))
        assert cert['status']=='PASS_FRESH_COMPLETE_SCALAR_HALF_LINE' and cert['source_sha256']==source_hash
        # Slacks are explicit certificate inputs, priced below and rechecked
        # by every scalar leaf.  This supports smaller repairs in tight regions.
        expected=cert['frozen_row']
        S.Q0_SLACK=F(expected['added_q0']);S.Q6_SLACK=F(expected['added_q6'])
        assert S.Q0_SLACK>=0 and S.Q6_SLACK>=0
        row=S.freeze(data,index)
        assert F(expected['r'])==row['r'] and F(expected['lambda_L1'])==row['mu'] and F(expected['eta'])==row['eta']
        assert list(map(F,expected['majorant']))==row['q']
        assert tuple(map(F,expected['phi_coefficients']))==row['coefficients']
        total_leaves+=compact_leaf_replay(row,cert['compact'])
        cutoff,tail=S.tail(row)
        assert cutoff==F(cert['tail']['domain'][0])==F(cert['compact']['domain'][1])
        assert cutoff>=max(b for b,_ in row['target'])
        assert cutoff>=max(b for b,_ in row['regression'])
        prices.append((row,S.prices(row)))
    return prices,total_leaves

def assemble(data,p,M,inflation,tail_change,prices,cap,coverage_interval=None,target_threshold=F(40,71)):
    t=[F(str(x))for x in data['t']];h=[F(str(x))for x in data['profile']['nodal_majorant']]
    slopes=[(h[i+1]-h[i])/(t[i+1]-t[i])for i in range(len(t)-1)]
    assert slopes[0]>=0 and all(a<=b for a,b in zip(slopes[:-1],slopes[1:]))
    assert F(str(data['profile']['tail_slope']))>=slopes[-1]
    rg=[F(str(x))for x in data['r_grid']];assert rg[0]==0 and all(a<b for a,b in zip(rg[:-1],rg[1:]))
    rho=F(849,1000);lo,hi=map(lambda z:F(str(z)),coverage_interval or data['profile']['amplitude_interval'])
    assert F(0)<=lo<=hi
    C=(2/arb.pi()).sqrt();v=(qa(lo)/C-1+qa(rho*rho))/qa(rho*rho)
    if rg[-1]>=rho:pass # The universal r<=rho bound covers every u.
    elif v<=0:assert rg[-1]>=rho
    else:assert qa(rg[-1])-qa(rho)*((-v.log()).sqrt()).erf()>0
    assert 194/arb(10).sqrt()-arb(633)/15>0
    eps,nu,delta=F(p['epsilon']),F(p['nu']),F(p['delta'])
    g=(nu-M[0][2])**2/M[2][2]-M[0][0]+eps-2*delta
    # After freezing, the exact scalar g itself supplies the sharp reserve.
    # The larger discovery bound is unnecessary for the finished theorem.
    L=max(F(0),-g)
    reserve=L*max(b-a for a,b in zip(rg[:-1],rg[1:]))**2/4
    endpoints=[(F(p['nodal_majorant'][0]),'r=0')]
    grouping={r:[]for r in rg[1:]}
    for row,records in prices:
        assert row['r'] in grouping and row['q'][3]>=0 and row['mu']>=0 and row['eta']>0
        grouping[row['r']].append(tuple(row['kappa_interval']))
        value=max(F(record['price_upper'])for record in records)
        endpoints.append((value+g*row['r']**2+inflation+tail_change*rho*rg[-1],row['index']))
    for intervals in grouping.values():
        intervals.sort();assert intervals and intervals[0][0]==0 and intervals[-1][1]==90
        assert all(a[1]==b[0]for a,b in zip(intervals[:-1],intervals[1:]))
    source_price,index=max(endpoints)
    ell=F(p['ell']);outer=max(lo*lo-2*ell*lo,hi*hi-2*ell*hi)
    total=source_price+reserve+outer
    assert total<cap<=target_threshold,(total,cap,target_threshold)
    return dict(status='PASS_COMPLETE_REGIONAL_BOUND',operator='P1-(18/25)P3',
                amplitude_interval=[str(lo),str(hi)],certified_rational_cap=str(cap),
                directed_total_upper=str(total),directed_total_display=float(total),
                margin_below_40_over_71=str(F(40,71)-total),active_price_row=index,
                target_threshold=str(target_threshold),margin_below_target=str(target_threshold-total),
                exact_gram_cost=str(g),negative_gram_bound=str(L),source_chord_reserve=str(reserve),
                local_height_inflation=str(inflation),tail_slope_inflation=str(tail_change),
                scope='All signs whose rotated agreement amplitude lies in the saved interval; every allowed cubic source law, shape, cumulant and source scale. This is a regional theorem and supplies no unrestricted operator cap or new Grothendieck endpoint by itself.')

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--data',type=Path,default=HERE.parent/'data')
    ap.add_argument('--output',type=Path);args=ap.parse_args();ctx.prec=512;ctx.cap=5;start=time.time()
    root=args.data;meta=json.load(open(root/'assembly.json'));cap=F(meta['regional_cap'])
    data=json.load(open(root/'source_profile.json'));local=json.load(open(root/'local_summary.json'))
    assert hashlib.sha256((root/'source_profile.json').read_bytes()).hexdigest()==local['source_sha256']
    p,M,inflation,tail_change,boxes=local_replay(root,data,local)
    prices,leaves=scalar_replay(root,data)
    result=assemble(data,p,M,inflation,tail_change,prices,cap,meta.get('coverage_interval'),F(meta.get('target_threshold','40/71')))
    result.update(local_boxes_replayed=boxes,scalar_leaves_replayed=leaves,precision_bits=ctx.prec,seconds=time.time()-start)
    print(json.dumps(result,indent=2),flush=True)
    if args.output:args.output.write_text(json.dumps(result,indent=2)+'\n')
if __name__=='__main__':main()
