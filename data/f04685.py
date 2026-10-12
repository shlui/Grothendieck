"""Directed (u,r,kappa) boxes from a fully accepted native fibre frame.

Every scalar inequality is pointwise and hence may be repriced on any
cumulant interval in [0,90].  A fixed row is selected per endpoint and
cumulant band.  Positive interpolation preserves convexity in kappa and
affinity in shape.  This exports complementary cells, not a global norm.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse,json,hashlib,sys,bisect
from flint import arb,ctx

def qa(x):
    x=F(x);return arb(x.numerator)/x.denominator

def upper(x):
    return F(str(x.upper().fmpq()))

def ceil_decimal(x,digits=12):
    scale=10**digits;return F(-((-x.numerator*scale)//x.denominator),scale)

def corner_price(S,row,k,shape,rho,l1_row=None):
    q=row['q'];eta=row['eta'];k=F(k);shape=F(shape)
    assert q[3]>=0 and row['eta']>0 and row['mu']>=0
    C=(2/arb.pi()).sqrt();K=194/arb(10).sqrt()
    sixth=15+15*qa(k)+(K-qa(F(633,15)*shape))*qa(k)*qa(k).sqrt()
    if l1_row is not None:
        pure_q=list(map(F,l1_row['majorant']));pure_eta=F(l1_row['eta'])
        assert pure_q[3]>=0 and pure_eta>0
        q=[a+row['mu']*b for a,b in zip(q,pure_q)]
        eta+=row['mu']*pure_eta
    value=row['target_mean']+qa(q[0]+q[1]+q[2]*(3+k))+qa(q[3])*sixth
    value+=qa(eta*k*(F(7,90)+shape/F(30)))
    if l1_row is None:value+=qa(row['mu'])*(qa(rho)-C)
    return upper(value)

def quadratic_max(a,b,g,lo,hi):
    points=[lo,hi]
    if g<0:
        vertex=-b/(2*g)
        if lo<vertex<hi:points.append(vertex)
    return max(a+b*x+g*x*x for x in points)

def bin_rho(ledger,lo,hi,global_rho):
    if ledger is None:return global_rho
    pieces=[p for p in ledger['leaves']
            if F(p['kappa_interval'][0])>=lo and F(p['kappa_interval'][1])<=hi]
    assert pieces and F(pieces[0]['kappa_interval'][0])==lo
    assert F(pieces[-1]['kappa_interval'][1])==hi
    assert all(F(a['kappa_interval'][1])==F(b['kappa_interval'][0])
               for a,b in zip(pieces[:-1],pieces[1:]))
    return min(global_rho,max(F(p['cap']) for p in pieces))

def main(args):
    ctx.prec=512
    package=Path(args.package);data_dir=package/'data'
    sys.path.insert(0,str(package/'code'));import scalar_halfline as S
    data=json.load(open(data_dir/'source_profile.json'))
    local=json.load(open(data_dir/'local_summary.json'));assert local['status']=='PASS'
    source_hash=hashlib.sha256((data_dir/'source_profile.json').read_bytes()).hexdigest()
    assert source_hash==local['source_sha256']
    p=local['profile'];M=[[F(z)for z in row]for row in p['M']]
    old_b=F(p['epsilon']);new_b=F(args.operator_b)
    assert 0<new_b<=old_b
    g=(F(p['nu'])-M[0][2])**2/M[2][2]-M[0][0]+new_b-2*F(p['delta'])
    shift=F(p['beta'])**2*(1/new_b-1/old_b)
    ell=F(p['ell']);rg=[F(str(z))for z in data['r_grid']]
    assert rg[0]==0 and all(a<b for a,b in zip(rg[:-1],rg[1:]))
    infl=F(p['nodal_majorant'][0])-F(str(data['profile']['nodal_majorant'][0]))
    tail=F(p['tail_slope'])-F(str(data['profile']['tail_slope']))
    global_rho=F(args.rho);assert 0<global_rho<=F(849,1000)
    by_r={r:[] for r in rg[1:]}
    for index in range(len(data['rows'])):
        cert=json.load(open(data_dir/'scalars'/f'row{index:02d}.json'))
        assert cert['status']=='PASS_FRESH_COMPLETE_SCALAR_HALF_LINE' and cert['source_sha256']==source_hash
        expected=cert['frozen_row'];S.Q0_SLACK=F(expected['added_q0']);S.Q6_SLACK=F(expected['added_q6'])
        assert S.Q0_SLACK>=0 and S.Q6_SLACK>=0
        row=S.freeze(data,index)
        assert row['q']==list(map(F,expected['majorant']))
        by_r[row['r']].append((row,infl,tail,str(package)))
    alternative_bindings=[]
    for name in args.alternative_package:
        alternate=Path(name);adir=alternate/'data'
        adata=json.load(open(adir/'source_profile.json'));alocal=json.load(open(adir/'local_summary.json'))
        assert alocal['status']=='PASS'
        ahash=hashlib.sha256((adir/'source_profile.json').read_bytes()).hexdigest()
        assert ahash==alocal['source_sha256']
        ap=alocal['profile']
        assert [[F(z)for z in row]for row in ap['M']]==M
        for key in ['ell','beta','delta','nu','epsilon','tail_slope']:
            assert F(ap[key])==F(p[key])
        assert list(map(F,ap['nodal_majorant']))==list(map(F,p['nodal_majorant']))
        assert list(map(lambda x:F(str(x)),adata['t']))==list(map(lambda x:F(str(x)),data['t']))
        ai=F(ap['nodal_majorant'][0])-F(str(adata['profile']['nodal_majorant'][0]))
        at=F(ap['tail_slope'])-F(str(adata['profile']['tail_slope']))
        assert ai>=0 and at>=0
        for index in range(len(adata['rows'])):
            cert=json.load(open(adir/'scalars'/f'row{index:02d}.json'))
            assert cert['status']=='PASS_FRESH_COMPLETE_SCALAR_HALF_LINE' and cert['source_sha256']==ahash
            expected=cert['frozen_row'];S.Q0_SLACK=F(expected['added_q0']);S.Q6_SLACK=F(expected['added_q6'])
            assert S.Q0_SLACK>=0 and S.Q6_SLACK>=0
            row=S.freeze(adata,index);assert row['q']==list(map(F,expected['majorant']))
            if row['r'] in by_r:by_r[row['r']].append((row,ai,at,str(alternate)))
        alternative_bindings.append(dict(package=str(alternate),source_sha256=ahash))
    ledger=json.load(open(args.l1_ledger)) if args.l1_ledger else None
    if ledger:assert ledger['status']=='PASS_CERTIFIED_KAPPA_SHAPE_CUBIC_L1_CAPS'
    l1_rows={}
    if args.l1_rows_file:
        l1_file=Path(args.l1_rows_file);l1_data=json.load(open(l1_file))
        assert l1_data['status']=='PASS_COMPLETE_CUBIC_L1_REFIT_PRICE_COVER'
        for name,record in l1_data['rows'].items():
            path=l1_file.parent/record['path'];source=json.load(open(path))
            assert hashlib.sha256(path.read_bytes()).hexdigest()==record['sha256']
            assert source['status'].startswith('PASS_DIRECTED_')
            assert list(map(F,source['majorant']))==list(map(F,record['majorant']))
            assert F(source['eta'])==F(record['eta'])
            l1_rows[name]=record
    kp=list(map(F,args.kappa_points));up=list(map(F,args.u_points))
    assert kp[0]==0 and kp[-1]==90 and all(a<b for a,b in zip(kp[:-1],kp[1:]))
    assert all(a<b for a,b in zip(up[:-1],up[1:]))
    # A frame is an upper bound at every u and every saved r.  Complete
    # all-source coverage is a separate union question, so an inward last
    # endpoint must never be silently promoted to full amplitude coverage.
    assert up[0]>=0
    C=(2/arb.pi()).sqrt();rho_source=F(849,1000)
    v=(qa(up[0])/C-1+qa(rho_source*rho_source))/qa(rho_source*rho_source)
    if rg[-1]>=rho_source:complete_source=True
    elif v<=0:complete_source=False
    else:complete_source=bool(qa(rg[-1])-qa(rho_source)*((-v.log()).sqrt()).erf()>0)
    step=F(args.r_step);assert step>0
    grid=set(rg);x=F(0)
    while x<rg[-1]:grid.add(x);x+=step
    grid=sorted(grid);threshold=F(args.threshold)
    cells=[];danger=[];endpoint_ledgers=[]
    for klo,khi in zip(kp[:-1],kp[1:]):
        rho=bin_rho(ledger,klo,khi,global_rho)
        corners=[(k,s) for k in (klo,khi) for s in (F(0),F(1))]
        endpoint_prices=[[F(p['nodal_majorant'][0])+shift]*4];selected=['r=0']
        for r in rg[1:]:
            candidates=[]
            for row,row_infl,row_tail,origin in by_r[r]:
                prices=[corner_price(S,row,k,s,rho)+row_infl+row_tail*rho_source*r+shift for k,s in corners]
                candidates.append((max(prices),origin,row['index'],'constant',prices))
                for name,l1_row in l1_rows.items():
                    pp=[corner_price(S,row,k,s,rho,l1_row)+row_infl+row_tail*rho_source*r+shift for k,s in corners]
                    candidates.append((max(pp),origin,row['index'],name,pp))
            _,origin,index,pure,prices=min(candidates)
            selected.append(dict(native_row=index,source_package=origin,l1_row=pure));endpoint_prices.append(prices)
        endpoint_ledgers.append(dict(kappa_interval=list(map(str,(klo,khi))),rho=str(rho),
            selected_rows=selected,prices=[[str(z)for z in row]for row in endpoint_prices]))
        for ulo,uhi in zip(up[:-1],up[1:]):
            outer=max(ulo*ulo-2*ell*ulo,uhi*uhi-2*ell*uhi)
            open_intervals=[]
            for rlo,rhi in zip(grid[:-1],grid[1:]):
                j=bisect.bisect_right(rg,rlo)-1
                if j==len(rg)-1:j-=1
                assert rg[j]<=rlo<rhi<=rg[j+1]
                bounds=[]
                for c in range(4):
                    slope=(endpoint_prices[j+1][c]-endpoint_prices[j][c])/(rg[j+1]-rg[j])
                    intercept=endpoint_prices[j][c]-slope*rg[j]
                    bounds.append(quadratic_max(intercept,slope,g,rlo,rhi)+outer)
                cap=ceil_decimal(max(bounds));accepted=cap<threshold
                cells.append(dict(u_interval=list(map(str,(ulo,uhi))),r_interval=list(map(str,(rlo,rhi))),
                    kappa_interval=list(map(str,(klo,khi))),cap=str(cap),accepted=accepted))
                if not accepted:
                    if open_intervals and open_intervals[-1][1]==rlo:open_intervals[-1][1]=rhi
                    else:open_intervals.append([rlo,rhi])
            for rlo,rhi in open_intervals:
                danger.append(dict(u_interval=list(map(str,(ulo,uhi))),r_interval=list(map(str,(rlo,rhi))),
                    kappa_interval=list(map(str,(klo,khi)))))
    out=dict(status='PASS_DIRECTED_NATIVE_CELL_PRICING_BASE_LEAF_CERTIFICATE_REQUIRED',
        package=str(package),u_grid=list(map(str,up)),kappa_grid=list(map(str,kp)),
        source_grid=list(map(str,rg)),beta_squared=str(F(p['beta'])**2),
        source_sha256=source_hash,operator=f'P1-({new_b})P3',threshold=str(threshold),
        cells=cells,unclosed_boxes=danger,endpoint_prices=endpoint_ledgers,
        exact_gram_cost=str(g),ell=str(ell),source_endpoint=str(rg[-1]),
        analytic_L1_rho=str(global_rho),source_bridge_rho='849/1000',
        original_operator_b=str(old_b),operator_b=str(new_b),local_constant_transport=str(shift),
        source_coverage_complete_at_lower_u=complete_source,
        alternative_scalar_bindings=alternative_bindings,
        scalar_rows_validity='Every accepted pointwise half-line comparison is valid for every normalized cubic with kappa in [0,90]; original discovery bands do not restrict repricing.',
        continuous_proof='Positive source chords preserve convex kappa and affine shape pricing. Exact quadratic maximization in r and convex outer-u endpoint bounds cover each whole rectangular cell.',
        scope='Complementary native upper-bound cells only. Unclosed cells are not counterexamples. No global norm cap or Grothendieck endpoint follows until native and R3 cells cover every physical source triple.',
        l1_ledger_sha256=hashlib.sha256(Path(args.l1_ledger).read_bytes()).hexdigest() if args.l1_ledger else None)
    out['l1_rows_file_sha256']=hashlib.sha256(Path(args.l1_rows_file).read_bytes()).hexdigest() if args.l1_rows_file else None
    Path(args.output).write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(dict(status=out['status'],cells=len(cells),accepted=sum(c['accepted']for c in cells),
                         unclosed_boxes=len(danger),output=args.output)))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--package',required=True);ap.add_argument('--output',required=True)
    ap.add_argument('--u-points',nargs='+',required=True)
    ap.add_argument('--kappa-points',nargs='+',default=['0','.1','.3','.53','.8','1','1.1','2','4','10','12','20','40','60','90'])
    ap.add_argument('--r-step',default='1/500');ap.add_argument('--threshold',default='40/71')
    ap.add_argument('--rho',default='849/1000');ap.add_argument('--l1-ledger')
    ap.add_argument('--l1-rows-file')
    ap.add_argument('--operator-b',default='18/25')
    ap.add_argument('--alternative-package',nargs='+',default=[])
    main(ap.parse_args())
