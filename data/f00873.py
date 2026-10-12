"""Fresh independent expanded native prices on the curved physical domain.

Imports only the frozen independent certificate reader and MPFR backend.
No new-round producer, pricer, or boundary helper is imported.
"""
import os, sys
sys.dont_write_bytecode = True
os.environ['GR_MPFR_BITS'] = '768'
from pathlib import Path
from fractions import Fraction as F
from functools import lru_cache
import argparse, json, hashlib, time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
BASE = ROOT/'verified_2026-10-01/tonight/operator/preferred_global'
sys.path.insert(0,str(BASE/'tonight/audit'))
from mpfr_interval import I, UP, DOWN, BITS, selftest
assert BITS == 768
import replay_native_R3_unions as N
ap = argparse.ArgumentParser()
ap.add_argument('--source',type=Path,default=ROOT/'new_round_2026-10-02/native_curved_endpoint_prices.json')
ap.add_argument('--output',type=Path,default=HERE/'CURVED_NATIVE_ENDPOINT_PRICE_REPLAY.json')
args = ap.parse_args()
assert args.output.resolve().is_relative_to(HERE)
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p): return json.loads(p.read_text())
def bound_record(p): return {'path':str(p.relative_to(ROOT)),'sha256':sha(p)}
source = args.source.resolve()
source_sha = sha(source)
d = read(source)
assert d['status'] == 'COMPLETE_CURVED_NATIVE_ENDPOINT_PRICE_PRODUCER_PENDING_INDEPENDENT_REPLAY'
assert d['physical_constraint'] == 'kappa*zeta<=36'
b = F(d['operator_b'])
assert b == F(711,1000) and d['precision_bits'] == 768
configpath = Path(d['input_config']['path'])
assert sha(configpath) == d['input_config']['sha256']
conf = read(configpath)
assert conf['status'] == 'PASS_ALL_FRESH_FINE_KAPPA_ENDPOINT_LEDGERS'
points = list(map(F,conf['kappa_points']))
assert len(points) == 68 and points[0] == 0 and points[-1] == 90
assert all(a<c for a,c in zip(points,points[1:]))
bands = list(zip(points,points[1:]))
proofpath = ROOT/'new_round_2026-10-02/audit/PHYSICAL_SHAPE_BOUNDARY_PRICING.md'
tensorpath = ROOT/'new_round_2026-10-02/audit/TENSOR_SHAPE_NECESSARY_CONSTRAINT.md'
l1pricepath = HERE/'ACCEPTED_L1_KAPPA_SPECIFIC_CAPS.json'
l1price = read(l1pricepath)
assert l1price['status'].startswith('PASS')

def candidates(ka,kb,za,zb):
    # Derived independently from straight-boundary endpoints and transitions.
    assert 0 <= ka < kb <= 90 and F(2,5) <= za < zb <= 1
    result = {(k,z) for k in (ka,kb) for z in (za,zb) if k*z<=36}
    for k in (ka,kb):
        if k>0:
            z=F(36)/k
            if za<=z<=zb: result.add((k,z))
    for z in (za,zb):
        k=F(36)/z
        if ka<=k<=kb: result.add((k,z))
    assert result
    return sorted(result)

frames = {}
premises = {str(configpath):sha(configpath),str(N.L1FILE):sha(N.L1FILE),str(l1pricepath):sha(l1pricepath)}
for bind in conf['frames']:
    mask = BASE/bind['mask']
    raw = read(mask)
    p = N.package(bind['package'])
    assert raw['source_sha256'] == p['source_sha256'] == bind['source_sha256']
    assert raw['status'] == 'PASS_DIRECTED_NATIVE_CELL_PRICING_BASE_LEAF_CERTIFICATE_REQUIRED'
    assert raw['l1_rows_file_sha256'] == sha(N.L1FILE)
    origins = {bind['package']}
    for alt in raw['alternative_scalar_bindings']:
        other = N.package(alt['package'])
        assert other['source_sha256'] == alt['source_sha256']
        N.sameframe(p,other)
        origins.add(alt['package'])
    M=p['M'];pr=p['profile'];old=F(pr['epsilon'])
    assert 0<b<=old
    gates=[]
    for aa,oo,cc in [(M[0][0],M[0][2],M[2][2]),(M[1][1],M[1][3],M[3][3]),(M[0][0]+1,M[0][2],M[2][2]-b),(M[1][1]+1,M[1][3],M[3][3]-b),(M[1][1]-b,M[1][3],M[3][3])]:
        det=aa*cc-oo*oo
        assert aa>=0 and cc>=0 and det>=0
        gates.append({'diagonal':[str(aa),str(cc)],'off_diagonal':str(oo),'determinant':str(det)})
    g=(F(pr['nu'])-M[0][2])**2/M[2][2]-M[0][0]+b-2*F(pr['delta'])
    shift=F(pr['beta'])**2*(1/b-1/old)
    assert F(raw['local_constant_transport']) == F(pr['beta'])**2*(1/F(raw['operator_b'])-1/old)
    groups={tuple(map(F,x['kappa_interval'])):x for x in raw['endpoint_prices']}
    assert sorted(groups) == bands
    frames[bind['mask']]={'p':p,'bind':bind,'raw':raw,'groups':groups,'origins':origins,'g':g,'shift':shift,'gates':gates,'mask':mask}
    premises[str(mask)] = sha(mask)
assert len(frames)==d['frames']==7

pi=I.pi();nu=(I(2)/pi).sqrt();K=I(194)/I(10).sqrt()
assert (K-F(633,15)).above(0)
@lru_cache(None)
def gaussian_mean(package,r):
    p=N.package(package)
    total=I(p['h'][0])
    increments=[p['slopes'][0]]+[right-left for left,right in zip(p['slopes'],p['slopes'][1:])]+[p['tau']-p['slopes'][-1]]
    assert len(increments)==len(p['t']) and all(v>=0 for v in increments)
    for knot,slope in zip(p['t'],increments):
        x=I(knot/r)
        total += I(r*slope)*(nu*(-x*x/2).exp()-x*(x/I(2).sqrt()).erfc())
    return total

scalar_bindings={};l1_bindings={};calls=0
@lru_cache(None)
def expanded_price(package,index,k,z,rho,l1):
    global calls
    row=N.scalar(package,index)
    local=N.package(package)
    mu=row['mu'];eta=row['eta'];q=list(row['q'])
    assert mu>=0 and eta>0 and q[3]>=0
    scalarpath=local['path']/f'data/scalars/row{index:02d}.json'
    scalar_bindings[(package,index)]={'package':package,'index':index,'certificate':bound_record(scalarpath),'r':str(row['r']),'q':list(map(str,q)),'lambda_L1':str(mu),'eta':str(eta)}
    if l1!='constant':
        qq,ee=N.pure(l1)
        assert qq[3]>=0 and ee>0
        q=[a+mu*c for a,c in zip(q,qq)]
        eta+=mu*ee
        purepath=N.L1FILE.parent/N.L1['rows'][l1]['path']
        l1_bindings[l1]={'certificate':bound_record(purepath),'majorant':list(map(str,qq)),'eta':str(ee)}
    s=(5*z-2)/3
    rational=q[0]+q[1]+(3+k)*q[2]+(15+15*k)*q[3]+eta*k*(F(7,90)+s/30)
    result=gaussian_mean(package,row['r'])+I(rational)+I(q[3])*(K-F(633,15)*s)*I(k)*I(k).sqrt()
    if l1=='constant': result+=I(mu)*(I(rho)-nu)
    # Nonnegative combined q6 gives convexity on every curved boundary.
    assert q[3]>=0 and eta>0
    calls+=1
    return result

seen=set();outtables=[];evaluation_count=0;start=time.time()
for table in d['tables']:
    name=table['frame'];assert name in frames
    frame=frames[name];p=frame['p'];pr=p['profile'];rg=p['rg']
    parentband=tuple(map(F,table['parent_kappa_interval']))
    ka,kb=map(F,table['kappa_interval']);za,zb=map(F,table['shape_interval'])
    assert parentband==(ka,kb) and (za,zb)==(F(2,5),F(1))
    key=(name,parentband);assert key not in seen;seen.add(key)
    group=frame['groups'][parentband];rho=F(group['rho'])
    assert 0<rho<=F(847,1000)
    candidate_points=candidates(ka,kb,za,zb)
    assert [list(map(str,v)) for v in candidate_points] == table['feasible_curved_vertices']
    assert list(map(F,table['r_grid']))==rg
    assert table['native_mask_sha256']==sha(frame['mask'])
    assert table['source_package']==frame['bind']['package'] and table['source_sha256']==p['source_sha256']
    assert F(table['ell'])==F(pr['ell']) and F(table['r_quadratic'])==frame['g']
    caps=list(map(F,table['endpoint_uniform_caps']))
    assert len(caps)==len(rg)==len(group['selected_rows'])
    assert caps[0]==F(pr['nodal_majorant'][0])+frame['shift']
    assert len(table['details'])==len(rg)-1
    fresh=[]
    for j,detail in enumerate(table['details'],1):
        choice=group['selected_rows'][j]
        package=choice['source_package'];assert package in frame['origins']
        row=N.scalar(package,choice['native_row']);local=N.package(package)
        assert detail['index']==j and row['r']==F(detail['r'])==rg[j]
        for field in ('source_package','native_row','l1_row'):assert detail[field]==choice[field]
        assert detail['native_scalar_sha256']==row['certificate_sha256']
        assert [(F(v['kappa']),F(v['shape'])) for v in detail['vertices']]==candidate_points
        saved=list(map(lambda v:F(v['upper']),detail['vertices']))
        assert max(saved)==F(detail['upper'])==caps[j]
        checked=[]
        for point,saved_upper in zip(candidate_points,saved):
            k,z=point
            value=expanded_price(package,choice['native_row'],k,z,rho,choice['l1_row'])
            value+=I(local['inflation']+local['tail']*F(849,1000)*rg[j]+frame['shift'])
            assert value.below(saved_upper),(name,parentband,j,point,value.record(),str(saved_upper))
            evaluation_count+=1
            checked.append({'kappa':str(k),'shape':str(z),'lower':value.lo.text(DOWN,90),'upper':value.hi.text(UP,90),'saved_upper':str(saved_upper)})
        fresh.append({'index':j,'certificate_sha256':row['certificate_sha256'],'vertices':checked,'saved_uniform_upper':str(caps[j])})
    outtables.append({'frame':name,'kappa_interval':[str(ka),str(kb)],'shape_interval':[str(za),str(zb)],'r_grid':list(map(str,rg)),'r0_exact':str(caps[0]),'details':fresh})
    if len(outtables)%67==0:print(json.dumps({'tables':len(outtables),'fresh_directed_endpoint_prices':evaluation_count,'seconds':round(time.time()-start,2)}),flush=True)
assert len(seen)==len(frames)*len(bands)==469
assert seen=={(name,band) for name in frames for band in bands}
for name,p in N.PACKAGES.items():
    assert d['accepted_packages'][name]=={k:p[k] for k in ('source_sha256','manifest_sha256','replay_sha256')}
    premises[str(p['path']/'SHA256.json')]=p['manifest_sha256']
    premises[str(p['path']/'replay.json')]=p['replay_sha256']
    premises[str(p['path']/'data/source_profile.json')]=p['source_sha256']
assert set(d['accepted_packages'])==set(N.PACKAGES)
assert sha(source)==source_sha and all(sha(Path(p))==digest for p,digest in premises.items())
out={'status':'PASS_INDEPENDENT_CURVED_NATIVE_SELECTED_ENDPOINT_PRICES',
 'source':bound_record(source),'input_sha256':source_sha,'input_config':bound_record(configpath),
 'operator_b':str(b),'precision_bits':BITS,'frames':len(frames),'kappa_bands':len(bands),'tables':outtables,
 'curved_uniform_tables':len(outtables),'selected_endpoint_prices_recomputed':evaluation_count,'distinct_expanded_prices':calls,
 'frame_parameters':{name:{'target_five_PSD_gates':x['gates'],'r_quadratic':str(x['g']),'beta_squared_transport':str(x['shift']),'reference_b':x['p']['profile']['epsilon'],'package':x['bind']['package'],'mask_sha256':sha(x['mask'])} for name,x in frames.items()},
 'selected_scalar_bindings':list(scalar_bindings.values()),'selected_L1_bindings':l1_bindings,
 'accepted_package_bindings':d['accepted_packages'],'verified_premise_files':premises,
 'boundary_proof':bound_record(proofpath),'physical_tensor_proof':bound_record(tensorpath),
 'accepted_L1_price_binding':bound_record(l1pricepath),'backend':selftest(),
 'independent_reader':bound_record(Path(N.__file__)),'script_sha256':sha(Path(__file__)),
 'producer_modules_imported':False,'seconds':time.time()-start,
 'formula':'Fresh Gaussian hinge mean plus combined exact native/L1 q0+q2+(3+k)q4+(15+15k)q6, q6(194/sqrt10-633s/15)k^(3/2), eta*k(7/90+s/30), and constant-L1 credit when selected. Source shape s=(5z-2)/3. Exact local inflation, accepted tail*.849*r and beta²(1/b-1/.72) transport are added.',
 'scope':'All469 native endpoint tables on the necessary physical domain kappa*z<=36, with exact selected scalar/L1 certificate hashes, nonnegative curved-convexity guards and target711 fivePSD gates. Existing complete scalar/local proofs are reused; no leaf proofs, radial chord/global partition, lower amplitude cover or global norm theorem is replayed or asserted.'}
args.output.write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:out[k] for k in ('status','operator_b','curved_uniform_tables','selected_endpoint_prices_recomputed','distinct_expanded_prices','seconds')},indent=2))
