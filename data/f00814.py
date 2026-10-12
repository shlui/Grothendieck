"""Producer-free directed whole-cell R3/B0 fixed-multiplier audit.

Only independent definitions from the frozen proof are imported; all writes
are restricted to this new numeric directory. LP is proposal-only: fixed
rational multipliers are freshly priced at all product vertices afterwards.
"""
import os,sys
sys.dont_write_bytecode=True
os.environ['GR_MPFR_BITS']='768'
from pathlib import Path
from fractions import Fraction as F
from itertools import product
import argparse,gzip,hashlib,json
ap=argparse.ArgumentParser();ap.add_argument('--base',type=Path,default=Path('verified_2026-10-01/tonight/operator/preferred_global'));ap.add_argument('--index',type=int,default=195896);ap.add_argument('--profile',type=Path);ap.add_argument('--reserve0',default='1/10000');ap.add_argument('--reserve6',default='1/100000000');ap.add_argument('--output',type=Path);ap.add_argument('--b',default='711/1000');ap.add_argument('--u',nargs=2);ap.add_argument('--r',nargs=2);ap.add_argument('--kappa',nargs=2);args=ap.parse_args()
BASE=args.base.resolve();OUT=Path(__file__).resolve().parent
sys.path.insert(0,str(BASE/'tonight/audit'))
from mpfr_interval import I,BITS,UP,DOWN
assert BITS==768
sys.path.insert(0,str(BASE/'tonight/audit/global_b7171'))
import replay_R3_corners as R
from scipy.optimize import linprog
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):
    with (gzip.open(p,'rt')if p.suffix=='.gz'else p.open())as f:return json.load(f)
def upper(v):return F(v.hi.text(UP,240))
def rec(v):return {'lower':v.lo.text(DOWN,90),'upper':v.hi.text(UP,90)}
def exact(x):return {k:list(map(F,v))if isinstance(v,list)else F(v)for k,v in x.items()}
source=BASE/'tonight/recovery/profile_cover/global_upper_b713_native_kappa68_twentyninebank_envelope.json.gz';d=read(source)
scanpath=BASE/'tonight/recovery/profile_cover/assigned_twenty-nine_coefficient_scan.json';scan=read(scanpath)
assert scan['parent_sha256']==sha(source)
cell=dict(d['cells'][args.index]);
if args.u or args.r or args.kappa:
    assert args.profile and args.u and args.r and args.kappa
    cell={'u_interval':args.u,'r_interval':args.r,'kappa_interval':args.kappa,'winner':'R3:lastface7575','requested_new_rectangle':True}
name=cell['winner'];assert name.startswith('R3:');name=name[3:]
binding=scan['profile_beta_bindings']['R3:'+name];profile=BASE/binding['profile_path'];assert sha(profile)==binding['profile_sha256']
b=F(args.b)
if args.profile:
    profile=args.profile.resolve();name=profile.stem
    ledgerpath=None;ledger={'reserves':[args.reserve0,args.reserve6]}
    whole_q_record={'status':'PENDING_INDEPENDENT_WHOLE_Q_ACCEPTANCE','profile_sha256':sha(profile),'reserves':[args.reserve0,args.reserve6]}
else:
    cachepath=BASE/'tonight/audit/global_b7171/full_upper57_twentysevenbank_R3_corner_replay.json';cache=read(cachepath);whole=cache['complete_whole_q_bindings'][name]
    assert whole['profile']=={'path':binding['profile_path'],'sha256':sha(profile)}
    ledgerpath=BASE/whole['accepted_whole_q_gate']['path'];assert sha(ledgerpath)==whole['accepted_whole_q_gate']['sha256'];ledger=read(ledgerpath)
    assert ledger['profile_sha256']==sha(profile)
    for item in ledger['reports'].values():assert sha(BASE/item['path'])==item['sha256']
    whole_q_record={'path':str(ledgerpath),'sha256':sha(ledgerpath)}
profile_start_sha=sha(profile);p=exact(read(profile)['profile'])
res0,res6=map(F,ledger['reserves']);row=dict(p=p,kind='rank',reserve0=res0,reserve6=res6)
spectra=R.spectral(row,b);c=R.coefficients(row,b);assert c['r_quadratic']==b
E00,E02,E22=R.even(row);theta=c['theta'];q0,q2,q4,q6=p['major'];q0+=res0;q6+=res6
A=q6+p['zeta']**2;D=p['G']-p['d0']-p['zeta']**2
assert A>0 and D>0 and (I(194)/I(10).sqrt()-F(633,15)).above(0)
def price(u,r,k,z):
    sixth=(I(194)/I(10).sqrt()-F(633,9)*(z-F(2,5)))*I(k)*I(k).sqrt()
    fourth=I(k*k)+F(9,4)*I(k*z)*I(k*z).sqrt()
    y=p['y']+E02*(r-p['r'])
    tv=k*(1+z)/9-2*theta*k/3+theta*theta*(k+2)
    cov=k/3-theta*(k+2)
    terms={'outer_u':I(u*u-2*p['ell']*u),'beta_completion':I(6*p['beta']**2/b),
      'r_y_transport':I((b-E00)*r*r-2*p['gamma']*r+y*y/E22),
      'major_q0_q2_q4':I(q0+q2+(3+k)*q4),'major_q6':I(q6)*(15+15*k+sixth),
      'centered_variance':I(p['eta']*tv/2),'positive_nuisance_covariance':I(p['nuisance2']*cov),
      'zeta_second_and_sixth':I(p['zeta']**2)*(6+9*k+sixth),'negative_cumulant_credit':I(-p['c']*k),
      'fourth_moment_price':I(D)*fourth}
    total=sum(terms.values(),I(0));accepted=R.price(row,c,k,z,u,r)
    assert total.lo<=accepted.hi and accepted.lo<=total.hi
    return total,terms
bankpath=BASE/d['linear_source_bank_binding']['path'];assert sha(bankpath)==d['linear_source_bank_binding']['sha256'];bank=read(bankpath)
acceptedbankpath=BASE/'tonight/audit/linear_contact_bank_independent.json';acceptedbank=read(acceptedbankpath)
assert acceptedbank['input_sha256']==sha(bankpath)and acceptedbank['status']=='PASS_INDEPENDENT_COMPLETE_LINEAR_CONTACT_SOURCE_BANK'
assert len(bank['rows'])==60

def b0_price(item,k,z):
    qs=list(map(F,item['polynomial_coefficients_q0_q2_q4_q6']));assert qs[3]>0
    return I(qs[0]+qs[1]+(3+k)*qs[2])+I(qs[3])*(15+15*k+(I(194)/I(10).sqrt()-F(633,9)*(z-F(2,5)))*I(k)*I(k).sqrt())
ui,ri,ki=(list(map(F,cell[x]))for x in('u_interval','r_interval','kappa_interval'));zs=[F(2,5),F(1)]
vertices=[];basevals=[];slopes=[];constraint_rows=[]
for u,r,k,z in product(ui,ri,ki,zs):
    v,terms=price(u,r,k,z);s=[];gates=[]
    for j,item in enumerate(bank['rows']):
        value=b0_price(item,k,z);gate=value-F(item['a'])*u-r
        s.append(upper(gate));gates.append({'row':j,'a':item['a'],'slack':rec(gate),'strictly_excluded':gate.below(0)})
    vertices.append({'u':str(u),'r':str(r),'kappa':str(k),'shape':str(z),'R3_price':rec(v),'terms':{key:rec(value)for key,value in terms.items()},'source_rows':gates})
    basevals.append(upper(v));slopes.append(s)
# Proposal uses float. Every resulting lambda is rounded to nonnegative rational
# and checked against every vertex in exact directed arithmetic below.
res=linprog([0.0]*60+[1.0],A_ub=[[float(x)for x in s]+[-1.0]for s in slopes],b_ub=[-float(v)for v in basevals],bounds=[(0,None)]*60+[(None,None)],method='highs')
assert res.success,res.message
scale=10**9;lambdas=[F(max(0,round(float(x)*scale)),scale)for x in res.x[:60]]
cert=[v+sum((lam*s for lam,s in zip(lambdas,ss)),F(0))for v,ss in zip(basevals,slopes)]
baseline=max(basevals);maximum=max(cert);winning=basevals.index(baseline);newwinning=cert.index(maximum)
for j,v in enumerate(vertices):v['coupled_fixed_rational_upper']=str(cert[j])
# The q0 reserve and q6 reserve contributions are exact/MPFR respectively.
u,r,k,z=(F(vertices[winning][x])for x in('u','r','kappa','shape'))
sixth=I(15+15*k)+(I(194)/I(10).sqrt()-F(633,9)*(z-F(2,5)))*I(k)*I(k).sqrt()
transport_uniform=None if args.profile else F(cell['cap'])+6*p['beta']**2*(1/b-1/F(713,1000))+(b-F(713,1000))*ri[0]**2
chord=(F(ui[1]-ui[0])**2+b*F(ri[1]-ri[0])**2)/4
out={'status':'PASS_INDEPENDENT_RESTRICTED_WHOLE_CELL_R3_B0_FIXED_MULTIPLIER_PRICE','scope':'Exactly one frozen cell and selected accepted whole-q R3 profile; not a complete global norm certificate. Numerical LP only chooses fixed rational nonnegative multipliers; their complete product-corner certificate is directed/exact.',
 'bits':BITS,'operator_b':str(b),'target_cap':'40/71','index':args.index,'profile':name,'source':{'path':str(source),'sha256':sha(source)},'scan':{'path':str(scanpath),'sha256':sha(scanpath)},
 'profile_binding':{'path':str(profile),'sha256':sha(profile)},'whole_q_gate':whole_q_record,'bank':{'path':str(bankpath),'sha256':sha(bankpath)},'accepted_bank':{'path':str(acceptedbankpath),'sha256':sha(acceptedbankpath)},
 'independent_engine':{'path':str(Path(R.__file__)),'sha256':sha(Path(R.__file__))},'script_sha256':sha(Path(__file__)),
 'cell':cell,'coefficients':{x:str(v)for x,v in c.items()},'six_PSD_gates':spectra,'vertices':vertices,
 'base_directed_upper':str(baseline),'base_winning_vertex':winning,'assigned_uniform_transport_upper':None if transport_uniform is None else str(transport_uniform),'uniform_transport_excess':None if transport_uniform is None else str(transport_uniform-baseline),
 'maximum_possible_u_r_chord_excess':str(chord),'reserve_contributions_at_base_winner':{'q0_exact':str(res0),'q6_directed':rec(I(res6)*sixth)},
 'coupled_B0_fixed_multipliers':[{'row':j,'a':bank['rows'][j]['a'],'lambda':str(lam)}for j,lam in enumerate(lambdas)if lam],
 'coupled_directed_upper':str(maximum),'coupled_winning_vertex':newwinning,'certified_improvement':str(baseline-maximum),
 'target_pass':maximum<F(40,71),'remaining_target_gap':str(maximum-F(40,71)),
 'proof':'Each accepted B0 row gives a*u+r<=E J_a(Q). Add its fixed nonnegative multiple to the R3 upper bound. Price all terms at the same kappa/shape corner before taking a maximum. The combined price is separately convex in u,r,kappa,shape: u² coefficient1, rank r² coefficientb>0, A>0,D>0, and every source q6>0 with 194/sqrt10-633/15>0. Product-corner interpolation thus covers the entire rectangle. Fixed rational multipliers do not depend on the actual point. Shape is z in[2/5,1]; source s=(5/3)(z-2/5). Frozen theta uses profile kappa, not actual kappa.'}
if args.profile:out['status']='PASS_INDEPENDENT_RESTRICTED_SOURCE_PRICE_CONDITIONAL_ON_NEW_WHOLE_Q_GATE';out['scope']='Source pricing at all16 product vertices is certified; independent whole-q gate acceptance remains pending. No new global bound is asserted.'
assert sha(profile)==profile_start_sha and sha(source)==scan['parent_sha256']and sha(bankpath)==d['linear_source_bank_binding']['sha256']
output=args.output.resolve()if args.output else OUT/'worst_cell_B0_coupling.json';assert output.is_relative_to(OUT);output.write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'status':out['status'],'profile':name,'base_upper':float(baseline),'base_winner':{x:vertices[winning][x]for x in('u','r','kappa','shape')},'excluded_B0_rows_at_winner':[x['row']for x in vertices[winning]['source_rows']if x['strictly_excluded']],'coupled_upper':float(maximum),'improvement':float(baseline-maximum),'target_pass':out['target_pass'],'multipliers':out['coupled_B0_fixed_multipliers']},indent=2))
