"""Producer-free fresh-cache replay; accepted frozen old gates plus new components.

No scalar/local-tree replay or geometry producer is imported. All profile,
reserve, exact target PSD and complete accepted whole-q bindings are checked.
"""
import os,sys
sys.dont_write_bytecode=True;os.environ['GR_MPFR_BITS']='768'
from pathlib import Path
from fractions import Fraction as F
from itertools import product
from functools import lru_cache
import argparse,gzip,hashlib,json
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2];BASE=ROOT/'verified_2026-10-01/tonight/operator/preferred_global'
ap=argparse.ArgumentParser();ap.add_argument('source',type=Path);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args();assert args.output.resolve().is_relative_to(HERE)
sys.path.insert(0,str(BASE/'tonight/audit'))
from mpfr_interval import I,BITS,UP,DOWN,selftest
assert BITS==768
sys.path.insert(0,str(BASE/'tonight/audit/global_b7171'))
import replay_R3_corners as R
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):
 with(gzip.open(p,'rt')if p.suffix=='.gz'else p.open())as h:return json.load(h)
def resolve(x,base):
 p=Path(x);return p if p.is_absolute()else base/p
verified={}
def check_records(v,base):
 if isinstance(v,dict):
  if 'path'in v and 'sha256'in v:
   p=resolve(v['path'],base);assert sha(p)==v['sha256'],p;verified[str(p)]=v['sha256']
  for child in v.values():check_records(child,base)
 elif isinstance(v,list):
  for child in v:check_records(child,base)
source=args.source.resolve();initial_sha=sha(source);d=read(source);b=F(d['operator_b']);assert b==F(711,1000)
assert d['status']=='COMPLETE_PARENT_LINKED_KAPPA_SHAPE_RANK_REFINEMENT_PENDING_INDEPENDENT_AUDIT'
parent=Path(d['parent']['path']);assert sha(parent)==d['parent']['sha256'];parent_source=read(parent)
assert d['parent']['sha256']=='b5fd7a4edb6a8275d2f00caa9d21edceac85fefc812842236e1e2e25aa76269b'
base_source=Path(d['base_global']['path']);assert sha(base_source)==d['base_global']['sha256'];base_source_data=read(base_source)
assert base_source_data['physical_constraint']=='4*kappa*zeta-kappa<=54'
u_parent=Path(d['u_parent']['path']);assert sha(u_parent)==d['u_parent']['sha256'] and parent_source['parent']==d['u_parent']
assert parent_source['base_global']==d['base_global']
base_report_path=HERE/'BANK43_LINEAR_TENSOR_R3_CACHE_REPLAY.json';extra_report_path=HERE/'BANK66_MONOTONE_R3_CACHE_REPLAY.json'
base_report=read(base_report_path);extra_report=read(extra_report_path)
for report in (base_report,extra_report):
 assert report['status']=='PASS_INDEPENDENT_NEW_ROUND_COMPLETE_R3_CACHE_PRICES' and report['bits']==768 and F(report['operator_b'])==b
 for path,digest in report['verified_premise_files'].items():assert sha(Path(path))==digest;verified[path]=digest
assert base_report['source']['sha256']==sha(base_source)
extra_source=ROOT/extra_report['source']['path'];assert sha(extra_source)==extra_report['source']['sha256'];extra_data=read(extra_source)
assert parent_source['cache']==extra_data['cache'] and parent_source['profiles']==extra_data['profiles'] and parent_source['base_global']==extra_data['base_global']
assert parent_source['parent']==extra_data['parent']
parent_profile_metadata=dict(base_source_data['profiles']);assert not set(parent_profile_metadata)&set(parent_source['profiles']);parent_profile_metadata.update(parent_source['profiles'])
accepted_bindings=dict(base_report['complete_whole_q_bindings']);assert not set(accepted_bindings)&set(extra_report['complete_whole_q_bindings']);accepted_bindings.update(extra_report['complete_whole_q_bindings'])
assert len(accepted_bindings)==66 and set(accepted_bindings)==set(parent_profile_metadata)
used_profiles={x['profile']for x in d['cache']};assert used_profiles<=set(accepted_bindings) and len(used_profiles)==63
rows={};bindings={};coefficients={}
for name in sorted(used_profiles):
 item=parent_profile_metadata[name];priorbinding=accepted_bindings[name]
 path=Path(item['profile_path']);assert sha(path)==item['profile_sha256']==priorbinding['profile']['sha256'];assert str(path)==priorbinding['profile']['path']
 declared=Path(item['acceptance_path']);assert sha(declared)==item['acceptance_sha256']
 raw=read(path);assert raw['kind']=='EXACT_RATIONAL_R3_PROFILE' and priorbinding['kind']=='rank';kind='rank';p=R.exact(raw['profile'])
 reserves=priorbinding['reserves'];binding=priorbinding['whole_q_premise']
 res0,res6=map(F,reserves);assert res0>=0 and res6>0
 if name=='highk_native76625':assert (res0,res6)==(F(3,100),F(1,5000000))
 row={'p':p,'kind':kind,'reserve0':res0,'reserve6':res6};spectral=R.spectral(row,b)
 E00,E02,E22=R.even(row);theta=p['kappa']/(3*(p['kappa']+2));yc=p['y']-E02*p['r'];q0,q2,q4,q6=p['major'];q0+=res0;q6+=res6
 A=q6+p['zeta']**2;D=p['G']-p['d0']-p['zeta']**2;assert A>0 and D>0
 c={'theta':theta,'A':A,'B':D,'constant':6*p['beta']**2/b+yc*yc/E22+q0+q2+3*q4+15*q6+p['eta']*theta*theta-2*p['nuisance2']*theta+6*p['zeta']**2,'q4':q4,'q6':q6,'yc':yc,'r_quadratic':b-E00+E02**2/E22,'r_linear':-2*p['gamma']+2*yc*E02/E22,'ell':p['ell'],'b':b}
 assert c==R.coefficients(row,b)
 for rr in(F(0),F(1,3),F(1)):
  y=p['y']+E02*(rr-p['r']);assert y-E02*rr==yc
  assert (b-E00)*rr*rr-2*p['gamma']*rr+y*y/E22==c['r_quadratic']*rr*rr+c['r_linear']*rr+yc*yc/E22
 row.update(c=c,q=[q0,q2,q4,q6],E22=E22);rows[name]=row
 bindings[name]={'profile':{'path':str(path),'sha256':sha(path)},'kind':kind,'reserves':list(map(str,(res0,res6))),'target_spectral_gates':spectral,'whole_q_premise':binding}
 assert spectral==priorbinding['target_spectral_gates']
 assert {k:str(v)for k,v in c.items()}==((base_report if name in base_report['expanded_coefficients']else extra_report)['expanded_coefficients'][name])
 coefficients[name]={k:str(v)for k,v in c.items()}
assert len(rows)==63
K=I(194)/I(10).sqrt();assert (K-F(633,15)).above(0)
def candidates(ka,kb,za,zb):
 result={(k,z)for k,z in product((ka,kb),(za,zb))if k*(4*z-1)<=54}
 for k in(ka,kb):
  if k>0 and za<=(k+54)/(4*k)<=zb:result.add((k,(k+54)/(4*k)))
 for z in(za,zb):
  if ka<=F(54)/(4*z-1)<=kb:result.add((F(54)/(4*z-1),z))
 assert result
 return sorted(result)
@lru_cache(None)
def point(name,k,z):
 row=rows[name];p=row['p'];c=row['c'];theta=c['theta'];q0,q2,q4,q6=row['q']
 tv=k*(1+z)/9-2*theta*k/3+theta*theta*(k+2);cov=k/3-theta*(k+2)
 rational=6*p['beta']**2/b+c['yc']**2/row['E22']+q0+q2+(3+k)*q4+(15+15*k)*q6+p['eta']*tv/2+p['nuisance2']*cov+p['zeta']**2*(6+9*k)-p['c']*k
 sixth=(K-F(633,9)*(z-F(2,5)))*I(k)*I(k).sqrt();fourth=I(k*k)+F(9,4)*I(k*z)*I(k*z).sqrt()
 value=I(rational)+I(c['A'])*sixth+I(c['B'])*fourth;alternative=R.price(row,c,k,z)
 assert value.lo<=alternative.hi and alternative.lo<=value.hi
 return value
prior_split_report_path=HERE/'BANK63_SPLIT_LINEAR_R3_CACHE_REPLAY.json';prior_split_report=read(prior_split_report_path)
assert prior_split_report['status']=='PASS_INDEPENDENT_NEW_ROUND_COMPLETE_R3_CACHE_PRICES'and prior_split_report['bits']==768 and F(prior_split_report['operator_b'])==b
prior_split_source=ROOT/prior_split_report['source']['path'];assert sha(prior_split_source)==prior_split_report['source']['sha256'];prior_split_data=read(prior_split_source)
def tablekey(entry):return(entry['profile'],*map(F,entry['kappa_interval']+entry['shape_interval']))
prior_raw={tablekey(entry):entry for entry in prior_split_data['cache']};prior_proofs={tablekey(entry):entry for entry in prior_split_report['tables']};assert set(prior_raw)==set(prior_proofs)
reused_tables=fresh_tables=fresh_directed_prices=0
seen=set();bandsets={name:set()for name in rows};outtables=[]
for entry in d['cache']:
 name=entry['profile'];lo,hi=map(F,entry['kappa_interval']);zlo,zhi=map(F,entry['shape_interval']);assert name in rows and 0<=lo<hi<=90 and F(2,5)<=zlo<zhi<=1
 box=(lo,hi,zlo,zhi);assert(name,*box)not in seen;seen.add((name,*box));bandsets[name].add(box);saved=list(map(F,entry['corner_uppers']))
 evaluation=candidates(lo,hi,zlo,zhi);assert entry['evaluation_points']==[list(map(str,x))for x in evaluation] and len(saved)==len(evaluation) and max(saved)==F(entry['upper'])
 key=tablekey(entry)
 if key in prior_raw and entry==prior_raw[key]:
  assert coefficients[name]==prior_split_report['expanded_coefficients'][name]
  assert bindings[name]['reserves']==prior_split_report['complete_whole_q_bindings'][name]['reserves'] and bindings[name]['target_spectral_gates']==prior_split_report['complete_whole_q_bindings'][name]['target_spectral_gates']
  assert bindings[name]['profile']['sha256']==prior_split_report['complete_whole_q_bindings'][name]['profile']['sha256']
  outtables.append(prior_proofs[key]);reused_tables+=1
  continue
 fresh_tables+=1;fresh_directed_prices+=len(evaluation)
 values=[]
 for j,(k,z)in enumerate(evaluation):
  value=point(name,k,z);assert value.below(saved[j]),(name,k,z,value.record(),saved[j])
  values.append({'kappa':str(k),'shape':str(z),'lower':value.lo.text(DOWN,90),'upper':value.hi.text(UP,90),'saved_upper':str(saved[j]),'signed_nuisance_covariance':str(rows[name]['p']['nuisance2']*(k/3-rows[name]['c']['theta']*(k+2)))})
 outtables.append({'profile':name,'kappa_interval':[str(lo),str(hi)],'shape_interval':[str(zlo),str(zhi)],'saved_upper':entry['upper'],'corners':values})
common=next(iter(bandsets.values()));assert all(v==common for v in bandsets.values())and len(outtables)==len(rows)*len(common)
used=set()
for item in d['refinements']:
 for child in item['children']:
  if child['status']in('EXACT_EMPTY_SOURCE_SUBCELL','EXACT_EMPTY_PHYSICAL_MOMENT_SUBCELL'):continue
  assert child['status']=='EXACT_DIRECTED_RANK_MOMENT_SUBCELL_PENDING_AUDIT'
  used.add(tuple(map(F,child['kappa_interval']+child['shape_interval'])))
assert used==common and len(common)==175 and len(outtables)==11025
assert {tuple(map(F,t['kappa_interval']+t['shape_interval']))for t in d['plane_cache']}==common
assert all(table['physical_constraint']=='4*kappa*zeta-kappa<=54'for table in d['plane_cache'])
assert reused_tables==10325 and fresh_tables==700
assert sha(source)==initial_sha
out={'status':'PASS_INDEPENDENT_NEW_ROUND_COMPLETE_R3_CACHE_PRICES','source':{'path':str(source.relative_to(ROOT)),'sha256':initial_sha},'input_sha256':initial_sha,'operator_b':str(b),'bits':BITS,'profiles':len(rows),'kappa_bands':len({(lo,hi)for lo,hi,zlo,zhi in common}),'moment_boxes':len(common),'corner_tables':len(outtables),'directed_corner_prices':sum(len(t['corners'])for t in outtables),'distinct_MPFR_prices':point.cache_info().currsize,'fresh_tables':fresh_tables,'fresh_directed_prices':fresh_directed_prices,'reused_tables':reused_tables,'accepted_identical_delta_cache':{'path':str(prior_split_report_path),'sha256':sha(prior_split_report_path),'source':prior_split_report['source'],'exact_raw_table_profile_reserve_coefficient_spectrum_equality':True},
 'parent':d['parent'],'complete_whole_q_bindings':bindings,'expanded_coefficients':coefficients,'expected_moment_boxes':[list(map(str,x))for x in sorted(common)],'tables':outtables,'verified_premise_file_count':len(verified),'verified_premise_files':verified,'backend':selftest(),'script_sha256':sha(Path(__file__)),'independent_R3_engine':{'path':str(Path(R.__file__)),'sha256':sha(Path(R.__file__))},
 'accepted_input_binding_reports':[{'path':str(x),'sha256':sha(x)}for x in(base_report_path,extra_report_path)],'exact_parent_cache_profile_mathematical_equality_to_accepted66':True,'profile_registry':'Each criterion uses its ancestor accepted name/profile SHA/exact reserve0/reserve6/targetb binding; profile SHA alone never selects a budget.','scope':'All11025 declared physical-curved rank cached moment constants and complete accepted whole-q profile/reserve bindings, at targetb=.711. Existing universal gates are reused; compact local leaves and old global geometry are not replayed. New geometry, selected cell offers and same-b lower branch must be separately composed for any global norm claim.',
 'band_coverage':'Every profile has exactly the same declared moment-box cache. For split sources this equals all nonempty child kappa/shape boxes; for nonsplit sources it covers every changed parent cumulant interval. Parent mixed intervals need not form one uniform grid. The separate geometry checker establishes exact partitions, physical clipping and per-cell choices.',
 'physical_boundary_proof':{'path':'new_round_2026-10-02/audit/LINEAR_TENSOR_SHAPE_CONSTRAINT.md','sha256':sha(ROOT/'new_round_2026-10-02/audit/LINEAR_TENSOR_SHAPE_CONSTRAINT.md')},'convexity':'Curved maxima use feasible rectangle corners and exact k*(4*z-1)=54 boundary intersections. A=q6+reserve6+zeta²>0,D=G-d0-zeta²>0 and194/sqrt10-633/15>0 give separate actual kappa/shape convexity. Theta is frozen toprofile kappa. Signed nuisance covariance and exact y(r) transport are included. Repaired highκ price uses reserve6=1/5000000; invalid original1e-8 premise is excluded.'}
args.output.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:out[k]for k in('status','source','profiles','kappa_bands','moment_boxes','corner_tables','directed_corner_prices','verified_premise_file_count','fresh_tables','fresh_directed_prices','reused_tables')},indent=2))
