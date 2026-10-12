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
assert d['status']=='COMPLETE_MONOTONE_RANK_CHILD_UPGRADE_PENDING_INDEPENDENT_AUDIT'
is_split='refinements' in d
assert not is_split
parent=Path(d['parent']['path']);assert sha(parent)==d['parent']['sha256'];parent_source=read(parent)
base_record=d['base_global']
base_source=Path(base_record['path']);assert sha(base_source)==base_record['sha256'];oldsource=read(base_source);assert oldsource['physical_constraint']=='4*kappa*zeta-kappa<=54'
assert parent_source['parent']==d['base_global']
cachepath=BASE/'tonight/audit/global_b7171/full_upper57_twentysevenbank_R3_corner_replay.json';prior=read(cachepath);assert prior['status']=='PASS_INDEPENDENT_GENERIC_COMPLETE_R3_KAPPA_CACHE'
oldbindings=prior['complete_whole_q_bindings'];scan=read(BASE/'tonight/recovery/profile_cover/assigned_twenty-nine_coefficient_scan.json')
componentpaths=[ROOT/'new_round_2026-10-02/audit'/part/'COMPONENT_ACCEPTANCE.json'for part in ('kappa13_q25_source_box','low_r_highk_source_box','moderate_crossover_source_box','highu_kappa10_source_box','highu_kappa4_source_box','kappa28_source_box','secondary_kappa066_source_box','diffuse_kappa03_source_box','kappa43_source_box')]
components={}
for path in componentpaths:
 data=read(path);assert data['status'].startswith('PASS_INDEPENDENT_')and F(data['operator_b'])==b
 check_records(data,ROOT)
 for item in data['components']:
  key=item['profile']['sha256'];assert key not in components;components[key]=(item,path,data)
rows={};bindings={};coefficients={}
for name,item in d['profiles'].items():
 path=Path(item['profile_path']);assert sha(path)==item['profile_sha256'];acceptedpath=Path(item['acceptance_path']);assert sha(acceptedpath)==item['acceptance_sha256'];accepted=read(acceptedpath)
 raw=read(path);p=R.exact(raw['profile']);kind='full'if raw['kind']=='EXACT_RATIONAL_GENERAL_R3_PROFILE'else'rank';assert raw['kind']in('EXACT_RATIONAL_GENERAL_R3_PROFILE','EXACT_RATIONAL_R3_PROFILE')
 if name in oldbindings:
  old=oldbindings[name];check_records(old,BASE)
  info=old.get('frozen_gate_binding',old);record=info['profile'];assert record['sha256']==sha(path)
  reserves=info.get('reserves')or info.get('target_reserves')or info.get('global_reserves');assert reserves
  binding={'accepted_prior_cache':{'path':str(cachepath),'sha256':sha(cachepath)},'complete_frozen_gate_binding':old}
 elif 'R3:'+name in scan['profile_beta_bindings']:
  assert accepted['status'].startswith('PASS_INDEPENDENT_COMPLETE_')and accepted['profile_sha256']==sha(path)
  reserves=accepted['reserves'];check_records(accepted,BASE)
  for part,rec in accepted['reports'].items():
   proof=read(BASE/rec['path']);assert proof['status']==rec['status']
   if part!='tail':assert proof['profile_sha256']==sha(path)
  assert accepted['q_cover'][0]==['0','8']and accepted['q_cover'][-1][-1]=='infinity'
  binding={'accepted_whole_q_gate':{'path':str(acceptedpath),'sha256':sha(acceptedpath)},'complete_accepted_gate':accepted}
 else:
  entry,componentpath,component=components[item['profile_sha256']]
  assert entry['profile']['sha256']==sha(path)and F(entry['operator_b'])==b
  reserves=[entry['reserve0'],entry['reserve6']]
  for rec in entry['whole_q_reports']:
   proof=read(ROOT/rec['path']);assert proof['status']==rec['status']
   if rec['part']!='tail':assert proof['profile_sha256']==sha(path)
  binding={'accepted_component':{'path':str(componentpath),'sha256':sha(componentpath)},'complete_component':entry,'producer_declared_premise_status':accepted['status']}
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
 coefficients[name]={k:str(v)for k,v in c.items()}
assert len(rows)==len(d['profiles'])==9
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
delta_parent_path=HERE/'BANK50_MONOTONE_R3_CACHE_REPLAY.json';delta_parent=read(delta_parent_path)
assert delta_parent['status']=='PASS_INDEPENDENT_NEW_ROUND_COMPLETE_R3_CACHE_PRICES' and delta_parent['bits']==768
old_delta_source_path=ROOT/delta_parent['source']['path'];assert sha(old_delta_source_path)==delta_parent['source']['sha256'];old_delta_source=read(old_delta_source_path)
assert old_delta_source['parent']==d['parent'] and old_delta_source['base_global']==d['base_global']
def delta_key(e):return(e['profile'],*map(F,e['kappa_interval']),*map(F,e['shape_interval']))
old_raw={delta_key(e):e for e in old_delta_source['cache']};old_proofs={delta_key(e):e for e in delta_parent['tables']};assert set(old_raw)==set(old_proofs)
metadata_promotions=[]
for name,item in old_delta_source['profiles'].items():
 current=d['profiles'][name]
 differences={k:(v,current.get(k))for k,v in item.items()if v!=current.get(k)}
 assert set(current)==set(item) and set(differences)<={'acceptance_path','acceptance_sha256','proof_status'}
 assert current['profile_path']==item['profile_path']and current['profile_sha256']==item['profile_sha256']
 if differences:
  assert current['proof_status']=='INDEPENDENT_WHOLE_Q_SOURCE_AND_TARGET_B711_ACCEPTED'
  verified_component=bindings[name]['whole_q_premise']['accepted_component']
  assert current['acceptance_path']==verified_component['path']and current['acceptance_sha256']==verified_component['sha256']
  metadata_promotions.append({'profile':name,'changes':differences,'mathematical_profile_reserve_and_cache_unchanged':True})
 assert coefficients[name]==delta_parent['expanded_coefficients'][name]
 assert bindings[name]['reserves']==delta_parent['complete_whole_q_bindings'][name]['reserves'] and bindings[name]['target_spectral_gates']==delta_parent['complete_whole_q_bindings'][name]['target_spectral_gates']
seen=set();bandsets={name:set()for name in rows};outtables=[];fresh_tables=reused_tables=fresh_directed_prices=0

for entry in d['cache']:
 name=entry['profile'];lo,hi=map(F,entry['kappa_interval']);zlo,zhi=map(F,entry['shape_interval']);assert name in rows and 0<=lo<hi<=90 and F(2,5)<=zlo<zhi<=1
 box=(lo,hi,zlo,zhi);assert(name,*box)not in seen;seen.add((name,*box));bandsets[name].add(box);saved=list(map(F,entry['corner_uppers']));assert entry.get('physical_kappa_shape_constraint',True) is True
 evaluation=candidates(lo,hi,zlo,zhi);assert entry['evaluation_points']==[list(map(str,x))for x in evaluation] and len(saved)==len(evaluation) and max(saved)==F(entry['upper'])
 if delta_key(entry)in old_raw:
  assert entry==old_raw[delta_key(entry)]
  outtables.append(old_proofs[delta_key(entry)]);reused_tables+=1
  continue
 fresh_tables+=1;fresh_directed_prices+=len(evaluation)
 values=[]
 for j,(k,z)in enumerate(evaluation):
  value=point(name,k,z);assert value.below(saved[j]),(name,k,z,value.record(),saved[j])
  values.append({'kappa':str(k),'shape':str(z),'lower':value.lo.text(DOWN,90),'upper':value.hi.text(UP,90),'saved_upper':str(saved[j]),'signed_nuisance_covariance':str(rows[name]['p']['nuisance2']*(k/3-rows[name]['c']['theta']*(k+2)))})
 outtables.append({'profile':name,'kappa_interval':[str(lo),str(hi)],'shape_interval':[str(zlo),str(zhi)],'saved_upper':entry['upper'],'corners':values})
assert set(old_raw)<=seen and reused_tables==len(old_raw)==476
common=next(iter(bandsets.values()));assert all(v==common for v in bandsets.values())and len(outtables)==len(rows)*len(common)
expected={(tuple(map(F,item['kappa_interval']))+(F(2,5),F(1)))for item in parent_source['refinements']if any(not child['status'].startswith('EXACT_EMPTY')and F(child['upper'])>F(40,71)for child in item['children'])}
assert common==expected and len(common)==68
parent_cache_path=HERE/'BANK43_LINEAR_TENSOR_R3_CACHE_REPLAY.json';accepted_parent=read(parent_cache_path)
assert accepted_parent['status']=='PASS_INDEPENDENT_NEW_ROUND_COMPLETE_R3_CACHE_PRICES'and accepted_parent['source']['sha256']==d['base_global']['sha256']
assert accepted_parent['profiles']==43 and accepted_parent['corner_tables']==3010 and accepted_parent['bits']==768
for path,digest in accepted_parent['verified_premise_files'].items():assert sha(Path(path))==digest
for name,binding in accepted_parent['complete_whole_q_bindings'].items():
 assert sha(Path(binding['profile']['path']))==binding['profile']['sha256']
assert not set(bindings)&set(accepted_parent['complete_whole_q_bindings'])
base_cache={ (e['profile'],*map(F,e['kappa_interval']),*map(F,e['shape_interval'])):e for e in oldsource['cache']}
accepted_tables={ (e['profile'],*map(F,e['kappa_interval']),*map(F,e['shape_interval'])):e for e in accepted_parent['tables']}
assert set(base_cache)==set(accepted_tables)
for key,entry in base_cache.items():
 table=accepted_tables[key]
 assert entry['upper']==table['saved_upper'] and list(map(F,entry['corner_uppers']))==[F(v['saved_upper'])for v in table['corners']]
assert sha(source)==initial_sha
out={'status':'PASS_INDEPENDENT_NEW_ROUND_COMPLETE_R3_CACHE_PRICES','source':{'path':str(source.relative_to(ROOT)),'sha256':initial_sha},'input_sha256':initial_sha,'operator_b':str(b),'bits':BITS,'profiles':len(rows),'kappa_bands':len({(lo,hi)for lo,hi,zlo,zhi in common}),'moment_boxes':len(common),'corner_tables':len(outtables),'directed_corner_prices':sum(len(t['corners'])for t in outtables),'distinct_MPFR_prices':point.cache_info().currsize,'accepted_metadata_promotions':metadata_promotions,
 'parent':d['parent'],'accepted_identical_delta_cache':{'path':str(delta_parent_path),'sha256':sha(delta_parent_path),'source':delta_parent['source'],'exact_raw_table_equality':True},'fresh_tables':fresh_tables,'fresh_directed_prices':fresh_directed_prices,'reused_tables':reused_tables,'accepted_base_cache':{'path':str(parent_cache_path),'sha256':sha(parent_cache_path),'source':accepted_parent['source'],'unchanged_cache_tables':3010,'input_equality_checked':True},'base_global':d['base_global'],'composed_profiles':52,'complete_whole_q_bindings':bindings,'expanded_coefficients':coefficients,'expected_moment_boxes':[list(map(str,x))for x in sorted(common)],'tables':outtables,'verified_premise_file_count':len(verified),'verified_premise_files':verified,'backend':selftest(),'script_sha256':sha(Path(__file__)),'independent_R3_engine':{'path':str(Path(R.__file__)),'sha256':sha(Path(R.__file__))},
 'component_acceptances':[{'path':str(path),'sha256':sha(path)}for path in componentpaths],'scope':'All nine newly added physical-curved profile cached moment constants and complete accepted whole-q profile/reserve bindings, at targetb=.711. Existing universal gates are reused; compact local leaves and old global geometry are not replayed. New geometry, selected cell offers and same-b lower branch must be separately composed for any global norm claim.',
 'band_coverage':'Every new profile has exactly all68 still-failing nonempty parent-child moment bands. Accepted3010 base43 tables are retained by exact ancestral source/input equality. Every profile has exactly the same declared moment-box cache. For split sources this equals all nonempty child kappa/shape boxes; for nonsplit sources it covers every changed parent cumulant interval. Parent mixed intervals need not form one uniform grid. The separate geometry checker establishes exact partitions, physical clipping and per-cell choices.',
 'physical_boundary_proof':{'path':'new_round_2026-10-02/audit/LINEAR_TENSOR_SHAPE_CONSTRAINT.md','sha256':sha(ROOT/'new_round_2026-10-02/audit/LINEAR_TENSOR_SHAPE_CONSTRAINT.md')},'convexity':'Curved maxima use feasible rectangle corners and exact k*(4*z-1)=54 boundary intersections. A=q6+reserve6+zeta²>0,D=G-d0-zeta²>0 and194/sqrt10-633/15>0 give separate actual kappa/shape convexity. Theta is frozen toprofile kappa. Signed nuisance covariance and exact y(r) transport are included. Repaired highκ price uses reserve6=1/5000000; invalid original1e-8 premise is excluded.'}
args.output.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:out[k]for k in('status','source','profiles','fresh_tables','fresh_directed_prices','reused_tables','corner_tables','directed_corner_prices','verified_premise_file_count')},indent=2))
