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
assert d['status']=='COMPLETE_PARENT_LINKED_AFFINE_HYBRID_PENDING_INDEPENDENT_AUDIT'
parent=Path(d['parent']['path']);assert sha(parent)==d['parent']['sha256']=='322eb4c03bace797cf4f75377cd915979005a88959395ed5300ad6558b4d0eb6';parent_source=read(parent)
upper=Path(d['upper_parent']['path']);assert sha(upper)==d['upper_parent']['sha256'];upper_source=read(upper);assert parent_source['parent']==d['upper_parent']
base_source=Path(parent_source['base_global']['path']);assert sha(base_source)==parent_source['base_global']['sha256'];base_data=read(base_source)
base_report_path=HERE/'BANK43_LINEAR_TENSOR_R3_CACHE_REPLAY.json';extra_report_path=HERE/'BANK70_MONOTONE_R3_CACHE_REPLAY.json'
base_report=read(base_report_path);extra_report=read(extra_report_path)
for report in(base_report,extra_report):
 assert report['status']=='PASS_INDEPENDENT_NEW_ROUND_COMPLETE_R3_CACHE_PRICES'and report['bits']==768 and F(report['operator_b'])==b
 for path,digest in report['verified_premise_files'].items():assert sha(Path(path))==digest;verified[path]=digest
assert base_report['source']['sha256']==sha(base_source)and extra_report['source']['sha256']==sha(upper)
metadata=dict(base_data['profiles']);assert not set(metadata)&set(upper_source['profiles']);metadata.update(upper_source['profiles'])
accepted_bindings=dict(base_report['complete_whole_q_bindings']);assert not set(accepted_bindings)&set(extra_report['complete_whole_q_bindings']);accepted_bindings.update(extra_report['complete_whole_q_bindings'])
assert len(accepted_bindings)==70
componentpaths=[]; manifestpaths=[]
for manifest_record in d['extra_bindings']:
 new_manifest_path=Path(manifest_record['path']);assert sha(new_manifest_path)==manifest_record['sha256'];new_manifest=read(new_manifest_path);assert new_manifest['status']==manifest_record['status']
 manifestpaths.append(new_manifest_path)
 for item in new_manifest['profiles']:
  name=item['name'];assert name not in metadata;metadata[name]=item
  componentpath=Path(item['acceptance_path']);assert sha(componentpath)==item['acceptance_sha256'];component=read(componentpath);assert component['status'].startswith('PASS_INDEPENDENT_')and F(component['operator_b'])==b
  check_records(component,ROOT);componentpaths.append(componentpath)
  matching=[entry for entry in component['components']if(entry['profile']['sha256'],F(entry['reserve0']),F(entry['reserve6']),F(entry['operator_b']))==(item['profile_sha256'],F(item['reserve0']),F(item['reserve6']),b)]
  assert len(matching)==1;entry=matching[0]
  for record in entry['whole_q_reports']:
   proof=read(ROOT/record['path']);assert proof['status']==record['status']
   if record['part']!='tail':assert proof['profile_sha256']==item['profile_sha256']
  accepted_bindings[name]={'profile':{'path':item['profile_path'],'sha256':item['profile_sha256']},'kind':'rank','reserves':[entry['reserve0'],entry['reserve6']],'whole_q_premise':{'accepted_component':{'path':str(componentpath),'sha256':sha(componentpath)},'complete_component':entry}}
assert len(accepted_bindings)==76
reuse_report_path=HERE/'BANK73_AFFINE_HYBRID_R3_CACHE_REPLAY.json';reuse_report=read(reuse_report_path)
assert reuse_report['status']=='PASS_INDEPENDENT_NEW_ROUND_COMPLETE_R3_CACHE_PRICES'and reuse_report['bits']==768 and F(reuse_report['operator_b'])==b
reuse_source=ROOT/reuse_report['source']['path'];assert sha(reuse_source)==reuse_report['source']['sha256']=='ce49c57b2444742f75fa225f4bdc720116a4fa1e04cba41b0ef85486ac5d880e';reuse_data=read(reuse_source)
assert reuse_data['parent']==d['parent']and reuse_data['upper_parent']==d['upper_parent']and reuse_data['native_binding']==d['native_binding']
for path,digest in reuse_report['verified_premise_files'].items():assert sha(Path(path))==digest;verified[path]=digest
def tablekey(t):return(t['profile'],*t['kappa_interval'],*t['shape_interval'])
reuse_raw={tablekey(t):t for t in reuse_data['cache']};reuse_tables={tablekey(t):t for t in reuse_report['tables']}
assert len(reuse_raw)==len(reuse_tables)==1440
for name,prior in reuse_report['complete_whole_q_bindings'].items():
 assert name in accepted_bindings and prior['profile']==accepted_bindings[name]['profile']and prior['kind']==accepted_bindings[name]['kind']and prior['reserves']==accepted_bindings[name]['reserves']and prior['whole_q_premise']==accepted_bindings[name]['whole_q_premise']
used_profiles={entry['profile']for entry in d['cache']};assert len(used_profiles)==9 and {'full68','full70','highk7575'}<=used_profiles
rows={};bindings={};coefficients={}
for name in sorted(used_profiles):
 item=metadata[name];priorbinding=accepted_bindings[name];path=Path(item['profile_path']);assert sha(path)==item['profile_sha256']==priorbinding['profile']['sha256'];assert str(path)==priorbinding['profile']['path']
 raw=read(path);kind='full'if raw['kind']=='EXACT_RATIONAL_GENERAL_R3_PROFILE'else'rank';assert raw['kind']in('EXACT_RATIONAL_GENERAL_R3_PROFILE','EXACT_RATIONAL_R3_PROFILE')and kind==priorbinding['kind']
 assert (kind=='full')==(name in{'full68','full70','highk7575'})
 p=R.exact(raw['profile']);reserves=priorbinding['reserves'];binding=priorbinding['whole_q_premise']
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
 if name in base_report['expanded_coefficients']:assert {k:str(v)for k,v in c.items()}==base_report['expanded_coefficients'][name]and spectral==priorbinding['target_spectral_gates']
 coefficients[name]={k:str(v)for k,v in c.items()}
assert len(rows)==9
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
fresh_tables=reused_tables=fresh_prices=0
seen=set();bandsets={name:set()for name in rows};outtables=[]
for entry in d['cache']:
 name=entry['profile'];lo,hi=map(F,entry['kappa_interval']);zlo,zhi=map(F,entry['shape_interval']);assert name in rows and 0<=lo<hi<=90 and F(2,5)<=zlo<zhi<=1
 box=(lo,hi,zlo,zhi);assert(name,*box)not in seen;seen.add((name,*box));bandsets[name].add(box);saved=list(map(F,entry['corner_uppers']))
 evaluation=candidates(lo,hi,zlo,zhi);assert entry['evaluation_points']==[list(map(str,x))for x in evaluation] and len(saved)==len(evaluation) and max(saved)==F(entry['upper'])
 key=tablekey(entry)
 if key in reuse_raw:
  assert entry==reuse_raw[key]
  assert bindings[name]==reuse_report['complete_whole_q_bindings'][name]
  assert coefficients[name]==reuse_report['expanded_coefficients'][name]
  prior=reuse_tables[key];assert prior['saved_upper']==entry['upper']and len(prior['corners'])==len(evaluation)
  for v,(k,z),saved_upper in zip(prior['corners'],evaluation,saved):
   assert v['kappa']==str(k)and v['shape']==str(z)and F(v['saved_upper'])==saved_upper and F(v['upper'])<=saved_upper
  outtables.append(prior);reused_tables+=1;continue
 fresh_tables+=1;fresh_prices+=len(evaluation)
 values=[]
 for j,(k,z)in enumerate(evaluation):
  value=point(name,k,z);assert value.below(saved[j]),(name,k,z,value.record(),saved[j])
  values.append({'kappa':str(k),'shape':str(z),'lower':value.lo.text(DOWN,90),'upper':value.hi.text(UP,90),'saved_upper':str(saved[j]),'signed_nuisance_covariance':str(rows[name]['p']['nuisance2']*(k/3-rows[name]['c']['theta']*(k+2)))})
 outtables.append({'profile':name,'kappa_interval':[str(lo),str(hi)],'shape_interval':[str(zlo),str(zhi)],'saved_upper':entry['upper'],'corners':values})
common=next(iter(bandsets.values()));assert all(v==common for v in bandsets.values())and len(outtables)==len(rows)*len(common)
selected=[(item,child)for item in parent_source['refinements']for child in item['children']if child['status']=='EXACT_DIRECTED_RANK_MOMENT_SUBCELL_PENDING_AUDIT'and F(child['upper'])>F(40,71)]
selected_boxes={tuple(map(F,child['kappa_interval']+child['shape_interval']))for item,child in selected}
assert common==selected_boxes and len(common)==240 and len(selected)==d['selected_moment_children']==4124 and len(outtables)==2160 and reused_tables==1440 and fresh_tables==720
native_report_path=HERE/'CURVED_NATIVE_ENDPOINT_PRICE_REPLAY.json';native_report=read(native_report_path);native_path=Path(d['native_binding']['path'])
assert sha(native_path)==d['native_binding']['sha256']==native_report['source']['sha256']
assert native_report['status']=='PASS_INDEPENDENT_CURVED_NATIVE_SELECTED_ENDPOINT_PRICES'and F(native_report['operator_b'])==b and native_report['precision_bits']==768 and native_report['curved_uniform_tables']==469
for path,digest in native_report['verified_premise_files'].items():assert sha(Path(path))==digest;verified[path]=digest
assert 'plane_cache'not in d and 'source_plane_cache'not in d
assert sha(source)==initial_sha
out={'status':'PASS_INDEPENDENT_NEW_ROUND_COMPLETE_R3_CACHE_PRICES','source':{'path':str(source.relative_to(ROOT)),'sha256':initial_sha},'input_sha256':initial_sha,'operator_b':str(b),'bits':BITS,'profiles':len(rows),'kappa_bands':len({(lo,hi)for lo,hi,zlo,zhi in common}),'moment_boxes':len(common),'corner_tables':len(outtables),'fresh_tables':fresh_tables,'fresh_directed_prices':fresh_prices,'directed_corner_prices':sum(len(t['corners'])for t in outtables),'reused_tables':reused_tables,'reuse_premise':{'source':reuse_report['source'],'numeric_report':{'path':str(reuse_report_path),'sha256':sha(reuse_report_path)},'complete_raw_table_profile_reserve_and_coefficient_equality':True},'distinct_MPFR_prices':point.cache_info().currsize,
'parent':d['parent'],'upper_parent':d['upper_parent'],'complete_whole_q_bindings':bindings,'expanded_coefficients':coefficients,'profile_kinds':{name:rows[name]['kind']for name in rows},'six_target_b711_spectra':{name:bindings[name]['target_spectral_gates']for name in rows},'expected_moment_boxes':[list(map(str,key))for key in sorted(common)],'exact_selected_parent_moment_boxes_equal_cache':True,'selected_parent_moment_children':len(selected),'tables':outtables,'verified_premise_file_count':len(verified),'verified_premise_files':verified,'backend':selftest(),'script_sha256':sha(Path(__file__)),'independent_R3_engine':{'path':str(Path(R.__file__)),'sha256':sha(Path(R.__file__))},'accepted_input_binding_reports':[{'path':str(path),'sha256':sha(path)}for path in(base_report_path,extra_report_path)],'new_hybrid_component_manifests':[{'path':str(path),'sha256':sha(path)}for path in manifestpaths],'component_acceptances':[{'path':str(path),'sha256':sha(path)}for path in componentpaths],
'unchanged_native469_premise':{'input':d['native_binding'],'accepted_numeric_report':{'path':str(native_report_path),'sha256':sha(native_report_path)},'uniform_tables':469,'selected_endpoint_prices':13936,'source_bytes_exact':True,'old_kappa_zeta36_domain_valid_on_stronger54_subset_for_kappa_le90':True},'parent_uniform_source_planes_unchanged':True,'uniform_source_plane_parent':d['parent'],'physical_boundary_proof':{'path':'new_round_2026-10-02/audit/LINEAR_TENSOR_SHAPE_CONSTRAINT.md','sha256':sha(ROOT/'new_round_2026-10-02/audit/LINEAR_TENSOR_SHAPE_CONSTRAINT.md')},
'convexity':'Strong54 physical candidate sets rebuilt exactly. Full E00/E02/E22 are used for full68/full70/highk7575; rank matrices are used only for the six newly accepted rank criteria. Frozen theta, signed nuisance covariance and affine y(r) transfer are preserved. A,D>0 and all six target711 spectra pass. Radial coefficients are independently exact-bound; full68/full70 have negative curvature, handled by tangents in the separate accepted affine geometry theorem.',
'scope':'All2160 nine-criterion moment tables/240 selected parentmoment domains;1440 exact73 tables reused and720 fresh tables, accepted whole-q/reserve and unchanged native469 premises. New affine envelope geometry, native interpolation/choice and complete global composition are separate independent proof components; no scalar/local tree replay or producer imports.'}
args.output.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({key:out[key]for key in('status','source','profiles','profile_kinds','moment_boxes','corner_tables','fresh_directed_prices','verified_premise_file_count')},indent=2))
