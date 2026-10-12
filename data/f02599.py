#!/usr/bin/env python3
"""Exact coverage ledger: complete amplitude theorems only.

Conditional R3/source cells and numerical discovery are deliberately excluded.
Each package below has a complete optimizer-free local/scalar/outer replay.
"""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json
from mpfr_interval import I
HERE=Path(__file__).resolve().parent;BASE=HERE.parents[1]
ROOT=BASE/'tonight/regional_certificates'
accepted_names=['u014_042','u021_030','u030_031','u030_0485','u040_054','u050_0638','u060_0662','u0662_0682','u077_c','u079_0797885']
rho=F(847,1000);b=F(18,25);target=F(40,71);fallback=F(451,800);ceiling=F(159577,200000)
assert (I(2)/I.pi()).sqrt().below(ceiling)
entries=[];pending=[]
for name in accepted_names:
 path=ROOT/name;rp=path/'replay.json';mp=path/'data/assembly.json'
 if not rp.exists():pending.append({'package':name,'reason':'No completed portable replay'});continue
 replay=json.loads(rp.read_text());meta=json.loads(mp.read_text())
 assert replay['status']=='PASS_COMPLETE_REGIONAL_BOUND' and meta['operator']=='P1-(18/25)P3'
 assert not any(k in meta for k in ('kappa_interval','source_interval','conditional_on','discovery_only'))
 lo,hi=map(F,meta.get('coverage_interval',meta['amplitude_interval']));cap=F(meta['regional_cap'])
 assert list(map(F,replay['amplitude_interval']))==[lo,hi]
 assert F(replay['directed_total_upper'])<cap<=F(meta.get('target_threshold','40/71'))
 entries.append({'kind':'complete_amplitude_package','name':name,'interval':[str(lo),str(hi)],'cap':str(cap),'local_boxes':replay['local_boxes_replayed'],'scalar_leaves':replay['scalar_leaves_replayed'],'replay_sha256':hashlib.sha256(rp.read_bytes()).hexdigest(),'manifest_sha256':hashlib.sha256((path/'SHA256.json').read_bytes()).hexdigest()})
energycap=F(21,100)**2+b*rho*rho;assert energycap<target
entries.append({'kind':'universal_energy_corollary','name':'energy_u0_021','interval':['0','21/100'],'cap':str(energycap),'rho':str(rho),'dependency':'Ten fresh whole-line comparisons and complete 2780-leaf global cubic L1 price cover; pure-third projection duality.'})
# The accepted lower package uses the same local/scalar data on the wider
# interval. The independently proved outer endpoint extension is retained
# separately under the weaker threshold, without changing package metadata.
extpath=HERE/'region_u6620_6850_fallback_independent_endpoint_ledger.json'
if extpath.exists():
 ext=json.loads(extpath.read_text());assert ext['status']=='PASS_INDEPENDENT_REGIONAL_ENDPOINT_ASSEMBLY' and F(ext['total_upper'])<fallback
 entries.append({'kind':'exact_outer_extension_corollary','name':'lower662_685_fallback','interval':['331/500','137/200'],'cap':str(fallback),'ledger_sha256':hashlib.sha256(extpath.read_bytes()).hexdigest(),'dependency':'Completed u0662_0682 package unchanged plus independent exact endpoint extension.'})
def coverage(threshold):
 selected=[(F(x['interval'][0]),F(x['interval'][1]),x['name'])for x in entries if F(x['cap'])<=threshold]
 selected.sort();merged=[]
 for lo,hi,name in selected:
  if merged and lo<=merged[-1][1]:merged[-1][1]=max(merged[-1][1],hi);merged[-1][2].append(name)
  else:merged.append([lo,hi,[name]])
 gaps=[];last=F(0)
 for lo,hi,names in merged:
  if lo>last:gaps.append([str(last),str(lo)])
  last=max(last,hi)
 if last<ceiling:gaps.append([str(last),str(ceiling)])
 return {'threshold':str(threshold),'complete':not gaps,'closed_covered_intervals':[[str(lo),'sqrt(2/pi)' if hi>=ceiling else str(hi)]for lo,hi,_ in merged],'open_uncovered_intervals':gaps,'participating_theorems':len(selected)}
result={'status':'PASS_EXACT_AMPLITUDE_COVERAGE_LEDGER','operator':'P1-(18/25)P3','physical_domain':['0','sqrt(2/pi)'],'rho_for_energy':str(rho),'entries':entries,'target_1775':coverage(target),'fallback_800_451':coverage(fallback),'pending_packages':pending,'excluded_scopes':['Conditional R3 cells restrict r and/or kappa and do not count as whole-amplitude coverage.','Conditional complement masks require a full multidimensional union proof.','Sampled or optimizer discovery profiles never count as acceptance.'],'retained_project_premises':{'lower':'1.773, conditional on unreplayed historical55-interval norm cover.','upper':'1.780498480511303, an unreplayed earlierprojectrounding premise.'},'global_conclusion':'Neither threshold has complete all-amplitude coverage; no new globalGrothendieck endpoint follows.'}
(HERE/'amplitude_coverage_ledger.json').write_text(json.dumps(result,indent=2)+'\n')
lines=['# Exact amplitude coverage ledger','',f'For T=P1-(18/25)P3, u ranges from 0 to sqrt(2/pi). The fresh global cubic moment cap used for the elementary energy region is rho=847/1000.','', '| Threshold | Complete covered amplitude intervals | Remaining open gaps |','|---|---|---|']
for key in ('target_1775','fallback_800_451'):
 x=result[key];lines.append('| '+x['threshold']+' | '+', '.join('['+a+', '+c+']'for a,c in x['closed_covered_intervals'])+' | '+', '.join('('+a+', '+c+')'for a,c in x['open_uncovered_intervals'])+' |')
lines += ['', 'Only complete amplitude theorems and exact analytic extensions are included. Conditional R3/source strips and complement masks are excluded from this union until their full multidimensional coverage is certified.','', 'The project values1.773 and1.780498480511303 remain unreplayed earlier premises. These regional proofs do not establish either a new global lower endpoint or a three-digit determination.']
(HERE/'AMPLITUDE_COVERAGE_LEDGER.md').write_text('\n'.join(lines)+'\n')
print(json.dumps({k:result[k]for k in ('status','target_1775','fallback_800_451')},indent=2))
