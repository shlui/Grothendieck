#!/usr/bin/env python3
"""Exact uncapping of the fully replayed frozen .711 source cover.
The companion fresh geometry replays validate every local bound used here.
This script independently checks the parent-chain bindings and retains the
actual rational affine maxima instead of replacing them by40/71.
"""
from pathlib import Path
from fractions import Fraction as F
import gzip,json,hashlib,sys
R=Path('/workspace/scratch/8242d848994f'); O=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):
 p=Path(p)
 with(gzip.open(p,'rt')if p.suffix=='.gz'else p.open())as f:return json.load(f)
def bind(rec):
 p=Path(rec['path']); assert sha(p)==rec['sha256'];return read(p)
def ceilf(v,digits=15):
 a=10**digits;return F(-((-v.numerator*a)//v.denominator),a)
final_path=R/'new_round_2026-10-02/cubic_search/b711_bank89_final_moment2_goal_proposal.json.gz'
final=read(final_path); affine=bind(final['parent']); moment=bind(affine['parent']); rank=bind(moment['parent']); usplit=bind(rank['parent']); strong=bind(usplit['parent']); source=bind(strong['parent']); scan=bind(strong['protected_transport_scan'])
assert all(F(d['operator_b'])==F(711,1000)for d in(final,affine,moment,rank,usplit,strong))
assert final['moment_ancestor']==affine['parent'] and affine['upper_parent']==moment['parent'] and moment['u_parent']==rank['parent'] and rank['base_global']==usplit['parent']
assert scan['parent_sha256']==sha(strong['parent']['path'])
base_root=Path(str(strong['protected_transport_scan']['path']).split('/tonight/recovery/')[0]); CB=F(bind({'path':str(base_root/source['contact_gate_binding']['path']),'sha256':source['contact_gate_binding']['sha256']})['parameters']['B'])
b,b0=F(711,1000),F(source['operator_b']);betas={}
for n,q in scan['profile_beta_bindings'].items():
 p=base_root/q['profile_path'];assert sha(p)==q['profile_sha256'];beta=F(read(p)['profile']['beta']);bsq=beta*beta*(6 if n.startswith('R3:')else 1);assert bsq==F(q['beta_norm_squared']);betas[n]=(bsq,beta)
def canonical(s):
 while s.startswith('inherited:'):s=s[len('inherited:'):]
 return s.replace('native-contact:','native:',1).replace('R3-contact:','R3:',1)
leaves={}; ids={}
for i,row in enumerate(source['cells']):
 bsq,beta=betas[canonical(row['winner'])];rlo=F(row['r_interval'][0]);shift=bsq*(1/b-1/b0)+(b-b0)*rlo*rlo
 if 'contact:'in row['winner']:
  assert canonical(row['winner']).startswith('R3:');lam=F(row['contact_credit']['lambda_value']);assert lam>=0;shift+=(1/b-1/b0)*(6*lam*beta*CB+6*lam*lam*CB*CB/4)
 leaves[i,]=ceilf(F(row['cap'])+shift);ids[i,]='transported_713_parent'
assert len(leaves)==272579
for rec in strong['changes']:
 key=rec['parent_index'],;assert leaves[key]==F(rec['old_b711_offer']);v=F(rec['new_offer']['cap']);assert v<=leaves[key];leaves[key]=v;ids[key]='strong43'
counts={'protected_original_parents':len(leaves),'strong43_changed':len(strong['changes'])}
for rec in usplit['refinements']:
 key=rec['parent_index'],;old=leaves.pop(key);ids.pop(key);assert old==F(rec['old_offer']);cs=[]
 for ci,c in enumerate(rec['children']):
  nk=key+(ci,);assert nk not in leaves
  v=F(0)if c['status'].startswith('EXACT_EMPTY')else F(c['upper']);assert v<=old
  leaves[nk]=v;ids[nk]='u_split';cs.append(v)
 assert max(cs)==F(rec['new_parent_upper'])
counts['u_split_leaf_count']=len(leaves)
for rec in rank['changes']:
 key=rec['parent_index'],rec['child_index'];assert leaves[key]==F(rec['old_upper']);v=F(rec['new_offer'].get('upper',rec['new_offer'].get('cap')));assert v<=leaves[key];leaves[key]=v;ids[key]='rank70'
counts['rank70_changed']=len(rank['changes'])
for rec in moment['refinements']:
 key=rec['parent_index'],rec['child_index'];old=leaves.pop(key);ids.pop(key);assert old==F(rec['old_upper']);cs=[]
 for mi,c in enumerate(rec['children']):
  nk=key+(mi,);assert nk not in leaves
  v=F(0)if c['status'].startswith('EXACT_EMPTY')else F(c['upper']);assert v<=old
  leaves[nk]=v;ids[nk]='moment70';cs.append(v)
 assert max(cs)==F(rec['new_parent_upper'])
counts['moment_split_leaf_count']=len(leaves)
def affine_upper(proof):
 pieces=[]
 for p in proof['pieces']:
  vals=[]
  for reg in p['dominance_proof']['regions']:
   # Companion producer-free replay recomputes coefficients and these vertices.
   assert len(reg['vertices'])==len(reg['vertex_values']) and reg['vertices'];v=max(map(F,reg['vertex_values']));assert v==F(reg['upper']);vals.append(v)
  v=max(vals);assert v==F(p['upper']);pieces.append(v)
 return max(pieces)
uncapped=0
for rec in affine['changes']:
 key=rec['parent_index'],rec['child_index'],rec['moment_child_index'];old=leaves[key];assert old==F(rec['old_upper']);raw=affine_upper(rec['hybrid_proof']);v=min(old,raw)
 # Frozen verifier deliberately accepts threshold whenever raw is below it.
 recorded=min(old,F(40,71)if raw<=F(40,71)else ceilf(raw));assert recorded==F(rec['new_upper']);assert v<=recorded
 uncapped+=recorded!=v;leaves[key]=v;ids[key]='affine89_uncapped'
counts['affine_uncapped']=uncapped
for rec in final['refinements']:
 key=rec['parent_index'],rec['child_index'],rec['moment_child_index'];old=leaves.pop(key);ids.pop(key);assert old<=F(rec['old_upper']);cs=[]
 for di,c in enumerate(rec['children']):
  nk=key+(di,);assert nk not in leaves
  if c['status'].startswith('EXACT_EMPTY'):v=F(0)
  else:
   raw=affine_upper(c['hybrid_proof']);assert raw==F(c['raw_affine_upper']);v=min(old,raw);assert v<=F(c['upper'])
  leaves[nk]=v;ids[nk]='final8_uncapped';cs.append(v)
 assert max(cs)<=F(rec['new_parent_upper'])
counts['final_leaf_count']=len(leaves)
maximum=max(leaves.values());worst=[(k,str(v),ids[k])for k,v in leaves.items()if v==maximum];threshold=F(40,71);assert maximum<threshold
B=ceilf(maximum,15);assert maximum<=B<threshold
answer={'status':'PASS_EXACT_UNCAPPED_COMPLETE_PARENT_ASSEMBLY','operator_b':str(b),'source_sha256':sha(final_path),'counts':counts,'exact_upper_maximum':str(maximum),'rational_uniform_cap':str(B),'rational_cap_decimal':format(float(B),'.17g'),'strict_margin':str(threshold-B),'strict_margin_decimal':float(threshold-B),'reciprocal_lower':str(1/B),'reciprocal_decimal':float(1/B),'worst':worst,'threshold_equal_leaf_count':sum(v==threshold for v in leaves.values()),'strict_all_leaves':all(v<threshold for v in leaves.values()),'scope':'Exact uncapping and complete parent assembly; valid after companion fresh producer-free geometry replays and inherited analytic/numeric premises.'}
(O/'STRICT_UNCAPPED_ASSEMBLY.json').write_text(json.dumps(answer,indent=2)+'\n');print(json.dumps(answer,indent=2))
