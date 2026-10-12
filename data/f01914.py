#!/usr/bin/env python3
"""Fresh exact same-coefficient outer glue at a smaller cubic coefficient.

The already accepted complete outer partition remains physical. Native
cost changes are beta_norm^2*(1/b-1/b0)+(b-b0)*r^2; the rank70 remainder
uses 6*beta_raw^2 in the constant. Spectral gates strengthen for b<=b0.
The uniform beta bound gives a simple complete outer cap, with no new
scalar integration or interpolation assumptions.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse,json,sys
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/'tonight/audit'))
import replay_native_R3_unions as N
sys.path.insert(0,str(Path(__file__).resolve().parent))
import replay_R3_corners as R
ap=argparse.ArgumentParser();ap.add_argument('--b',type=F,required=True);ap.add_argument('--output',type=Path,required=True)
args=ap.parse_args();b=args.b;b0=F(7171,10000);assert 0<b<=b0
base_cap=F(28169,50000)
bindings={};native_profiles={};spectral={}
def accepted(path):
 d=json.loads((ROOT/path).read_text());assert d['status'].startswith('PASS')
 bindings[path]=dict(sha256=N.sha(ROOT/path),status=d['status'])
 return d
old=accepted('tonight/audit/global_b7171/STRENGTHENED_GLOBAL_B7171_ACCEPTANCE.json')
assert old['operator_b']==str(b0) and old['historical_1773_or_1780498_premises_used'] is False
first=accepted('tonight/audit/native_b7171_union/first_four_transfers_independent.json')
transferpath='tonight/recovery/profile_cover/native_b7171_regional_transfer_assembly.json'
assert N.sha(ROOT/transferpath)==first['input_sha256']
transfers=json.loads((ROOT/transferpath).read_text())['rows'][:4]
def profile(name):
 if name in native_profiles:return
 d=N.package(name);p=d['profile'];M=d['M'];beta=F(p['beta'])
 assert beta*beta<=F(1,16)
 for a,o,c in ((M[0][0]+1,M[0][2],M[2][2]-b),
               (M[1][1]+1,M[1][3],M[3][3]-b),
               (M[1][1]-b,M[1][3],M[3][3])):N.psd(a,o,c)
 native_profiles[name]=dict(beta_normalized=str(beta),beta_squared=str(beta*beta),
                            manifest_sha256=d['manifest_sha256'],
                            source_sha256=d['source_sha256'],replay_sha256=d['replay_sha256'])
for row,report in zip(transfers,first['transfers']):
 assert row['package']==report['package'] and F(row['directed_total_upper'])==F(report['total'])<base_cap
 profile(row['package'])
regions=[
 ('tonight/audit/native_b7171_union/u600_682_independent_replay.json',F(3,5),F(341,500)),
 ('tonight/audit/native_b7171_union/u682_6905_independent_replay.json',F(341,500),F(1381,2000)),
 ('tonight/audit/native_b7171_union/u775_c_independent_replay.json',F(31,40),F(159577,200000))]
region_results=[];total_native=total_R3=0
for path,lo,hi in regions:
 report=accepted(path);assert F(report['b'])==b0 and list(map(F,report['interval']))==[lo,hi]
 assert F(report['maximum_rounded_cap'])<base_cap
 source=ROOT/report['path'];assert N.sha(source)==report['input_sha256']
 data=json.loads(source.read_text());frames={v['mask']:v for v in data['bindings']};count=other=0
 for cell in data['cells']:
  ulo,uhi=map(F,cell['u_interval'])
  if not lo<=ulo<uhi<=hi:continue
  if cell['accepted']:
   assert F(cell['cap'])<base_cap
   profile(frames[cell['frame']]['package']);count+=1
  else:
   assert path==regions[1][0];other+=1
 assert count==report['native_cells'] and other==report['R3_cells']
 total_native+=count;total_R3+=other
 region_results.append(dict(interval=list(map(str,(lo,hi))),native_cells=count,R3_cells=other,original_accepted_report=path))
assert total_native==45025 and total_R3==41
rp=ROOT/'tonight/operator/projection/lower_amp_certificate/profile_u.70.json'
box=accepted('tonight/audit/diffuse_b7171_R3/independent_extended_remainder_box.json')
assert N.sha(rp)==box['profile_sha256'] and F(box['cap'])<base_cap
p=R.exact(json.loads(rp.read_text())['profile']);rawbeta=p['beta'];assert 6*rawbeta*rawbeta<=F(1,16)
row=dict(kind='rank',p=p)
rank_spectral=R.spectral(row,b)
energy_end=F(21,100);energy=energy_end**2+b*F(847,1000)**2
assert first['coverage']==['0','319/500'] and energy<base_cap
assert list(map(F,first['transfers'][-1]['amplitude_interval']))[1]==F(319,500)
shift=(1/b-1/b0)/16;cap=base_cap+shift
out=dict(status='PASS_INDEPENDENT_COMPLETE_SMALLER_B_OUTER_SOURCE_TRANSPORT',
 operator_b=str(b),baseline_b=str(b0),outer_norm_upper_bound=str(cap),
 outer_norm_upper_display=float(cap),baseline_uniform_outer_cap=str(base_cap),uniform_beta_squared_bound='1/16',
 exact_uniform_completion_shift=str(shift),
 complete_outer_amplitude_coverage=[['0','1381/2000'],['31/40','sqrt(2/pi)']],
 elementary_low_u_end=str(energy_end),elementary_energy_cap=str(energy),
 accepted_native_cells_transported=total_native,accepted_R3_cells_transported=total_R3,
 first_four_complete_regions_transported=4,native_profiles=native_profiles,
 rank70_profile=dict(path=str(rp.relative_to(ROOT)),sha256=N.sha(rp),beta_raw=str(rawbeta),
                     normalized_beta_squared=str(6*rawbeta*rawbeta),target_b_spectral_gates=rank_spectral),
 region_results=region_results,bindings=bindings,
 analytic_review='The same frozen local potentials remain valid. Reducing b strengthens every spectral inequality. Native beta_norm squared and rank70 six times beta_raw squared are each at most1/16. Their exact constant shift is at most(1/b-1/b0)/16; the additional(b-b0)*r^2 is nonpositive and may be dropped uniformly. Every old complete outer-cell upper bound and first-four regional total is independently checked below .56338, so the new complete outer cap is .56338 plus the recorded shift. Physical source coverage does not depend on the operator coefficient.',
 scope='Complete outer-source theorem only, at the displayed new coefficient. A global bound additionally needs a complete middle proof at EXACTLY that coefficient and a common cap exceeding this outer cap.')
args.output.write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:v for k,v in out.items()if k not in('bindings','native_profiles','rank70_profile')},indent=2))
