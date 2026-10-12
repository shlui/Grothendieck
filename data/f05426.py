"""Copy exact accepted R3 branches into a path-independent replay bundle.

Accepted original files are read only. Copied independent source auditors are
changed only to resolve relative hash-ledger keys against the extracted root.
"""
from pathlib import Path
import hashlib,json,shutil,zipfile
ROOT=Path(__file__).resolve().parents[3]
DEST=Path(__file__).resolve().parent/'replay_accepted_R3'

def copy(rel):
 p=ROOT/rel;q=DEST/rel;q.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,q)

def main():
 DEST.mkdir(exist_ok=True)
 certificate='tonight/operator/projection/certificate/'
 for name in ['mpfr_interval.py','gaussian_max_certificate.py','certify_local_point.py','certify_local_bernstein_region.py','rational_phi_bounds.py','bernstein_integer_bounds.py','price_frozen_profile.py','certify_local_region.py']:
  copy(certificate+name)
 lower='tonight/operator/projection/lower_amp_certificate/'
 names=['general_R3_cover.py','general_R3_price.py','price_finer_support_cells.py','price_source_cells.py','price_general_source_polygons.py','price_b7171_native_remainder.py','price_b7171_native_remainder_extended.py','lower_source_domain_finer.json','lower_source_cells_tight.json','general_source_polygons.json','b7171_native_remainder_box.json','b7171_native_remainder_extended_box.json']
 for lab in ['68','70']:
  names += [f'profile_u.{lab}.json',f'tail_u{lab}_tight.json',f'general_profile_u{lab}_high_exact.json',f'general_u{lab}_high_price_and_tail.json']
  for part in ['body','exterior']:
   names += [f'{part}_u{lab}_finer.json.gz',f'{part}_u{lab}_tight.json.gz',f'general_u{lab}_high_{part}.json.gz']
 for name in names:copy(lower+name)
 diffuse='tonight/operator/projection/diffuse_certificate/'
 for name in ['price_critical_source_mask.py','price_critical_source_mask_b7171.py','critical_source_upper_masks.json','critical_source_upper_masks_b7171_newworst.json']:
  copy(diffuse+name)
 for label in ['interval_u7478_k0_01','interval_b7171_newworst']:
  for name in [f'profile_{label}.json',f'body_{label}.json.gz',f'exterior_{label}.json.gz',f'tail_{label}.json']:
   copy(diffuse+name)
 for rel in ['tonight/operator/cubic_L1_refinement/kappa_shape_caps.json','tonight/audit/mpfr_interval.py','tonight/audit/replay_R3_physical_q_trees.py','tonight/recovery/profile_cover/directed_cubic_integral.py','tonight/audit/lower_amp_R3/independent_source_price.py','tonight/audit/lower_amp_R3/COMPLETE_RESTRICTED_DOMAIN_REPLAY.json','tonight/audit/lower_amp_R3/PHYSICAL_Q_DOMAIN_REVIEW.md','tonight/audit/diffuse_R3/independent_diffuse_source_price.py','tonight/audit/refined_l1/replay_refined_prices.py','tonight/audit/refined_l1/independent_refined_l1_price_ledger.json','tonight/theory/replay_general_r3_cover.py','tonight/theory/audit_general_r3_parameters_tail.py','tonight/theory/audit_general_r3_source_bands.py','tonight/theory/general_r3_u70_whole_q_audit.json','tonight/theory/general_r3_u68_whole_q_audit.json','tonight/theory/general_r3_source_bands_audit.json','tonight/theory/GENERAL_R3_EXTENSION_AUDIT.md','tonight/theory/general_r3_extension.tex']:
  copy(rel)
 energy=ROOT/'reduced_search/energy'
 for p in [energy/'cubic_L1_certificate.json',*energy.glob('stein_L1_row_*.json')]:copy(str(p.relative_to(ROOT)))
 # Copied ledger keys become relative. The referenced immutable bytes do not
 # change. This is a storage-path adaptation, with no mathematical changes.
 for f in DEST.rglob('*.json'):
  d=json.loads(f.read_text())
  if 'sha256' in d and isinstance(d['sha256'],dict):
   values={}
   for key,value in d['sha256'].items():
    p=Path(key)
    relative=str(p.relative_to(ROOT))if p.is_absolute()else key
    assert hashlib.sha256((DEST/relative).read_bytes()).hexdigest()==value,(f,relative)
    values[relative]=value
   d['sha256']=values;f.write_text(json.dumps(d,indent=2)+'\n')
 # Relative-only adapters in copied independent auditors. Original producer
 # and scalar/local certificate code are byte-for-byte unchanged.
 f=DEST/'tonight/audit/lower_amp_R3/independent_source_price.py'
 s=f.read_text().replace('path=Path(filename)','path=ROOT/filename').replace("source['sha256'][str(profile_path)]","source['sha256'][str(profile_path.relative_to(ROOT))]")
 f.write_text(s)
 f=DEST/'tonight/audit/diffuse_R3/independent_diffuse_source_price.py'
 s=f.read_text().replace('file=Path(filename)','file=ROOT/filename');f.write_text(s)
 print(DEST)
if __name__=='__main__':main()
