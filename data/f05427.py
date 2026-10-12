"""Build a new relative-path package; never alter earlier published ZIPs."""
from pathlib import Path
import json, shutil, hashlib
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
DEST=HERE/'replay_accepted_R3_expanded'
OLD=HERE/'replay_accepted_R3_with_targeted'
SPECS=[
 ('interval_robust_moderatek_box','replay_R3_robust_current_q_trees.py','robustcurrent','robust_moderatek_R3_recheck','moderatek',20),
 ('interval_primary_u7425_k4_10','replay_R3_highk_q_trees.py','highk','primary_u7425_R3_recheck','primary',32),
 ('interval_lowerface_u7175_k053_08','replay_R3_lowerface_q_trees.py','lowerface','lowerface_u7175_R3_recheck','lowerface',32),
 ('interval_nextface_u7375_k4_10','replay_R3_nextface_q_trees.py','nextface','nextface_u7375_R3_recheck','nextface',32),
 ('interval_lowerdiffuse_u72125_k0_01','replay_R3_lowerface_q_trees.py','lowerface','lowerdiffuse_u72125_R3_recheck','lowerdiffuse',32),
]

def main():
 if DEST.exists(): shutil.rmtree(DEST)
 shutil.copytree(OLD,DEST,ignore=shutil.ignore_patterns('__pycache__','replay_outputs','replay_outputs_targeted'))
 def copy(rel):
  source=ROOT/rel; target=DEST/rel; target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(source,target)
 target='tonight/operator/projection/targeted_certificate/'
 for lab,checker,branch,audit,ledger,R in SPECS:
  for stem,suffix in [('profile_','.json'),('body_','.json.gz'),('exterior_','.json.gz'),('tail_','.json')]: copy(target+stem+lab+suffix)
  for f in (ROOT/'tonight/audit'/audit).glob('*.json'):copy(str(f.relative_to(ROOT)))
  copy('tonight/audit/'+checker)
  copy(target+ledger+'_source_box.json')
  copy(target+'price_'+ledger+'_source_box.py')
 shutil.copy2(HERE/'verify_expanded_R3_bundle.py',DEST/'verify_expanded.py')
 (DEST/'EXPANDED_SPECS.json').write_text(json.dumps([dict(label=lab,checker=checker,branch=branch,audit=audit,ledger=ledger,R=R)for lab,checker,branch,audit,ledger,R in SPECS],indent=2))
 with (DEST/'README.md').open('a')as f:f.write('''\nAdditional accepted gates: moderate15, primary18, lower19, next20, and lower-diffuse21.\nRun `python3 verify_expanded.py` to freshly replay all five complete q partitions\nand exact tails, and re-evaluate all eighty restricted-source corners with the\nindependently derived expanded price formula. Run `python3 verify_expanded.py --all`\nto replay the earlier targeted gates too. All paths are relative; no optimizer,\narchive, network, or original workspace is used. These restricted source boxes\nalone do not assert complete global source coverage or a new lower endpoint.\n''')
 excluded={'MANIFEST.json','COMPLETE_PORTABLE_REPLAY.json','VALIDATION.json','tonight/audit/global_b7171/targeted_32_outer_vertices.json','tonight/audit/global_b7171/robust_32_outer_vertices.json'}
 hashes={str(f.relative_to(DEST)):hashlib.sha256(f.read_bytes()).hexdigest()for f in DEST.rglob('*')if f.is_file()and'__pycache__'not in f.parts and str(f.relative_to(DEST))not in excluded}
 (DEST/'MANIFEST.json').write_text(json.dumps(dict(kind='RELATIVE_SHA256_ACCEPTED_RESTRICTED_R3_EXPANDED_INPUTS',sha256=hashes),indent=2))
 print(DEST,len(hashes))
if __name__=='__main__':main()
