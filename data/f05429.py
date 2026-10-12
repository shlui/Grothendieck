"""Extend the frozen portable R3 package without changing its original ZIP."""
from pathlib import Path
import json,shutil,hashlib
ROOT=Path(__file__).resolve().parents[3];HERE=Path(__file__).resolve().parent
DEST=HERE/'replay_accepted_R3_with_targeted'
OLD=HERE/'replay_accepted_R3'

def copy(rel):
 p=ROOT/rel;q=DEST/rel;q.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,q)

def main():
 if DEST.exists():shutil.rmtree(DEST)
 shutil.copytree(OLD,DEST,ignore=shutil.ignore_patterns('__pycache__','replay_outputs'))
 target='tonight/operator/projection/targeted_certificate/'
 labels=['interval_lowr_u76_r16_k0_2','interval_middle_u7475_r267_k08_1','interval_middlehigh_u7475_r275_k11_2','interval_robust_middlehigh_box','interval_robust_currentworst_box','interval_robust_broadcurrent_box']
 for lab in labels:
  for stem,suffix in [('profile_','.json'),('body_','.json.gz'),('exterior_','.json.gz'),('tail_','.json')]:copy(target+stem+lab+suffix)
 for name in ['general_profile_highk_u7575_r259_k2_4.json','general_highk_body.json.gz','general_highk_exterior16_loose.json.gz','general_highk_tail16_loose.json','targeted_source_boxes.json','robust_source_boxes.json','price_targeted_source_boxes.py','price_robust_source_boxes.py','broadcurrent_source_box.json','price_broadcurrent_source_box.py']:
  copy(target+name)
 for name in ['replay_R3_generous_q_trees.py','replay_targeted_rank_q_trees.py','replay_R3_robust_q_trees.py','replay_R3_robust_current_q_trees.py']:
  copy('tonight/audit/'+name)
 for folder in ['targeted_lowr_R3_recheck','targeted_middle_R3_recheck','targeted_middlehigh_R3_recheck','robust_middlehigh_R3_recheck','robust_currentworst_R3_recheck','robust_broadcurrent_R3_recheck']:
  for p in (ROOT/'tonight/audit'/folder).glob('*.json'):copy(str(p.relative_to(ROOT)))
 for name in ['replay_R3_corners.py','replay_middlebank_R3_corners.py','replay_targeted_outer_boxes.py','replay_robust_outer_boxes.py']:
  copy('tonight/audit/global_b7171/'+name)
 for name in ['replay_targeted_highk.py','targeted_highk_whole_q_audit.json']:
  copy('tonight/theory/'+name)
 shutil.copy2(HERE/'verify_targeted_R3_bundle.py',DEST/'verify_targeted.py')
 readme=(DEST/'README.md').read_text()+'''\nThe additional targeted profiles are independently accepted whole-q gates. Run\n`python3 verify_targeted.py` for their fresh complete partitions, exact tails,\nand independently expanded restricted source-box prices. `--quick` checks tails\nand source prices using the frozen accepted body/exterior records. The full-E\nwhole-q auditor runs in a temporary copy so accepted input reports remain immutable.\nThe added source boxes also do not assert a global norm bound.\n'''
 (DEST/'README.md').write_text(readme)
 generated={'MANIFEST.json','COMPLETE_PORTABLE_REPLAY.json','VALIDATION.json','tonight/audit/global_b7171/targeted_32_outer_vertices.json','tonight/audit/global_b7171/robust_32_outer_vertices.json'}
 manifest={str(p.relative_to(DEST)):hashlib.sha256(p.read_bytes()).hexdigest()for p in DEST.rglob('*')if p.is_file()and'__pycache__'not in p.parts and str(p.relative_to(DEST))not in generated}
 (DEST/'MANIFEST.json').write_text(json.dumps(dict(kind='RELATIVE_SHA256_ACCEPTED_RESTRICTED_R3_INPUTS',sha256=manifest),indent=2))
 print(DEST,len(manifest))
if __name__=='__main__':main()
