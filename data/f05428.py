from pathlib import Path
import hashlib,json,shutil
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
OLD=HERE/'replay_accepted_R3_expanded';DEST=HERE/'replay_accepted_R3_final'
NEW=[('interval_finalface_u7475_r269_k4_10','replay_R3_q50_q_trees.py','q50','finalface_u7475_R3_recheck','finalface'),('interval_lastface_u7575_r238_k35_12','replay_R3_lowerface_q_trees.py','lowerface','lastface_u7575_R3_recheck','lastface'),('interval_refinedface_u75125_r262_k10_12','replay_R3_lowerface_q_trees.py','lowerface','refinedface_u75125_R3_recheck','refinedface')]
if DEST.exists():shutil.rmtree(DEST)
shutil.copytree(OLD,DEST,ignore=shutil.ignore_patterns('__pycache__','replay_outputs*'))
def copy(rel):
 f=ROOT/rel;g=DEST/rel;g.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(f,g)
specs=json.loads((DEST/'EXPANDED_SPECS.json').read_text())
for lab,checker,branch,audit,ledger in NEW:
 for stem,suffix in [('profile_','.json'),('body_','.json.gz'),('exterior_','.json.gz'),('tail_','.json')]:copy('tonight/operator/projection/targeted_certificate/'+stem+lab+suffix)
 for f in (ROOT/'tonight/audit'/audit).glob('*.json'):copy(str(f.relative_to(ROOT)))
 copy('tonight/audit/'+checker)
 for f in [ledger+'_source_box.json','price_'+ledger+'_source_box.py']:copy('tonight/operator/projection/targeted_certificate/'+f)
 specs.append(dict(label=lab,checker=checker,branch=branch,audit=audit,ledger=ledger,R=32))
(DEST/'EXPANDED_SPECS.json').write_text(json.dumps(specs,indent=2))
p=DEST/'verify_expanded.py';s=p.read_text().replace('fresh_partitions=5,exact_tails=5,fresh_source_corners=80','fresh_partitions=len(specs),exact_tails=len(specs),fresh_source_corners=16*len(specs)');p.write_text(s)
with (DEST/'README.md').open('a')as f:f.write('\nThe final package also includes the three accepted late operator gates finalface7475, lastface7575, and refinedface75125. The same command freshly replays all eight new complete partitions and128 source corners; `--all` includes the older targeted gates too. These are restricted-source results; a separate complete global envelope is needed for any new constant bound.\n')
excluded={'MANIFEST.json','COMPLETE_PORTABLE_REPLAY.json','VALIDATION.json','tonight/audit/global_b7171/targeted_32_outer_vertices.json','tonight/audit/global_b7171/robust_32_outer_vertices.json'}
sha={str(f.relative_to(DEST)):hashlib.sha256(f.read_bytes()).hexdigest()for f in DEST.rglob('*')if f.is_file()and'__pycache__'not in f.parts and str(f.relative_to(DEST))not in excluded}
(DEST/'MANIFEST.json').write_text(json.dumps(dict(kind='RELATIVE_SHA256_FINAL_ACCEPTED_RESTRICTED_R3_INPUTS',sha256=sha),indent=2));print(DEST,len(sha))
