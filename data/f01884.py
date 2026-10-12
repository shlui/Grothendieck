#!/usr/bin/env python3
"""Optimizer-free complete targeted R3 replay; all input paths are relative."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import argparse,hashlib,json,os,subprocess,sys,tempfile,shutil,time
ROOT=Path(__file__).resolve().parent;OUT=ROOT/'replay_outputs_targeted';OUT.mkdir(exist_ok=True)

def run(tag,args,root=ROOT):
 env=dict(os.environ);env['PYTHONPATH']=os.pathsep.join(str(root/p)for p in ['tonight/audit','tonight/theory','tonight/operator/projection/certificate']);env['PYTHONINTMAXSTRDIGITS']='0'
 log=OUT/(tag+'.log')
 with log.open('w')as f:r=subprocess.run([sys.executable,*map(str,args)],cwd=root,env=env,stdout=f,stderr=subprocess.STDOUT)
 if r.returncode:raise RuntimeError(f'{tag} failed; see {log}')
 return tag

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--quick',action='store_true');ap.add_argument('--workers',type=int,default=4);a=ap.parse_args();start=time.time()
 m=json.loads((ROOT/'MANIFEST.json').read_text())
 for p,h in m['sha256'].items():assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h,p
 specs=[('interval_lowr_u76_r16_k0_2','replay_R3_generous_q_trees.py','generous',64),('interval_middle_u7475_r267_k08_1','replay_targeted_rank_q_trees.py','middle',16),('interval_middlehigh_u7475_r275_k11_2','replay_targeted_rank_q_trees.py','middle',16),('interval_robust_middlehigh_box','replay_R3_robust_q_trees.py','robust',20),('interval_robust_currentworst_box','replay_R3_robust_current_q_trees.py','robustcurrent',20),('interval_robust_broadcurrent_box','replay_R3_robust_current_q_trees.py','robustcurrent',20)]
 jobs=[]
 for label,checker,branch,R in specs:
  for part in (['tail']if a.quick else ['body','exterior','tail']):
   args=[ROOT/'tonight/audit'/checker,'--profile',label,'--branch',branch,'--part',part,'--data-dir',ROOT/'tonight/operator/projection/targeted_certificate','--output-dir',OUT]
   if checker!='replay_R3_generous_q_trees.py':args+=['--tail-R',str(R)]
   jobs.append((label+'_'+part,args))
 with ThreadPoolExecutor(a.workers)as pool:
  for tag in pool.map(lambda x:run(*x),jobs):print('PASS',tag,flush=True)
 if not a.quick:
  with tempfile.TemporaryDirectory(prefix='r3-targeted-full-E-')as t:
   copied=Path(t)/'bundle';shutil.copytree(ROOT,copied,ignore=shutil.ignore_patterns('__pycache__','replay_outputs','replay_outputs_targeted'))
   run('targeted_full_E_whole_q',[copied/'tonight/theory/replay_targeted_highk.py','--workers',str(a.workers)],copied)
   shutil.copy2(copied/'tonight/theory/targeted_highk_whole_q_audit.json',OUT/'fresh_targeted_highk_whole_q.json')
 # Price checkers that have historical default output locations run in a
 # temporary copy, so optional fresh verification never alters immutable data.
 with tempfile.TemporaryDirectory(prefix='r3-targeted-prices-')as t:
  copied=Path(t)/'bundle';shutil.copytree(ROOT,copied,ignore=shutil.ignore_patterns('__pycache__','replay_outputs','replay_outputs_targeted'))
  for tag,checker,result_name in [('targeted_32_prices','replay_targeted_outer_boxes.py','targeted_32_outer_vertices.json'),('robust_32_prices','replay_robust_outer_boxes.py','robust_32_outer_vertices.json')]:
   run(tag,[copied/'tonight/audit/global_b7171'/checker],copied)
   shutil.copy2(copied/'tonight/audit/global_b7171'/result_name,OUT/('fresh_'+result_name))
   print('PASS',tag,flush=True)
  run('broadcurrent_16_prices',[copied/'tonight/audit/global_b7171/replay_robust_outer_boxes.py','--source',copied/'tonight/operator/projection/targeted_certificate/broadcurrent_source_box.json','--output',OUT/'broadcurrent_16_prices.json','--threshold','451/800'],copied)
  print('PASS broadcurrent_16_prices',flush=True)
 for p,h in m['sha256'].items():assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h,p
 d=dict(status='PASS_PORTABLE_TARGETED_RESTRICTED_R3_REPLAY',mode='quick_tails_prices_accepted_q_premises'if a.quick else'complete_partitions_tails_prices',seconds=time.time()-start,scope='Only the listed source boxes; no global operator norm or new Grothendieck lower endpoint follows without complete source-domain coverage.')
 (OUT/'result.json').write_text(json.dumps(d,indent=2));print(json.dumps(d,indent=2))
if __name__=='__main__':main()
