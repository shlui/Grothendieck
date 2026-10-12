#!/usr/bin/env python3
"""Path-independent, optimizer-free replay of restricted R3 certificates.

Default runs complete frozen-q trees, exact tails and source-price proofs.
--quick runs hashes, all tails, and source prices; it uses the included
independent all-q acceptance records as premises instead of rerunning trees.
"""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import argparse,hashlib,json,os,subprocess,sys,time
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'replay_outputs';OUT.mkdir(exist_ok=True)

def run(tag,args):
 env=dict(os.environ);env['PYTHONPATH']=os.pathsep.join(str(ROOT/x)for x in ['tonight/audit','tonight/theory','tonight/operator/projection/certificate'])
 command=[sys.executable,*[str(x)for x in args]]
 path=OUT/(tag+'.log')
 with path.open('w')as f:r=subprocess.run(command,cwd=ROOT,env=env,stdout=f,stderr=subprocess.STDOUT)
 if r.returncode:raise RuntimeError(f'{tag} failed; see {path}')
 return tag

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--quick',action='store_true');ap.add_argument('--workers',type=int,default=4);args=ap.parse_args();start=time.time()
 manifest=json.loads((ROOT/'MANIFEST.json').read_text())
 for rel,digest in manifest['sha256'].items():assert hashlib.sha256((ROOT/rel).read_bytes()).hexdigest()==digest,rel
 print(json.dumps({'stage':'RELATIVE_MANIFEST_PASS','files':len(manifest['sha256'])}),flush=True)
 qjobs=[]
 for label,folder in [('68','lower_amp_certificate'),('70','lower_amp_certificate'),('interval_u7478_k0_01','diffuse_certificate'),('interval_b7171_newworst','diffuse_certificate')]:
  for part in (['tail']if args.quick else ['body','exterior','tail']):
   qjobs.append((label+'_'+part,[ROOT/'tonight/audit/replay_R3_physical_q_trees.py','--profile',label,'--branch','finer','--part',part,'--data-dir',ROOT/'tonight/operator/projection'/folder,'--output-dir',OUT]))
 if not args.quick:
  for label in ['68','70']:
   for part in ['body','exterior','tail']:
    qjobs.append((label+'_tight_'+part,[ROOT/'tonight/audit/replay_R3_physical_q_trees.py','--profile',label,'--branch','tight','--part',part,'--output-dir',OUT]))
 with ThreadPoolExecutor(args.workers)as pool:
  for name in pool.map(lambda x:run(*x),qjobs):print(json.dumps({'stage':'Q_REPLAY_PASS','job':name}),flush=True)
 for tag,command in [('kappa_L1_prices',[ROOT/'tonight/audit/refined_l1/replay_refined_prices.py']),('lower_source_prices',[ROOT/'tonight/audit/lower_amp_R3/independent_source_price.py']),('diffuse_old_source_prices',[ROOT/'tonight/audit/diffuse_R3/independent_diffuse_source_price.py']),('diffuse_new_b7171_source_prices',[ROOT/'tonight/audit/diffuse_b7171_R3/independent_source_pricing.py']),('free_b_native_remainder',[ROOT/'verify_b_transport_box.py']),('general_exact_gates_tails',[ROOT/'tonight/theory/audit_general_r3_parameters_tail.py'])]:
  run(tag,command);print(json.dumps({'stage':'SOURCE_REPLAY_PASS','job':tag}),flush=True)
 if not args.quick:
  for label in ['70','68']:
   run('general_q_'+label,[ROOT/'tonight/theory/replay_general_r3_cover.py','--label',label,'--workers',str(args.workers)]);print(json.dumps({'stage':'GENERAL_Q_PASS','label':label}),flush=True)
 run('general_source_bands',[ROOT/'tonight/theory/audit_general_r3_source_bands.py'])
 result={'status':'PASS_PORTABLE_ACCEPTED_RESTRICTED_R3_BRANCHES','mode':'quick_tails_prices_and_accepted_q_premises'if args.quick else 'complete_q_trees_tails_and_source_prices','seconds':time.time()-start,'scope':'Only listed physical source domains. This package does not establish an all-source operator norm or any new global Grothendieck lower bound. Global 633/moment and accepted global cubic-L1 scalar inequalities are external analytic premises documented in README.'}
 (OUT/'result.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
