#!/usr/bin/env python3
import os,json,subprocess,hashlib
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
r=Path('/workspace/scratch/8242d848994f');out=Path(__file__).resolve().parent;num=r/'new_round_2026-10-02/audit/numeric/strict_all_price_replay';num.mkdir(exist_ok=True)
config=r/'new_round_2026-10-02/audit/portable_replay/FINAL_GOAL_REPLAY_CONFIG.json';c=json.loads(config.read_text());env=dict(os.environ,GR_FROZEN_BASELINE=str(r/'verified_2026-10-01/tonight/operator/preferred_global'),GR_MPFR_BITS='768',PYTHONDONTWRITEBYTECODE='1')
jobs=[j for j in c['jobs']if '/audit/numeric/'in j['script']]
def run(j):
 def arg(x):return x.replace('{numeric}',str(num)).replace('{output}',str(out))if'{'in x else str(r/x)if x.startswith(('new_round','verified_'))else x
 dst=Path(arg(j['result']));log=out/(j['name']+'.log')
 with log.open('w')as f:p=subprocess.run(['python',str(r/j['script']),*map(arg,j['args'])],env=env,cwd=r,stdout=f,stderr=subprocess.STDOUT)
 assert p.returncode==0,(j['name'],p.returncode,str(log));d=json.loads(dst.read_text());assert d['status'].startswith('PASS');copy=out/dst.name;copy.write_bytes(dst.read_bytes());result={'name':j['name'],'status':d['status'],'fresh_report':str(copy),'fresh_sha256':hashlib.sha256(copy.read_bytes()).hexdigest(),'checker_sha256':hashlib.sha256((r/j['script']).read_bytes()).hexdigest()};print(json.dumps(result),flush=True);return result
with ThreadPoolExecutor(4)as pool:results=list(pool.map(run,jobs))
(out/'ALL_DEFAULT_NUMERIC_REPLAYS.json').write_text(json.dumps({'status':'PASS_FRESH_ALL_DEFAULT_NUMERIC_PRICE_REPLAYS','jobs':len(results),'bits':768,'config_sha256':hashlib.sha256(config.read_bytes()).hexdigest(),'results':results,'local_whole_line_trees_rerun':False,'scope':'All current directed price/plane/native jobs repriced; inherited accepted whole-line criteria remain analytic inputs.'},indent=2)+'\n')
