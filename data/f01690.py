#!/usr/bin/env python3
"""Collect a complete accepted row set, including mixed exact allowances."""
from pathlib import Path
from fractions import Fraction as F
import argparse,hashlib,json
p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--output-dir',type=Path,required=True);a=p.parse_args()
raw=json.loads(a.source.read_text());digest=hashlib.sha256(a.source.read_bytes()).hexdigest();rows=[];allow=[];branches={}
for i,rr in enumerate(raw['rows']):
 path=a.output_dir/f'continuous_scalar_row_{i:02d}.json';r=json.loads(path.read_text());fr=r['frozen_row']
 assert r['row_index']==i and r['status']=='PASS_FRESH_COMPLETE_SCALAR_HALF_LINE' and r['source_sha256']==digest
 assert F(fr['r'])==F(str(rr['r'])) and F(fr['lambda_L1'])==F(str(rr['lambda_L1'])) and F(fr['eta'])==F(str(rr['eta']))
 q0,q6=F(fr['added_q0']),F(fr['added_q6']);assert q0>=0 and q6>=0
 assert list(map(F,fr['majorant']))==[F(str(x))+(q0 if j==0 else q6 if j==3 else 0)for j,x in enumerate(rr['majorant'])]
 assert tuple(map(F,fr['kappa_interval']))==tuple(F(str(x))for x in rr['kappa_interval'])
 last=max(max(F(str(t))for t in raw['t'])/F(fr['r']),max(map(F,fr['regression_hinges'])))
 assert F(r['tail']['domain'][0])>=last and F(r['tail']['domain'][0])==F(r['compact']['domain'][1])
 b=r['tail']['branch'];branches[b]=branches.get(b,0)+1
 rows.append({'row':i,'status':r['status'],'leaves':r['compact']['leaf_count'],'tail_cutoff':r['tail']['domain'][0],'largest_corner':max(x['price_display']for x in r['prices']),'seconds':r['seconds'],'file':path.name})
 allow.append({'row':i,'q0':str(q0),'q6':str(q6)})
out={'status':'PASS','source':str(a.source.resolve()),'rows':rows,'aggregate_note':'All source rows have complete accepted half-lines, exact source binding, explicit row-specific allowances and valid all-active tail domains.','allowances':allow,'tail_branches':branches,'total_leaves':sum(x['leaves']for x in rows)}
(a.output_dir/'continuous_scalar_summary.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'rows':len(rows),'leaves':out['total_leaves'],'tail_branches':branches,'nondefault_q6':[(x['row'],x['q6'])for x in allow if x['q6']!='1/1000000000'],'large_cutoffs':[(x['row'],x['tail_cutoff'])for x in rows if x['tail_cutoff']!='64']},indent=2))
