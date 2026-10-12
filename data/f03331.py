"""Floating cross-grid check of every coarse local optimizer, not only top six."""
import json
from pathlib import Path
import numpy as np
from search import Engine, optimize

HERE = Path(__file__).resolve().parent
sources=['results.json','diverse_results.json']
data={'status':'DISCOVERY_ONLY','all_coarse':[],'refined':[]}
e=Engine(96)
for source in sources:
    D=json.loads((HERE/source).read_text())
    for row in D['coarse']:
        checked=e.evaluate(np.array(row['parameters']),row['b'],row['symmetric'],False)
        data['all_coarse'].append({'source':source, 'start':row['start'],
            'b':row['b'],'symmetric':row['symmetric'],'parameters':row['parameters'],
            'value_96':checked['value'],'coarse_value':row['discovery_value']})
    print(json.dumps({'source':source,'checked':len(D['coarse']),
        'best96':max(r['value_96'] for r in data['all_coarse'] if r['source']==source)}),flush=True)
    (HERE/'all_coarse_crosscheck.json').write_text(json.dumps(data,indent=2)+'\n')
for b in [.72,.7171]:
    rows=sorted((r for r in data['all_coarse'] if r['b']==b),
                key=lambda r:r['value_96'],reverse=True)[:3]
    for r in rows:
        row=optimize(e,np.array(r['parameters']),b,r['symmetric'],400)
        row['source']=r['source'];row['coarse_start']=r['start'];row['checks']=[]
        data['refined'].append(row)
del e
for n in [128,160]:
    e=Engine(n)
    for row in data['refined']:
        q=e.evaluate(np.array(row['parameters']),row['b'],row['symmetric'],False)
        q['nodes']=n;row['checks'].append(q)
    del e
    print(json.dumps({'nodes':n,'best_by_b':{str(b):max(r['checks'][-1]['value'] for r in data['refined'] if r['b']==b) for b in [.72,.7171]}}),flush=True)
    (HERE/'all_coarse_crosscheck.json').write_text(json.dumps(data,indent=2)+'\n')
data['finished']=True
(HERE/'all_coarse_crosscheck.json').write_text(json.dumps(data,indent=2)+'\n')
