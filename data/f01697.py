#!/usr/bin/env python3
"""Exact independent domain check for saved all-active scalar tails."""
from pathlib import Path
from fractions import Fraction as F
import json
HERE=Path(__file__).resolve().parent
checks=[]
for summary in sorted(HERE.rglob('continuous_scalar_summary.json')):
 d=json.loads(summary.read_text());src=Path(d['source']);raw=json.loads(src.read_text())
 for item in d['rows']:
  row=item['row'];cert=json.loads((summary.parent/item['file']).read_text())
  r=F(str(raw['rows'][row]['r']));last=max(max(F(str(x))for x in raw['t'])/r,F(12));cutoff=F(cert['tail']['domain'][0])
  assert cutoff>=last,(summary,row,last,cutoff)
  assert cutoff==F(cert['compact']['domain'][1])
  checks.append({'directory':str(summary.parent.relative_to(HERE)),'row':row,'last_hinge':str(last),'tail_cutoff':str(cutoff)})
out={'status':'PASS_EXACT_ALL_ACTIVE_TAIL_DOMAIN_AUDIT','rows_checked':len(checks),'checks':checks}
(HERE/'tail_domain_audit.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'status':out['status'],'rows_checked':len(checks)}))
