#!/usr/bin/env python3
"""Fresh directed price comparison for the improved global cubic L1 cap.

This reads only exact row data and imports the independently derived MPFR
price/tail formulas. Accepted compact scalar inequalities and the fresh
633 moment theorem are separately identified premises.
"""
from fractions import Fraction as F
from pathlib import Path
import hashlib,json
from mpfr_interval import I,UP,selftest
from replay_stein_prices_and_tails import prices,tail
HERE=Path(__file__).resolve().parent
SOURCE=HERE.parents[1]/'reduced_search/energy'
cap=F(848521,10**6)
certificate=json.loads((SOURCE/'cubic_L1_certificate.json').read_text())
intervals=[(item['row_path'],*map(F,item['kappa_interval'])) for item in certificate['rows']]
assert len(intervals)==8 and intervals[0][1]==0 and intervals[-1][2]==90
assert all(lo<hi for _,lo,hi in intervals)
assert all(a[2]==b[1] for a,b in zip(intervals[:-1],intervals[1:]))
assert (I(194)/I(10).sqrt()-F(633,15)).above(0)
rows=[];maximum=F(-1000);active=None
for filename,lo,hi in intervals:
 path=SOURCE/filename;row=json.loads(path.read_text());q0,q2,q4,q6=map(F,row['majorant']);eta=F(row['eta'])
 assert q6>0 and eta>0
 values=prices(row,lo,hi,cap)
 for v in values:
  upper=F(v['price_interval']['upper']);assert upper<cap
  if upper>maximum:maximum=upper;active={'row_path':filename,'kappa':v['kappa'],'sixth_coefficient_offset':v['sixth_coefficient_offset']}
 rows.append({'row_path':filename,'row_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'kappa_interval':[str(lo),str(hi)],'majorant':list(map(str,(q0,q2,q4,q6))),'eta':str(eta),'independent_directed_corner_prices':values,'exact_tail_checks':tail(row)})
result={'status':'PASS_IMPROVED_GLOBAL_CUBIC_L1_PRICE_CAP','engine':selftest(),'rho_star':str(cap),'rho_star_squared':str(cap*cap),'corner_count':32,'maximum_price_upper':str(maximum),'maximum_price_display':float(maximum),'active_corner':active,'strict_margin':str(cap-maximum),'rows':rows,'convexity':'For both shape endpoints the kappa^(3/2) coefficient is q6*(194/sqrt(10)-offset)>0. The price is convex in kappa; the shared variance price is affine in shape. Gap-free8-intervalcover therefore reduces to these32 corners.','dependencies':['Eight accepted global whole-line scalar comparisons, with compact intervals independently replayed earlier.','Fresh sharp633 cubic moment theorem and0<=kappa<=90,7kappa/45<=V<=2kappa/9.'],'scope':'Universal E|Q|<848521/1e6 for normalized pure thirdGaussianchaos. No globalGrothendieckendpoint asserted.'}
(HERE/'cubic_L1_improved_rho_audit.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:result[k] for k in ('status','rho_star','corner_count','maximum_price_display','active_corner')},indent=2))
