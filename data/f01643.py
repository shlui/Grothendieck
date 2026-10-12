#!/usr/bin/env python3
"""Independent rational extrapolation of frozen whole-frame transfer costs.
This is a failure of one sufficient proof route, not a scalar norm lower bound.
"""
from fractions import Fraction as F
from pathlib import Path
import hashlib,json
HERE=Path(__file__).resolve().parent
OLD=HERE/'frozen_easy_tail_audit.json'
x=json.loads(OLD.read_text())
lo=next(t for t in x['old_five_price_transfers'] if F(t['epsilon'])==F(71,100))
hi=next(t for t in x['old_five_price_transfers'] if F(t['epsilon'])==F(717,1000))
e=F(18,25); target=F(40,71)
rows=[]
for r0,r1 in zip(lo['rows'],hi['rows']):
 assert r0['interval']==r1['interval'] and r0['index']==r1['index']
 q0=list(map(F,r0['new_quadratic'])); q1=list(map(F,r1['new_quadratic']))
 assert q0[:2]==q1[:2]
 c0,c1=F(r0['new_cap']),F(r1['new_cap'])
 # fixed profile cost: mass*max(epsilon-epsilon_profile,0)+beta^2/epsilon.
 # It is convex on (0,infinity), so a preceding secant bounds it below.
 lower=c1+(e-F(717,1000))*(c1-c0)/(F(717,1000)-F(71,100))
 assert lower>target
 rows.append({'interval':r1['interval'],'profile_index':r1['index'],
  'epsilon720_transfer_cap_lower':str(lower),'lower_decimal':float(lower),
  'excess_above_40_71':str(lower-target)})
A=F(2403891349687,2000000000000); B=F(12706580435619,10000000000000)
C=F(797884561,10**9); l=F(19947,25000); erf=F(56419,50000); sq=F(2391,10**6)
high=C*C+e*erf*sq-e*(3*l-2*C)**2/6
assert high<target
result={'status':'PASS_EXACT_FROZEN_TRANSFER_ROUTE_AUDIT',
 'epsilon':str(e),'target':str(target),
 'inherited_premises':'Only arithmetic is regenerated: frozen transfer records and the proved mixed-energy/high-amplitude estimates are inputs.',
 'low_u_energy_cap_at_zero':str(e/B),'low_u_energy_cap_decimal':float(e/B),
 'low_u_energy_excess':str(e/B-target),'needed_energy_B_at_zero':str(e/target),
 'needed_squared_cubic_L1_cap':str(target/e),
 'high_u_interval':[str(l),'sqrt(2/pi)'],'high_u_cap':str(high),'high_u_cap_decimal':float(high),
 'whole_frame_transfers_obstructed':len(rows),'minimum_transfer_cap_lower':min(r['lower_decimal'] for r in rows),
 'rows':rows,'frozen_input_sha256':hashlib.sha256(OLD.read_bytes()).hexdigest(),
 'scope':'32 sufficient whole-frame costs exceed target at epsilon .72. This does not prove the true operator norm exceeds target and does not rule out profile refits or clipped source cells.'}
(HERE/'epsilon720_transfer_audit.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:result[k] for k in ['status','whole_frame_transfers_obstructed','minimum_transfer_cap_lower','low_u_energy_cap_decimal','high_u_cap_decimal']},indent=2))
