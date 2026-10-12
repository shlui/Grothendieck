#!/usr/bin/env python3
"""Independent exact/MPFR source-Lipschitz extension of the high-u profile."""
from pathlib import Path
from fractions import Fraction as F
import json
from mpfr_interval import I,UP
HERE=Path(__file__).resolve().parent
base=HERE.parents[1]
data=json.loads((base/'tonight/operator/high_u_literature/source_profile.json').read_text())
local=json.loads((base/'tonight/recovery/profile_cover/high_u_literature_summary.json').read_text())
accepted=json.loads((HERE/'high_u_independent_endpoint_ledger.json').read_text())
assert accepted['status']=='PASS_INDEPENDENT_HIGH_U_ENDPOINT_ASSEMBLY'
p=local['profile'];rho=F(849,1000);lo=F(3941,5000);hi=F(159577,200000);ell=F(p['ell']);tau=F(p['tail_slope']);g=F(accepted['exact_gram_cost'])
assert ell>=hi and g<0 and tau+2*g>0
rm=F(str(data['r_grid'][-1]));reserve=F(accepted['source_chord_reserve'])
old=max([F(p['nodal_majorant'][0])]+[F(x['independent_price_upper']) for x in accepted['endpoint_ledger']])+reserve
terminal=max(F(x['independent_price_upper']) for x in accepted['endpoint_ledger'] if F(x['r'])==rm)
c=(I(2)/I.pi()).sqrt();v=(I(lo)/c-1+rho*rho)/(rho*rho)
assert v.above(0) and v.below(1)
R=I(rho)*(-v.log()).sqrt().erf();assert R.above(rm) and R.below(rho)
beyond=I(terminal)+I(tau*rho)*(R-rm)+I(g)*(R*R-rm*rm)
beyond_upper=F(beyond.hi.text(UP,180));source=max(old,beyond_upper)
total=source+lo*lo-2*ell*lo;assert total<F(11267,20000)<F(40,71)
out={'status':'PASS_INDEPENDENT_HIGH_U_LIPSCHITZ_EXTENSION','amplitude_interval':[str(lo),str(hi)],'total_upper':str(total),'total_display':float(total),'rational_cap':'11267/20000','old_source_bound':str(old),'terminal_source_bound':str(terminal),'extended_source_bound':str(beyond_upper),'source_cap_interval':R.record(),'exact_derivative_margin':str(tau+2*g),'basis':'Independent MPFR all48-corner endpoint ledger; accepted repaired global Lipschitz target; fresh rho and source bridge. Regional theorem only.'}
(HERE/'high_u_lipschitz_independent_ledger.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:out[k]for k in ('status','total_display','rational_cap','amplitude_interval')},indent=2))
