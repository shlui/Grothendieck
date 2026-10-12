"""Directed reflected-sign tests; ceilings for the cubic construction only."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,runpy,sys
ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'reduced_search/operator/cubic720/certify_cubic720_norm_witness.py'
sys.path.insert(0,str(BASE.parent))
module=runpy.run_path(str(BASE));compute=module['compute'];g=compute.__globals__
rows=json.loads((ROOT/'tonight/optimize_cubic_witness_coefficients.json').read_text())
out=[]
for row in rows:
 if row['b'] not in ('141/200','713/1000','7171/10000'):continue
 g.update(A=F(row['a']),B=F(row['c']),EPS=F(row['b']))
 r=compute(320,10,5,check_fixed_witness=False)
 r['scope']='Directed norm LOWER witness and pure-cubic reciprocal-norm construction CEILING at the recorded exact coefficient. No Grothendieck lower bound is inferred.'
 lo=F(r['scalar_norm_lower_witness']['lower'])
 hi=F(r['scalar_norm_lower_witness']['upper'])
 r['excludes_1775_for_this_coefficient']=lo>F(40,71)
 r['formula_source']={'path':str(BASE.relative_to(ROOT)),'sha256':hashlib.sha256(BASE.read_bytes()).hexdigest()}
 out.append(r)
 print(json.dumps({'b':row['b'],'norm_lower':float(lo),'norm_upper':float(hi),'method_ceiling_upper':float(1/lo),'excludes_1775':r['excludes_1775_for_this_coefficient']}),flush=True)
p=ROOT/'tonight/selected_witness_coefficient_certificates.json';p.write_text(json.dumps(out,indent=2)+'\n')
