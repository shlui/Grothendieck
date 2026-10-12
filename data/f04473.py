#!/usr/bin/env python3
"""Complete directed price cover using eight accepted plus two new rows."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,sys
from flint import ctx
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
from certify_kappa_shape_caps import load_rows,corner_cap


def run():
    ctx.prec=512
    rows=load_rows()
    for name,file in [('68','row_68.json'),('96','row_96.json')]:
        path=HERE/file;r=json.loads(path.read_text())
        assert r['status']=='PASS_DIRECTED_SHARP_GENERIC_ODD_CHEBYSHEV_RATIONAL_HALF_LINE'
        q=list(map(F,r['majorant']));eta=F(r['eta']);assert q[3]>0 and eta>0
        rows[name]=dict(q=q,eta=eta,sha256=hashlib.sha256(path.read_bytes()).hexdigest())
    prior=json.loads((HERE.parent/'kappa_shape_caps.json').read_text())
    leaves=[]
    for source in prior['leaves']:
        lo,hi=map(F,source['kappa_interval'])
        cap,name=min((corner_cap(row,lo,hi)[0],name)for name,row in rows.items())
        leaves.append(dict(kappa_interval=[str(lo),str(hi)],row=name,cap=str(cap),
                           corners=corner_cap(rows[name],lo,hi)[1]))
    maximum=max(F(r['cap'])for r in leaves)
    assert maximum<F(847,1000)
    table=[]
    for source in prior['table']:
        lo,hi=map(F,source['kappa_interval']);pieces=[r for r in leaves if F(r['kappa_interval'][0])>=lo and F(r['kappa_interval'][1])<=hi]
        cap=max(F(r['cap'])for r in pieces)
        table.append(dict(kappa_interval=[str(lo),str(hi)],cap=str(cap),display=float(cap),pieces=len(pieces)))
    out=dict(status='PASS_COMPLETE_CUBIC_L1_REFIT_PRICE_COVER',target='847/1000',
        maximum_price=str(maximum),maximum_display=float(maximum),precision_bits=ctx.prec,
        rows={n:dict(majorant=list(map(str,r['q'])),eta=str(r['eta']),sha256=r['sha256'],
            path=(f'row_{n}.json' if n in ['68','96'] else f'../../../../reduced_search/energy/stein_L1_row_{n}.json'))for n,r in rows.items()},
        leaves=leaves,table=table,
        proof_scope='Complete kappa0..90/shape0..1 corner cover. Two NEW complete-halfline scalar rows and eight already accepted complete-halfline rows; same fresh633/Stein analytic inputs. New independent halfline replay tracked separately.',
        source_cover_sha256=hashlib.sha256((HERE.parent/'kappa_shape_caps.json').read_bytes()).hexdigest())
    (HERE/'complete_price_cover.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:out[k]for k in ['status','target','maximum_price','maximum_display','table']},indent=2))


if __name__=='__main__':run()
