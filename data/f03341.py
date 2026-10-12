#!/usr/bin/env python3
"""Directed kappa/shape repricing of eight already accepted Stein rows.

No new scalar halfline is asserted. Every row's accepted scalar inequality
holds for all real x, so its analytic moment price is valid on all kappa.
Rational subdivisions and convex corner prices prove each interval cap.
"""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json
from flint import arb,ctx

HERE=Path(__file__).resolve().parent
ENERGY=HERE.parents[2]/'reduced_search/energy'
NAMES=['0','2','3','4','k8','5','6','7']


def aa(x):
    x=F(x)
    return arb(x.numerator)/x.denominator


def upper(x,digits=12):
    v=F(str(x.upper().fmpq()));scale=10**digits
    return F(-((-v.numerator*scale)//v.denominator),scale)


def load_rows():
    accepted=json.loads((ENERGY/'cubic_L1_certificate.json').read_text())
    assert accepted['status']=='PASS_DIRECTED_COMPLETE_CUBIC_L1_BOUND'
    rows={}
    for rec in accepted['rows']:
        path=ENERGY/rec['row_path']
        assert hashlib.sha256(path.read_bytes()).hexdigest()==rec['row_sha256']
        row=json.loads(path.read_text())
        assert row['status']=='PASS_DIRECTED_SHARP_GENERIC_ODD_CHEBYSHEV_RATIONAL_HALF_LINE'
        q=list(map(F,row['majorant']));eta=F(row['eta'])
        assert q[3]>0 and eta>0
        name=path.name[len('stein_L1_row_'):-len('.json')]
        rows[name]=dict(q=q,eta=eta,sha256=rec['row_sha256'])
    assert set(rows)==set(NAMES)
    return rows


def price(row,k,shape):
    q,eta=row['q'],row['eta'];k,shape=F(k),F(shape)
    assert 0<=k<=90 and 0<=shape<=1
    c=(2/arb.pi()).sqrt();K=194/arb(10).sqrt()
    assert K-aa(F(633,15))>0
    return c+aa(q[0]+q[1]+(3+k)*q[2])+aa(q[3])*(
        15+15*aa(k)+(K-aa(F(633,15)*shape))*aa(k)*aa(k).sqrt()
        )+aa(eta*k*(F(7,90)+shape/F(30)))


def corner_cap(row,lo,hi,slo=F(0),shi=F(1)):
    corners=[upper(price(row,k,s)) for k in (lo,hi) for s in (slo,shi)]
    return max(corners),list(map(str,corners))


def run():
    ctx.prec=512
    rows=load_rows()
    bins=[('coarse0_1',F(0),F(1)),('coarse1_10',F(1),F(10)),
          ('coarse10_90',F(10),F(90)),('fine0_03',F(0),F(3,10)),
          ('fine03_11',F(3,10),F(11,10)),('fine11_4',F(11,10),F(4)),
          ('fine4_10',F(4),F(10)),('high10_20',F(10),F(20)),
          ('high20_40',F(20),F(40)),('high40_60',F(40),F(60)),
          ('high60_90',F(60),F(90))]
    grid={F(i,100) for i in range(1001)}|{F(i,20) for i in range(200,1801)}
    grid|={F(i,1000) for i in range(9500,9701)}
    grid|={a for _,a,_ in bins}|{b for _,_,b in bins}
    grid=sorted(grid)
    leaves=[]
    for lo,hi in zip(grid[:-1],grid[1:]):
        caps=[(corner_cap(rows[n],lo,hi)[0],n) for n in NAMES]
        cap,name=min(caps)
        leaves.append(dict(kappa_interval=[str(lo),str(hi)],row=name,
                           cap=str(cap),corners=corner_cap(rows[name],lo,hi)[1]))
    tables=[]
    for name,lo,hi in bins:
        pieces=[(i,p)for i,p in enumerate(leaves)
                if F(p['kappa_interval'][0])>=lo and F(p['kappa_interval'][1])<=hi]
        assert pieces and F(pieces[0][1]['kappa_interval'][0])==lo
        assert F(pieces[-1][1]['kappa_interval'][1])==hi
        maximum,active=max((F(p['cap']),i) for i,p in pieces)
        tables.append(dict(name=name,kappa_interval=[str(lo),str(hi)],rho=str(maximum),
                           rho_display=float(maximum),pieces=len(pieces),active_leaf=active))
    points=[]
    for k in map(F,['0','.3','.53','1','1.1','4','6','8','10','12','20','28','40','60','90']):
        cap,name=min((corner_cap(rows[n],k,k)[0],n)for n in NAMES)
        points.append(dict(kappa=str(k),rho=str(cap),rho_display=float(cap),row=name,
                           shape_corners=corner_cap(rows[name],k,k)[1]))
    maximum=max(F(p['cap'])for p in leaves)
    assert maximum<F(849,1000)
    out=dict(status='PASS_CERTIFIED_KAPPA_SHAPE_CUBIC_L1_CAPS',precision_bits=ctx.prec,
        shape_parameter='E U^2 = kappa*(7/45+shape/15), 0<=shape<=1',
        formula='L_i(k,s)=c+q0_i+q2_i+(3+k)q4_i+[15+15k+(K-633*s/15)*k^(3/2)]q6_i+eta_i*k*(7/90+s/30)',
        constants=dict(c='sqrt(2/pi)',K='194/sqrt(10)'),
        price_validity='Each accepted whole-line comparison is valid for every normalized cubic; discovery/old pricing intervals do not restrict its new moment price. For fixed row, price is affine in shape and convex in kappa.',
        rows={n:dict(majorant=list(map(str,r['q'])),eta=str(r['eta']),sha256=r['sha256'])for n,r in rows.items()},
        table=tables,points=points,leaves=leaves,
        global_rho=str(maximum),global_rho_display=float(maximum),
        accepted_energy_certificate_sha256=hashlib.sha256((ENERGY/'cubic_L1_certificate.json').read_bytes()).hexdigest(),
        scope='Fresh directed repricing and exact rational interval cover using previously accepted complete scalar rows and fresh633 theorem. No new scalar potential and no global Grothendieck operator cap.')
    (HERE/'kappa_shape_caps.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:out[k]for k in ['status','global_rho','global_rho_display','table','points']},indent=2))


if __name__=='__main__':
    run()
