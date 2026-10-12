"""Native and frozen-contact prices on exact physical source polygons.

Source planes are globally proved on each complete cumulant/shape band.
Native source chords are split at every saved radial knot before applying
the exact boundary maximization. Floating routines only nominate a frame.
"""
import bisect
from fractions import Fraction as F
import math
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'theory'))
from source_polygon_prices import clip_source_polygon,max_source_quadratic
from price_native_complement import ceil_decimal
from contact_credit import contact_band

def native_polygon_bound(d,ulo,uhi,klo,khi,rlo,rhi,planes,lam=F(0),constant=F(0),a=F(0)):
    rg=d['_r']
    if rhi>rg[-1]:return None
    endpoints=d['_endpoints'][(klo,khi)]
    cuts=sorted(set([rlo,rhi]+[r for r in rg if rlo<r<rhi]))
    intervals=list(zip(cuts[:-1],cuts[1:])) if rlo<rhi else [(rlo,rhi)]
    values=[]
    for lo,hi in intervals:
        polygon=clip_source_polygon(ulo,uhi,lo,hi,planes)
        if not polygon:continue
        j=min(bisect.bisect_right(rg,lo)-1,len(rg)-2)
        assert rg[j]<=lo<=hi<=rg[j+1]
        for c in range(4):
            slope=(endpoints[j+1][c]-endpoints[j][c])/(rg[j+1]-rg[j])
            intercept=endpoints[j][c]-slope*rg[j]+constant
            value,_=max_source_quadratic(polygon,-2*d['_ell']-lam*a,d['_g'],slope-lam,intercept)
            values.append(value)
    return ceil_decimal(max(values)) if values else None

def source_polygon_bound(t,corner,polygon,lam=F(0),constant=F(0),a=F(0)):
    result=max_source_quadratic(polygon,-2*t['p']['ell']-lam*a,t['_g'],t['_lin']-lam,corner+constant)
    return None if result is None else ceil_decimal(result[0])

def polygon_contact(t_or_d,corner,ulo,uhi,klo,khi,rlo,rhi,planes,polygon,contact,b,kind):
    Aq,index,w=contact_band(contact,klo,khi)
    if kind=='native':beta=t_or_d['_contact_beta_term']
    else:beta=6*t_or_d['p']['beta']*contact['B']/b
    S=Aq-contact['a']*ulo-rlo+beta
    if S>=0:return None
    ideal=-2*b*S/contact['d2'];scale=10**10
    lam=F((ideal.numerator*scale)//ideal.denominator,scale)
    assert lam>=0
    constant=lam*(Aq+beta)+lam*lam*contact['d2']/(4*b)
    if kind=='native':cap=native_polygon_bound(t_or_d,ulo,uhi,klo,khi,rlo,rhi,planes,lam,constant,contact['a'])
    else:cap=source_polygon_bound(t_or_d,corner,polygon,lam,constant,contact['a'])
    cert=dict(contact_row=index,contact_w=w,gate_price_upper=str(Aq),
              beta_d_over_b_upper=str(beta),S_upper=str(S),lambda_value=str(lam),polygon_pricing=True)
    return (cap,cert) if cap is not None else None
