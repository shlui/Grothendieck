#!/usr/bin/env python3
"""Exact physical-source bank from cubic-contact homogeneity, on actual gaps.

All base rows use the accepted global analytic contact family at B=0. Exact
positive rational scaling generates the bank without new all-line tests.
Fixed-row convex-kappa/affine-shape endpoint prices certify each frozen box.
"""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import sys

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
FAMILY=HERE.parent/'cubic_head_contact'
sys.path.insert(0,str(FAMILY))
import certify_direct_contact_majorant as C

BASE_A=F(147,50)
CONTACTS=[F(1,2),F(1),F(2),F(11,2),F(10),F(20)]
HEADS=[F(7,4),F(2),F(9,4),F(5,2),F(11,4),F(147,50),F(13,4),F(7,2),F(4),F(5)]
LEDGERS=['tonight/recovery/profile_cover/native_union_b7171_target.json',
         'tonight/recovery/profile_cover/native_upper_b7171_target_union.json']


def main():
    C.A0=BASE_A
    C.B0=F(0)
    bases=[C.row(w) for w in CONTACTS]
    rows=[]
    for a in HEADS:
        scale=a/BASE_A
        for i,base in enumerate(bases):
            q=list(map(F,base['polynomial_coefficients_q0_q2_q4_q6']))
            qs=[q[j]*scale**(1-2*j) for j in range(4)]
            assert qs[3]>0
            rows.append({'a':str(a),'base_row':i,'scale':str(scale),
                        'w':base['w'],'contact_s':str(F(base['w'])*a*a),
                        'polynomial_coefficients_q0_q2_q4_q6':list(map(str,qs))})
    cache={}
    def endpoint_price(row,k,shape):
        key=(row,str(k),shape)
        if key not in cache:
            q=rows[row]['polynomial_coefficients_q0_q2_q4_q6']
            cache[key]=C.upper(C.price(q,F(k),shape),10**20)
        return cache[key]
    ledger_results=[]
    source_bindings=[]
    for relative in LEDGERS:
        source=ROOT/relative
        data=json.loads(source.read_text())
        source_bindings.append({'path':relative,'sha256':hashlib.sha256(source.read_bytes()).hexdigest()})
        boxes=[]
        for box in data['unclosed_boxes']:
            ulo,uhi=map(F,box['u_interval'])
            rlo,rhi=map(F,box['r_interval'])
            klo,khi=map(F,box['kappa_interval'])
            candidates=[]
            for row in range(len(rows)):
                corners=[endpoint_price(row,k,shape) for k in [klo,khi]for shape in [0,1]]
                cap=max(corners)-F(rows[row]['a'])*ulo
                candidates.append((cap,row,corners))
            cap,row,corners=min(candidates,key=lambda t:t[0])
            result=dict(box)
            result.update({'selected_row':row,'uniform_source_upper':str(cap),
                           'moment_price_corners':list(map(str,corners)),
                           'excluded':cap<rlo,'clipped_r_upper':str(min(rhi,cap)),
                           'source_clip_gain':str(max(F(0),rhi-cap))})
            boxes.append(result)
        ledger_results.append({'source':relative,'boxes':boxes,
                              'box_count':len(boxes),
                              'excluded_count':sum(b['excluded']for b in boxes),
                              'source_reduced_count':sum(F(b['clipped_r_upper'])<F(b['r_interval'][1])for b in boxes)})
    prices=[{'row':row,'kappa':k,'shape':shape,'upper':str(bound)}
            for (row,k,shape),bound in sorted(cache.items())]
    result={'status':'PASS_EXACT_CONTACT_PHYSICAL_SOURCE_BANK_AND_FROZEN_BOX_EXCLUSIONS',
            'bits':C.ctx.prec,'base_A':str(BASE_A),'base_B':'0',
            'base_rows':bases,'rows':rows,'all_endpoint_prices':prices,
            'source_ledger_bindings':source_bindings,'ledgers':ledger_results,
            'homogeneity_proof':'J_(scale*A)(x)=scale*J_A(x/scale), scale>0, '
                  'so an accepted polynomial majorant coefficient q_(2j) '
                  'scales by scale^(1-2j). Every scaled row is global.',
            'continuous_box_proof':'One fixed row is selected per entire box. '
                  'Its q6>0 moment price is convex in kappa and affine in '
                  'shape, so four corners bound the whole moment rectangle. '
                  'Its physical r-plane is price-a*u, decreasing in u. '
                  'The source cap at u_lo is therefore uniform; only boxes '
                  'with source cap strictly below r_lo are excluded.',
            'scope':'Additional physical source constraints/exclusions only; '
                  'the remaining global norm assembly is separate.',
            'analytic_family_sha256':hashlib.sha256((FAMILY/'certify_direct_contact_majorant.py').read_bytes()).hexdigest(),
            'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (HERE/'source_bank.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':result['status'],'base_count':len(bases),
             'row_count':len(rows),'directed_endpoint_price_count':len(prices),
             'ledgers':[{k:x[k]for k in ['source','box_count','excluded_count','source_reduced_count']}for x in ledger_results]},indent=2))


if __name__=='__main__':
    main()
