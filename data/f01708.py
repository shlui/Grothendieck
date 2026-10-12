#!/usr/bin/env python3
"""Independent directed prices and exact single-b assembly of four old gates."""
from fractions import Fraction as F
from pathlib import Path
import json,time
import replay_native_R3_unions as N
from mpfr_interval import I,BITS

ledger=N.COVER/'native_b7171_regional_transfer_assembly.json'
d=json.loads(ledger.read_text());assert d['operator_b']=='7171/10000'
b=F(7171,10000);out=[];nprices=0
for row in d['rows'][:4]:
    name=row['package'];p=N.package(name);old=json.loads((p['path']/'replay.json').read_text());M=p['M'];prof=p['profile'];rg=p['rg']
    assert old['status']=='PASS_COMPLETE_REGIONAL_BOUND'
    lo,hi=map(F,row['amplitude_interval']);assert [str(lo),str(hi)]==old['amplitude_interval']
    assert rg[-1]>=N.source_upper(lo,F(847,1000))
    beta,ell,nu,delta=map(F,(prof['beta'],prof['ell'],prof['nu'],prof['delta']))
    g=(nu-M[0][2])**2/M[2][2]-M[0][0]+b-2*delta
    shift=beta**2*(1/b-1/F(prof['epsilon']))
    values=[(F(prof['nodal_majorant'][0])+shift,'r=0')]
    assert len(row['bindings']['scalars'])==len(p['source']['rows'])
    for index in range(len(p['source']['rows'])):
        path=p['path']/f'data/scalars/row{index:02d}.json';cert=json.loads(path.read_text());rr=N.scalar(name,index)
        assert N.sha(path)==row['bindings']['scalars'][path.name]
        ki=list(map(F,cert['frozen_row']['kappa_interval']));expected={(k,s)for k in ki for s in(F(0),F(1))}
        assert {(F(c['kappa']),F(c['shape']))for c in cert['prices']}==expected
        for c in cert['prices']:
            value=N.price(name,index,F(c['kappa']),F(c['shape']),F(849,1000),'constant')
            assert value.below(F(c['price_upper'])),(name,index,c)
            nprices+=1
        r=rr['r'];value=max(F(c['price_upper'])for c in cert['prices'])+g*r*r+p['inflation']+p['tail']*F(849,1000)*rg[-1]+shift
        values.append((value,path.name))
    # The row bands are separate at each positive r, and cover all kappa.
    for r in rg[1:]:
        bands=sorted(tuple(map(F,x['kappa_interval']))for x in p['source']['rows']if F(str(x['r']))==r)
        assert bands[0][0]==0 and bands[-1][1]==90 and all(x[1]==y[0]for x,y in zip(bands,bands[1:]))
    base,active=max(values);reserve=max(F(0),-g)*max(y-x for x,y in zip(rg,rg[1:]))**2/4
    total=base+reserve+max(lo*lo-2*ell*lo,hi*hi-2*ell*hi)
    assert total==F(row['directed_total_upper'])<F(40,71)
    assert g==F(row['exact_gram_cost']) and shift==F(row['local_constant_transport']) and reserve==F(row['source_chord_reserve'])
    assert active==row['active_row'] and N.sha(p['path']/'replay.json')==row['bindings']['replay_sha256']
    out.append({'package':name,'amplitude_interval':list(map(str,(lo,hi))),'total':str(total),'display':float(total),'active_row':active,'g':str(g),'shift':str(shift),'reserve':str(reserve),'manifest_sha256':p['manifest_sha256']})
# These overlap through .638; elementary low-u energy covers through .21.
assert all(F(out[i+1]['amplitude_interval'][0])<=F(out[i]['amplitude_interval'][1])for i in range(3))
energy=F(21,100)**2+b*F(847,1000)**2;assert energy<F(40,71)
result={'status':'PASS_INDEPENDENT_DIRECTED_FIRST_FOUR_B7171_TRANSFERS','operator_b':str(b),'coverage':['0','319/500'],'threshold':'40/71','transfers':out,'low_energy_end':'21/100','low_energy_cap':str(energy),'endpoint_prices_recomputed':nprices,'precision_bits':BITS,'source_rho_in_transferred_gates':'849/1000','source_rho_in_energy':'847/1000','input_sha256':N.sha(ledger),'producer_modules_imported':False}
path=Path(__file__).parent/'native_b7171_union/first_four_transfers_independent.json';path.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items()if k!='transfers'},indent=2))
