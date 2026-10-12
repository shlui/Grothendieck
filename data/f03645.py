#!/usr/bin/env python3
"""Exact complete low-amplitude target transport, retaining radial loss.

Accepted baseline local/integration/partition replays are hash-bound. This
checker performs no producer imports and exhausts every accepted lower cell.
"""
from pathlib import Path
from fractions import Fraction as F
from collections import defaultdict
import argparse,json,os,sys
os.environ['GR_MPFR_BITS']='768'
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/'tonight/audit'))
import replay_native_R3_unions as N
sys.path.insert(0,str(Path(__file__).resolve().parent))
import replay_R3_corners as R

def main(b,output):
    b0=F(7171,10000);assert 0<b<=b0;T=F(40,71);delta=1/b-1/b0
    bindings={};profiles={}
    def accepted(path):
        p=ROOT/path;d=json.loads(p.read_text());assert d['status'].startswith('PASS')
        bindings[path]={'sha256':N.sha(p),'status':d['status']};return d
    def native(name):
        if name in profiles:return F(profiles[name]['beta_squared'])
        d=N.package(name);p=d['profile'];M=d['M'];beta=F(p['beta'])
        for a,o,c in [(M[0][0]+1,M[0][2],M[2][2]-b),
                      (M[1][1]+1,M[1][3],M[3][3]-b),(M[1][1]-b,M[1][3],M[3][3])]:N.psd(a,o,c)
        profiles[name]={'beta_squared':str(beta*beta),'manifest_sha256':d['manifest_sha256'],
                        'source_sha256':d['source_sha256'],'replay_sha256':d['replay_sha256']}
        return beta*beta
    first=accepted('tonight/audit/native_b7171_union/first_four_transfers_independent.json')
    assert first['coverage']==['0','319/500'] and F(first['operator_b'])==b0
    transferpath=ROOT/'tonight/recovery/profile_cover/native_b7171_regional_transfer_assembly.json'
    assert N.sha(transferpath)==first['input_sha256']
    raw=json.loads(transferpath.read_text())['rows'][:4];early=[]
    for x,y in zip(raw,first['transfers']):
        assert x['package']==y['package'] and F(x['directed_total_upper'])==F(y['total'])
        native(y['package']);pkg=N.package(y['package']);prof=pkg['profile'];M=pkg['M'];rg=pkg['rg']
        beta,ell,nu,dd=map(F,(prof['beta'],prof['ell'],prof['nu'],prof['delta']))
        g=(nu-M[0][2])**2/M[2][2]-M[0][0]+b-2*dd
        shift=beta**2*(1/b-1/F(prof['epsilon']))
        values=[(F(prof['nodal_majorant'][0])+shift,'r=0')];fresh=0
        for index in range(len(pkg['source']['rows'])):
            path=pkg['path']/f'data/scalars/row{index:02d}.json';cert=json.loads(path.read_text());rr=N.scalar(y['package'],index)
            assert N.sha(path)==x['bindings']['scalars'][path.name]
            ki=list(map(F,cert['frozen_row']['kappa_interval']))
            assert {(F(v['kappa']),F(v['shape']))for v in cert['prices']}=={(k,z)for k in ki for z in(F(0),F(1))}
            for v in cert['prices']:
                actual=N.price(y['package'],index,F(v['kappa']),F(v['shape']),F(849,1000),'constant')
                assert actual.below(F(v['price_upper']));fresh+=1
            r=rr['r'];value=max(F(v['price_upper'])for v in cert['prices'])+g*r*r+pkg['inflation']+pkg['tail']*F(849,1000)*rg[-1]+shift
            values.append((value,path.name))
        for r in rg[1:]:
            bands=sorted(tuple(map(F,v['kappa_interval']))for v in pkg['source']['rows']if F(str(v['r']))==r)
            assert bands[0][0]==0 and bands[-1][1]==90 and all(a[1]==c[0]for a,c in zip(bands,bands[1:]))
        base,active=max(values);reserve=max(F(0),-g)*max(c-a for a,c in zip(rg,rg[1:]))**2/4
        lo,hi=map(F,y['amplitude_interval']);assert rg[-1]>=N.source_upper(lo,F(847,1000))
        cap=base+reserve+max(lo*lo-2*ell*lo,hi*hi-2*ell*hi)
        assert cap<T
        early.append({'package':y['package'],'amplitude_interval':y['amplitude_interval'],
                      'upper':str(cap),'radial_negative_transport_discarded':False,
                      'gram_cost':str(g),'local_completion_shift':str(shift),'chord_reserve':str(reserve),
                      'active_row':active,'fresh_directed_corner_prices':fresh})
    assert len(early)==4
    intervals=list(map(lambda x:tuple(map(F,x['amplitude_interval'])),early))
    assert intervals[0][0]<=F(21,100) and intervals[-1][1]==F(319,500)
    assert all(c[0]<=a[1]for a,c in zip(intervals,intervals[1:]))
    rank=accepted('tonight/audit/diffuse_b7171_R3/independent_extended_remainder_box.json')
    rp=ROOT/'tonight/operator/projection/lower_amp_certificate/profile_u.70.json'
    assert R.sha(rp)==rank['profile_sha256'] and F(rank['operator_b'])==b0
    for v in rank['checked_hash_ledger']:assert N.sha(ROOT/v['path'])==v['sha256']
    p=R.exact(json.loads(rp.read_text())['profile']);rankbeta=6*p['beta']**2
    rankspectral=R.spectral(dict(kind='rank',p=p),b)
    region_paths=[('tonight/audit/native_b7171_union/u600_682_independent_replay.json',F(3,5),F(341,500)),
                  ('tonight/audit/native_b7171_union/u682_6905_independent_replay.json',F(341,500),F(1381,2000))]
    regions=[];grandmax=max(F(x['upper'])for x in early);count=other=0
    for reportpath,lo,hi in region_paths:
        a=accepted(reportpath);assert F(a['b'])==b0 and tuple(map(F,a['interval']))==(lo,hi)
        source=ROOT/a['path'];assert N.sha(source)==a['input_sha256'];d=json.loads(source.read_text())
        frames={v['mask']:v for v in d['bindings']}
        for v in frames.values():
            assert N.sha(ROOT/v['mask'])==v['mask_sha256']
            native(v['package']);assert profiles[v['package']]['source_sha256']==v['source_sha256']
        groups=defaultdict(list);localmax=F(0);nc=rc=0;active=None
        for cell in d['cells']:
            u=tuple(map(F,cell['u_interval']));k=tuple(map(F,cell['kappa_interval']));r=tuple(map(F,cell['r_interval']))
            if not lo<=u[0]<u[1]<=hi:continue
            assert 0<=r[0]<r[1] and 0<=k[0]<k[1]<=90
            groups[u,k].append(r)
            if cell['accepted']:
                oldcap=F(cell['cap']);assert oldcap<T
                cap=oldcap+native(frames[cell['frame']]['package'])*delta+(b-b0)*r[0]**2;nc+=1
            else:
                assert lo==F(341,500)
                for axis,v in [('u',u),('r',r),('kappa',k)]:
                    a0,a1=map(F,rank[axis]);assert a0<=v[0]<v[1]<=a1
                cap=F(rank['cap'])+rankbeta*delta+(b-b0)*r[0]**2;rc+=1
            assert cap<T,('TARGET_TRANSPORT_FAILURE',b,u,k,r,str(cap))
            if cap>localmax:localmax=cap;active=cell
        up=sorted({u for u,k in groups});ks=sorted({k for u,k in groups})
        assert up[0][0]==lo and up[-1][1]==hi and all(x[1]==y[0]for x,y in zip(up,up[1:]))
        assert ks[0][0]==0 and ks[-1][1]==90 and all(x[1]==y[0]for x,y in zip(ks,ks[1:]))
        assert len(groups)==len(up)*len(ks)==a['source_band_groups']
        for (u,k),parts in groups.items():
            parts.sort();assert parts[0][0]==0 and all(x[1]==y[0]for x,y in zip(parts,parts[1:]))
            assert parts[-1][1]==N.source_upper(u[0],F(847,1000))
        assert nc==a['native_cells'] and rc==a['R3_cells']
        count+=nc;other+=rc;grandmax=max(grandmax,localmax)
        regions.append({'amplitude_interval':list(map(str,(lo,hi))),'source':R.file_record(source),
                        'complete_groups':len(groups),'native_cells':nc,'R3_cells':rc,
                        'upper':str(localmax),'active_cell':active})
    assert count==39929 and other==41
    energy=F(21,100)**2+b*F(847,1000)**2;assert energy<T
    grandmax=max(grandmax,energy)
    ceil=F(-((-grandmax.numerator*10**12)//grandmax.denominator),10**12);assert ceil<T
    out={'status':'PASS_INDEPENDENT_COMPLETE_SMALLER_B_LOWER_TARGET_BRANCH',
         'operator_b':str(b),'baseline_b':str(b0),'threshold':'40/71',
         'complete_amplitude_coverage':['0','1381/2000'],'clean_rational_upper':str(ceil),
         'exact_upper':str(grandmax),'upper_display':float(grandmax),
         'source_rho':'847/1000','energy_end':'21/100','energy_cap':str(energy),
         'first_four_regions':early,'regions':regions,'native_cells':count,'R3_cells':other,
         'profiles':profiles,'rank70_profile':R.file_record(rp),'rank70_spectral':rankspectral,
         'bindings':bindings,'producer_modules_imported':False,
         'analytic_review':'The same physical partition and frozen local inequalities remain valid at the single smaller coefficient. Constant transport is beta_norm^2*(1/b-1/b0); radial transport is (b-b0)r^2, whose maximum on each nonnegative source interval occurs at its lower endpoint. Adding this maximum to an already independently accepted baseline cell cap is a valid conservative bound. All cells, all source ranges, all cumulant groups and every R3 remainder containment are checked; early complete regional profiles are freshly repriced at every source endpoint with their actual negative radial coefficient and the exact adjusted chord reserve. No coefficient choices are mixed.',
         'scope':'Complete low-amplitude branch only; the global norm also requires a complete same-coefficient proof on u≥.6905.'}
    output.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items()if k not in ('bindings','profiles','rank70_spectral','regions','first_four_regions')},indent=2))

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--b',required=True,type=F);a.add_argument('--output',required=True,type=Path)
    x=a.parse_args();main(x.b,x.output)
