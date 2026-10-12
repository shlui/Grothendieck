#!/usr/bin/env python3
"""Producer-free source-box prices for operator's accepted targeted gates."""
from pathlib import Path
from fractions import Fraction as F
from itertools import product
import argparse,json,sys
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path.insert(0,str(HERE/'global_b7171'))
import replay_R3_corners as R

def main(source,gatepath,output):
    d=R.load(source);g=R.load(gatepath)
    assert g['status'].startswith('PASS_INDEPENDENT_COMPLETE_')
    profile=ROOT/g['profile_path'];p=R.exact(R.load(profile)['profile'])
    assert R.sha(profile)==g['profile_sha256']
    rows=[]
    for box in d['boxes']:
        if box['profile']!=profile.name:continue
        assert box['kind']=='rank' and (F(box['reserve0']),F(box['reserve6']))==tuple(map(F,g['reserves']))
        for path,digest in box['sha256'].items():assert R.sha(ROOT/path)==digest
        bounds={k:tuple(map(F,box[k]))for k in ('u','r','kappa','shape')}
        assert all(len(v)==2 and v[0]<=v[1]for v in bounds.values())
        row=dict(kind='rank',p=p,reserve0=F(box['reserve0']),reserve6=F(box['reserve6']))
        old={tuple(F(v[k])for k in ('u','r','kappa','shape')):F(v['upper'])for v in box['vertices']}
        assert len(old)==16
        results=[]
        for b in (F(7171,10000),F(713,1000),F(141,200)):
            spectral=R.spectral(row,b);c=R.coefficients(row,b);assert c['A']>0 and c['B']>0
            corners=[];maximum=F(0)
            for u,r,k,z in product(*(bounds[x]for x in ('u','r','kappa','shape'))):
                val=R.price(row,c,k,z,u,r)
                if b==F(box['operator_b']):
                    assert val.below(old[u,r,k,z]) and val.below(F(box['cap']))
                rec=val.record();maximum=max(maximum,F(rec['upper']))
                corners.append(dict(u=str(u),r=str(r),kappa=str(k),shape=str(z),price=rec))
            ceil=F(-((-maximum.numerator*10**9)//maximum.denominator),10**9)
            results.append(dict(operator_b=str(b),spectral=spectral,clean_rational_upper=str(ceil),corners=corners))
        rows.append(dict(name=box['name'],source_rectangle={k:list(map(str,v))for k,v in bounds.items()},results=results))
    assert len(rows)==1
    for part,v in g['reports'].items():assert R.sha(ROOT/v['path'])==v['sha256']
    out=dict(status='PASS_INDEPENDENT_TARGETED_RANK_REGIONAL_BOX_PRICES',source=R.file_record(source),
             accepted_whole_q_gate=R.file_record(gatepath),profile=R.file_record(profile),engine=R.selftest(),
             boxes=rows,directed_prices=48,
             convexity='All target-b PSD gates hold exactly. Positive amplitude/source curvatures and A>0,B>0 cumulant/shape coefficients prove the sixteen-vertex enclosure at each common operator coefficient.',
             scope='Only the displayed source boxes; exhaustive global amplitude/source coverage is separate.')
    output.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'status':out['status'],'directed_prices':48,'caps':[{k:v for k,v in x.items()if k not in ('spectral','corners')}for x in rows[0]['results']]},indent=2))
if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('source',type=Path);a.add_argument('gate',type=Path);a.add_argument('--output',type=Path,required=True)
    x=a.parse_args();main(x.source.resolve(),x.gate.resolve(),x.output.resolve())
