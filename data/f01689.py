#!/usr/bin/env python3
"""Hash-bind freshly replayed rank q-gates without producer imports."""
from pathlib import Path
from fractions import Fraction as F
import argparse, hashlib, json

ROOT=Path(__file__).resolve().parents[2]
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def rel(p): return str(p.resolve().relative_to(ROOT))

def main(profile, reports, label, adapter, branch, r0, r6, end, aliases):
    source=json.loads(profile.read_text())
    assert source['kind']=='EXACT_RATIONAL_R3_PROFILE'
    out={'status':'PASS_INDEPENDENT_COMPLETE_'+label+'_R3_SCALAR_GATE',
         'profile_path':rel(profile),'profile_sha256':sha(profile),
         'reserves':[str(r0),str(r6)],
         'q_cover':[['0','8'],['8',str(end)],[str(end),'infinity']],
         'reports':{},'adapter_sha256':sha(adapter),
         'adapter_path':rel(adapter),
         'adapter_review':'Complete compact trees and physical exclusions freshly checked; exact rational tail regenerated. This adapter changes only literal reserves relative to the independently reviewed rank replay.',
         'scope':'Complete whole-q physical local gate only. Outer pricing and exhaustive amplitude/source coverage require separate acceptance.'}
    for part in ('body','exterior','tail'):
        matches=list(reports.glob('*_'+branch+'_'+part+'.json'))
        assert len(matches)==1
        p=matches[0]; d=json.loads(p.read_text())
        assert d['status']==('PASS_INDEPENDENT_EXACT_R3_TAIL' if part=='tail' else 'PASS_INDEPENDENT_PHYSICAL_Q_TREE')
        if part!='tail': assert d['profile_sha256']==sha(profile)
        else: assert F(d['R'])==end and all(F(x)>=0 for x in d['shifted_coefficients'])
        out['reports'][part]={'path':rel(p),'sha256':sha(p),'status':d['status'],
            'accepted':d.get('accepted'),'excluded':d.get('excluded')}
    if aliases:
        records=json.loads(aliases.read_text())
        for alias,v in records.items():
            assert sha(ROOT/alias)==sha(ROOT/v['original'])==v['sha256']
        out['byte_identical_input_aliases']=records
    (reports/'ACCEPTANCE.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--profile',type=Path,required=True)
    a.add_argument('--reports',type=Path,required=True);a.add_argument('--label',required=True)
    a.add_argument('--adapter',type=Path,required=True);a.add_argument('--branch',required=True)
    a.add_argument('--reserve0',type=F,required=True);a.add_argument('--reserve6',type=F,required=True)
    a.add_argument('--end',type=F,required=True);a.add_argument('--aliases',type=Path)
    x=a.parse_args();main(x.profile.resolve(),x.reports.resolve(),x.label,x.adapter.resolve(),x.branch,x.reserve0,x.reserve6,x.end,x.aliases)
