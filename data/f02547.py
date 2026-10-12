#!/usr/bin/env python3
"""Fresh readonly replay of saved physical-moment box trees.

Box reconstruction, exclusions and quadratic radii are independent of the
cover producer. Directed cubic point integration uses its separately audited
Arb root backend. This script makes no stationary-point assumption.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse
import gzip
import hashlib
import importlib.util
import json
from flint import arb,ctx

HERE=Path(__file__).resolve().parent
COVER=HERE.parents[1]/'tonight/recovery/profile_cover'


def exact(x):
    x=F(x)
    return arb(x.numerator)/x.denominator


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--prefix',default='local_physical_d256_inflation1_125')
    args=parser.parse_args()
    ctx.prec=384
    spec=importlib.util.spec_from_file_location('audited_cubic_backend',COVER/'directed_cubic_integral.py')
    backend=importlib.util.module_from_spec(spec);spec.loader.exec_module(backend)
    ctx.prec=384
    summary=json.loads((COVER/f'{args.prefix}_summary.json').read_text())
    assert summary['status']=='PASS'
    p=summary['profile'];m=[[F(x) for x in row] for row in p['M']]
    initial=[tuple(map(F, pair)) for pair in summary['initial_box']]
    c=(2/arb.pi()).sqrt()
    caps=[arb(1),c,2*c*(-arb(1)/2).exp()/arb(2).sqrt(),
          c*(1+4*(-arb(3)/2).exp())/arb(6).sqrt()]
    assert all(-a==b and bound<=exact(b) for (a,b),bound in zip(initial,caps))
    output=[]
    for index in (0,25,50):
        path=COVER/f'{args.prefix}_node{index:03d}.json.gz'
        with gzip.open(path,'rt') as stream:record=json.load(stream)
        assert record['status']=='PASS' and record['pending_boxes']==0
        assert record['q']==summary['t'][index]
        assert record['target']==p['nodal_majorant'][index]
        stack=[initial];accepted=excluded=splits=0;largest=None
        for event in record['events']:
            assert stack
            box=stack.pop()
            minimum=sum((a*a if a>0 else b*b if b<0 else F(0)
                         for a,b in box),F(0))
            if event=='N':
                assert minimum>1
                excluded+=1
            elif event=='C':
                center=[(a+b)/2 for a,b in box]
                widths=[(b-a)/2 for a,b in box]
                radius=sum((abs(m[i][j])*widths[i]*widths[j]
                            for i in range(4) for j in range(4)),F(0))
                value=backend.local_point(p,record['q'],[str(x) for x in center])+exact(radius)
                gap=value-exact(F(record['target']))
                assert gap<=0,('BOX_FAIL',index,box,gap)
                high=F(str(gap.upper().fmpq()))
                largest=high if largest is None else max(largest,high)
                accepted+=1
            else:
                assert event.startswith('S') and len(event)==2
                axis=int(event[1]);assert 0<=axis<4
                a,b=box[axis];mid=(a+b)/2
                left,right=list(box),list(box)
                left[axis]=(a,mid);right[axis]=(mid,b)
                stack.extend([right,left]);splits+=1
        assert not stack
        assert accepted==record['accepted'] and excluded==record['excluded']
        assert accepted+excluded+splits==record['nodes']==len(record['events'])
        output.append({'node':index,'q':record['q'],'nodes':record['nodes'],
                       'accepted':accepted,'excluded':excluded,'splits':splits,
                       'largest_directed_gap':str(largest),
                       'tree_sha256':hashlib.sha256(path.read_bytes()).hexdigest()})
    result={'status':'PASS_FRESH_SELECTED_PHYSICAL_TREE_REPLAY','precision_bits':ctx.prec,
            'nodes':output,'total_boxes':sum(x['nodes'] for x in output),
            'backend_sha256':hashlib.sha256((COVER/'directed_cubic_integral.py').read_bytes()).hexdigest(),
            'scope':'Three complete trees freshly replayed; all 51 producer trees separately passed.'}
    output_name='selected_physical_nodes_replay.json' if args.prefix=='local_physical_d256_inflation1_125' else f'{args.prefix}_selected_replay.json'
    (HERE/output_name).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':result['status'],'total_boxes':result['total_boxes'],
                      'indices':[x['node'] for x in output]},indent=2))


if __name__=='__main__':
    main()
