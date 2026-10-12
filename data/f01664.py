#!/usr/bin/env python3
"""Exact maximum epsilon for fixed 2x2 profile blocks and fixed delta.

Reads either {"profile":{"M":4x4,"delta":...}} or
{"E":2x2,"R":2x2,"delta":...}. Decimal inputs are exact rationals.
This checks matrix gates only, never the scalar potential or amplitude cover.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse, hashlib, json

def rational(x):
    return x if isinstance(x,F) else F(str(x))

def check(data,epsilon):
    q=data.get('profile',data)
    delta=rational(q['delta'])
    if 'M' in q:
        A=[[rational(x) for x in row] for row in q['M']]
        assert len(A)==4 and all(len(row)==4 for row in A)
        assert all(A[i][j]==A[j][i] for i in range(4) for j in range(4))
        assert all(A[i][j]==0 for i in (0,2) for j in (1,3))
        E=[[A[0][0],A[0][2]],[A[2][0],A[2][2]]]
        R=[[A[1][1],A[1][3]],[A[3][1],A[3][3]]]
    else:
        E=[[rational(x) for x in row] for row in q['E']]
        R=[[rational(x) for x in row] for row in q['R']]
    for B in (E,R):
        assert B[0][1]==B[1][0]
        assert B[0][0]>0 and B[1][1]>0
        assert B[0][0]*B[1][1]-B[0][1]**2>0
    assert delta>0
    caps={
        'E_minus_diag_minus1_epsilon':E[1][1]-E[0][1]**2/(E[0][0]+1),
        'E_minus_diag_epsilon_minus_delta_0':delta+E[0][0]-E[0][1]**2/E[1][1],
        'R_minus_diag_minus1_epsilon':R[1][1]-R[0][1]**2/(R[0][0]+1),
        'R_minus_diag_epsilon_0':R[0][0]-R[0][1]**2/R[1][1]}
    maximum=min(caps.values())
    return {'status':'PASS_FIXED_GATE_ARITHMETIC_ONLY',
            'epsilon_tested':str(epsilon),
            'epsilon_maximum_exact':str(maximum),
            'epsilon_maximum_float_display':float(maximum),
            'all_four_gates_pass':epsilon<=maximum,
            'gate_margins_exact':{k:str(v-epsilon) for k,v in caps.items()},
            'limiting_gate':min(caps,key=caps.get),
            'scope':'Fixed matrix gates only. No potential or amplitude certificate is implied.'}

def main():
    p=argparse.ArgumentParser();p.add_argument('input');p.add_argument('--epsilon',default='.72');p.add_argument('--output')
    a=p.parse_args();path=Path(a.input)
    raw=path.read_bytes();data=json.loads(raw,parse_float=F)
    result=check(data,F(a.epsilon));result['input_sha256']=hashlib.sha256(raw).hexdigest()
    encoded=json.dumps(result,indent=2)+'\n'
    if a.output:Path(a.output).write_text(encoded)
    print(encoded)

if __name__=='__main__':main()
