#!/usr/bin/env python3
"""Exact, dimension-free verification of cubic-chaos cumulant identities.

Connected Wick diagrams are enumerated as loopless multigraphs of
degree3. The right sides are expanded as contractions of a symmetric
cubic tensor. Equality of coefficients of every graph is checked.
Only standard-library exact rational/integer arithmetic is used.
"""
from fractions import Fraction as F
from itertools import permutations, product
from math import factorial
from collections import defaultdict
from pathlib import Path
import json

def graph_key(matrix):
    n=len(matrix)
    return min(tuple(matrix[p[i]][p[j]] for i in range(n) for j in range(i+1,n))
               for p in permutations(range(n)))

def connected(matrix):
    seen={0};queue=[0]
    while queue:
        i=queue.pop()
        for j,v in enumerate(matrix[i]):
            if v and j not in seen:seen.add(j);queue.append(j)
    return len(seen)==len(matrix)

def wick_graphs(n):
    matrix=[[0]*n for _ in range(n)];remaining=[3]*n;out=defaultdict(F)
    labeled=0
    def distribute(i,js,total,allocated):
        if not js:
            if total==0:yield allocated
            return
        j=js[0]
        for v in range(min(total,remaining[j])+1):
            yield from distribute(i,js[1:],total-v,allocated+[(j,v)])
    def run(i):
        nonlocal labeled
        if i==n:
            if any(remaining) or not connected(matrix):return
            den=1
            for a in range(n):
                for b in range(a+1,n):den*=factorial(matrix[a][b])
            out[graph_key(matrix)]+=F(6**n,den);labeled+=1
            return
        if remaining[i]==0:run(i+1);return
        for row in distribute(i,list(range(i+1,n)),remaining[i],[]):
            old=remaining[i];remaining[i]=0
            for j,v in row:matrix[i][j]=matrix[j][i]=v;remaining[j]-=v
            run(i+1)
            for j,v in row:matrix[i][j]=matrix[j][i]=0;remaining[j]+=v
            remaining[i]=old
    run(0)
    return dict(out),labeled

def contraction_key(vertices):
    n=len(vertices);indices=defaultdict(list)
    for i,slots in enumerate(vertices):
        assert len(slots)==3
        for index in slots:indices[index].append(i)
    matrix=[[0]*n for _ in range(n)]
    for index,ends in indices.items():
        assert len(ends)==2,(index,ends)
        i,j=ends;assert i!=j
        matrix[i][j]+=1;matrix[j][i]+=1
    assert connected(matrix)
    return graph_key(matrix)

def block_A(indices,name,pairing=0):
    a,b,c,d=indices
    if pairing==1:b,c=c,b
    if pairing==2:b,c,d=d,b,c
    return [(a,b,name),(c,d,name)]

def trace_product(kinds):
    # Matrix pair-space indices: (a_i,b_i), cyclically identified.
    n=len(kinds);out=defaultdict(F)
    for choices in product(*(range(3) if k=='S' else range(1) for k in kinds)):
        vs=[];weight=F(1)
        for i,(k,choice) in enumerate(zip(kinds,choices)):
            ext=(f'a{i}',f'b{i}',f'a{(i+1)%n}',f'b{(i+1)%n}')
            vs+=block_A(ext,f'e{i}',choice)
            if k=='S':weight/=3
        out[contraction_key(vs)]+=weight
    return dict(out)

def sgg():
    out=defaultdict(F)
    for choice in range(3):
        vs=block_A(('a','b','c','d'),'e',choice)
        # <S,G tensor G>, with G_ab = sum_uv f_auv f_buv.
        vs += [('a','u','v'),('b','u','v'),('c','w','z'),('d','w','z')]
        out[contraction_key(vs)]+=F(1,3)
    return dict(out)

def sumterms(terms):
    out=defaultdict(F)
    for price,t in terms:
        for k,v in t.items():out[k]+=price*v
    return {k:v for k,v in out.items() if v}

def verify():
    w4,l4=wick_graphs(4)
    # s=tr A^2, h=tr S^2; S is the full symmetrization of A.
    rhs4=sumterms([(F(1296),trace_product('AA')),(F(1944),trace_product('SS'))])
    assert w4==rhs4,('fourth',w4,rhs4)
    w6,l6=wick_graphs(6)
    rhs6=sumterms([(F(46656*135,2),trace_product('SSS')),
                   (F(46656*45),trace_product('ASS')),
                   (F(46656*165,2),trace_product('AAS')),
                   (F(-46656*20),trace_product('AAA')),
                   (F(46656*135,2),sgg())])
    assert w6==rhs6,('sixth',w6,rhs6)
    result={'status':'PASS','scope':'dimension-free connected Wick-diagram identities',
      'fourth':{'labeled_connected_degree3_graphs':l4,'isomorphism_classes':len(w4),
                'sum_coefficients':str(sum(w4.values())),
                'coefficients':[{'canonical_edges':list(k),'coefficient':str(v)} for k,v in sorted(w4.items())]},
      'sixth':{'labeled_connected_degree3_graphs':l6,'isomorphism_classes':len(w6),
               'sum_coefficients':str(sum(w6.values())),
               'coefficients':[{'canonical_edges':list(k),'coefficient':str(v)} for k,v in sorted(w6.items())]}}
    Path(__file__).with_name('wick_identity_certificate.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v if k not in ('fourth','sixth') else {a:b for a,b in v.items() if a!='coefficients'} for k,v in result.items()}))

if __name__=='__main__':verify()
