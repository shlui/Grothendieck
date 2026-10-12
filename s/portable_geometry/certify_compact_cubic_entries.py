#!/usr/bin/env python3
"""Independent exact coefficient expansion for the compact cubic formulas."""
from fractions import Fraction as F
from itertools import combinations
import hashlib
import json
from pathlib import Path


def polynomial_add(poly, monomial, value):
    poly[monomial] = poly.get(monomial, F(0)) + value
    if poly[monomial] == 0:
        del poly[monomial]


def product(poly1, poly2):
    out = {}
    for m1, c1 in poly1.items():
        for m2, c2 in poly2.items():
            polynomial_add(out, tuple(a + b for a, b in zip(m1, m2)), c1*c2)
    return out


def scaled(poly, coefficient):
    return {m: c*coefficient for m, c in poly.items() if c*coefficient}


def add(*polys):
    out = {}
    for p in polys:
        for m, c in p.items():
            polynomial_add(out, m, c)
    return out


def variable(d, j):
    exponent = [0]*(2*d)
    exponent[j] = 1
    return {tuple(exponent): F(1)}


def pow_poly(p, n, d):
    out = {(0,)*(2*d): F(1)}
    for _ in range(n):
        out = product(out, p)
    return out


def check_kernel(d):
    x = [variable(d, j) for j in range(d)]
    y = [variable(d, j+d) for j in range(d)]
    dot = add(*(product(xj,yj) for xj,yj in zip(x,y)))
    rx = add(*(pow_poly(v,2,d) for v in x))
    ry = add(*(pow_poly(v,2,d) for v in y))
    compact = scaled(add(pow_poly(dot,3,d),
                         scaled(product(dot,add(rx,ry)),-3),
                         scaled(dot,3*(d+2))), F(1,6))
    hermite = {}
    for j in range(d):
        hx = add(pow_poly(x[j],3,d), scaled(x[j],-3))
        hy = add(pow_poly(y[j],3,d), scaled(y[j],-3))
        hermite = add(hermite,scaled(product(hx,hy),F(1,6)))
    one = {(0,)*(2*d): F(1)}
    for i in range(d):
        for j in range(d):
            if i == j:
                continue
            hx = product(add(pow_poly(x[i],2,d),scaled(one,-1)),x[j])
            hy = product(add(pow_poly(y[i],2,d),scaled(one,-1)),y[j])
            hermite = add(hermite,scaled(product(hx,hy),F(1,2)))
    for indices in combinations(range(d),3):
        px, py = one, one
        for j in indices:
            px, py = product(px,x[j]), product(py,y[j])
        hermite = add(hermite,product(px,py))
    assert compact == hermite
    return len(compact)


def check_conditioning(d):
    # Algebraic conditional moments may be arbitrary rational numbers.
    a = [[F((7*j+3*k+1)%19-9,k+2) for k in range(3)] for j in range(d)]
    b = [[F((11*j+5*k+2)%23-11,k+3) for k in range(3)] for j in range(d)]
    u = [aa[0]*bb[0] for aa,bb in zip(a,b)]
    e = [(aa[1]-1)*(bb[1]-1) for aa,bb in zip(a,b)]
    t = [(aa[2]-3*aa[0])*(bb[2]-3*bb[0]) for aa,bb in zip(a,b)]
    s = sum(u,F(0))
    compact = (sum(t,F(0))+3*s*sum(e,F(0))
               -3*sum((ei*ui for ei,ui in zip(e,u)),F(0))
               +s**3-3*s*sum((ui**2 for ui in u),F(0))
               +2*sum((ui**3 for ui in u),F(0)))/6
    direct = sum(t,F(0))/6
    direct += sum((e[i]*u[j]/2 for i in range(d) for j in range(d) if i!=j),F(0))
    direct += sum((u[i]*u[j]*u[k] for i,j,k in combinations(range(d),3)),F(0))
    assert compact == direct


def main():
    here = Path(__file__).resolve().parent
    fragment = here/'COMPACT_CUBIC_CELL_ENTRIES.tex'
    result = {
        'status':'PASS_EXACT_CUBIC_KERNEL_AND_PRODUCT_CELL_ALGEBRA',
        'symbolic_coefficient_expansions':
            [{'dimension':d,'nonzero_compact_coefficients':check_kernel(d)}
             for d in range(1,8)],
        'conditional_pattern_expansions': [check_conditioning(d) or d for d in range(1,33)],
        'arithmetic':'fractions.Fraction; exact polynomial coefficient dictionaries',
        'producer_imports':False,
        'fragment_sha256':hashlib.sha256(fragment.read_bytes()).hexdigest(),
        'entrypoint_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'scope':'Exact algebra; the displayed general proof supplies arbitrary-dimensional validity.'
    }
    output=here/'COMPACT_CUBIC_CELL_ENTRY_CERTIFICATE.json'
    output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
