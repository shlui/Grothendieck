"""Portable exact Bernstein bounds, using a common integer denominator.

No numerical sampling or floating input is used.  The Bernstein coefficients
are constructed by polynomial Horner evaluation in the Bernstein basis.
"""
from fractions import Fraction
from functools import lru_cache
from math import lcm


def exact(x):
    if not isinstance(x,(int,str,Fraction)):
        raise TypeError('Exact integer, decimal string, or Fraction required')
    return Fraction(x)


@lru_cache(maxsize=256)
def integer_polynomial(coefficients):
    p=list(coefficients)
    while len(p)>1 and not p[-1]:p.pop()
    common=1
    for x in p:common=lcm(common,x.denominator)
    return tuple(x.numerator*(common//x.denominator) for x in p),common


def bernstein_integer_coefficients(polynomial,a,b):
    """Return (integer Bernstein coefficients, positive common denominator).

    For B of degree j-1, multiplication by q=a(1-t)+bt gives degree-j
    coefficients a*(j-k)/j*B[k]+b*k/j*B[k-1].  Adding a constant adds
    it to every coefficient.  A common denominator S*d**j*j! keeps
    this entire recurrence in integers.
    """
    p=tuple(exact(x) for x in polynomial)
    if not p:raise ValueError('Empty polynomial')
    a,b=exact(a),exact(b)
    if a>b:raise ValueError('Reversed interval')
    ns,S=integer_polynomial(p)
    d=lcm(a.denominator,b.denominator)
    A=a.numerator*(d//a.denominator)
    B=b.numerator*(d//b.denominator)
    v=[ns[-1]];factor=1
    for j,c in enumerate(ns[-2::-1],1):
        factor*=d*j
        constant=c*factor
        v=[constant+(A*(j-k)*v[k] if k<j else 0)
                    +(B*k*v[k-1] if k else 0) for k in range(j+1)]
    return tuple(v),S*factor


def bernstein_bounds(polynomial,a,b):
    ns,den=bernstein_integer_coefficients(polynomial,a,b)
    return Fraction(min(ns),den),Fraction(max(ns),den)


if __name__=='__main__':
    # A simple identity: q^2 on [0,1] has Bernstein coefficients 0,0,1.
    assert bernstein_integer_coefficients([0,0,1],0,1)==((0,0,2),2)
    assert bernstein_bounds([0,0,1],0,1)==(0,1)
    print('PASS_EXACT_INTEGER_BERNSTEIN_SMOKE_CHECK')
