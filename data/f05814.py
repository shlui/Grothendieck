"""Minimal directed real intervals backed by system MPFR (LP64 platforms).
Every public Interval input is an integer, exact decimal string, or Fraction.
No floating-point input is accepted. MPFR errors/NaNs are checked in conversion.
"""
from __future__ import annotations
import ctypes as C
import ctypes.util
from fractions import Fraction
import os

if C.sizeof(C.c_long) != 8 or C.sizeof(C.c_void_p)!=8:
    raise RuntimeError('This MPFR ABI adapter supports LP64 systems only.')
class _Raw(C.Structure):
    _fields_=[('prec',C.c_long),('sign',C.c_int),('exp',C.c_long),('limbs',C.c_void_p)]
_P=C.POINTER(_Raw)
_lib=C.CDLL(ctypes.util.find_library('mpfr') or 'libmpfr.so.6')
def _bind(name,args,restype=C.c_int):
    fn=getattr(_lib,'mpfr_'+name);fn.argtypes=args;fn.restype=restype;return fn
_init=_bind('init2',[_P,C.c_long],None);_clear=_bind('clear',[_P],None)
_set=_bind('set_str',[_P,C.c_char_p,C.c_int,C.c_int]);_cmp=_bind('cmp',[_P,_P]);_nan=_bind('nan_p',[_P])
_ops={n:_bind(n,[_P,_P,_P,C.c_int]) for n in ['add','sub','mul','div']}
_un={n:_bind(n,[_P,_P,C.c_int]) for n in ['neg','sqrt','exp','log','erf','erfc','atan']}
_pi=_bind('const_pi',[_P,C.c_int])
_get=_bind('get_str',[C.c_void_p,C.POINTER(C.c_long),C.c_int,C.c_size_t,_P,C.c_int],C.c_void_p)
_free=_bind('free_str',[C.c_void_p],None)
_getd=_bind('get_d',[_P,C.c_int],C.c_double)
_version=_bind('get_version',[],C.c_char_p)
BITS=int(os.environ.get('GR_MPFR_BITS','256'))
DOWN,UP=3,2
class Real:
    __slots__=('raw',)
    def __init__(self, value=None, rnd=DOWN):
        self.raw=_Raw();_init(C.byref(self.raw),BITS)
        if value is not None:
            if not isinstance(value,(str,int)):raise TypeError('Exact inputs only')
            rc=_set(C.byref(self.raw),str(value).encode(),10,rnd)
            if rc or _nan(C.byref(self.raw)):raise ValueError('Invalid MPFR input')
    def __del__(self):
        if hasattr(self,'raw'):_clear(C.byref(self.raw))
    def __lt__(self,o):return _cmp(C.byref(self.raw),C.byref(o.raw))<0
    def __le__(self,o):return _cmp(C.byref(self.raw),C.byref(o.raw))<=0
    def op(self,o,name,rnd):
        z=Real();_ops[name](C.byref(z.raw),C.byref(self.raw),C.byref(o.raw),rnd);return z
    def un(self,name,rnd):
        z=Real();_un[name](C.byref(z.raw),C.byref(self.raw),rnd);return z
    def text(self,rnd,digits=65):
        exp=C.c_long();ptr=_get(None,C.byref(exp),10,digits,C.byref(self.raw),rnd)
        raw=C.string_at(ptr).decode();_free(ptr)
        if '@' in raw:raise ArithmeticError('Nonfinite MPFR result: '+raw)
        sign=''
        if raw.startswith('-'):sign='-';raw=raw[1:]
        return sign+raw[0]+'.'+raw[1:]+'e'+str(exp.value-1)
    def display_float(self):return _getd(C.byref(self.raw),0)
class I:
    __slots__=('lo','hi')
    def __init__(self,value=0,hi=None):
        if isinstance(value,I):self.lo,self.hi=value.lo,value.hi;return
        if isinstance(value,Real):
            if not isinstance(hi,Real):raise TypeError
            self.lo,self.hi=value,hi
        elif isinstance(value,Fraction):
            t=I(value.numerator)/I(value.denominator);self.lo,self.hi=t.lo,t.hi
        elif isinstance(value,(int,str)):
            self.lo=Real(value,DOWN);self.hi=Real(value if hi is None else hi,UP)
        else:raise TypeError('Use exact integer, decimal string, Fraction, or Interval')
        if not self.lo<=self.hi:raise ValueError('Reversed interval')
    @staticmethod
    def pi():
        l,r=Real(),Real();_pi(C.byref(l.raw),DOWN);_pi(C.byref(r.raw),UP);return I(l,r)
    def __add__(self,o):
        o=I(o);return I(self.lo.op(o.lo,'add',DOWN),self.hi.op(o.hi,'add',UP))
    __radd__=__add__
    def __neg__(self):return I(self.hi.un('neg',DOWN),self.lo.un('neg',UP))
    def __sub__(self,o):return self+-I(o)
    def __rsub__(self,o):return I(o)+-self
    def __mul__(self,o):
        o=I(o)
        if I(0).lo<=self.lo and I(0).lo<=o.lo:
            return I(self.lo.op(o.lo,'mul',DOWN),self.hi.op(o.hi,'mul',UP))
        pairs=[(a,b) for a in (self.lo,self.hi) for b in (o.lo,o.hi)]
        return I(min(a.op(b,'mul',DOWN) for a,b in pairs),max(a.op(b,'mul',UP) for a,b in pairs))
    __rmul__=__mul__
    def inv(self):
        z=Real(0)
        if self.lo<=z and z<=self.hi:raise ZeroDivisionError
        one=Real(1)
        return I(one.op(self.hi,'div',DOWN),one.op(self.lo,'div',UP))
    def __truediv__(self,o):return self*I(o).inv()
    def __rtruediv__(self,o):return I(o)*self.inv()
    def square(self):
        z=Real(0)
        if z<=self.lo:return I(self.lo.op(self.lo,'mul',DOWN),self.hi.op(self.hi,'mul',UP))
        if self.hi<=z:return (-self).square()
        return I(z,max(self.lo.op(self.lo,'mul',UP),self.hi.op(self.hi,'mul',UP)))
    def __pow__(self,n):
        if not isinstance(n,int) or n<0:raise ValueError
        if n==0:return I(1)
        if n==1:return self
        t=self**(n//2);z=t.square()
        return z*self if n%2 else z
    def increasing(self,n):return I(self.lo.un(n,DOWN),self.hi.un(n,UP))
    def exp(self):return self.increasing('exp')
    def sqrt(self):return self.increasing('sqrt')
    def log(self):return self.increasing('log')
    def erf(self):return self.increasing('erf')
    def atan(self):return self.increasing('atan')
    def erfc(self):return I(self.hi.un('erfc',DOWN),self.lo.un('erfc',UP))
    def below(self,o):return self.hi<I(o).lo
    def above(self,o):return I(o).hi<self.lo
    def contains(self,o):
        o=I(o);return self.lo<=o.lo and o.hi<=self.hi
    def overlaps(self,o):
        o=I(o);return self.lo<=o.hi and o.lo<=self.hi
    def record(self):return {'lower':self.lo.text(DOWN),'upper':self.hi.text(UP)}
    def __repr__(self):return str(self.record())

def selftest():
    for a in [Fraction(-7,3),Fraction(0),Fraction(2,3),Fraction(17,5)]:
        for b in [Fraction(-11,7),Fraction(1,9),Fraction(19)]:
            assert (I(a)+I(b)).overlaps(I(a+b))
            assert (I(a)-I(b)).overlaps(I(a-b))
            assert (I(a)*I(b)).overlaps(I(a*b))
            assert (I(a)/I(b)).overlaps(I(a/b))
    assert I(0).exp().contains(1)
    assert I(0).erf().contains(0)
    assert I(4).sqrt().contains(2)
    assert I(2).sqrt().square().contains(2)
    assert I.pi().above('3.14159265358979323846264338327950288419716939937510')
    assert I.pi().below('3.14159265358979323846264338327950288419716939937511')
    assert (16*I(Fraction(1,5)).atan()-4*I(Fraction(1,239)).atan()).overlaps(I.pi())
    return {'MPFR':_version().decode(),'bits':BITS,'status':'PASS'}
if __name__=='__main__':print(selftest())
