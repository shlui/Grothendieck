"""Complete curved-boundary candidates for 4*kappa*zeta-kappa <= 54.

For the guarded current moment/source families, pricing is separately
convex and convex along zeta=1/4+27/(2*kappa). These exact rational boundary
points therefore suffice to maximize every such price on a moment box.
"""
from fractions import Fraction as F

def points(klo,khi,zlo=F(2,5),zhi=F(1)):
    klo,khi,zlo,zhi=map(F,(klo,khi,zlo,zhi))
    assert 0<=klo<=khi<=90 and F(2,5)<=zlo<=zhi<=1
    c,d=F(1,4),F(27,2)
    if klo*(zlo-c)>d:return []
    result={(k,z)for k in (klo,khi)for z in (zlo,zhi)if k*(z-c)<=d}
    for k in (klo,khi):
        if k:
            z=c+d/k
            if zlo<=z<=zhi:result.add((k,z))
    for z in (zlo,zhi):
        k=d/(z-c)
        if klo<=k<=khi:result.add((k,z))
    return sorted(result)

if __name__=='__main__':
    assert points(60,90)==[(F(60),F(2,5)),(F(60),F(19,40)),(F(90),F(2,5))]
    assert len(points(0,1))==4
    assert not points(80,90,F(1,2),1)
    assert points(90,90)==[(F(90),F(2,5))]
    print('PASS_STRONG_PHYSICAL_SHAPE_BOUNDARY_CANDIDATES')
