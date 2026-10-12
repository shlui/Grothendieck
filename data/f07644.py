#!/usr/bin/env python3
"""Exact physical source polygon clipping and quadratic upper prices.

All inputs are rational. Each plane a*u+r<=cap must already have a proved
uniform cap over the whole cumulant/shape cell. Intersection preserves all
physical points. For a positive u-square coefficient, the maximum of a
separable quadratic on the polygon lies on its boundary; each edge is an
exact one-variable quadratic, with its interior maximum checked explicitly.
"""
from fractions import Fraction as F


def clip_source_polygon(ulo, uhi, rlo, rhi, planes=()):
    ulo, uhi, rlo, rhi = map(F, (ulo, uhi, rlo, rhi))
    assert ulo <= uhi and rlo <= rhi
    polygon = [(ulo, rlo), (uhi, rlo), (uhi, rhi), (ulo, rhi)]
    for a, cap in planes:
        a, cap = F(a), F(cap)
        assert a > 0
        if not polygon:
            break
        clipped = []
        previous = polygon[-1]
        fp = a * previous[0] + previous[1] - cap
        for current in polygon:
            fc = a * current[0] + current[1] - cap
            if (fp <= 0) != (fc <= 0):
                t = fp / (fp - fc)
                assert 0 <= t <= 1
                clipped.append(tuple(x + t * (y - x)
                                     for x, y in zip(previous, current)))
            if fc <= 0:
                clipped.append(current)
            previous, fp = current, fc
        polygon = []
        for point in clipped:
            if not polygon or point != polygon[-1]:
                polygon.append(point)
        if len(polygon) > 1 and polygon[0] == polygon[-1]:
            polygon.pop()
    return polygon


def max_source_quadratic(polygon, u_linear, r_square, r_linear,
                         constant=F(0), u_square=F(1)):
    """Return exact maximum and an attaining boundary point, or None.

    The price is u_square*u^2+u_linear*u+r_square*r^2+r_linear*r+constant.
    The strictly positive u_square proves that checking the boundary is
    sufficient even when r_square is negative.
    """
    u_linear, r_square, r_linear, constant, u_square = map(
        F, (u_linear, r_square, r_linear, constant, u_square))
    assert u_square > 0
    if not polygon:
        return None
    polygon = [tuple(map(F, point)) for point in polygon]

    def price(point):
        u, r = point
        return (u_square*u*u + u_linear*u + r_square*r*r
                + r_linear*r + constant)

    candidates = [(price(point), point) for point in polygon]
    for start, end in zip(polygon, polygon[1:] + polygon[:1]):
        du, dr = (end[0] - start[0], end[1] - start[1])
        square = u_square*du*du + r_square*dr*dr
        linear = ((2*u_square*start[0] + u_linear)*du
                  + (2*r_square*start[1] + r_linear)*dr)
        if square < 0:
            t = -linear / (2*square)
            if 0 < t < 1:
                point = (start[0] + t*du, start[1] + t*dr)
                candidates.append((price(point), point))
    return max(candidates)


if __name__ == '__main__':
    # The diagonal constraint r+u<=1 and price u^2-2r^2+3r give 5/4 at
    # u=r=1/2, exceeding every polygon vertex. This exercises the essential
    # edge-interior maximum rather than merely checking rectangle corners.
    polygon = clip_source_polygon(0, 1, 0, 1, [(1, 1)])
    result = max_source_quadratic(polygon, 0, -2, 3)
    assert result == (F(5, 4), (F(1, 2), F(1, 2)))
    assert clip_source_polygon(1, 2, 1, 2, [(1, 1)]) == []
    print('PASS_EXACT_POLYGON_EDGE_INTERIOR_MAXIMUM_AND_EMPTY_DOMAIN')
