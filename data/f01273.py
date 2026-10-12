"""Exact O(d) evaluation of folded Gaussian cubic matrix entries.

Arithmetic is delegated to the input scalar type: Fraction works exactly for
formal rational moment inputs; directed intervals can be used for Gaussian
cell probabilities and moments. No numerical integration is performed here.
"""


def projected_cubic_kernel(x, y):
    """The kernel of the orthonormal total-degree-three Hermite projector."""
    if len(x) != len(y):
        raise ValueError("The two vectors must have the same dimension")
    s = sum(a * b for a, b in zip(x, y))
    nx = sum(a * a for a in x)
    ny = sum(b * b for b in y)
    return (s ** 3 - 3 * s * (nx + ny) + 3 * (len(x) + 2) * s) / 6


def conditional_cubic_kernel(x_moments, y_moments):
    """E K3(X,Y) for independent product laws.

    Each argument is a sequence of coordinate triples (E X_i, E X_i^2,
    E X_i^3). Coordinates are independent within each law and X,Y are
    independent. The implementation uses O(d) scalar operations and O(1)
    additional space; no list of Hermite multi-indices is needed.
    """
    if len(x_moments) != len(y_moments):
        raise ValueError("The two moment arrays must have the same dimension")
    s = e_sum = eu_sum = u2_sum = u3_sum = t_sum = 0
    for (a1, a2, a3), (b1, b2, b3) in zip(x_moments, y_moments):
        u = a1 * b1
        e = (a2 - 1) * (b2 - 1)
        t = (a3 - 3 * a1) * (b3 - 3 * b1)
        s += u
        e_sum += e
        eu_sum += e * u
        u2_sum += u * u
        u3_sum += u * u * u
        t_sum += t
    return (t_sum + 3 * s * e_sum - 3 * eu_sum + s ** 3
            - 3 * s * u2_sum + 2 * u3_sum) / 6


def folded_cubic_entry(x_moments, y_moments, probability_x,
                       probability_y, cubic_coefficient):
    """A_CD = 4 gamma(C) gamma(D) E[K1(X,Y)-b K3(X,Y)].

    C,D are representatives of antipodal cell pairs. Conditional moments
    refer to the representative cells, before multiplying by their masses.
    """
    linear = sum(a[0] * b[0] for a, b in zip(x_moments, y_moments))
    cubic = conditional_cubic_kernel(x_moments, y_moments)
    return 4 * probability_x * probability_y * (
        linear - cubic_coefficient * cubic)
