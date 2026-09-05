"""Line/Euclidean-ball histories with rational inner and outer certificates."""
from fractions import Fraction as F
from math import isqrt


def sqrt_bracket(value, bits=40):
    value = F(value)
    if value < 0 or type(bits) is not int or not 0 <= bits <= 128:
        raise ValueError('nonnegative radicand and bounded precision required')
    n, d = isqrt(value.numerator), isqrt(value.denominator)
    if n*n == value.numerator and d*d == value.denominator:
        return F(n, d), F(n, d)
    scale = 2**bits
    low = isqrt(value.numerator*scale*scale // value.denominator)
    return F(low, scale), F(low+1, scale)


def intersect(left, right):
    if left is None or right is None:
        return None
    lo, hi = max(left[0], right[0]), min(left[1], right[1])
    return (lo, hi) if lo <= hi else None


def propagate(intervals, times, rates):
    reachable = []
    a, c = rates
    for m, interval in enumerate(intervals):
        if m:
            dt = times[m]-times[m-1]
            interval = intersect(interval, (reachable[-1][0]+a*dt,
                                             reachable[-1][1]+c*dt))
        if interval is None:
            return None
        reachable.append(interval)
    return reachable


def witness(reachable, times, rates):
    values = [sum(reachable[-1])/2]
    a, c = rates
    for m in range(len(times)-2, -1, -1):
        dt = times[m+1]-times[m]
        allowed = intersect(reachable[m], (values[-1]-c*dt, values[-1]-a*dt))
        assert allowed is not None
        values.append(sum(allowed)/2)
    return tuple(reversed(values))


def verify(values, times, positions, radius, direction, rates):
    a, c = rates
    return (len(values) == len(times)
            and all(sum((u*s-y)**2 for u, y in zip(direction, point)) <= radius**2
                    for s, point in zip(values, positions))
            and all(a*(v-u) <= r-q <= c*(v-u)
                    for u, v, q, r in zip(times, times[1:], values, values[1:])))


def decode(times, positions, radius, directions, rates=(F(1), F(6, 5)), bits=40):
    """No directional pruning or true state is used by this solver.

    For a unit line direction u, a ball's scalar preimage is centered at
    u.y and has radius sqrt(b^2-|y|^2+(u.y)^2). Irrational roots are enclosed.
    """
    times = tuple(map(F, times))
    positions = tuple(tuple(map(F, p)) for p in positions)
    directions = tuple(tuple(map(F, u)) for u in directions)
    radius, rates = F(radius), tuple(map(F, rates))
    if (not times or len(times) != len(positions) or radius < 0
            or len(rates) != 2 or not 0 < rates[0] <= rates[1]
            or any(v <= u for u, v in zip(times, times[1:]))
            or not directions or any(sum(u*u for u in tau) != 1 for tau in directions)
            or any(len(p) != len(tau) for p in positions for tau in directions)):
        raise ValueError('valid times, unit directions, balls and positive rates required')
    sqrt_bracket(F(0), bits)
    branches = []
    for tau in directions:
        outers, inners = [], []
        for point in positions:
            center = sum(u*y for u, y in zip(tau, point))
            radicand = radius**2-sum(y*y for y in point)+center**2
            if radicand < 0:
                outers.append(None)
                inners.append(None)
            else:
                low, high = sqrt_bracket(radicand, bits)
                outers.append((center-high, center+high))
                inners.append((center-low, center+low))
        out, inside = propagate(outers, times, rates), propagate(inners, times, rates)
        q = witness(inside, times, rates) if inside is not None else None
        if q is not None:
            assert verify(q, times, positions, radius, tau, rates)
        branches.append(dict(outer=out[-1] if out is not None else None,
                             witness=q, sample_outers=outers))
    retained = [i for i, b in enumerate(branches) if b['outer'] is not None]
    if not retained:
        status, selected = 'INCOMPATIBLE', None
    elif len(retained) == 1 and branches[retained[0]]['witness'] is not None:
        status, selected = 'UNIQUE', retained[0]
    elif len(retained) > 1 and all(branches[i]['witness'] is not None for i in retained):
        status, selected = 'AMBIGUOUS', None
    else:
        status, selected = 'UNRESOLVED', None
    return dict(status=status, selected=selected, branches=branches)
