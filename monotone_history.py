"""Branch recovery on a cubic path using rational interval bounds.

Fixed inverse precision can leave boundary cases unresolved. The decoder
requires the chart coverage and positive parameter-rate bounds checked below.
"""
from fractions import Fraction as F

CHARTS = {'negative': (F(-5, 4), F(-3, 4)),
          'positive': (F(3, 4), F(5, 4))}
RATES = (F(9, 10), F(11, 10))


def intersect(a, b):
    if a is None or b is None:
        return None
    low, high = max(a[0], b[0]), min(a[1], b[1])
    return (low, high) if low <= high else None


def inverse_bracket(f, value, domain, steps):
    """f must be strictly increasing; caller has checked the image range."""
    low, high = domain
    if f(low) == value:
        return low, low
    if f(high) == value:
        return high, high
    for _ in range(steps):
        mid = (low + high) / 2
        y = f(mid)
        if y == value:
            return mid, mid
        if y < value:
            low = mid
        else:
            high = mid
    return low, high


def preimage(f, values, domain, steps):
    """Return (outer, inner) for a closed interval under increasing f."""
    low, high = domain
    a, b = values
    f_low, f_high = f(low), f(high)
    if b < f_low or a > f_high:
        return None, None
    left = (low, low) if a <= f_low else inverse_bracket(f, a, domain, steps)
    right = (high, high) if b >= f_high else inverse_bracket(f, b, domain, steps)
    outer = intersect((left[0], right[1]), domain)
    inner = intersect((left[1], right[0]), domain)
    return outer, inner


def forward(intervals, times):
    """Exact reachable prefix intervals for a scalar two-sided rate band."""
    result = []
    for k, interval in enumerate(intervals):
        if k:
            dt = times[k] - times[k-1]
            previous = result[-1]
            interval = intersect(interval, (previous[0] + RATES[0]*dt,
                                            previous[1] + RATES[1]*dt))
        if interval is None:
            return None
        result.append(interval)
    return tuple(result)


def backward(reachable, times):
    q = [(reachable[-1][0] + reachable[-1][1]) / 2]
    for k in range(len(times)-2, -1, -1):
        dt = times[k+1] - times[k]
        allowed = intersect(reachable[k], (q[-1] - RATES[1]*dt,
                                           q[-1] - RATES[0]*dt))
        assert allowed is not None
        q.append((allowed[0] + allowed[1]) / 2)
    return tuple(reversed(q))


def curve(q, lam):
    return q*q-1, lam*q*(q*q-1), F(0)


def check_witness(q, history, lam, domain):
    t, y, e = history
    return (len(q) == len(t)
            and all(domain[0] <= z <= domain[1] for z in q)
            and all(abs(a-b) <= error for z, point, error in zip(q, y, e)
                    for a, b in zip(curve(z, lam), point))
            and all(RATES[0]*(t[k]-t[k-1]) <= q[k]-q[k-1]
                    <= RATES[1]*(t[k]-t[k-1]) for k in range(1, len(t))))


def decode(raw, lam, *, directional=True, steps=24):
    lam = F(lam)
    if lam <= 0 or type(steps) is not int or not 0 <= steps <= 64:
        raise ValueError('positive scale and bounded integer precision required')
    t = tuple(F(v) for v in raw['times'])
    y = tuple(tuple(F(v) for v in row) for row in raw['positions'])
    e = tuple(F(v) for v in raw['linf_errors'])
    if not t or len(t) != len(y) or len(t) != len(e):
        raise ValueError('nonempty equal-length history required')
    if any(len(row) != 3 for row in y) or any(v < 0 for v in e):
        raise ValueError('three coordinates and nonnegative error required')
    if any(b <= a for a, b in zip(t, t[1:])):
        raise ValueError('strict timestamps required')
    duration = t[-1] - t[0]
    covered = (all(e[k] <= F(1, 16) and abs(y[k][0]) <= F(3, 8)
                   for k in (0, len(t)-1)) and RATES[1]*duration < F(3, 2))
    if not covered:
        return {'status': 'OUTSIDE_COVERAGE', 'selected': None, 'branches': {},
                'evaluations': 0, 'pruned': 0}
    calls, pruned, branches = 0, 0, {}
    for name, domain in CHARTS.items():
        if directional and len(t) > 1:
            m = (F(-11, 4), F(-27, 20)) if name == 'negative' else (F(27, 20), F(11, 4))
            z, error = y[-1][0]-y[0][0], e[0]+e[-1]
            if not duration*m[0]-error <= z <= duration*m[1]+error:
                branches[name] = {'state': 'EMPTY', 'outer': None, 'witness': None}
                pruned += 1
                continue
        outers, inners = [], []
        for point, error in zip(y, e):
            outer = inner = domain
            if abs(point[2]) > error:
                outer = inner = None
            for axis in (0, 1):
                if outer is None:
                    break
                sign = -1 if axis == 0 and name == 'negative' else 1

                def f(q):
                    nonlocal calls
                    calls += 1
                    return sign*curve(q, lam)[axis]

                values = sorted((sign*(point[axis]-error), sign*(point[axis]+error)))
                o, i = preimage(f, values, domain, steps)
                outer, inner = intersect(outer, o), intersect(inner, i)
            outers.append(outer)
            inners.append(inner)
        reachable = forward(outers, t)
        if reachable is None:
            branches[name] = {'state': 'EMPTY', 'outer': None, 'witness': None}
            continue
        interior = forward(inners, t)
        witness = backward(interior, t) if interior is not None else None
        if witness is not None:
            assert check_witness(witness, (t, y, e), lam, domain)
        branches[name] = {'state': 'WITNESSED' if witness is not None else 'POSSIBLE',
                          'outer': reachable[-1], 'witness': witness}
    retained = [name for name, b in branches.items() if b['state'] != 'EMPTY']
    if not retained:
        status, selected = 'INCOMPATIBLE', None
    elif len(retained) == 1 and branches[retained[0]]['state'] == 'WITNESSED':
        status, selected = 'UNIQUE', retained[0]
    elif all(branches[n]['state'] == 'WITNESSED' for n in retained):
        status, selected = 'AMBIGUOUS', None
    else:
        status, selected = 'UNRESOLVED', None
    return {'status': status, 'selected': selected, 'branches': branches,
            'evaluations': calls, 'pruned': pruned}
