"""Branch recovery on quartic and quintic paths.

The three methods share the same inverse intervals. DC checks all pairwise
time constraints; MI propagates intervals; DMI first excludes by direction.
"""
from fractions import Fraction as F
from monotone_history import CHARTS, RATES, intersect, preimage, forward, backward


def polynomial(q, family, scale):
    # Expanded expressions, independent of the factored witness checker.
    if family == 'quartic':
        second = q**4/F(5) + q**3 - q*q/F(5) - q
    elif family == 'quintic':
        second = q**5-q
    else:
        raise ValueError('undeclared family')
    return q*q-1, scale*second, F(0)


def verify_witness(q, raw, family, scale, domain):
    """Check original factored equations; no calls to polynomial/preimage."""
    t = tuple(map(F, raw['times']))
    if len(q) != len(t):
        return False
    for value, point, error in zip(q, raw['positions'], raw['linf_errors']):
        if not domain[0] <= value <= domain[1]:
            return False
        x = (value-1)*(value+1)
        factor = value+value*value/F(5) if family == 'quartic' else value*(value*value+1)
        actual = (x, F(scale)*x*factor, F(0))
        if any(abs(a-F(b)) > F(error) for a, b in zip(actual, point)):
            return False
    return all(RATES[0]*(t[k]-t[k-1]) <= q[k]-q[k-1]
               <= RATES[1]*(t[k]-t[k-1]) for k in range(1, len(t)))


def pairwise(intervals, times):
    """Exact difference-constraint feasibility/projection, without forward().

    Return final projection and a greatest feasible vector. All pair bounds
    are checked, including nonadjacent pairs. This is the closed form of the
    shortest-path constraints for a scalar chain with a constant rate band.
    """
    if any(v is None for v in intervals):
        return None
    a, c = RATES
    for j in range(len(times)):
        for i in range(j+1):
            dt = times[j]-times[i]
            if (intervals[i][0]+a*dt > intervals[j][1]
                    or intervals[j][0] > intervals[i][1]+c*dt):
                return None
    q = tuple(min(intervals[j][1] + (c if t >= times[j] else a)*(t-times[j])
                  for j in range(len(times))) for t in times)
    final = (max(v[0]+a*(times[-1]-t) for v, t in zip(intervals, times)),
             min(v[1]+c*(times[-1]-t) for v, t in zip(intervals, times)))
    return final, q


def decode(raw, family, scale, method='DMI', steps=24):
    scale = F(scale)
    if family not in ('quartic', 'quintic') or scale <= 0:
        raise ValueError('declared family and positive scale required')
    if method not in ('DC', 'MI', 'DMI') or type(steps) is not int or not 0 <= steps <= 64:
        raise ValueError('declared method and bounded integer precision required')
    t = tuple(map(F, raw['times']))
    y = tuple(tuple(map(F, row)) for row in raw['positions'])
    errors = tuple(map(F, raw['linf_errors']))
    if (not t or len(t) != len(y) or len(t) != len(errors)
            or any(len(row) != 3 for row in y) or any(e < 0 for e in errors)
            or any(u >= v for u, v in zip(t, t[1:]))):
        raise ValueError('invalid finite history')
    duration = t[-1]-t[0]
    if (any(errors[k] > F(1, 16) or abs(y[k][0]) > F(3, 8)
            for k in (0, len(t)-1)) or RATES[1]*duration >= F(3, 2)):
        return dict(status='OUTSIDE_COVERAGE', selected=None, branches={}, evaluations=0, pruned=0)
    calls, pruned, branches = 0, 0, {}
    for name, domain in CHARTS.items():
        if method == 'DMI' and len(t) > 1:
            low, high = ((F(-11, 4), F(-27, 20)) if name == 'negative'
                         else (F(27, 20), F(11, 4)))
            error = errors[0]+errors[-1]
            if not duration*low-error <= y[-1][0]-y[0][0] <= duration*high+error:
                branches[name] = dict(state='EMPTY', outer=None, witness=None)
                pruned += 1
                continue
        outers, inners = [], []
        for point, error in zip(y, errors):
            outer = inner = domain if abs(point[2]) <= error else None
            for axis in (0, 1):
                if outer is None:
                    break
                sign = -1 if axis == 0 and name == 'negative' else 1

                def evaluate(q):
                    nonlocal calls
                    calls += 1
                    return sign*polynomial(q, family, scale)[axis]

                values = sorted((sign*(point[axis]-error), sign*(point[axis]+error)))
                o, i = preimage(evaluate, values, domain, steps)
                outer, inner = intersect(outer, o), intersect(inner, i)
            outers.append(outer)
            inners.append(inner)
        if method == 'DC':
            out, inside = pairwise(outers, t), pairwise(inners, t)
            final = out[0] if out is not None else None
            witness = inside[1] if inside is not None else None
        else:
            out, inside = forward(outers, t), forward(inners, t)
            final = out[-1] if out is not None else None
            witness = backward(inside, t) if inside is not None else None
        if witness is not None:
            assert verify_witness(witness, raw, family, scale, domain)
        branches[name] = dict(state=('EMPTY' if final is None else
                                    'WITNESSED' if witness is not None else 'POSSIBLE'),
                              outer=final, witness=witness)
    retained = [n for n, b in branches.items() if b['state'] != 'EMPTY']
    if not retained:
        status, selected = 'INCOMPATIBLE', None
    elif len(retained) == 1 and branches[retained[0]]['state'] == 'WITNESSED':
        status, selected = 'UNIQUE', retained[0]
    elif all(branches[n]['state'] == 'WITNESSED' for n in retained):
        status, selected = 'AMBIGUOUS', None
    else:
        status, selected = 'UNRESOLVED', None
    return dict(status=status, selected=selected, branches=branches, evaluations=calls, pruned=pruned)
