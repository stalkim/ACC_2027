"""Small, fixed checks of the geometric bound and informative histories.

Run once into a new directory. Inputs, witnesses, source hashes and timing
samples are retained. Only Python's standard library is required.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
import json
from math import factorial
from pathlib import Path
import platform
import shutil
from statistics import median
from time import perf_counter

from extension_history import decode, verify_witness
from monotone_history import CHARTS
import straight_history as straight

METHODS = ('DC', 'MI', 'DMI')
INFORMATIVE_EXAMPLE = dict(times=['-0.005', '-0.0025', '0.005'], positions=[
    ['0.000025', '-0.0049875624375', '0'],
    ['-0.01579375', '-0.00298751405859375', '0'],
    ['0.000025', '0.0050125625625', '0']], linf_errors=['0.012']*3)


def point(q):
    x = (q-1)*(q+1)
    return x, x*(q+q*q/5)/2, F(0)


def endpoints(raw):
    return {k: [v[0], v[-1]] for k, v in raw.items()}


def informative_cases():
    cases = [('informative-example', INFORMATIVE_EXAMPLE,
              tuple(F(1)+F(t) for t in INFORMATIVE_EXAMPLE['times']))]
    for middle, shift in product((F(-3, 1000), F(-1, 400), F(-1, 500)),
                                 (F(-1, 5000), F(0), F(1, 5000))):
        times = (F(-1, 200), middle, F(1, 200))
        positive, negative = tuple(1+t+shift for t in times), tuple(-1+t+shift for t in times)
        observations = [tuple((a+b)/2 for a, b in zip(point(p), point(n)))
                        for p, n in zip(positive, negative)]
        observations[1] = tuple(a+b for a, b in zip(point(positive[1]), (F(-27, 2500), 0, 0)))
        raw = dict(times=times, positions=observations, linf_errors=[F(3, 250)]*3)
        name = f'interior={middle},shift={shift}'
        for branch, q in (('positive', positive), ('negative', negative)):
            assert verify_witness((q[0], q[-1]), endpoints(raw), 'quartic', F(1, 2), CHARTS[branch])
        assert verify_witness(positive, raw, 'quartic', F(1, 2), CHARTS['positive'])
        # Analytic contradiction: the middle x upper bound forces q_1>-1.
        assert observations[1][0]+F(3, 250) < 0
        q2_lower = -1+F(9, 10)*(times[2]-times[1])
        assert q2_lower*q2_lower-1 < observations[2][0]-F(3, 250)
        cases.append((name, raw, positive))
    return cases


def densified(count):
    # Keep all original observations; add noiseless samples of q=1+t.
    rows = {F(t): tuple(map(F, p))
            for t, p in zip(INFORMATIVE_EXAMPLE['times'], INFORMATIVE_EXAMPLE['positions'])}
    needed = count-len(rows)
    for j in range(1, needed+2):
        if len(rows) == count:
            break
        t = F(-1, 200)+F(1, 100)*F(j, needed+2)
        if t not in rows:
            rows[t] = point(1+t)
    assert len(rows) == count
    times = sorted(rows)
    return dict(times=times, positions=[rows[t] for t in times], linf_errors=[F(3, 250)]*count)


def trig_bounds(x, sine):
    """Adjacent alternating partial sums; valid here for 0<=x<=1."""
    assert 0 <= x <= 1
    terms = [(-1)**j*x**(2*j+int(sine))/factorial(2*j+int(sine)) for j in range(13)]
    left, right = sum(terms[:12]), sum(terms)
    return min(left, right), max(left, right)


def curvature_rows():
    h, k, b, reach = F(3, 5), F(4, 5), F(1, 100), F(3, 10)
    rows = []
    for curvature, rho in product((F(1, 2), F(1), F(2)), (F(0), F(1, 20), F(1, 10))):
        cos_l, _ = trig_bounds(curvature*reach, False)
        _, sin_u = trig_bounds(curvature*reach, True)
        assert h*cos_l-k*sin_u > 0  # x' positive on the entire chart
        assert k*cos_l-h*sin_u > 0  # x' strictly increasing
        for j in range(1, 51):
            span = F(j, 200)
            covered = rho+F(6, 5)*span < reach
            row = dict(curvature=curvature, rho=rho, span=span, covered=covered)
            if covered:
                s0, s1 = trig_bounds(curvature*rho, True), trig_bounds(curvature*(rho+span), True)
                c0, c1 = trig_bounds(curvature*rho, False), trig_bounds(curvature*(rho+span), False)
                low = (k*(c1[0]-c0[1])+h*(s1[0]-s0[1]))/curvature
                high = (k*(c1[1]-c0[0])+h*(s1[1]-s0[0]))/curvature
                general = h*span-curvature*(F(6, 5)*span*rho+F(18, 25)*span**2)
                assert general <= low
                exact_pass, general_pass = low > 2*b, general > 2*b
                assert not general_pass or exact_pass
                row.update(chord_outer=(low, high), general_chord_lower=general,
                           exact_projected_pass=exact_pass, general_pass=general_pass,
                           threshold_unresolved=low <= 2*b < high)
            rows.append(row)
    return rows


def encode(value):
    if isinstance(value, F):
        return str(value)
    if isinstance(value, dict):
        return {str(k): encode(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [encode(v) for v in value]
    return value


def write_json(path, value):
    path.write_text(json.dumps(encode(value), indent=2, sort_keys=True)+'\n')


def straight_cases():
    for (h, k), alpha, count in product(
            ((F(5, 13), F(12, 13)), (F(3, 5), F(4, 5)), (F(4, 5), F(3, 5)), (F(1), F(0))),
            (F(4, 5), F(99, 100), F(1), F(101, 100), F(6, 5)), (2, 3, 9, 33)):
        radius, span = F(1, 100), F(1, 50)*alpha/h
        times = tuple(-span/2+span*F(m, count-1) for m in range(count))
        directions = ((h, k, F(0)), (-h, k, F(0)))
        base = dict(h=h, k=k, alpha=alpha, count=count, times=times,
                    radius=radius, directions=directions)
        yield dict(base, kind='midpoint', positions=[(F(0), k*t, F(0)) for t in times]), None
        for branch, variable, anchor in product((0, 1), (False, True), (F(-1, 3), F(-2, 3))):
            q = [anchor*span]
            for m in range(1, count):
                rate = F(6, 5) if variable and m % 2 == 0 else F(1)
                q.append(q[-1]+rate*(times[m]-times[m-1]))
            observations = [tuple(u*s+(F(9, 10)*radius*(-1)**m if axis == 0 else 0)
                                  for axis, u in enumerate(directions[branch])) for m, s in enumerate(q)]
            yield dict(base, kind='valid', positions=observations), dict(branch=branch, q=q)
        yield dict(base, kind='off_plane', positions=[(F(0), F(0), 2*radius)]*count), None


def run(output, protocol=None):
    output.mkdir(parents=True, exist_ok=False)
    source = Path(__file__).resolve().parent
    snapshot = output/'source'
    snapshot.mkdir()
    for name in ('geometric_checks.py', 'straight_history.py', 'extension_history.py',
                 'monotone_history.py', 'test_geometric_checks.py'):
        shutil.copy2(source/name, snapshot/name)
    if protocol is not None:
        shutil.copy2(protocol, snapshot/'protocol.md')
    write_json(output/'manifest.json', dict(python=platform.python_version(), platform=platform.platform(),
               source_sha256={p.name: sha256(p.read_bytes()).hexdigest() for p in snapshot.iterdir()}))
    # Write all generated inputs and truth before solver execution.
    informative = informative_cases()
    lines = list(straight_cases())
    dense = [(n, bits, densified(n)) for n, bits in product((3, 9, 33), (18, 24, 32))]
    write_json(output/'inputs.json', dict(straight=[r for r, _ in lines],
               informative=[dict(id=i, raw=r) for i, r, _ in informative],
               dense=[dict(count=n, bits=b, raw=r) for n, b, r in dense]))
    write_json(output/'truth.json', dict(straight=[t for _, t in lines],
               informative=[dict(id=i, q=q) for i, _, q in informative]))
    log = (output/'results.jsonl').open('w')
    def record(row):
        log.write(json.dumps(encode(row), sort_keys=True)+'\n')
        log.flush()
    start = perf_counter()
    outcomes = Counter()
    for index, (raw, truth) in enumerate(lines):
        result = straight.decode(raw['times'], raw['positions'], raw['radius'], raw['directions'])
        record(dict(group='straight', id=index, result=result))
        outcomes[raw['kind']+':'+result['status']] += 1
        if raw['kind'] == 'midpoint':
            assert result['status'] == ('AMBIGUOUS' if raw['alpha'] <= 1 else 'INCOMPATIBLE')
            if raw['alpha'] <= 1:
                for tau in raw['directions']:
                    assert straight.verify(raw['times'], raw['times'], raw['positions'], raw['radius'], tau, (F(1), F(6, 5)))
        elif raw['kind'] == 'off_plane':
            assert result['status'] == 'INCOMPATIBLE'
        else:
            branch, q = truth['branch'], truth['q']
            assert straight.verify(q, raw['times'], raw['positions'], raw['radius'], raw['directions'][branch], (F(1), F(6, 5)))
            outer = result['branches'][branch]['outer']
            assert outer is not None and outer[0] <= q[-1] <= outer[1]
            if raw['alpha'] > 1:
                assert (result['status'], result['selected']) == ('UNIQUE', branch)
    print('Straight checks passed:', dict(outcomes), flush=True)

    informative_rows = []
    for name, raw, truth in informative:
        for mode, data in (('endpoints', endpoints(raw)), ('full', raw)):
            results = []
            for method in METHODS:
                result = decode(data, 'quartic', F(1, 2), method)
                row = dict(group='informative', id=name, mode=mode, method=method, result=result)
                record(row)
                informative_rows.append(row)
                results.append(result)
                assert result['status'] == ('AMBIGUOUS' if mode == 'endpoints' else 'UNIQUE')
                assert result['pruned'] == 0
                if mode == 'full':
                    assert result['selected'] == 'positive'
                    out = result['branches']['positive']['outer']
                    assert out[0] <= truth[-1] <= out[1]
            assert results[0]['branches']['positive']['outer'] == results[1]['branches']['positive']['outer']
    print('Informative checks passed:', len(informative), 'histories', flush=True)

    dense_rows = []
    for count, bits, raw in dense:
        results = []
        for method in METHODS:
            result = decode(raw, 'quartic', F(1, 2), method, bits)
            results.append(result)
            record(dict(group='scaling', count=count, bits=bits, method=method, result=result))
            assert (result['status'], result['selected'], result['pruned']) == ('UNIQUE', 'positive', 0)
            dense_rows.append(dict(count=count, bits=bits, method=method, evaluations=result['evaluations']))
        assert results[0]['branches']['positive']['outer'] == results[1]['branches']['positive']['outer']
        assert results[1]['evaluations'] == results[2]['evaluations']

    curved = curvature_rows()
    write_json(output/'curvature.json', curved)
    curve_summary = dict(covered=sum(r['covered'] for r in curved),
                        general_pass=sum(r.get('general_pass', False) for r in curved),
                        exact_projected_pass=sum(r.get('exact_projected_pass', False) for r in curved),
                        unresolved=sum(r.get('threshold_unresolved', False) for r in curved))
    print('Curvature checks passed:', curve_summary, flush=True)
    timings, batches = [], {m: [] for m in ('MI', 'DMI')}
    for repetition in range(3):
        for method in (('MI', 'DMI') if repetition % 2 == 0 else ('DMI', 'MI')):
            durations = []
            for name, raw, _ in informative[1:]:
                tic = perf_counter()
                result = decode(raw, 'quartic', F(1, 2), method)
                elapsed = perf_counter()-tic
                assert (result['status'], result['selected'], result['pruned']) == ('UNIQUE', 'positive', 0)
                durations.append(elapsed)
                timings.append(dict(repetition=repetition, method=method, id=name, seconds=elapsed))
            batches[method].append(sum(durations))
    log.close()
    write_json(output/'timings.json', timings)
    summary = dict(straight_outcomes=dict(outcomes), straight_calls=len(lines),
                   informative_histories=len(informative), informative_calls=len(informative_rows),
                   scaling=dense_rows, curvature=curve_summary, timing_batches=batches,
                   timing_medians={m: median(v) for m, v in batches.items()},
                   elapsed_seconds=perf_counter()-start, all_checks_passed=True)
    write_json(output/'summary.json', summary)
    print(json.dumps(encode(summary), sort_keys=True), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--protocol', type=Path, help='optional frozen experiment plan to archive')
    args = parser.parse_args()
    run(args.output, args.protocol)
