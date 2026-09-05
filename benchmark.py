"""Compare three branch-recovery methods on 60 deterministic histories."""
import argparse
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import platform
import signal
from statistics import median
import sys
from time import perf_counter

from extension_history import decode, verify_witness
from monotone_history import CHARTS


def factored(q, family, scale):
    x = (q-1)*(q+1)
    factor = q+q*q/F(5) if family == 'quartic' else q*(q*q+1)
    return x, scale*x*factor, F(0)


def prepare():
    cases, truth = [], {}
    for family, scale in product(('quartic', 'quintic'), (F(1, 2), F(1), F(2))):
        prefix = family+'_'+str(scale).replace('/', '_')
        for anchor, span, error in product((F(-1), F(1)), (F(1, 50), F(2, 25)), (F(1, 2000), F(1, 100))):
            times = tuple(span*t for t in (F(-1), F(-3, 4), F(-1, 2), F(-1, 5), F(0)))
            rates = (F(9, 10), F(11, 10), F(19, 20), F(21, 20))
            q = [anchor]
            for k in range(3, -1, -1):
                q.append(q[-1]-rates[k]*(times[k+1]-times[k]))
            q = tuple(reversed(q))
            patterns = ((8, -6, 2, -9, 7), (-7, 3, 9, -4, -8), (1, -2, 3, -4, 5))
            positions = [tuple(v+error*patterns[axis][k]/10 for axis, v in enumerate(factored(z, family, scale)))
                         for k, z in enumerate(q)]
            key = prefix+'_nominal_'+str(len(cases))
            cases.append(dict(case=key, family=family, scale=scale, kind='nominal',
                              history=dict(times=times, positions=positions, linf_errors=[error]*5)))
            truth[key] = dict(branch='negative' if anchor < 0 else 'positive', histories=[q], expected='UNIQUE')
        for error, kind in ((F(3, 50), 'ambiguous'), (F(1, 2000), 'incompatible')):
            key = prefix+'_'+kind
            cases.append(dict(case=key, family=family, scale=scale, kind=kind,
                              history=dict(times=[F(-1, 200), F(1, 200)],
                                           positions=[(0, 0, 0)]*2, linf_errors=[error]*2)))
            histories = [(a-F(1, 200), a+F(1, 200)) for a in (F(-1), F(1))] if kind == 'ambiguous' else []
            truth[key] = dict(branch=None, histories=histories, expected=kind.upper())
    assert len(cases) == 60 and Counter(c['kind'] for c in cases) == dict(nominal=48, ambiguous=6, incompatible=6)
    return cases, truth


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def write_json(path, value):
    path.write_text(json.dumps(value, default=str, indent=2, sort_keys=True)+'\n')


def timeout(*_):
    raise TimeoutError('fixed five-second call limit')


def run(output, *, check_only=False):
    if not __debug__:
        raise RuntimeError('Run without -O; the benchmark checks its results with assertions.')
    if not hasattr(signal, 'setitimer'):
        raise RuntimeError('The benchmark requires Linux, macOS or WSL for per-call timeouts.')
    output.mkdir(parents=True, exist_ok=False)
    directory = Path(__file__).resolve().parent
    paths = [directory/name for name in (
        'benchmark.py', 'extension_history.py', 'monotone_history.py',
        'test_extension_history.py', 'test_monotone_history.py',
        'test_benchmark.py', 'check_precision.py', 'demo.py', 'README.md')]
    pins = {p.name: digest(p) for p in paths}
    cases, truth = prepare()
    write_json(output/'inputs.json', cases)
    write_json(output/'truth.json', truth)
    manifest = dict(experiment='polynomial-branch-comparison', check_only=check_only,
                    pins=pins, python=sys.version,
                    platform=platform.platform(), precision=24, threads=1,
                    inputs_sha256=digest(output/'inputs.json'), truth_sha256=digest(output/'truth.json'))
    write_json(output/'manifest.json', manifest)
    # Check declared truth in the independent original-equation scorer.
    for case in cases:
        for q in truth[case['case']]['histories']:
            domain = CHARTS['negative' if q[-1] < 0 else 'positive']
            assert verify_witness(q, case['history'], case['family'], case['scale'], domain)
    rows, baseline, totals = [], {}, {m:[] for m in ('DC', 'MI', 'DMI')}
    start = perf_counter()
    signal.signal(signal.SIGALRM, timeout)
    for repetition in range(1 if check_only else 6):
        batch = {m:0.0 for m in totals}
        for index, case in enumerate(cases):
            order = ('DC', 'MI', 'DMI')
            shift = (index+repetition)%3
            order = order[shift:]+order[:shift]
            outputs = {}
            for method in order:
                if perf_counter()-start >= 180:
                    raise TimeoutError('fixed 180-second suite limit')
                signal.setitimer(signal.ITIMER_REAL, 5)
                tick = perf_counter()
                try:
                    result = decode(case['history'], case['family'], case['scale'], method)
                    elapsed = perf_counter()-tick
                finally:
                    signal.setitimer(signal.ITIMER_REAL, 0)
                for name, branch in result['branches'].items():
                    if branch['witness'] is not None:
                        assert verify_witness(branch['witness'], case['history'], case['family'], case['scale'], CHARTS[name])
                        assert branch['outer'][0] <= branch['witness'][-1] <= branch['outer'][1]
                tr = truth[case['case']]
                assert result['status'] == tr['expected'], (case['case'], method, result)
                assert result['selected'] == tr['branch']
                for q in tr['histories']:
                    name = 'negative' if q[-1] < 0 else 'positive'
                    outer = result['branches'][name]['outer']
                    assert outer is not None and outer[0] <= q[-1] <= outer[1]
                outputs[method] = result
                batch[method] += elapsed
                row = dict(case=case['case'], kind=case['kind'], family=case['family'], method=method,
                           repetition=repetition, seconds=elapsed, result=result)
                rows.append(row)
                with (output/'calls.jsonl').open('a') as handle:
                    handle.write(json.dumps(row, default=str, sort_keys=True)+'\n')
                if repetition == 0:
                    baseline[(case['case'], method)] = result
                else:
                    assert baseline[(case['case'], method)] == result
            for name in CHARTS:
                assert outputs['MI']['branches'][name]['outer'] == outputs['DC']['branches'][name]['outer']
            assert len({(r['status'], r['selected']) for r in outputs.values()}) == 1
        if repetition > 0:
            for method, seconds in batch.items():
                totals[method].append(seconds)
    for p in paths:
        assert pins[p.name] == digest(p)
    summary = dict(status='COMPLETE', check_only=check_only, seconds=perf_counter()-start, cases=60, calls=len(rows), violations=[],
                   source_pins_unchanged=True, methods={})
    for method, seconds in totals.items():
        chosen = [r for r in rows if r['method'] == method and r['repetition'] == 0]
        summary['methods'][method] = dict(statuses=dict(Counter(r['result']['status'] for r in chosen)),
                                         evaluations=sum(r['result']['evaluations'] for r in chosen),
                                         pruned=sum(r['result']['pruned'] for r in chosen),
                                         batch_seconds=seconds, median_seconds=median(seconds) if seconds else None,
                                         min_seconds=min(seconds) if seconds else None,
                                         max_seconds=max(seconds) if seconds else None)
    write_json(output/'summary.json', summary)
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check-only', action='store_true',
                        help='check all 60 histories once, without timing repetitions')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    run(args.output.resolve(), check_only=args.check_only)
