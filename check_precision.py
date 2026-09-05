"""Exact checks of the analytical margin application; no decoder calls."""
import argparse
from fractions import Fraction as F
import json
from pathlib import Path


def check(directory):
    # Bounds on the full union of charts, lambda <= 2.
    width, qmax = F(1, 2), F(5, 4)
    quartic_bound = 2*(3*qmax*qmax+1+(4*qmax**3+2*qmax)/5)
    quintic_bound = 2*(5*qmax**4-1)  # derivative positive: >= lambda*149/256
    assert max(2*qmax, quartic_bound, quintic_bound) < 23
    assert width/F(2**18) < F(1, 2000)/230
    cases = json.loads((directory/'inputs.json').read_text())
    truth = json.loads((directory/'truth.json').read_text())
    checked = 0
    for case in cases:
        if case['kind'] != 'nominal':
            continue
        raw = case['history']
        phases = tuple(map(F, truth[case['case']]['histories'][0]))
        for q, point, error in zip(phases, raw['positions'], raw['linf_errors']):
            error, scale = F(error), F(case['scale'])
            domain = (F(-5, 4), F(-3, 4)) if q < 0 else (F(3, 4), F(5, 4))
            assert min(q-domain[0], domain[1]-q) >= F(81, 500)
            x = (q-1)*(q+1)
            factor = q+q*q/5 if case['family'] == 'quartic' else q*(q*q+1)
            position = (x, scale*x*factor, F(0))
            assert all(abs(a-F(b)) <= F(9, 10)*error for a, b in zip(position, point))
            assert width/F(2**18) < error/230 < F(81, 500)
        checked += 1
    assert checked == 48
    print('48 nominal histories: derivative/noise/chart-margin premises checked; B=18 sufficient for inner feasibility.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('results', type=Path)
    check(parser.parse_args().results)
