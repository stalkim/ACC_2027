"""Tests for inverse bounds, feasible histories and branch decisions."""
import unittest
from fractions import Fraction as F
from monotone_history import (CHARTS, RATES, backward, check_witness, curve,
                              decode, forward, inverse_bracket, preimage)


def data(phases, lam=F(1, 2), error=F(1, 2000)):
    return {'times': [str(F(k, 50)) for k in range(len(phases))],
            'positions': [[str(v) for v in curve(q, lam)] for q in phases],
            'linf_errors': [str(error)]*len(phases)}


class Tests(unittest.TestCase):
    def test_affine_inverse(self):
        for target in (F(0), F(1, 3), F(1)):
            lo, hi = inverse_bracket(lambda q: 2*q+1, 2*target+1, (F(0), F(1)), 24)
            self.assertLessEqual(lo, target)
            self.assertGreaterEqual(hi, target)
            self.assertLessEqual(hi-lo, F(1, 2**24))

    def test_endpoint_equalities(self):
        f = lambda q: q*q
        for value in (F(1), F(4)):
            outer, inner = preimage(f, (value, value), (F(1), F(2)), 0)
            self.assertEqual(outer, inner)
            self.assertEqual(outer[0], outer[1])

    def test_empty_preimage(self):
        self.assertEqual(preimage(lambda q: q, (F(2), F(3)), (F(0), F(1)), 24), (None, None))

    def test_outer_is_not_existence(self):
        outer, inner = preimage(lambda q: q*q, (F(2), F(2)), (F(1), F(2)), 2)
        self.assertIsNotNone(outer)
        self.assertIsNone(inner)

    def test_backward_shared_history(self):
        times = (F(0), F(1), F(2))
        reach = forward(((F(0), F(1)), (F(1), F(2)), (F(2), F(3))), times)
        q = backward(reach, times)
        for a, b in zip(q, q[1:]):
            self.assertLessEqual(RATES[0], b-a)
            self.assertLessEqual(b-a, RATES[1])

    def test_rate_infeasibility(self):
        self.assertIsNone(forward(((F(0), F(0)), (F(2), F(2))), (F(0), F(1))))

    def test_positive_and_negative_polynomial_branches(self):
        for lam in (F(1, 2), F(1), F(2)):
            for name, sign in (('negative', -1), ('positive', 1)):
                raw = data((F(sign)-F(1, 100), F(sign)+F(1, 100)), lam)
                for directional in (False, True):
                    r = decode(raw, lam, directional=directional)
                    self.assertEqual(r['selected'], name)
                    q = r['branches'][name]['witness']
                    history = (tuple(map(F,raw['times'])), tuple(tuple(map(F,y)) for y in raw['positions']), tuple(map(F,raw['linf_errors'])))
                    self.assertTrue(check_witness(q, history, lam, CHARTS[name]))
                    lo, hi = r['branches'][name]['outer']
                    self.assertLessEqual(lo, F(sign)+F(1,100))
                    self.assertGreaterEqual(hi, F(sign)+F(1,100))

    def test_two_exact_histories(self):
        raw = {'times':['0','1/50'], 'positions':[['0','0','0']]*2, 'linf_errors':['3/50']*2}
        for lam in (F(1,2),F(1),F(2)):
            self.assertEqual(decode(raw,lam)['status'],'AMBIGUOUS')

    def test_zero_small_error_incompatible(self):
        raw = {'times':['0','1/50'], 'positions':[['0','0','0']]*2, 'linf_errors':['1/2000']*2}
        for direction in (False,True):
            self.assertEqual(decode(raw,1,directional=direction)['status'],'INCOMPATIBLE')

    def test_one_zero_is_ambiguous(self):
        raw = {'times':['0'], 'positions':[['0','0','0']], 'linf_errors':['0']}
        self.assertEqual(decode(raw,1)['status'],'AMBIGUOUS')

    def test_coverage_failure(self):
        raw = data((F(-2), F(-2)+F(1,50)))
        self.assertEqual(decode(raw,1)['status'],'OUTSIDE_COVERAGE')

    def test_nonmonotone_time_rejected(self):
        raw = data((F(-1), F(-1)))
        raw['times']=['0','0']
        with self.assertRaises(ValueError):
            decode(raw,1)


if __name__ == '__main__':
    unittest.main()
