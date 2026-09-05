"""Tests for polynomial paths and the two temporal solvers."""
from fractions import Fraction as F
from itertools import product
import unittest
from monotone_history import CHARTS, RATES, forward, backward, preimage
from extension_history import decode, pairwise, polynomial, verify_witness


class ExtensionTests(unittest.TestCase):
    def test_pairwise_exhaustive_small_chains(self):
        # Exact endpoints, rate equality, nonlocal inconsistency and singletons.
        choices = [(F(a, 2), F(b, 2)) for a in range(-1, 3) for b in range(a, 3)]
        times = (F(0), F(1, 3), F(1))
        for intervals in product(choices, repeat=3):
            rec, oracle = forward(intervals, times), pairwise(intervals, times)
            self.assertEqual(rec is None, oracle is None)
            if rec is not None:
                self.assertEqual(rec[-1], oracle[0])
                q = oracle[1]
                self.assertTrue(all(lo <= z <= hi for (lo, hi), z in zip(intervals, q)))
                self.assertTrue(all(RATES[0]*(times[k]-times[k-1]) <= q[k]-q[k-1]
                                    <= RATES[1]*(times[k]-times[k-1]) for k in (1, 2)))

    def test_nonlocal_contradiction(self):
        chain = [(F(0), F(0)), (F(-10), F(10)), (F(3), F(3))]
        self.assertIsNone(pairwise(chain, (F(0), F(1), F(2))))

    def test_precision_positive_margin(self):
        # q=t lies at distance mu=1/8 from every endpoint; h=1/16 < mu.
        times = (F(0), F(1, 3), F(2, 3))
        mu, h = F(1, 8), F(1, 16)
        inner = [(t-mu+h, t+mu-h) for t in times]
        self.assertIsNotNone(forward(inner, times))

    def test_precision_exclusion_and_equality(self):
        # The two singleton parameters violate the upper rate by g=1/10.
        times = (F(0), F(1))
        for h, survives in ((F(1, 10), True), (F(1, 20), True), (F(1, 32), False)):
            outer = [(-h, h), (F(6, 5)-h, F(6, 5)+h)]
            self.assertEqual(forward(outer, times) is not None, survives)
        # A feasible boundary-only chain need not have an inner witness.
        self.assertIsNotNone(pairwise([(F(0), F(0)), (F(9, 10), F(9, 10))], times))
        self.assertIsNone(forward([None, None], times))

    def test_nonlinear_inverse_inclusions(self):
        for steps in (0, 3, 11):
            out, inside = preimage(lambda q: q**3, (F(1, 3), F(2, 3)), (F(0), F(1)), steps)
            self.assertLessEqual(out[0]**3, F(1, 3))
            self.assertGreaterEqual(out[1]**3, F(2, 3))
            if inside is not None:
                self.assertGreaterEqual(inside[0]**3, F(1, 3))
                self.assertLessEqual(inside[1]**3, F(2, 3))

    def test_factored_witness_checker(self):
        for family in ('quartic', 'quintic'):
            for q in (F(-6, 5), F(-1), F(4, 5), F(11, 10)):
                raw = dict(times=[F(0)], positions=[polynomial(q, family, F(3, 4))], linf_errors=[F(0)])
                domain = CHARTS['negative' if q < 0 else 'positive']
                self.assertTrue(verify_witness((q,), raw, family, F(3, 4), domain))
                self.assertFalse(verify_witness((q+F(1, 100),), raw, family, F(3, 4), domain))

    def test_decoder_separate_development_input(self):
        for family in ('quartic', 'quintic'):
            q = (F(19, 20), F(99, 100), F(103, 100))
            raw = dict(times=[F(0), F(1, 25), F(2, 25)],
                       positions=[polynomial(z, family, F(3, 4)) for z in q],
                       linf_errors=[F(1, 1000)]*3)
            results = [decode(raw, family, F(3, 4), m) for m in ('DC', 'MI', 'DMI')]
            self.assertEqual({r['status'] for r in results}, {'UNIQUE'})
            self.assertEqual({r['selected'] for r in results}, {'positive'})
            self.assertEqual(results[0]['branches']['positive']['outer'], results[1]['branches']['positive']['outer'])

    def test_closed_pruning_endpoint(self):
        # An admitted chord exactly at the pruning boundary must not be removed.
        raw = dict(times=[F(0), F(1, 100)], positions=[(0, 0, 0), (F(27, 2000), 0, 0)], linf_errors=[0, 0])
        result = decode(raw, 'quartic', 1)
        self.assertEqual(result['pruned'], 1)

    def test_input_validation(self):
        with self.assertRaises(ValueError):
            decode(dict(times=[], positions=[], linf_errors=[]), 'quartic', 1)
        with self.assertRaises(ValueError):
            decode({}, 'undeclared', 1)


if __name__ == '__main__':
    unittest.main()
