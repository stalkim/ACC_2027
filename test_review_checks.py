"""Boundary, enclosure and informative-history regression tests."""
from fractions import Fraction as F
from itertools import product
import unittest
from extension_history import decode, pairwise
from monotone_history import forward
import straight_history as straight
from review_checks import REVIEW, endpoints, curvature_rows, densified


class ReviewTests(unittest.TestCase):
    def test_exact_sqrt(self):
        for q in (F(0), F(1, 3), F(7, 13), F(19, 7)):
            self.assertEqual(straight.sqrt_bracket(q*q), (q, q))

    def test_irrational_sqrt(self):
        for value, bits in product((F(2), F(2, 3), F(7, 11)), (0, 5, 40)):
            low, high = straight.sqrt_bracket(value, bits)
            self.assertLessEqual(low*low, value)
            self.assertGreaterEqual(high*high, value)
            self.assertLessEqual(high-low, F(1, 2**bits))

    def test_straight_exact_threshold(self):
        tau = ((F(3, 5), F(4, 5)), (F(-3, 5), F(4, 5)))
        for alpha in (F(99, 100), F(1), F(101, 100)):
            span = F(1, 30)*alpha
            times = (-span/2, F(0), span/2)
            y = [(F(0), F(4, 5)*t) for t in times]
            r = straight.decode(times, y, F(1, 100), tau)
            self.assertEqual(r['status'], 'AMBIGUOUS' if alpha <= 1 else 'INCOMPATIBLE')

    def test_straight_invalid(self):
        for times, y, b, tau in (([], [], 1, [(1, 0)]), ([0], [(0, 0)], 1, [(2, 0)]),
                                 ([0, 0], [(0, 0)]*2, 1, [(1, 0)]), ([0], [(0, 0)], -1, [(1, 0)])):
            with self.assertRaises(ValueError):
                straight.decode(times, y, b, tau)
        with self.assertRaises(ValueError):
            straight.sqrt_bracket(F(-1))

    def test_review_example(self):
        for method in ('DC', 'MI', 'DMI'):
            for raw, status in ((endpoints(REVIEW), 'AMBIGUOUS'), (REVIEW, 'UNIQUE')):
                result = decode(raw, 'quartic', '1/2', method)
                self.assertEqual(result['status'], status)
                self.assertEqual(result['pruned'], 0)

    def test_dense_preserves_original(self):
        for n in (3, 9, 33):
            raw = densified(n)
            self.assertEqual(len(raw['times']), n)
            for t, y in zip(REVIEW['times'], REVIEW['positions']):
                index = raw['times'].index(F(t))
                self.assertEqual(raw['positions'][index], tuple(map(F, y)))

    def test_gap_equivalence(self):
        times = (F(0), F(1, 3), F(1))
        options = [(F(a, 2), F(b, 2)) for a in range(-1, 3) for b in range(a, 3)]
        for intervals in product(options, repeat=3):
            gap = max(max(intervals[j][0]+F(9, 10)*(times[k]-times[j])-intervals[k][1],
                          intervals[k][0]-intervals[j][1]-F(11, 10)*(times[k]-times[j]))
                      for j in range(3) for k in range(j, 3))
            self.assertEqual(gap > 0, forward(intervals, times) is None)
            self.assertEqual(gap > 0, pairwise(intervals, times) is None)

    def test_curvature_example(self):
        row = next(r for r in curvature_rows() if r['curvature'] == 2
                   and r['rho'] == F(1, 10) and r['span'] == F(3, 50))
        self.assertTrue(row['covered'])
        self.assertFalse(row['general_pass'])
        self.assertTrue(row['exact_projected_pass'])


if __name__ == '__main__':
    unittest.main()
