#!/usr/bin/env python3
"""Independent smooth-Voronoi geometry tests (atomic units)."""
import math
import unittest
import numpy as np


def smooth_voronoi(point, atoms, r_cut, smoothing_length=1.0):
    """Return unscaled W; all lengths are bohr.

    Ratios are formed from weights rescaled by their maximum.  With a > 0 the
    largest raw weight is finite and positive, so this cannot divide by zero.
    """
    p = np.asarray(point, dtype=float)
    axyz = np.asarray(atoms, dtype=float)
    if axyz.ndim != 2 or axyz.shape[1] != 3 or len(axyz) == 0:
        raise ValueError("atoms must be a nonempty (n,3) array")
    if smoothing_length <= 0 or r_cut < 0:
        raise ValueError("smoothing_length must be > 0 and r_cut >= 0")
    d2 = np.einsum("ij,ij->i", axyz-p, axyz-p)
    delta = d2 - d2.min() + smoothing_length**2
    scaled = (delta.min()/delta)**2
    rwa2 = np.dot(scaled, d2)/scaled.sum()
    rwa = math.sqrt(max(0.0, rwa2))
    return max(0.0, rwa-r_cut)**2


class SmoothVoronoiTests(unittest.TestCase):
    def setUp(self):
        self.atoms = np.array([[-1.2, 0., 0.], [1.2, 0., 0.]])

    def test_nonnegative_and_one_center(self):
        for r in np.linspace(0, 20, 101):
            got = smooth_voronoi([r, 0, 0], [[0, 0, 0]], 3.0)
            self.assertGreaterEqual(got, 0.0)
            self.assertAlmostEqual(got, max(0.0, r-3.0)**2, places=13)

    def test_atom_order_translation_rotation(self):
        p = np.array([2.1, -0.3, 0.7]); ref = smooth_voronoi(p, self.atoms, 1.7)
        self.assertEqual(ref, smooth_voronoi(p, self.atoms[::-1], 1.7))
        t = np.array([1e5, -2e5, 3e5])
        self.assertAlmostEqual(ref, smooth_voronoi(p+t, self.atoms+t, 1.7), places=10)
        q, _ = np.linalg.qr(np.array([[1.,2.,3.],[-2.,1.,1.],[1.,-4.,2.]]))
        self.assertAlmostEqual(ref, smooth_voronoi(q@p, self.atoms@q.T, 1.7), places=13)

    def test_diatomic_point_group_symmetry(self):
        for p in ([2.3,.4,.8], [-2.3,.4,.8], [2.3,-.4,.8], [2.3,.4,-.8]):
            self.assertAlmostEqual(smooth_voronoi([2.3,.4,.8], self.atoms, 1),
                                   smooth_voronoi(p, self.atoms, 1), places=13)

    def test_voronoi_boundary_continuity(self):
        eps = np.logspace(-2, -10, 9)
        jumps = [abs(smooth_voronoi([e,.3,.2], self.atoms, .5)-
                     smooth_voronoi([-e,.3,.2], self.atoms, .5)) for e in eps]
        self.assertLess(max(jumps), 1e-13)

    def test_onset_is_c1(self):
        h = 1e-6
        self.assertEqual(smooth_voronoi([3,0,0], [[0,0,0]], 3), 0)
        right_derivative = smooth_voronoi([3+h,0,0], [[0,0,0]], 3)/h
        self.assertLess(abs(right_derivative), 2*h)

    def test_extreme_points_are_finite_and_order_invariant(self):
        cases = ([0,0,0], [0,.1,0], [1.2,0,0], [1e150,-1e150,1e150])
        for p in cases:
            x = smooth_voronoi(p, self.atoms, 1)
            y = smooth_voronoi(p, self.atoms[::-1], 1)
            self.assertTrue(math.isfinite(x)); self.assertEqual(x, y)


if __name__ == "__main__":
    unittest.main()
