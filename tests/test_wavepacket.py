"""Physical invariants for the independent radial packet validation."""
import unittest
import numpy as np
from cap_wavepacket import radial_packet, radial_grid, summarize


class RadialPacketTests(unittest.TestCase):
    def test_exact_half_line_normalization_and_boundary(self):
        r = np.arange(0, 800, .02)
        for energy in (.05, .2, .5):
            for time in (0, 20, 100):
                u = radial_packet(r, time, energy)
                self.assertEqual(u[0], 0)
                self.assertAlmostEqual(float(np.trapezoid(abs(u)**2, r)), 1, places=11)

    def test_free_reference_and_absorber_sign(self):
        free = radial_grid(.2, 0, dx=.1, dt=.02, duration=10)
        cap = radial_grid(.2, .001, dx=.1, dt=.02, duration=10)
        self.assertLess(summarize(free)['maximum_inner_squared_error'], 1e-5)
        self.assertAlmostEqual(free['norms'][0], free['norms'][-1], places=11)
        self.assertLess(cap['norms'][-1], .999)
        self.assertLess(summarize(cap)['maximum_sampled_norm_increase'], 1e-12)
        with self.assertRaises(ValueError):
            radial_grid(.2, -.001)


if __name__ == '__main__':
    unittest.main()
