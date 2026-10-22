"""Resonance sign, width, c-product, and false-stationarity guards."""
import unittest
import numpy as np
from cap_resonance import eigenstate,trajectory,resonance_parameters,prepare


class ResonanceTests(unittest.TestCase):
    def test_invalid_trajectory_inputs(self):
        h=np.diag([.1,.3]);s=np.eye(2);w=np.eye(2)
        etas=np.geomspace(1e-5,.1,31)
        for seed in (0.,-1.,np.nan,np.inf):
            with self.subTest(seed=seed), self.assertRaises(ValueError):
                trajectory(h,s,w,etas,.1-.001j,seed_eta=seed)
        for invalid in (1.,np.ones((7,2)),[],[1.]*7):
            with self.subTest(etas=invalid), self.assertRaises(ValueError):
                trajectory(h,s,w,invalid,.1-.001j)
        with self.assertRaises(ValueError):
            prepare(np.empty((0,0)),np.empty((0,0)),np.empty((0,0)))

    def test_width_and_lifetime_convention(self):
        result=resonance_parameters(.2-.01j)
        self.assertAlmostEqual(result['width_hartree'],.02)
        self.assertAlmostEqual(result['lifetime_au'],50.)
        self.assertIsNone(resonance_parameters(.2+0j)['lifetime_au'])

    def test_c_product_derivative_matches_finite_difference(self):
        h=np.array([[.1,.03],[.03,.25]]);w=np.array([[.2,.1],[.1,2.]])
        eta=.08
        e,d,u,_,_=eigenstate(h,w,eta,target=.1-.01j)
        step=1e-6
        plus=eigenstate(h,w,eta+step,previous=u)[0]
        minus=eigenstate(h,w,eta-step,previous=u)[0]
        self.assertLess(abs(d-(plus-minus)/(2*step)),1e-9)
        incorrect=-1j*(u.conj()@w@u)/(u.conj()@u)
        self.assertGreater(abs(d-incorrect),1e-3)

    def test_linear_trajectory_has_no_false_endpoint_optimum(self):
        result=trajectory(np.diag([.1,.3]),np.eye(2),np.diag([1.,2.]),
                          np.geomspace(1e-5,.1,31),.1-.001j)
        self.assertEqual(result['status'],'no_interior_stationary_point')
        self.assertIsNone(result['selected_candidate'])
        self.assertFalse(result['resonance_validated'])

    def test_ao_scaling_preserves_complex_spectrum(self):
        h=np.array([[.1,.03],[.03,.25]]);s=np.array([[1.,.1],[.1,1.]])
        w=np.array([[.2,.1],[.1,2.]])
        scale=np.diag([3.,.2])
        h1,w1,_=prepare(h,s,w)
        h2,w2,_=prepare(scale@h@scale,scale@s@scale,scale@w@scale)
        np.testing.assert_allclose(np.sort_complex(np.linalg.eigvals(h1-.01j*w1)),
                                   np.sort_complex(np.linalg.eigvals(h2-.01j*w2)),atol=1e-13)
        with self.assertRaises(ValueError):
            prepare(h,s,-np.eye(2))


if __name__=='__main__':
    unittest.main()
