#!/usr/bin/env python3
"""Validate every distributed Lebedev-Laikov grid."""
from pathlib import Path
import math
import unittest
import numpy as np
from scipy.spatial import cKDTree

ORDERS = dict(zip(
    [6,14,26,38,50,74,86,110,146,170,194,230,266,302,350,434,590,
     770,974,1202,1454,1730,2030,2354,2702,3074,3470,3890,4334,4802,5294,5810],
    [3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,35,41,47,53,59,65,
     71,77,83,89,95,101,107,113,119,125,131]))

def load(path):
    raw = np.loadtxt(path)
    if raw.ndim != 2 or raw.shape[1] != 3:
        raise ValueError(f"malformed grid {path}")
    az, polar, w = raw.T
    az = np.deg2rad(az); polar = np.deg2rad(polar)
    xyz = np.column_stack((np.cos(az)*np.sin(polar),
                           np.sin(az)*np.sin(polar), np.cos(polar)))
    return xyz, w

def sphere_monomial(i,j,k):
    if i%2 or j%2 or k%2: return 0.0
    return (2*math.gamma((i+1)/2)*math.gamma((j+1)/2)*math.gamma((k+1)/2)
            / math.gamma((i+j+k+3)/2))

class LebedevTests(unittest.TestCase):
    def test_all_grids(self):
        root = Path(__file__).resolve().parents[1]/"lebedev"
        files = {int(p.stem.split('_')[1]):p for p in root.glob('lebedev_*.txt')}
        self.assertEqual(set(files), set(ORDERS))
        for n, degree in ORDERS.items():
            with self.subTest(points=n, degree=degree):
                xyz,w = load(files[n])
                self.assertEqual(len(w), n)
                self.assertTrue(np.isfinite(xyz).all() and np.isfinite(w).all())
                self.assertLess(np.max(abs(np.linalg.norm(xyz,axis=1)-1)), 3e-15)
                self.assertLess(abs(w.sum()-1), 5e-13)
                self.assertFalse(np.any(w == 0))
                # Inversion must preserve the weighted point multiset.
                tree = cKDTree(xyz)
                inversion_error, opposite = tree.query(-xyz)
                self.assertLess(float(inversion_error.max()), 2e-12)
                self.assertLess(float(np.max(abs(w-w[opposite]))), 2e-15)
                # A nearest-neighbour distance of zero (apart from self) would
                # expose duplicate points without relying on decimal rounding.
                nearest = tree.query(xyz, k=2)[0][:,1]
                self.assertGreater(float(nearest.min()), 1e-8)
                # Representative monomials at every total degree through the
                # advertised order; includes mixed even and all odd axes.
                for d in range(degree+1):
                    probes={(d,0,0),(0,d,0),(0,0,d)}
                    if d>=2: probes.add((d-2,2,0))
                    if d>=4: probes.add((d-4,2,2))
                    for i,j,k in probes:
                        got=4*math.pi*np.dot(w,xyz[:,0]**i*xyz[:,1]**j*xyz[:,2]**k)
                        want=sphere_monomial(i,j,k)
                        # Files carry only 15 decimal digits; the worst measured
                        # error is 4.51e-12 and is the constant/weight sum.
                        self.assertLess(abs(got-want), 5e-12)

if __name__ == '__main__': unittest.main()
