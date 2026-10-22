#!/usr/bin/env python3
"""Regression checks on locally generated Gaussian/reference logs."""
from pathlib import Path
import re
import unittest

ROOT = Path(__file__).resolve().parent / "gaussian"

def scalar(path, label):
    text = path.read_text(errors="replace")
    m = re.search(re.escape(label) + r"\s*=\s*([-+0-9.Ee]+)", text)
    if not m: raise AssertionError(f"missing {label} in {path}")
    return float(m.group(1))

class GaussianReferences(unittest.TestCase):
    def test_shell_encodings_and_route_semantics(self):
        pure = (ROOT/"one_center_spdfghi_pure.fchk").read_text()
        cart = (ROOT/"one_center_spdfghi_cartesian.fchk").read_text()
        self.assertIn("Number of basis functions                  I               50", pure)
        self.assertRegex(pure, r"Shell types\s+I\s+N=\s+7\s+0\s+-1\s+-2\s+-3\s+-4\s+-5\s+-6")
        self.assertIn("Number of basis functions                  I               85", cart)
        self.assertRegex(cart, r"Shell types\s+I\s+N=\s+7\s+0\s+-1\s+2\s+3\s+4\s+5\s+6")

    def test_all_shell_w1_overlap(self):
        for kind in ("pure", "cartesian"):
            p = ROOT/f"one_center_spdfghi_{kind}_w1.log"
            self.assertLess(scalar(p, "max |S_num - S_analytical|"), 2e-12)

    def test_n2_radial_convergence(self):
        errors = [scalar(ROOT/f"n2_pure_w1_gauss{n}.log",
                         "max |S_num - S_analytical|") for n in (100,200,300)]
        self.assertGreater(errors[0], errors[1])
        self.assertGreater(errors[1], errors[2])
        self.assertLess(errors[2], 2e-12)

if __name__ == "__main__": unittest.main()
