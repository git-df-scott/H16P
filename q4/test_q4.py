#!/usr/bin/env python3
import unittest

import numpy as np

from q4_integrals import basis_float, basis_orbit_float, q4_coefficients
from q4_search import track_parameter_path, zhao_five_filter


class Q4Tests(unittest.TestCase):
    def test_original_normal_form(self):
        data = q4_coefficients(1.0)
        self.assertAlmostEqual(data["b"] ** 2 + data["c"] ** 2, 4.0)
        self.assertAlmostEqual(data["kappa"], 2.0)

    def test_independent_integrals(self):
        area = basis_float(4.0, 2.5, order=192)
        orbit = basis_orbit_float(4.0, 2.5)
        self.assertLess(np.max(np.abs(area - orbit)), 1.0e-10)

    def test_zhao_strip_excludes_kappa_four(self):
        keep, _reason = zhao_five_filter(4.0, (1, 2, 3, 4))
        self.assertFalse(keep)

    def test_continuation_interface(self):
        records = track_parameter_path([(4.0, (1, 0, 0, 0))], grid_points=31)
        self.assertEqual(records[0]["roots"], [])


if __name__ == "__main__":
    unittest.main()
