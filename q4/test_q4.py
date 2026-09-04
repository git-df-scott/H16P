#!/usr/bin/env python3
import unittest
from fractions import Fraction as Q

import numpy as np

from q4_integrals import basis_float, basis_orbit_float, q4_coefficients
from q4_search import track_parameter_path, zhao_reduced_filter


class Q4Tests(unittest.TestCase):
    def test_original_normal_form(self):
        data = q4_coefficients(1.0)
        self.assertAlmostEqual(data["b"] ** 2 + data["c"] ** 2, 4.0)
        self.assertAlmostEqual(data["kappa"], 2.0)

    def test_independent_integrals(self):
        area = basis_float(4.0, 2.5, order=192)
        orbit = basis_orbit_float(4.0, 2.5)
        self.assertLess(np.max(np.abs(area - orbit)), 1.0e-10)

    @staticmethod
    def reduced_coefficients(kappa, beta0, p_at_beta):
        # P2 is linear and satisfies P2(kappa) + kappa - beta0 = 0.
        a1 = p_at_beta / (beta0 - kappa) - 1
        a0 = beta0 - kappa - kappa * a1
        return (a0, a1, Q(0), beta0, Q(1))

    def test_corrected_strip_has_no_kappa_four_cutoff(self):
        coefficients = self.reduced_coefficients(Q(4), Q(0), Q(1))
        self.assertEqual(zhao_reduced_filter(Q(4), coefficients), (True, "survives"))

    def test_strip_boundaries_exact(self):
        kappa = Q(2)
        lower = Q(8, 31)
        for beta in (lower - Q(1, 1000), lower, Q(1)):
            coefficients = self.reduced_coefficients(kappa, beta, Q(1))
            self.assertEqual(zhao_reduced_filter(kappa, coefficients), (False, "beta0_strip"))
        coefficients = self.reduced_coefficients(kappa, lower + Q(1, 1000), Q(1))
        self.assertEqual(zhao_reduced_filter(kappa, coefficients), (True, "survives"))

    def test_cubic_comment_bound_exact(self):
        kappa, beta = Q(2), Q(1, 2)
        threshold = Q(25, 3456)
        for p_value in (threshold - Q(1, 100000), threshold):
            coefficients = self.reduced_coefficients(kappa, beta, p_value)
            self.assertEqual(zhao_reduced_filter(kappa, coefficients), (False, "zhao_comment_bound"))
        # 1/50 is above the valid cubic threshold and below the old linear one.
        coefficients = self.reduced_coefficients(kappa, beta, Q(1, 50))
        self.assertEqual(zhao_reduced_filter(kappa, coefficients), (True, "survives"))
        for factor in (Q(-7), Q(1, 10**20)):
            self.assertEqual(zhao_reduced_filter(kappa, tuple(factor * c for c in coefficients)), (True, "survives"))

    def test_beta1_zero_exact(self):
        self.assertEqual(zhao_reduced_filter(Q(2), (Q(1), Q(1), Q(1), Q(1), Q(0))), (False, "beta1_zero"))

    def test_continuation_interface(self):
        records = track_parameter_path([(4.0, (1, 0, 0, 0))], grid_points=31)
        self.assertEqual(records[0]["roots"], [])


if __name__ == "__main__":
    unittest.main()
