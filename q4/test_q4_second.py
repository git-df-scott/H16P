#!/usr/bin/env python3
"""Independent bounded second-strike regression checks; no search."""
import os
import resource
for name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
             "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[name] = "1"
resource.setrlimit(resource.RLIMIT_CPU, (10, 10))

import unittest
from fractions import Fraction as Q
import mpmath as mp
import numpy as np

from q4_integrals import alpha_beta_from_mu, basis_mp
from q4_lobe_certificate import PARAMETERS, primitive_enclosure
from q4_lobe_anchors import from_primitive_anchors
from q4_reconstruction import reconstruct, original_values, mu_from_universal


class SecondStrikeTests(unittest.TestCase):
    @staticmethod
    def direct_primitive(t, coefficients):
        A, B, eta = coefficients
        def integrand(u):
            F = mp.hyp2f1(mp.mpf(1)/6, mp.mpf(5)/6, 1, u)
            Fp = mp.mpf(5)/36*mp.hyp2f1(mp.mpf(7)/6, mp.mpf(11)/6, 2, u)
            M = 1-6*(1-u)*Fp/F
            return u*F*(A+B*u-1+(u-eta)*M)
        return mp.quad(integrand, [0, t])

    @staticmethod
    def to_mp(q):
        return mp.mpf(q.numerator)/q.denominator

    def test_exact_enclosure_against_independent_quadrature(self):
        # This cross-check is numerical; the certificate's tail proof is exact.
        with mp.workdps(80):
            coefficients = tuple(Q(PARAMETERS[key]) for key in ("A", "B", "eta"))
            for t in (Q(1, 8), Q(7, 8)):
                low, high, _ = primitive_enclosure(*coefficients, t)
                independent = self.direct_primitive(self.to_mp(t),
                                      tuple(map(self.to_mp, coefficients)))
                # At 1/8 the analytic enclosure is narrower than 80-digit
                # quadrature precision. The outward 30-place bounds remain
                # meaningful, while the 7/8 exact tail is resolved directly.
                if t == Q(1, 8):
                    from q4_lobe_certificate import outward_decimal_interval
                    low, high = map(mp.mpf, outward_decimal_interval(low, high))
                else:
                    low, high = self.to_mp(low), self.to_mp(high)
                self.assertLess(low, independent)
                self.assertLess(independent, high)

    def test_primitive_anchors_against_independent_quadrature(self):
        with mp.workdps(60):
            anchors = (mp.mpf(1)/4, mp.mpf(1)/2, mp.mpf(3)/4)
            coefficients = from_primitive_anchors(anchors)
            self.assertLess(abs(self.direct_primitive(anchors[1], coefficients)),
                            mp.mpf("1e-50"))

    def test_original_coefficient_transport_round_trip(self):
        A, B, eta = (float(Q(PARAMETERS[key])) for key in ("A", "B", "eta"))
        for k in (2., 4.):
            alpha0, alpha1, alpha2, beta0, beta1 = alpha_beta_from_mu(
                                   k, mu_from_universal(k, A, B, eta))
            np.testing.assert_allclose((-(alpha1+2*k*alpha2), (k-1)*alpha2,
                               (k-beta0)/(k-1), beta1), (A, B, eta, 1),
                               rtol=2e-12, atol=2e-12)
            self.assertLess(abs(alpha0+k*alpha1+k*k*alpha2+k-beta0), 2e-12)

    def test_corrected_pf_against_independent_original_area(self):
        A, B, eta = (float(Q(PARAMETERS[key])) for key in ("A", "B", "eta"))
        k, a, t = 4., .75, .5
        mu = mu_from_universal(k, A, B, eta)
        independent = float(sum(m*v for m, v in zip(mu, basis_mp(k, k-(k-1)*t,
                                                                               dps=40))))
        value = float(original_values(a, reconstruct(a, A, B, eta), t))
        self.assertLess(abs(independent-value), 2e-11)


if __name__ == "__main__":
    unittest.main()
