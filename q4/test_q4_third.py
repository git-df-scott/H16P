#!/usr/bin/env python3
"""Four bounded third-strike regressions; no scan or repeated quadrature."""
import os
import resource
for _name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
              "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_name] = "1"
resource.setrlimit(resource.RLIMIT_CPU,(10,10))

from fractions import Fraction as Q
from pathlib import Path
import hashlib
import json
import unittest

import mpmath as mp

from q4_lobe_anchors import primitive_basis
from q4_threshold_path import (primitive_basis_closed,threshold_anchors,
    coefficients_from_r,primitive_value_closed,certify_frozen)


class ThirdStrikeTests(unittest.TestCase):
    def test_exact_escape_threshold_identities(self):
        ratio_bound=Q(601,136136)
        a_max=Q(2,3)*(1/(192*ratio_bound)-1)
        self.assertEqual(a_max,Q(2593,21636))
        self.assertEqual(1/(1-a_max),Q(21636,19043))
        self.assertEqual(Q(1,216)-ratio_bound,Q(395,1837836))
        # This exact positive-power inequality is equivalent to the bound
        # on the first Green-peak integral through t=5/11.
        self.assertLess(Q(11,6)**7,Q(55,27)**6)
        self.assertEqual(Q(6,7)*(Q(55,27)-1),Q(8,9))
        self.assertEqual(11*3**11,1948617)
        self.assertEqual(2**7*5**6,2000000)

    def test_closed_moments_against_positive_period_series(self):
        with mp.workdps(60):
            t=mp.mpf(1)/2
            direct_series=primitive_basis(t)
            closed=primitive_basis_closed(t)
            self.assertLess(max(abs(x-y) for x,y in zip(direct_series,closed)),
                            mp.mpf("1e-52"))

    def test_analytic_threshold_path_anchor_residuals(self):
        with mp.workdps(60):
            for r in (mp.mpf(2)/5,mp.mpf(5)/11,mp.mpf(3)/4):
                anchors=threshold_anchors(r)
                self.assertEqual(anchors[0],r)
                self.assertTrue(0<anchors[0]<anchors[1]<anchors[2]<1)
                coefficients=coefficients_from_r(r)
                self.assertTrue(1<coefficients[2]<mp.mpf(54)/31)
                for t in anchors:
                    self.assertLess(abs(primitive_value_closed(t,coefficients)),
                                    mp.mpf("1e-52"))

    def test_frozen_rational_certificate_hash_and_box(self):
        directory=Path(__file__).parent
        frozen=json.loads((directory/"data/third_threshold_certificate.json").read_text())
        self.assertEqual(hashlib.sha256((directory/"q4_threshold_path.py").read_bytes()).hexdigest(),
                         frozen["script_sha256"])
        self.assertEqual(hashlib.sha256((directory/"q4_lobe_certificate.py").read_bytes()).hexdigest(),
                         frozen["series_helper_sha256"])
        # Replay exact rational signs and analytic tails, not merely JSON flags.
        self.assertEqual(certify_frozen(),frozen)
        box=frozen["certified_parameter_box"]
        radius=Q(box["closed_linf_radius"])
        last=Q(frozen["witnesses"][-1]["t"])
        perturbation=radius*(last**2/(1-last)+last**3/(3*(1-last)))
        self.assertEqual(perturbation,Q(box["uniform_witness_perturbation_bound"]))
        self.assertLess(perturbation,Q(1,10**6))
        self.assertGreater(Q(frozen["primitive_root_intervals"][0]["endpoints"][0]),Q(5,11))
        self.assertFalse(frozen["original_Abelian_integral_five_zero_claim"])


if __name__ == "__main__":
    unittest.main()
