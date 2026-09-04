#!/usr/bin/env python3
"""Claude hostile regression wrapper: runs each independent checker as a
separate process (the repository's own scripts install per-process CPU
fuses, so these are deliberately not imported into one interpreter).
Run from anywhere:  python3 audit/test_claude_hostile.py"""
import os, subprocess, sys, unittest
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
SCRIPTS = {
    "strike1_identities": ("claude_check_strike1.py", "ALL STRIKE-1 IDENTITY CHECKS PASSED"),
    "threshold_certificate": ("claude_check_threshold.py", "THRESHOLD CERTIFICATE INDEPENDENTLY VERIFIED"),
    "strike2_ingredients": ("claude_check_strike2.py", "STRIKE-2 INGREDIENT CHECKS PASSED"),
    "boundary_controls": ("claude_check_boundaries.py", "BOUNDARY CONTROLS PASSED"),
    "endpoint_identities": ("claude_check_endpoint_identities.py", "Variation-of-parameters decomposition Y=Phi*y+P*y2: OK"),
    "reconstruction_random": ("claude_check_reconstruction.py", "PASS"),
}
class ClaudeHostileTests(unittest.TestCase):
    pass
def _make(name, script, marker):
    def test(self):
        proc = subprocess.run([sys.executable, os.path.join(HERE, script)], cwd=ROOT,
                              capture_output=True, text=True, timeout=1800)
        self.assertEqual(proc.returncode, 0, proc.stdout[-2000:]+proc.stderr[-2000:])
        self.assertIn(marker, proc.stdout, proc.stdout[-2000:])
    test.__name__ = "test_"+name
    return test
for name, (script, marker) in SCRIPTS.items():
    setattr(ClaudeHostileTests, "test_"+name, _make(name, script, marker))
if __name__ == "__main__":
    unittest.main(verbosity=2)
