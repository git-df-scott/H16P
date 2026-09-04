#!/usr/bin/env python3
"""Two independent original-area checks for the frozen third-strike shots."""
import os
for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
            "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[key] = "1"
import hashlib
import json
import resource
import time
from pathlib import Path
from q4_integrals import basis_mp, alpha_beta_from_mu
from q4_reconstruction import reconstruct, original_values, mu_from_universal
from q4_search import zhao_reduced_filter


def main():
    resource.setrlimit(resource.RLIMIT_CPU, (10, 10))
    os.nice(10)
    started = time.process_time()
    folder = Path(__file__).with_name("data")
    results = []
    for file, index, t in (("third_tuned_shoot.json", 1, .7),
                           ("third_shape_shoot.json", 2, .95)):
        row = json.loads((folder/file).read_text())["rows"][index]
        k, a = row["kappa"], row["a"]
        coefficients = tuple(map(float, row["A_B_eta"]))
        mu = mu_from_universal(k, *coefficients)
        sol = reconstruct(a, *coefficients, t_end=max(t, .99))
        value = float(original_values(a, sol, t))
        s = k-(k-1)*t
        independent = float(sum(m*v for m, v in zip(mu, basis_mp(k, s, dps=40))))
        difference = abs(independent-value)
        assert difference < 2e-11
        passed, reason = zhao_reduced_filter(k, alpha_beta_from_mu(k, mu))
        assert passed
        results.append({"record_file": file, "record_index": index,
            "kappa": k, "t": t, "s": s, "original_mu": mu.tolist(),
            "scalar_PF_I": value, "independent_original_area_I": independent,
            "absolute_difference": difference,
            "corrected_Zhao_filter": [passed, reason]})
        print(file, "kappa", k, "independent discrepancy", difference)
    record = {"status": "NUMERICAL_ONLY", "checks": results,
        "cpu_seconds": time.process_time()-started,
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (folder/"third_independent_checks.json").write_text(json.dumps(record, indent=2)+"\n")
    print("Both original-area comparisons and corrected filters passed.")


if __name__ == "__main__":
    main()
