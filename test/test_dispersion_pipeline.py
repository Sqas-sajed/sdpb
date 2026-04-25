import json
import os
import tempfile
import unittest
from fractions import Fraction

from eft_bounds.dispersion_pipeline import (
    build_eq11_constraints_json,
    compute_D_closed_form_example,
    compute_D_coefficient,
    d4_D_spin_polynomial,
    d4_eq11_kernel_polynomial,
    enumerate_eq11_null_constraints,
    run_reliability_checks,
    write_pipeline_outputs,
)


def eval_poly(coeffs, x):
    total = Fraction(0)
    power = Fraction(1)
    for coeff in coeffs:
        total += coeff * power
        power *= x
    return total


class DispersionPipelineTests(unittest.TestCase):
    def test_closed_form_example_matches_paper(self):
        for ell in [2, 4, 6]:
            self.assertEqual(
                compute_D_coefficient(1, 2, ell, Fraction(1, 2)),
                compute_D_closed_form_example(ell, Fraction(1, 2)),
            )

    def test_d4_D_polynomial_reconstructs_direct_values(self):
        coeffs = d4_D_spin_polynomial(1, 3)
        self.assertEqual(eval_poly(coeffs, Fraction(5)), compute_D_coefficient(1, 3, 5))

    def test_d4_eq11_kernel_polynomial_reconstructs_direct_values(self):
        coeffs = d4_eq11_kernel_polynomial(2, 3)
        direct = (2 * 4 + 1) * compute_D_coefficient(2, 3, 4, Fraction(1, 2))
        self.assertEqual(eval_poly(coeffs, Fraction(4)), direct)

    def test_enumeration_respects_cutoff(self):
        self.assertEqual(
            enumerate_eq11_null_constraints(8),
            [(1, 2), (1, 3), (1, 4), (1, 5), (1, 6), (2, 3), (2, 4)],
        )

    def test_json_payload_and_reliability_checks(self):
        payload = build_eq11_constraints_json(8, precision=40)
        self.assertTrue(payload["constraints"])
        checks = run_reliability_checks(precision=40)
        self.assertTrue(checks["passed"])
        self.assertGreaterEqual(len(checks["examples"]), 5)

    def test_full_pipeline_writes_outputs(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            result = write_pipeline_outputs(tmpdir, max_order=6, precision=30, max_spin=4)
            self.assertTrue(result["checks"]["passed"])
            for path in result["files"].values():
                self.assertTrue(os.path.exists(path))
            with open(result["files"]["eq11_constraints"], "r", encoding="utf-8") as handle:
                payload = json.load(handle)
            self.assertIn("constraints", payload)


if __name__ == "__main__":
    unittest.main()
