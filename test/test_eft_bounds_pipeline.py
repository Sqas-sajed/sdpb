import unittest
from fractions import Fraction

from eft_bounds.crossing_pipeline import (
    basis_functional_to_ordinary_terms,
    build_constraint_pmp,
    crossing_to_extremal_transform,
    derive_null_constraints,
    run_reliability_checks,
    symmetric_indices,
)


class CrossingPipelineTests(unittest.TestCase):
    def test_reliability_checks_cover_at_least_five_examples(self):
        result = run_reliability_checks(max_degree=6, m_sq=Fraction(1))
        self.assertTrue(result["passed"])
        self.assertGreaterEqual(len(result["examples"]), 5)

    def test_massless_crossing_to_extremal_is_diagonal_scaling(self):
        crossing_basis, extremal_basis, transform = crossing_to_extremal_transform(
            max_degree=8,
            m_sq=Fraction(0),
        )
        self.assertEqual(crossing_basis, extremal_basis)
        for row_index, basis_index in enumerate(crossing_basis):
            p, q = basis_index
            expected = Fraction((-1) ** q, 2 ** p)
            for col_index, value in enumerate(transform[row_index]):
                if col_index == row_index:
                    self.assertEqual(value, expected)
                else:
                    self.assertEqual(value, 0)

    def test_extremal_basis_functional_recovers_g3_from_stu(self):
        terms = basis_functional_to_ordinary_terms(
            max_degree=6,
            basis_index=(0, 1),
            basis="extremal",
            m_sq=Fraction(0),
        )
        # One exact left-inverse choice is enough to recover g3 on the constrained space.
        self.assertEqual(terms, {(1, 2): Fraction(-1)})

    def test_constraint_pmp_contains_paired_blocks(self):
        pmp, metadata = build_constraint_pmp(max_degree=6, m_sq=Fraction(1), precision=40)
        self.assertIn("PositiveMatrixWithPrefactorArray", pmp)
        self.assertEqual(
            len(pmp["PositiveMatrixWithPrefactorArray"]),
            2 * len(metadata["null_constraints"]),
        )

    def test_constraints_exist_beyond_trivial_cutoff(self):
        constraints = derive_null_constraints(max_degree=6, m_sq=Fraction(1))
        self.assertTrue(constraints)
        symmetric = symmetric_indices(6)
        self.assertTrue(symmetric)


if __name__ == "__main__":
    unittest.main()
