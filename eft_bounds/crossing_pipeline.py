"""
crossing_pipeline.py — Exact crossing-basis derivations and SDPB helpers.

This module connects three representations of an identical-scalar EFT
amplitude up to a chosen Mandelstam cutoff degree:

1. Ordinary polynomial basis in (s, t) with u = 4m² - s - t
2. Crossing-symmetric basis from Sinha-Zahed eq. (3):
      M = Σ W_{p,q} x^p y^q
   with x = -(st + tu + us), y = -stu
3. Extremal-EFT symmetric basis from Caron-Huot/Duong eq. (2.3):
      M = Σ G_{r,q} E2^r E3^q
   with E2 = s² + t² + u², E3 = stu

The code computes exact rational change-of-basis matrices, derives null
constraints on the ordinary (s, t) coefficients from crossing symmetry,
and generates JSON outputs (including optional PMP JSON encodings of the
constraints) for a chosen cutoff.
"""

from __future__ import annotations

import argparse
import json
import os
from collections import defaultdict
from fractions import Fraction
from typing import Dict, Iterable, List, Sequence, Tuple

from .physics import fraction_to_str
from .pmp_generator import generate_pmp_for_linear_functional_bound, write_pmp_json

BivariatePolynomial = Dict[Tuple[int, int], Fraction]
Index = Tuple[int, int]


def ordinary_indices(max_degree: int) -> List[Index]:
    """Return ordinary (s^a t^b) indices with a + b <= max_degree."""
    indices: List[Index] = []
    for total in range(max_degree + 1):
        for a in range(total + 1):
            indices.append((a, total - a))
    return indices


def symmetric_indices(max_degree: int) -> List[Index]:
    """Return symmetric-basis indices (p, q) with 2p + 3q <= max_degree."""
    indices: List[Index] = []
    for q in range(max_degree // 3 + 1):
        for p in range((max_degree - 3 * q) // 2 + 1):
            indices.append((p, q))
    indices.sort(key=lambda item: (2 * item[0] + 3 * item[1], item[1], item[0]))
    return indices


def _poly_add(lhs: BivariatePolynomial, rhs: BivariatePolynomial) -> BivariatePolynomial:
    result = defaultdict(Fraction, lhs)
    for key, value in rhs.items():
        result[key] += value
        if result[key] == 0:
            del result[key]
    return dict(result)


def _poly_scale(poly: BivariatePolynomial, factor: Fraction) -> BivariatePolynomial:
    if factor == 0:
        return {}
    return {key: factor * value for key, value in poly.items() if factor * value != 0}


def _poly_multiply(lhs: BivariatePolynomial, rhs: BivariatePolynomial) -> BivariatePolynomial:
    result: Dict[Index, Fraction] = defaultdict(Fraction)
    for (a1, b1), c1 in lhs.items():
        for (a2, b2), c2 in rhs.items():
            result[(a1 + a2, b1 + b2)] += c1 * c2
    return {key: value for key, value in result.items() if value != 0}


def _poly_pow(poly: BivariatePolynomial, exponent: int) -> BivariatePolynomial:
    result: BivariatePolynomial = {(0, 0): Fraction(1)}
    if exponent == 0:
        return result
    base = dict(poly)
    power = exponent
    while power > 0:
        if power & 1:
            result = _poly_multiply(result, base)
        power >>= 1
        if power:
            base = _poly_multiply(base, base)
    return result


def x_polynomial(m_sq: Fraction = Fraction(1)) -> BivariatePolynomial:
    """
    Return x = -(st + tu + us) as a polynomial in s and t.

    With u = 4m² - s - t, this equals:
        x = s² + st + t² - 4m²(s + t)
    """
    mu = 4 * m_sq
    return {
        (2, 0): Fraction(1),
        (1, 1): Fraction(1),
        (0, 2): Fraction(1),
        (1, 0): -mu,
        (0, 1): -mu,
    }


def y_polynomial(m_sq: Fraction = Fraction(1)) -> BivariatePolynomial:
    """
    Return y = -stu as a polynomial in s and t.

    With u = 4m² - s - t, this equals:
        y = s²t + st² - 4m² st
    """
    mu = 4 * m_sq
    return {
        (2, 1): Fraction(1),
        (1, 2): Fraction(1),
        (1, 1): -mu,
    }


def e2_polynomial(m_sq: Fraction = Fraction(0)) -> BivariatePolynomial:
    """Return E2 = s² + t² + u² as a polynomial in s and t."""
    mu = 4 * m_sq
    return {
        (2, 0): Fraction(2),
        (1, 1): Fraction(2),
        (0, 2): Fraction(2),
        (1, 0): -2 * mu,
        (0, 1): -2 * mu,
        (0, 0): mu * mu,
    }


def e3_polynomial(m_sq: Fraction = Fraction(0)) -> BivariatePolynomial:
    """Return E3 = stu as a polynomial in s and t."""
    return _poly_scale(y_polynomial(m_sq), Fraction(-1))


def expand_crossing_monomial(index: Index, m_sq: Fraction = Fraction(1)) -> BivariatePolynomial:
    """Expand x^p y^q into ordinary s,t monomials."""
    p, q = index
    return _poly_multiply(_poly_pow(x_polynomial(m_sq), p), _poly_pow(y_polynomial(m_sq), q))


def expand_extremal_monomial(index: Index, m_sq: Fraction = Fraction(0)) -> BivariatePolynomial:
    """Expand E2^r E3^q into ordinary s,t monomials."""
    r, q = index
    return _poly_multiply(_poly_pow(e2_polynomial(m_sq), r), _poly_pow(e3_polynomial(m_sq), q))


def _matrix_from_expansions(
    ordinary_basis: Sequence[Index],
    symmetric_basis: Sequence[Index],
    expansion_fn,
    m_sq: Fraction,
) -> List[List[Fraction]]:
    matrix: List[List[Fraction]] = []
    expanded_columns = [expansion_fn(index, m_sq) for index in symmetric_basis]
    for ordinary_index in ordinary_basis:
        matrix.append([
            column.get(ordinary_index, Fraction(0)) for column in expanded_columns
        ])
    return matrix


def _transpose(matrix: Sequence[Sequence[Fraction]]) -> List[List[Fraction]]:
    if not matrix:
        return []
    return [list(column) for column in zip(*matrix)]


def _rref(matrix: Sequence[Sequence[Fraction]]) -> Tuple[List[List[Fraction]], List[int]]:
    work = [list(row) for row in matrix]
    if not work:
        return [], []
    row_count = len(work)
    col_count = len(work[0])
    pivot_columns: List[int] = []
    pivot_row = 0

    for col in range(col_count):
        pivot = None
        for row in range(pivot_row, row_count):
            if work[row][col] != 0:
                pivot = row
                break
        if pivot is None:
            continue
        work[pivot_row], work[pivot] = work[pivot], work[pivot_row]
        pivot_value = work[pivot_row][col]
        work[pivot_row] = [entry / pivot_value for entry in work[pivot_row]]

        for row in range(row_count):
            if row == pivot_row or work[row][col] == 0:
                continue
            factor = work[row][col]
            work[row] = [
                current - factor * pivot_entry
                for current, pivot_entry in zip(work[row], work[pivot_row])
            ]

        pivot_columns.append(col)
        pivot_row += 1
        if pivot_row == row_count:
            break

    return work, pivot_columns


def _nullspace(matrix: Sequence[Sequence[Fraction]]) -> List[List[Fraction]]:
    rref_matrix, pivot_columns = _rref(matrix)
    if not rref_matrix:
        return []
    col_count = len(rref_matrix[0])
    free_columns = [col for col in range(col_count) if col not in pivot_columns]
    basis: List[List[Fraction]] = []

    for free_col in free_columns:
        vector = [Fraction(0)] * col_count
        vector[free_col] = Fraction(1)
        for row, pivot_col in enumerate(pivot_columns):
            vector[pivot_col] = -rref_matrix[row][free_col]
        basis.append(vector)
    return basis


def _invert_square(matrix: Sequence[Sequence[Fraction]]) -> List[List[Fraction]]:
    size = len(matrix)
    if size == 0:
        return []
    work = [
        list(row) + [Fraction(int(i == j)) for j in range(size)]
        for i, row in enumerate(matrix)
    ]

    for col in range(size):
        pivot = None
        for row in range(col, size):
            if work[row][col] != 0:
                pivot = row
                break
        if pivot is None:
            raise ValueError("Matrix is singular and cannot be inverted exactly.")
        work[col], work[pivot] = work[pivot], work[col]
        pivot_value = work[col][col]
        work[col] = [entry / pivot_value for entry in work[col]]
        for row in range(size):
            if row == col or work[row][col] == 0:
                continue
            factor = work[row][col]
            work[row] = [
                current - factor * pivot_entry
                for current, pivot_entry in zip(work[row], work[col])
            ]
    return [row[size:] for row in work]


def _left_inverse_from_full_rank_matrix(matrix: Sequence[Sequence[Fraction]]) -> List[List[Fraction]]:
    """
    Return L such that L * matrix = I for a full-column-rank matrix.
    """
    transposed = _transpose(matrix)
    _, pivot_columns = _rref(transposed)
    if len(pivot_columns) != len(transposed):
        raise ValueError("Expected a full-column-rank matrix.")
    selected_rows = [list(matrix[row]) for row in pivot_columns]
    inverse = _invert_square(selected_rows)

    row_count = len(matrix)
    basis_count = len(transposed)
    left_inverse = [[Fraction(0) for _ in range(row_count)] for _ in range(basis_count)]
    for basis_index in range(basis_count):
        for local_row, global_row in enumerate(pivot_columns):
            left_inverse[basis_index][global_row] = inverse[basis_index][local_row]
    return left_inverse


def derive_null_constraints(max_degree: int, m_sq: Fraction = Fraction(1)) -> List[Dict[Index, Fraction]]:
    """
    Derive all exact full-crossing null constraints on ordinary coefficients.

    The result is the left nullspace of the map from symmetric-basis coefficients
    to ordinary (s,t) coefficients.
    """
    ordinary_basis = ordinary_indices(max_degree)
    sym_basis = symmetric_indices(max_degree)
    matrix = _matrix_from_expansions(ordinary_basis, sym_basis, expand_crossing_monomial, m_sq)
    left_nullspace = _nullspace(_transpose(matrix))

    constraints: List[Dict[Index, Fraction]] = []
    for vector in left_nullspace:
        constraint = {
            ordinary_basis[i]: value for i, value in enumerate(vector) if value != 0
        }
        if constraint:
            constraints.append(_normalize_constraint(constraint))
    return constraints


def _normalize_constraint(constraint: Dict[Index, Fraction]) -> Dict[Index, Fraction]:
    ordered = sorted(constraint.items())
    pivot_coeff = ordered[0][1]
    normalized = {key: value / pivot_coeff for key, value in constraint.items()}
    if next(iter(sorted(normalized.items())))[1] < 0:
        normalized = {key: -value for key, value in normalized.items()}
    return normalized


def _json_ready_terms(data: Dict[Index, Fraction], precision: int) -> List[Dict[str, object]]:
    return [
        {
            "index": [index[0], index[1]],
            "coefficient": fraction_to_str(coeff, precision),
        }
        for index, coeff in sorted(data.items())
    ]


def build_crossing_to_ordinary_map(
    max_degree: int,
    m_sq: Fraction = Fraction(1),
    precision: int = 50,
) -> Dict[str, object]:
    ordinary_basis = ordinary_indices(max_degree)
    sym_basis = symmetric_indices(max_degree)
    matrix = _matrix_from_expansions(ordinary_basis, sym_basis, expand_crossing_monomial, m_sq)
    expansions = []
    for column_index, sym_index in enumerate(sym_basis):
        expansions.append({
            "crossing_index": [sym_index[0], sym_index[1]],
            "ordinary_terms": _json_ready_terms(
                {
                    ordinary_basis[row]: matrix[row][column_index]
                    for row in range(len(ordinary_basis))
                    if matrix[row][column_index] != 0
                },
                precision,
            ),
        })
    return {
        "max_degree": max_degree,
        "m_sq": fraction_to_str(m_sq, precision),
        "crossing_basis": [list(index) for index in sym_basis],
        "ordinary_basis": [list(index) for index in ordinary_basis],
        "expansions": expansions,
    }


def crossing_to_extremal_transform(
    max_degree: int,
    m_sq: Fraction = Fraction(0),
) -> Tuple[List[Index], List[Index], List[List[Fraction]]]:
    ordinary_basis = ordinary_indices(max_degree)
    crossing_basis = symmetric_indices(max_degree)
    extremal_basis = symmetric_indices(max_degree)
    crossing_matrix = _matrix_from_expansions(
        ordinary_basis, crossing_basis, expand_crossing_monomial, m_sq
    )
    extremal_matrix = _matrix_from_expansions(
        ordinary_basis, extremal_basis, expand_extremal_monomial, m_sq
    )
    left_inverse = _left_inverse_from_full_rank_matrix(extremal_matrix)
    transform: List[List[Fraction]] = []
    for left_row in left_inverse:
        row: List[Fraction] = []
        for cross_col in range(len(crossing_basis)):
            entry = Fraction(0)
            for ordinary_row, coeff in enumerate(left_row):
                entry += coeff * crossing_matrix[ordinary_row][cross_col]
            row.append(entry)
        transform.append(row)
    return crossing_basis, extremal_basis, transform


def build_translation_metadata(
    max_degree: int,
    m_sq: Fraction = Fraction(0),
    precision: int = 50,
) -> Dict[str, object]:
    crossing_basis, extremal_basis, transform = crossing_to_extremal_transform(max_degree, m_sq)
    rows = []
    for extremal_index, row in zip(extremal_basis, transform):
        rows.append({
            "extremal_index": [extremal_index[0], extremal_index[1]],
            "extremal_label": extremal_label(extremal_index),
            "crossing_terms": _json_ready_terms(
                {
                    crossing_basis[col]: row[col]
                    for col in range(len(crossing_basis))
                    if row[col] != 0
                },
                precision,
            ),
        })
    return {
        "max_degree": max_degree,
        "m_sq": fraction_to_str(m_sq, precision),
        "notes": [
            "Extremal basis means powers of E2 = s^2 + t^2 + u^2 and E3 = stu.",
            "Crossing basis means powers of x = -(st+tu+us) and y = -stu.",
        ],
        "rows": rows,
    }


def extremal_label(index: Index) -> str:
    r, q = index
    degree = 2 * r + 3 * q
    if q == 0:
        if degree == 0:
            return "constant"
        return f"g{degree}"
    if q == 1:
        return f"g{degree}"
    return f"g{degree}_power_{q}_of_stu"


def basis_functional_to_ordinary_terms(
    max_degree: int,
    basis_index: Index,
    basis: str = "crossing",
    m_sq: Fraction = Fraction(0),
) -> Dict[Index, Fraction]:
    """
    Return the exact ordinary-coefficient linear functional that extracts one
    symmetric-basis coefficient.
    """
    ordinary_basis = ordinary_indices(max_degree)
    sym_basis = symmetric_indices(max_degree)
    if basis == "crossing":
        matrix = _matrix_from_expansions(ordinary_basis, sym_basis, expand_crossing_monomial, m_sq)
    elif basis == "extremal":
        matrix = _matrix_from_expansions(ordinary_basis, sym_basis, expand_extremal_monomial, m_sq)
    else:
        raise ValueError(f"Unknown basis '{basis}'.")

    if basis_index not in sym_basis:
        raise ValueError(f"Basis index {basis_index} is outside the cutoff degree {max_degree}.")

    left_inverse = _left_inverse_from_full_rank_matrix(matrix)
    row = left_inverse[sym_basis.index(basis_index)]
    return {
        ordinary_basis[i]: coeff for i, coeff in enumerate(row) if coeff != 0
    }


def build_constraint_pmp(
    max_degree: int,
    m_sq: Fraction = Fraction(1),
    precision: int = 80,
) -> Tuple[Dict[str, object], Dict[str, object]]:
    """
    Encode all null constraints as paired 1x1 PMP inequality blocks.

    The PMP is meant as a machine-readable translation artifact. It is not a
    standalone optimization problem unless the caller adds a normalization and
    objective of interest.
    """
    ordinary_basis = ordinary_indices(max_degree)
    constraints = derive_null_constraints(max_degree, m_sq)
    zero_vector = ["0" for _ in ordinary_basis]
    blocks = []

    for sign in (Fraction(1), Fraction(-1)):
        for constraint in constraints:
            coeffs = []
            for ordinary_index in ordinary_basis:
                coeff = sign * constraint.get(ordinary_index, Fraction(0))
                coeffs.append([fraction_to_str(coeff, precision)])
            blocks.append({
                "DampedRational": {
                    "base": "1",
                    "constant": "1",
                    "poles": [],
                },
                "polynomials": [[coeffs]],
            })

    return (
        {
            "objective": list(zero_vector),
            "normalization": list(zero_vector),
            "PositiveMatrixWithPrefactorArray": blocks,
        },
        {
            "ordinary_basis": [list(index) for index in ordinary_basis],
            "null_constraints": [
                {"terms": _json_ready_terms(constraint, precision)} for constraint in constraints
            ],
        },
    )


def generate_two_sided_bound_pmps(
    output_dir: str,
    max_degree: int,
    target_index: Index,
    normalization_index: Index,
    basis: str = "extremal",
    m_sq: Fraction = Fraction(0),
    max_spin: int = 10,
    precision: int = 80,
) -> Dict[str, str]:
    """
    Generate lower/upper PMP JSON files for one symmetric-basis ratio.
    """
    os.makedirs(output_dir, exist_ok=True)
    objective_terms = basis_functional_to_ordinary_terms(
        max_degree=max_degree,
        basis_index=target_index,
        basis=basis,
        m_sq=m_sq,
    )
    normalization_terms = basis_functional_to_ordinary_terms(
        max_degree=max_degree,
        basis_index=normalization_index,
        basis=basis,
        m_sq=m_sq,
    )

    lower = generate_pmp_for_linear_functional_bound(
        objective_terms=objective_terms,
        normalization_terms=normalization_terms,
        max_spin=max_spin,
        max_order=max_degree,
        use_null_constraints=True,
        crossing_type="crossing_basis",
        m_sq=m_sq,
        precision=precision,
        bound_direction="lower",
    )
    upper = generate_pmp_for_linear_functional_bound(
        objective_terms=objective_terms,
        normalization_terms=normalization_terms,
        max_spin=max_spin,
        max_order=max_degree,
        use_null_constraints=True,
        crossing_type="crossing_basis",
        m_sq=m_sq,
        precision=precision,
        bound_direction="upper",
    )

    lower_path = os.path.join(output_dir, f"lower_{basis}_{target_index[0]}_{target_index[1]}_over_{normalization_index[0]}_{normalization_index[1]}.json")
    upper_path = os.path.join(output_dir, f"upper_{basis}_{target_index[0]}_{target_index[1]}_over_{normalization_index[0]}_{normalization_index[1]}.json")
    write_pmp_json(lower, lower_path)
    write_pmp_json(upper, upper_path)
    return {"lower": lower_path, "upper": upper_path}


def _evaluate_constraints(
    coefficients: Dict[Index, Fraction],
    constraints: Iterable[Dict[Index, Fraction]],
) -> List[Fraction]:
    residuals = []
    for constraint in constraints:
        total = Fraction(0)
        for index, coeff in constraint.items():
            total += coeff * coefficients.get(index, Fraction(0))
        residuals.append(total)
    return residuals


def run_reliability_checks(max_degree: int = 6, m_sq: Fraction = Fraction(1)) -> Dict[str, object]:
    """
    Validate the exact pipeline on at least five examples.
    """
    constraints = derive_null_constraints(max_degree, m_sq)
    ordinary_basis = ordinary_indices(max_degree)
    sym_basis = symmetric_indices(max_degree)
    matrix = _matrix_from_expansions(ordinary_basis, sym_basis, expand_crossing_monomial, m_sq)

    examples = {
        "constant": {(0, 0): Fraction(1)},
        "linear_x": {(1, 0): Fraction(1)},
        "linear_y": {(0, 1): Fraction(1)},
        "quadratic_x2_plus_xy": {(2, 0): Fraction(1), (1, 1): Fraction(1)},
        "mixed_x_plus_y": {(1, 0): Fraction(2), (0, 1): Fraction(-3)},
        "cubic_combo": {(3, 0): Fraction(1), (0, 2): Fraction(2), (1, 1): Fraction(1)},
    }

    results = []
    for name, symmetric_coeffs in examples.items():
        ordinary_coeffs: Dict[Index, Fraction] = defaultdict(Fraction)
        for column_index, sym_index in enumerate(sym_basis):
            weight = symmetric_coeffs.get(sym_index, Fraction(0))
            if weight == 0:
                continue
            for row_index, ordinary_index in enumerate(ordinary_basis):
                ordinary_coeffs[ordinary_index] += weight * matrix[row_index][column_index]
        residuals = _evaluate_constraints(ordinary_coeffs, constraints)
        max_residual = max((abs(value) for value in residuals), default=Fraction(0))
        results.append({
            "name": name,
            "symmetric_coefficients": _json_ready_terms(symmetric_coeffs, 30),
            "max_residual": str(max_residual),
            "passed": max_residual == 0,
        })

    all_passed = all(result["passed"] for result in results)
    return {"passed": all_passed, "examples": results}


def write_pipeline_outputs(
    output_dir: str,
    max_degree: int,
    m_sq: Fraction,
    precision: int,
    max_spin: int,
) -> Dict[str, object]:
    os.makedirs(output_dir, exist_ok=True)

    crossing_map = build_crossing_to_ordinary_map(max_degree, m_sq, precision)
    translation = build_translation_metadata(max_degree, m_sq, precision)
    constraints = derive_null_constraints(max_degree, m_sq)
    constraints_payload = {
        "max_degree": max_degree,
        "m_sq": fraction_to_str(m_sq, precision),
        "constraints": [
            {"terms": _json_ready_terms(constraint, precision)} for constraint in constraints
        ],
    }
    constraint_pmp, constraint_metadata = build_constraint_pmp(max_degree, m_sq, precision)
    checks = run_reliability_checks(max_degree=max_degree, m_sq=m_sq)

    files = {
        "crossing_map": os.path.join(output_dir, "crossing_to_ordinary.json"),
        "translation": os.path.join(output_dir, "eq3_to_eq23_translation.json"),
        "constraints": os.path.join(output_dir, "null_constraints.json"),
        "constraint_pmp": os.path.join(output_dir, "null_constraints_pmp.json"),
        "constraint_pmp_metadata": os.path.join(output_dir, "null_constraints_pmp_metadata.json"),
        "checks": os.path.join(output_dir, "reliability_checks.json"),
    }

    for path, payload in [
        (files["crossing_map"], crossing_map),
        (files["translation"], translation),
        (files["constraints"], constraints_payload),
        (files["constraint_pmp_metadata"], constraint_metadata),
        (files["checks"], checks),
    ]:
        with open(path, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2)

    write_pmp_json(constraint_pmp, files["constraint_pmp"])

    files.update(
        generate_two_sided_bound_pmps(
            output_dir=output_dir,
            max_degree=max_degree,
            target_index=(1, 0),
            normalization_index=(0, 0),
            basis="extremal",
            m_sq=m_sq,
            max_spin=max_spin,
            precision=precision,
        )
    )
    return {"files": files, "checks": checks}


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Exact crossing-basis derivations and SDPB JSON generation."
    )
    parser.add_argument("--output-dir", required=True, help="Directory for generated JSON files.")
    parser.add_argument("--max-degree", type=int, default=6, help="Maximum Mandelstam degree.")
    parser.add_argument("--mass-squared", type=str, default="1", help="External mass squared as a rational number.")
    parser.add_argument("--precision", type=int, default=80, help="Decimal precision for JSON strings.")
    parser.add_argument("--max-spin", type=int, default=10, help="Maximum spin for generated PMP files.")
    args = parser.parse_args()

    m_sq = Fraction(args.mass_squared)
    result = write_pipeline_outputs(
        output_dir=args.output_dir,
        max_degree=args.max_degree,
        m_sq=m_sq,
        precision=args.precision,
        max_spin=args.max_spin,
    )
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
