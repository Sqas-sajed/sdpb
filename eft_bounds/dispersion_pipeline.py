"""
dispersion_pipeline.py — Paper-aligned positivity/null-constraint pipeline.

This module follows the notation of:

  - Sinha & Zahed, "Crossing Symmetric Dispersion Relations in QFTs"
  - Caron-Huot & Duong, "Extremal Effective Field Theories"

It focuses on the spin-and-s1/mass formulation of the null constraints:

  * Sinha-Zahed eq. (4): W_{n-m,m} from partial waves via B^{(ell)}_{n,m}(s1)
  * Sinha-Zahed eq. (11): null constraints for m > n in terms of
      D^{(n,m)}_{ell,alpha} / s1^{2n+m}
  * Extremal-EFT eqs. (3.16) and (3.19): lower/upper two-sided bounds

The implementation exports:

  * exact D-coefficients
  * exact d=4 eq.(11) spin polynomials in ell divided by a power of s1
  * JSON artifacts for a chosen EFT cutoff
  * translated eq.(3) -> eq.(2.3) metadata
  * example lower/upper SDPB PMP files
"""

from __future__ import annotations

import argparse
import json
import os
from fractions import Fraction
from math import factorial
from typing import Dict, List, Sequence, Tuple

from .crossing_pipeline import build_translation_metadata, generate_two_sided_bound_pmps
from .physics import fraction_to_str

Index = Tuple[int, int]


def _pochhammer(a: Fraction, n: int) -> Fraction:
    result = Fraction(1)
    for i in range(n):
        result *= a + i
    return result


def _gamma_positive_integer(n: int) -> Fraction:
    if n <= 0:
        raise ValueError(f"Gamma is singular for non-positive integer argument {n}.")
    return Fraction(factorial(n - 1))


def _evaluate_polynomial(coeffs: Sequence[Fraction], x: Fraction) -> Fraction:
    total = Fraction(0)
    x_power = Fraction(1)
    for coeff in coeffs:
        total += coeff * x_power
        x_power *= x
    return total


def _interpolate_polynomial(xs: Sequence[Fraction], ys: Sequence[Fraction]) -> List[Fraction]:
    if len(xs) != len(ys):
        raise ValueError("x/y sample sizes must match.")
    size = len(xs)
    matrix = []
    for x, y in zip(xs, ys):
        row = [Fraction(1)]
        for _ in range(1, size):
            row.append(row[-1] * x)
        row.append(y)
        matrix.append(row)

    for col in range(size):
        pivot = None
        for row in range(col, size):
            if matrix[row][col] != 0:
                pivot = row
                break
        if pivot is None:
            raise ValueError("Interpolation matrix is singular.")
        matrix[col], matrix[pivot] = matrix[pivot], matrix[col]
        pivot_value = matrix[col][col]
        matrix[col] = [entry / pivot_value for entry in matrix[col]]
        for row in range(size):
            if row == col or matrix[row][col] == 0:
                continue
            factor = matrix[row][col]
            matrix[row] = [
                current - factor * pivot_entry
                for current, pivot_entry in zip(matrix[row], matrix[col])
            ]

    coeffs = [row[-1] for row in matrix]
    while coeffs and coeffs[-1] == 0:
        coeffs.pop()
    return coeffs or [Fraction(0)]


def compute_D_coefficient(
    n: int,
    m: int,
    ell: int,
    alpha: Fraction = Fraction(1, 2),
) -> Fraction:
    """
    Compute D^{(n,m)}_{ell,alpha} from Sinha-Zahed eq. (11).
    """
    if n < 1 or m <= n:
        raise ValueError("eq.(11) applies only for m > n >= 1.")

    total = Fraction(0)
    gamma_factor = _gamma_positive_integer(m - n)
    for j in range(n, m + 1):
        numerator = (
            Fraction((-4) ** j)
            * _pochhammer(Fraction(-ell, 2), j)
            * _pochhammer(alpha + Fraction(ell, 2), j)
            * Fraction(3 * j - m - 2 * n)
            * gamma_factor
        )
        denominator = (
            Fraction(factorial(j))
            * Fraction(factorial(m - j))
            * Fraction(factorial(j - n))
            * _pochhammer(alpha + Fraction(1, 2), j)
        )
        total += numerator / denominator
    return total


def compute_D_closed_form_example(ell: int, alpha: Fraction = Fraction(1, 2)) -> Fraction:
    """
    Closed-form example quoted below eq. (11) for (n,m)=(1,2).
    """
    numerator = (
        Fraction(2)
        * ell
        * (ell + 2 * alpha)
        * (-11 - 10 * alpha + 2 * ell * (ell + 2 * alpha))
    )
    denominator = (2 * alpha + 1) * (2 * alpha + 3)
    return numerator / denominator


def enumerate_eq11_null_constraints(max_order: int) -> List[Index]:
    """
    Enumerate (n,m) with m > n >= 1 and 2n + m <= max_order.
    """
    constraints: List[Index] = []
    for n in range(1, max_order + 1):
        for m in range(n + 1, max_order + 1):
            if 2 * n + m <= max_order:
                constraints.append((n, m))
    return constraints


def d4_D_spin_polynomial(n: int, m: int) -> List[Fraction]:
    """
    Fit D^{(n,m)}_{ell,1/2} as an exact polynomial in ell.
    """
    degree = 2 * m
    xs = [Fraction(k) for k in range(degree + 1)]
    ys = [compute_D_coefficient(n, m, int(x), Fraction(1, 2)) for x in xs]
    return _interpolate_polynomial(xs, ys)


def d4_eq11_kernel_polynomial(n: int, m: int) -> List[Fraction]:
    """
    Return the d=4 eq.(11) integrand kernel polynomial in ell, excluding 1/(pi s1^{2n+m}).

    In d=4, alpha=1/2 and C_ell^{(1/2)}(1)=1, so the kernel is:
        (2 ell + 1) D^{(n,m)}_{ell,1/2}.
    """
    d_poly = d4_D_spin_polynomial(n, m)
    result = [Fraction(0)] * (len(d_poly) + 1)
    for power, coeff in enumerate(d_poly):
        result[power] += coeff
        result[power + 1] += 2 * coeff
    while result and result[-1] == 0:
        result.pop()
    return result or [Fraction(0)]


def build_eq11_constraints_json(
    max_order: int,
    precision: int = 80,
) -> Dict[str, object]:
    constraints = []
    for n, m in enumerate_eq11_null_constraints(max_order):
        constraints.append({
            "n": n,
            "m": m,
            "s1_denominator_power": 2 * n + m,
            "overall_prefactor": "1/pi",
            "d4_D_polynomial_in_spin_ell": [
                fraction_to_str(coeff, precision) for coeff in d4_D_spin_polynomial(n, m)
            ],
            "d4_eq11_kernel_polynomial_in_spin_ell": [
                fraction_to_str(coeff, precision) for coeff in d4_eq11_kernel_polynomial(n, m)
            ],
            "formula": (
                "integrand = sum_ell (2*ell+1) a_ell(s1) D_{ell,1/2}^{(n,m)}"
                f" / (pi * s1^{2 * n + m})"
            ),
        })
    return {
        "dimension": 4,
        "alpha": "1/2",
        "max_order": max_order,
        "constraints": constraints,
    }


def extremal_d4_functional_polynomials() -> Dict[str, Dict[str, object]]:
    """
    d=4 polynomial data from Extremal EFT eqs. (3.23)-(3.24), after m^2 = M^2(1+x).
    """
    return {
        "g2": {
            "denominator_power_in_(1+x)": 2,
            "spin_polynomial_in_J": ["1"],
            "x_prefactor_after_clearing_denominator": ["1", "2", "1"],
        },
        "M2_g3": {
            "denominator_power_in_(1+x)": 3,
            "spin_polynomial_in_J": ["3", "-2"],
            "spin_polynomial_variable": "J(J+1)",
            "x_prefactor_after_clearing_denominator": ["1", "1"],
        },
        "M4_g4": {
            "denominator_power_in_(1+x)": 4,
            "spin_polynomial_in_J": ["1/2"],
            "x_prefactor_after_clearing_denominator": ["1"],
        },
        "M4_n4": {
            "denominator_power_in_(1+x)": 4,
            "spin_polynomial_in_J": ["0", "-16", "2"],
            "spin_polynomial_variable": "J(J+1)",
            "x_prefactor_after_clearing_denominator": ["1"],
        },
    }


def run_reliability_checks(precision: int = 80) -> Dict[str, object]:
    """
    Run paper-focused checks on the dispersion/null-constraint rewrite.
    """
    checks = []

    for ell in [2, 4, 6]:
        direct = compute_D_coefficient(1, 2, ell, Fraction(1, 2))
        closed = compute_D_closed_form_example(ell, Fraction(1, 2))
        checks.append({
            "name": f"D_(1,2)_closed_form_ell_{ell}",
            "passed": direct == closed,
            "direct": fraction_to_str(direct, precision),
            "closed_form": fraction_to_str(closed, precision),
        })

    for (n, m, ell) in [(1, 3, 5), (2, 3, 4)]:
        poly = d4_D_spin_polynomial(n, m)
        reconstructed = _evaluate_polynomial(poly, Fraction(ell))
        direct = compute_D_coefficient(n, m, ell, Fraction(1, 2))
        checks.append({
            "name": f"D_polynomial_reconstruction_{n}_{m}_ell_{ell}",
            "passed": reconstructed == direct,
            "direct": fraction_to_str(direct, precision),
            "reconstructed": fraction_to_str(reconstructed, precision),
        })

    all_passed = all(check["passed"] for check in checks)
    return {"passed": all_passed, "examples": checks}


def write_pipeline_outputs(
    output_dir: str,
    max_order: int,
    precision: int = 80,
    max_spin: int = 10,
) -> Dict[str, object]:
    os.makedirs(output_dir, exist_ok=True)

    eq11_payload = build_eq11_constraints_json(max_order=max_order, precision=precision)
    translation_payload = build_translation_metadata(
        max_degree=max_order,
        m_sq=Fraction(0),
        precision=precision,
    )
    extremal_payload = extremal_d4_functional_polynomials()
    checks_payload = run_reliability_checks(precision=precision)

    files = {
        "eq11_constraints": os.path.join(output_dir, "eq11_null_constraints.json"),
        "eq3_to_eq23_translation": os.path.join(output_dir, "eq3_to_eq23_translation.json"),
        "extremal_d4_functionals": os.path.join(output_dir, "extremal_d4_functionals.json"),
        "reliability_checks": os.path.join(output_dir, "dispersion_reliability_checks.json"),
    }

    for path, payload in [
        (files["eq11_constraints"], eq11_payload),
        (files["eq3_to_eq23_translation"], translation_payload),
        (files["extremal_d4_functionals"], extremal_payload),
        (files["reliability_checks"], checks_payload),
    ]:
        with open(path, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2)

    files.update(
        generate_two_sided_bound_pmps(
            output_dir=output_dir,
            max_degree=max_order,
            target_index=(0, 1),
            normalization_index=(1, 0),
            basis="extremal",
            m_sq=Fraction(0),
            max_spin=max_spin,
            precision=precision,
        )
    )

    return {"files": files, "checks": checks_payload}


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate the paper-aligned dispersion/null-constraint EFT pipeline."
    )
    parser.add_argument("--output-dir", required=True, help="Directory for generated JSON/PMP files.")
    parser.add_argument("--max-order", type=int, default=8, help="EFT cutoff order in 2n+m.")
    parser.add_argument("--precision", type=int, default=80, help="Decimal precision for JSON strings.")
    parser.add_argument("--max-spin", type=int, default=10, help="Maximum spin for generated PMP files.")
    args = parser.parse_args()

    result = write_pipeline_outputs(
        output_dir=args.output_dir,
        max_order=args.max_order,
        precision=args.precision,
        max_spin=args.max_spin,
    )
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
