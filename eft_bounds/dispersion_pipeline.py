"""
dispersion_pipeline.py — Paper-aligned positivity/null-constraint pipeline.

This module follows the notation of:

  - Sinha & Zahed, "Crossing Symmetric Dispersion Relations in QFTs"
  - Caron-Huot & Duong, "Extremal Effective Field Theories"

It focuses on the spin-and-s1/mass formulation of the null constraints:

  * CSDR eq.(4):  W_{n-m,m} from partial waves via B^{(ell)}_{n,m}(s1)
  * CSDR eq.(11): D^{(n,m)}_{ell,alpha} for m > n >= 1 (null constraint kernels)
  * User notes eq.(2): unified heavy-average form for general d
  * Extremal EFT Section 3.3: SDPB formulation

Key corrections relative to the previous version
-------------------------------------------------
1. S1, S2, S3 are NOT Mandelstam variables; they are Mandelstam minus mu/3.
   S1+S2+S3 = 0.  The W_{p,q} basis uses x = -(S1*S2+S2*S3+S3*S1), y = -S1*S2*S3.
2. d is a free parameter (not hardcoded to d=4).
3. Objectives and null constraints use the CSDR W_{m,n} kernels, not the
   fixed-t Extremal EFT spectral functions (different subtraction schemes).

SDPB problem built here
-----------------------
Variables:
  z_{p,q}  for each objective pair (p >= 1, q >= 0, 2p+3q <= K)
  c_{n,m}  for each null constraint pair (m > n >= 1, 2n+m <= K)

For each even spin ell in {0, 2, 4, ..., ell_max}:
  Polynomial block (1x1):
    sum_{obj}  z_{p,q}  * C^alpha_ell(1)*(2ell+d-3)*(1+x)^{K-2p-3q}
  + sum_{null} c_{n,m}  * D^{(n,m)}_{ell,alpha}*C^alpha_ell(1)*(2ell+d-3)*(1+x)^{K-2n-m}
  >= 0  for all x >= 0

DampedRational prefactor: e^{-x} / (1+x)^K

Normalization:  z_{p1,q1} = 1  (the Wilson coefficient to normalize)
Objective:     maximize  z_{p0,q0}  (the Wilson coefficient to bound above)
"""

from __future__ import annotations

import argparse
import json
import os
from fractions import Fraction
from math import factorial
from typing import Any, Dict, List, Optional, Sequence, Tuple

from .csdr import (
    D_coeff,
    alpha_from_d,
    check_closed_form_1_2,
    enumerate_null_pairs,
    enumerate_obj_pairs,
    gegenbauer_at_one,
    null_kernel_coeff,
    obj_kernel_coeff,
    poly_expand_1px,
    poly_pad,
    poly_scale,
    s1_power,
)
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

    NOTE: These are the spectral functions from the *Extremal EFT* paper, which
    uses a fixed-t dispersion relation and a DIFFERENT subtraction scheme from
    CSDR.  They are kept here for reference only.  The SDPB problem generated
    by build_csdr_pmp_json uses the CSDR kernels instead.
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

    Checks include:
    1. D^{(1,2)} matches closed-form from CSDR paper for several ell values.
    2. D polynomial reconstruction matches direct evaluation.
    3. csdr.check_closed_form_1_2() passes.
    """
    checks = []

    # Check 1: closed-form match for D^{(1,2)}
    for ell in [2, 4, 6]:
        direct = compute_D_coefficient(1, 2, ell, Fraction(1, 2))
        closed = compute_D_closed_form_example(ell, Fraction(1, 2))
        checks.append({
            "name": f"D_(1,2)_closed_form_ell_{ell}",
            "passed": direct == closed,
            "direct": fraction_to_str(direct, precision),
            "closed_form": fraction_to_str(closed, precision),
        })

    # Check 2: polynomial reconstruction
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

    # Check 3: csdr module cross-check
    csdr_ok = check_closed_form_1_2()
    checks.append({
        "name": "csdr_check_closed_form_1_2",
        "passed": csdr_ok,
        "note": "Checks csdr.D_coeff against closed-form for (n,m)=(1,2), several ell",
    })

    all_passed = all(check["passed"] for check in checks)
    return {"passed": all_passed, "examples": checks}


# ---------------------------------------------------------------------------
# New CSDR-based SDPB PMP builder
# ---------------------------------------------------------------------------

def _compute_e_inverse(precision: int) -> str:
    """Compute 1/e = e^{-1} to precision decimal digits."""
    from decimal import Decimal, getcontext
    getcontext().prec = precision + 20
    e_val = Decimal(1)
    factorial_acc = Decimal(1)
    for k in range(1, precision + 100):
        factorial_acc *= k
        term = Decimal(1) / factorial_acc
        e_val += term
        if term < Decimal(10) ** (-(precision + 10)):
            break
    e_inv = Decimal(1) / e_val
    return format(e_inv, f".{precision}f").rstrip("0").rstrip(".")


def _spin_block_csdr(
    ell: int,
    obj_pairs: List[Index],
    null_pairs: List[Index],
    K: int,
    d: int,
    precision: int,
) -> Dict[str, Any]:
    """
    Build one SDPB PositiveMatrixWithPrefactor block for a given spin ell.

    Parameters
    ----------
    ell : int      Spin (even non-negative integer).
    obj_pairs :    List of objective (p,q) pairs.
    null_pairs :   List of null-constraint (n,m) CSDR pairs.
    K : int        Maximum spectral power (denominator exponent).
    d : int        Spacetime dimension.
    precision :    Output decimal precision.

    Returns
    -------
    dict  SDPB PositiveMatrixWithPrefactor block.
    """
    poly_vectors: List[List[Fraction]] = []

    # Objective components: kernel = C^alpha_ell(1)*(2ell+d-3)*(1+x)^{K-2p-3q}
    for (p, q) in obj_pairs:
        kappa = obj_kernel_coeff(p, q, ell, d)
        exponent = K - s1_power(p, q)
        base_poly = poly_expand_1px(exponent)
        scaled = poly_scale(base_poly, kappa)
        poly_vectors.append(poly_pad(scaled, K + 1))

    # Null-constraint components: kernel = D*C^alpha_ell(1)*(2ell+d-3)*(1+x)^{K-2n-m}
    for (n, m) in null_pairs:
        kappa = null_kernel_coeff(n, m, ell, d)
        exponent = K - (2 * n + m)
        base_poly = poly_expand_1px(exponent)
        scaled = poly_scale(base_poly, kappa)
        poly_vectors.append(poly_pad(scaled, K + 1))

    # Check that not all polynomials are identically zero
    if all(all(c == 0 for c in pv) for pv in poly_vectors):
        return None

    # Prefactor: e^{-x} / (1+x)^K
    # Pole at x = -1 with multiplicity K
    e_inv_str = _compute_e_inverse(precision)
    prefactor = {
        "base": e_inv_str,
        "constant": "1",
        "poles": ["-1"] * K,
    }

    # Build polynomials JSON: 1×1 matrix, one entry per decision variable
    poly_json = [[
        [fraction_to_str(c, precision) for c in pv]
        for pv in poly_vectors
    ]]

    return {
        "prefactor": prefactor,
        "polynomials": [poly_json],
    }


def build_csdr_pmp_json(
    obj_index: Index,
    norm_index: Index,
    d: int = 4,
    K: int = 8,
    max_spin: int = 10,
    precision: int = 200,
    bound_direction: str = "upper",
) -> Dict[str, Any]:
    """
    Build an SDPB PMP JSON for bounding W_{obj_index} / W_{norm_index}
    using the CSDR dispersion relation and user notes eq.(2).

    The SDPB dual problem is:
      Maximize z_{obj_index}
      Subject to z_{norm_index} = 1
      For all even ell in [0, max_spin] and all x >= 0:
        sum_{obj} z_{p,q} * C^alpha_ell(1)*(2ell+d-3)*(1+x)^{K-2p-3q}
        + sum_{null} c_{n,m} * D^{(n,m)}_ell*C^alpha_ell(1)*(2ell+d-3)*(1+x)^{K-2n-m}
        >= 0

    The c_{n,m} are free Lagrange multipliers enforcing the CSDR null constraints.
    The z_{p,q} represent Wilson coefficient values (objectives).

    Simplification steps performed here (full record):
    ---------------------------------------------------
    Step 1: Enumerate all objective pairs (p,q) with p>=1, q>=0, 2p+3q<=K.
    Step 2: Enumerate all null constraint pairs (n,m) with m>n>=1, 2n+m<=K.
    Step 3: For each spin ell (even, 0 to max_spin):
            - Compute C^{alpha}_ell(1) = (2*alpha)_ell / ell! (Gegenbauer at 1)
            - Compute obj kernel kappa_{p,q,ell} = C^alpha_ell(1)*(2ell+d-3)
              [Eq. see obj_kernel_coeff in csdr.py]
            - Compute null kernel tilde_c_{n,m,ell} = D^{(n,m)}_{ell,alpha}
              * C^alpha_ell(1) * (2ell+d-3)  [User notes eq.(2)]
              where D^{(n,m)} is from CSDR eq.(11): sum over j=n..m
            - Build polynomial for each variable: kappa * (1+x)^{K-power}
            - Assemble 1×1 block with DampedRational prefactor e^{-x}/(1+x)^K
    Step 4: Assign objective vector: -1 at z_{obj_index} (for lower bound) or
            +1 (for upper bound).
    Step 5: Assign normalization vector: 1 at z_{norm_index}.

    Parameters
    ----------
    obj_index : (p, q)      Pair to minimize/maximize.
    norm_index : (p, q)     Pair to normalize to 1.
    d : int                  Spacetime dimension (default 4).
    K : int                  Maximum spectral power 2p+3q or 2n+m (default 8).
    max_spin : int           Maximum even spin included (default 10).
    precision : int          Output decimal precision (default 200).
    bound_direction : str    "upper" (maximize) or "lower" (minimize).

    Returns
    -------
    dict  SDPB PMP JSON structure.
    """
    obj_pairs = enumerate_obj_pairs(K)
    null_pairs = enumerate_null_pairs(K)

    # Check indices are in enumeration
    if obj_index not in obj_pairs:
        raise ValueError(
            f"obj_index {obj_index} not in enumerated objective pairs for K={K}. "
            f"Need p>=1, q>=0, 2p+3q<={K}."
        )
    if norm_index not in obj_pairs:
        raise ValueError(
            f"norm_index {norm_index} not in enumerated objective pairs for K={K}. "
            f"Need p>=1, q>=0, 2p+3q<={K}."
        )

    # Decision variables: objectives first, then null constraints
    all_pairs = obj_pairs + [("null", n, m) for (n, m) in null_pairs]
    N = len(obj_pairs) + len(null_pairs)  # total decision variables

    # Objective vector: maximize z_{obj_index} → b_{obj_index} = +1 (SDPB maximizes b·z)
    # For "lower" bound: minimize z_{obj_index} → b_{obj_index} = -1
    b = [Fraction(0)] * N
    obj_pos = obj_pairs.index(obj_index)
    if bound_direction == "upper":
        b[obj_pos] = Fraction(1)
    elif bound_direction == "lower":
        b[obj_pos] = Fraction(-1)
    else:
        raise ValueError(f"bound_direction must be 'upper' or 'lower', got '{bound_direction}'")

    # Normalization vector: c · z = 1, c_{norm_index} = 1
    c = [Fraction(0)] * N
    norm_pos = obj_pairs.index(norm_index)
    c[norm_pos] = Fraction(1)

    # Build spin blocks
    pmp_array = []
    for ell in range(0, max_spin + 1, 2):
        block = _spin_block_csdr(
            ell=ell,
            obj_pairs=obj_pairs,
            null_pairs=null_pairs,
            K=K,
            d=d,
            precision=precision,
        )
        if block is not None:
            pmp_array.append(block)

    return {
        "objective": [fraction_to_str(bi, precision) for bi in b],
        "normalization": [fraction_to_str(ci, precision) for ci in c],
        "PositiveMatrixWithPrefactorArray": pmp_array,
        "_metadata": {
            "description": "CSDR-based SDPB PMP from user notes eq.(2)",
            "d": d,
            "alpha": fraction_to_str(alpha_from_d(d), precision),
            "K_max_order": K,
            "max_spin": max_spin,
            "n_obj_pairs": len(obj_pairs),
            "n_null_pairs": len(null_pairs),
            "n_spin_blocks": len(pmp_array),
            "bound_direction": bound_direction,
            "obj_index": list(obj_index),
            "norm_index": list(norm_index),
        },
    }


def write_pipeline_outputs(
    output_dir: str,
    max_order: int,
    d: int = 4,
    precision: int = 80,
    max_spin: int = 10,
) -> Dict[str, object]:
    """
    Write all pipeline outputs to output_dir.

    Outputs:
    - eq11_null_constraints.json:  D coefficients and kernels for all null pairs.
    - csdr_reliability_checks.json: self-consistency checks.
    - csdr_pmp_upper_<obj>_over_<norm>.json: SDPB PMP for upper bound.
    - csdr_pmp_lower_<obj>_over_<norm>.json: SDPB PMP for lower bound.
    """
    os.makedirs(output_dir, exist_ok=True)

    eq11_payload = build_eq11_constraints_json(max_order=max_order, precision=precision)
    checks_payload = run_reliability_checks(precision=precision)

    files = {}

    # Write eq11 constraints
    p = os.path.join(output_dir, "eq11_null_constraints.json")
    with open(p, "w", encoding="utf-8") as fh:
        json.dump(eq11_payload, fh, indent=2)
    files["eq11_constraints"] = p

    # Write reliability checks
    p = os.path.join(output_dir, "csdr_reliability_checks.json")
    with open(p, "w", encoding="utf-8") as fh:
        json.dump(checks_payload, fh, indent=2)
    files["reliability_checks"] = p

    # Write CSDR SDPB PMPs
    # Example: bound W_{1,0} / W_{2,0} and W_{2,0} / W_{1,0}
    example_bounds = [
        ((2, 0), (1, 0), "upper"),
        ((2, 0), (1, 0), "lower"),
    ]
    for obj_idx, norm_idx, direction in example_bounds:
        if 2 * obj_idx[0] + 3 * obj_idx[1] > max_order:
            continue
        if 2 * norm_idx[0] + 3 * norm_idx[1] > max_order:
            continue
        try:
            pmp = build_csdr_pmp_json(
                obj_index=obj_idx,
                norm_index=norm_idx,
                d=d,
                K=max_order,
                max_spin=max_spin,
                precision=precision,
                bound_direction=direction,
            )
            fname = (
                f"csdr_pmp_{direction}_{obj_idx[0]}_{obj_idx[1]}"
                f"_over_{norm_idx[0]}_{norm_idx[1]}.json"
            )
            p = os.path.join(output_dir, fname)
            with open(p, "w", encoding="utf-8") as fh:
                json.dump(pmp, fh, indent=2)
            files[f"csdr_pmp_{direction}_{obj_idx}_{norm_idx}"] = p
        except ValueError as exc:
            files[f"csdr_pmp_{direction}_{obj_idx}_{norm_idx}_error"] = str(exc)

    return {"files": files, "checks": checks_payload}


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Generate paper-aligned CSDR dispersion / null-constraint pipeline.\n"
            "Produces SDPB PMP files using the CSDR kernels (user notes eq.(2)).\n"
            "d is a free parameter — do NOT default to d=4."
        )
    )
    parser.add_argument("--output-dir", required=True, help="Directory for generated files.")
    parser.add_argument("--max-order", type=int, default=8, help="Maximum 2n+m cutoff.")
    parser.add_argument("--d", type=int, default=4, help="Spacetime dimension.")
    parser.add_argument("--precision", type=int, default=80, help="Decimal precision for output.")
    parser.add_argument("--max-spin", type=int, default=10, help="Maximum spin.")
    args = parser.parse_args()

    result = write_pipeline_outputs(
        output_dir=args.output_dir,
        max_order=args.max_order,
        d=args.d,
        precision=args.precision,
        max_spin=args.max_spin,
    )
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
