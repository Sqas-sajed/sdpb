"""
pmp_generator.py — Generate PMP JSON input files for SDPB.

This module translates the EFT dual optimization problem into SDPB's
Polynomial Matrix Program (PMP) JSON format.

The PMP problem solved by SDPB is:
  maximize  a · z  over  z ∈ R^{N+1}
  such that Σ_n z_n W_j^n(x) ≽ 0  for all x ≥ 0, j=1,...,J
  and       n · z = 1

With normalization n = (1, 0, ..., 0), this reduces to:
  maximize  b₀ + b · y  over  y ∈ R^N
  such that M_j^0(x) + Σ_{n=1}^N y_n M_j^n(x) ≽ 0  for all x ≥ 0

The physics mapping:
  - y_n ↔ Wilson coefficients (after eliminating dependent ones via null constraints)
  - j ↔ spin ℓ in the partial-wave expansion
  - x ↔ μ - μ_threshold (shifted mass squared)
  - M_j^n(x) ↔ spectral function v_ℓ^{(n)}(x)
"""

import json
import math
from collections import defaultdict
from fractions import Fraction
from typing import Any, Dict, List, Optional, Tuple

from .physics import (
    build_positivity_polynomials,
    fraction_to_str,
    gegenbauer_coefficients,
    _taylor_coeff_at_one,
    _binomial,
    _expand_binomial_power,
    _poly_multiply,
)
from .null_constraints import (
    eliminate_variables,
    get_null_constraints,
    get_crossing_symmetric_null_constraints,
)


def generate_pmp_json(
    objective_index: Tuple[int, int],
    max_spin: int = 10,
    max_order: int = 6,
    use_null_constraints: bool = True,
    crossing_type: str = "su",
    m_sq: Fraction = Fraction(1),
    d: int = 4,
    precision: int = 200,
) -> Dict[str, Any]:
    """
    Generate a PMP JSON structure for bounding a specific Wilson coefficient.

    Parameters
    ----------
    objective_index : (int, int)
        The (p, q) index of the Wilson coefficient to minimize (lower bound).
    max_spin : int
        Maximum spin ℓ to include in the partial-wave expansion.
    max_order : int
        Maximum total order p+q for Wilson coefficients.
    use_null_constraints : bool
        Whether to incorporate Sinha-Zahed null constraints.
    crossing_type : str
        "su" for s↔u crossing only, "full" for full S₃ crossing.
    m_sq : Fraction
        Mass squared of external scalars.
    d : int
        Spacetime dimension.
    precision : int
        Number of decimal digits for numerical output.

    Returns
    -------
    dict
        PMP JSON structure ready to be written to file.
    """
    mu_threshold = 4 * m_sq

    # 1. Enumerate all Wilson coefficient indices (p, q) up to max_order
    all_indices = []
    for total in range(max_order + 1):
        for p in range(total + 1):
            q = total - p
            all_indices.append((p, q))

    # 2. Apply null constraints if requested
    if use_null_constraints:
        if crossing_type == "full":
            constraints = get_crossing_symmetric_null_constraints(max_order, m_sq)
        else:
            constraints = get_null_constraints(max_order, m_sq)
        free_indices, substitution = eliminate_variables(constraints, all_indices)
    else:
        free_indices = list(all_indices)
        substitution = {}

    # 3. Normalization: fix one Wilson coefficient (typically g_{0,0} or g_{2,0}).
    #    We use the convention that the normalization index has value 1.
    #    Choose (0, 0) as the normalization index (first coefficient = 1).
    norm_index = (0, 0)

    # If norm_index was eliminated, find an alternative
    if norm_index not in free_indices:
        # Try (1, 0), (2, 0), etc.
        for candidate in [(1, 0), (2, 0), (0, 1)]:
            if candidate in free_indices:
                norm_index = candidate
                break

    # Build the z-vector mapping:
    #   z_0 = coefficient for normalization (fixed by n·z = 1)
    #   z_1, ..., z_N = remaining free Wilson coefficients
    #
    # With normalization n = (1, 0, ..., 0), z_0 is determined as:
    #   z_0 = (1 - Σ_{n≥1} n_n z_n) / n_0 = 1 (since n_n = 0 for n ≥ 1)
    #
    # Actually, in SDPB's formulation with normalization, we set:
    #   normalization[k] = 1 for the norm_index position, 0 elsewhere
    #   This forces W_{norm_index} = 1.

    # Map free indices to z-vector positions
    # Position 0 in z-vector: norm_index (fixed to 1 by normalization)
    # Positions 1..N: other free indices
    z_indices = [norm_index] + [idx for idx in free_indices if idx != norm_index]
    N = len(z_indices) - 1  # Number of free decision variables

    # 4. Build objective vector
    # We want to minimize W_{objective_index}, which means
    # maximize -W_{objective_index}.
    #
    # If objective_index is free:
    #   objective[k] = -1 if z_indices[k] == objective_index, else 0
    # If objective_index was eliminated by null constraints:
    #   W_{obj} = Σ coeff_i × W_{free_i}
    #   objective[k] = -coeff_i where z_indices[k] = free_i

    objective = [Fraction(0)] * (N + 1)

    if objective_index in substitution:
        # Eliminated variable: express in terms of free variables
        for free_idx, coeff in substitution[objective_index].items():
            if free_idx in z_indices:
                k = z_indices.index(free_idx)
                objective[k] = -coeff
    elif objective_index in z_indices:
        k = z_indices.index(objective_index)
        objective[k] = Fraction(-1)
    else:
        raise ValueError(
            f"Objective index {objective_index} not found in free or eliminated indices"
        )

    # 5. Build normalization vector
    normalization = [Fraction(0)] * (N + 1)
    normalization[0] = Fraction(1)  # norm_index is at position 0

    # 6. Build PositiveMatrixWithPrefactorArray
    # For each even spin ℓ = 0, 2, 4, ..., max_spin:
    #   Create a 1×1 matrix polynomial constraint
    #   The constraint is: Σ_n z_n v_ℓ^{(n)}(x) ≥ 0 for all x ≥ 0
    #   where v_ℓ^{(n)}(x) encodes the spectral function.

    pmp_array = []

    for ell in range(0, max_spin + 1, 2):
        block = _build_spin_block(
            ell=ell,
            z_indices=z_indices,
            substitution=substitution,
            m_sq=m_sq,
            d=d,
            max_order=max_order,
            precision=precision,
        )
        if block is not None:
            pmp_array.append(block)

    # 7. If using Method B for remaining null constraints (as cross-check),
    #    add 1×1 constant blocks for each constraint.
    #    (This is optional; Method A already eliminates them.)

    # 8. Assemble the PMP JSON
    pmp = {
        "objective": [fraction_to_str(a, precision) for a in objective],
        "normalization": [fraction_to_str(n, precision) for n in normalization],
        "PositiveMatrixWithPrefactorArray": pmp_array,
    }

    return pmp


def _enumerate_indices(max_order: int) -> List[Tuple[int, int]]:
    indices = []
    for total in range(max_order + 1):
        for p in range(total + 1):
            indices.append((p, total - p))
    return indices


def _get_constraints(
    max_order: int,
    m_sq: Fraction,
    crossing_type: str,
) -> List[Dict[Tuple[int, int], Fraction]]:
    if crossing_type == "full":
        return get_crossing_symmetric_null_constraints(max_order, m_sq)
    if crossing_type == "crossing_basis":
        from .crossing_pipeline import (
            derive_null_constraints as derive_crossing_basis_null_constraints,
        )
        return derive_crossing_basis_null_constraints(max_order, m_sq)
    return get_null_constraints(max_order, m_sq)


def _reduce_linear_functional(
    terms: Dict[Tuple[int, int], Fraction],
    z_indices: List[Tuple[int, int]],
    substitution: Dict[Tuple[int, int], Dict[Tuple[int, int], Fraction]],
) -> List[Fraction]:
    reduced: Dict[Tuple[int, int], Fraction] = defaultdict(Fraction)
    for index, coeff in terms.items():
        if index in substitution:
            for free_index, sub_coeff in substitution[index].items():
                reduced[free_index] += coeff * sub_coeff
        else:
            reduced[index] += coeff
    return [reduced.get(index, Fraction(0)) for index in z_indices]


def generate_pmp_for_linear_functional_bound(
    objective_terms: Dict[Tuple[int, int], Fraction],
    normalization_terms: Dict[Tuple[int, int], Fraction],
    max_spin: int = 10,
    max_order: int = 6,
    use_null_constraints: bool = True,
    crossing_type: str = "su",
    m_sq: Fraction = Fraction(1),
    d: int = 4,
    precision: int = 200,
    bound_direction: str = "lower",
) -> Dict[str, Any]:
    """
    Generate a PMP for a generic linear-functional lower or upper bound.

    The optimization variable is the ordinary EFT coefficient vector in the
    s,t polynomial basis. `objective_terms` and `normalization_terms` specify
    exact linear functionals on that basis.
    """
    all_indices = _enumerate_indices(max_order)
    if use_null_constraints:
        constraints = _get_constraints(max_order=max_order, m_sq=m_sq, crossing_type=crossing_type)
        free_indices, substitution = eliminate_variables(constraints, all_indices)
    else:
        free_indices = list(all_indices)
        substitution = {}

    z_indices = list(free_indices)
    objective = _reduce_linear_functional(objective_terms, z_indices, substitution)
    if bound_direction == "lower":
        objective = [-coeff for coeff in objective]
    elif bound_direction != "upper":
        raise ValueError(f"Unknown bound direction '{bound_direction}'.")
    normalization = _reduce_linear_functional(normalization_terms, z_indices, substitution)

    if all(coeff == 0 for coeff in normalization):
        raise ValueError("Normalization functional vanishes after applying constraints.")

    pmp_array = []
    for ell in range(0, max_spin + 1, 2):
        block = _build_spin_block(
            ell=ell,
            z_indices=z_indices,
            substitution=substitution,
            m_sq=m_sq,
            d=d,
            max_order=max_order,
            precision=precision,
        )
        if block is not None:
            pmp_array.append(block)

    return {
        "objective": [fraction_to_str(value, precision) for value in objective],
        "normalization": [fraction_to_str(value, precision) for value in normalization],
        "PositiveMatrixWithPrefactorArray": pmp_array,
    }


def _build_spin_block(
    ell: int,
    z_indices: List[Tuple[int, int]],
    substitution: Dict[Tuple[int, int], Dict[Tuple[int, int], Fraction]],
    m_sq: Fraction,
    d: int,
    max_order: int,
    precision: int,
) -> Optional[Dict[str, Any]]:
    """
    Build a PositiveMatrixWithPrefactor block for spin ℓ.

    The spectral positivity condition for spin ℓ is:
      Σ_n z_n v_ℓ^{(n)}(x) ≥ 0  for all x ≥ 0

    where v_ℓ^{(n)}(x) is the spectral function for the n-th Wilson coefficient.

    After variable elimination (null constraints), each z_n corresponds to
    a free Wilson coefficient, and v_ℓ^{(n)}(x) may be a linear combination
    of the original spectral functions.

    For the 1×1 case, the polynomial matrix is:
      polynomials = [[[v_ℓ^{(0)}(x), v_ℓ^{(1)}(x), ..., v_ℓ^{(N)}(x)]]]

    Parameters
    ----------
    ell : int
        Spin quantum number.
    z_indices : list of (int, int)
        Wilson coefficient indices for the z-vector.
    substitution : dict
        Eliminated variable substitutions.
    m_sq : Fraction
        Mass squared.
    d : int
        Spacetime dimension.
    max_order : int
        Maximum total order.
    precision : int
        Output precision.

    Returns
    -------
    dict or None
        PositiveMatrixWithPrefactor block, or None if trivially zero.
    """
    mu_threshold = 4 * m_sq
    N_plus_1 = len(z_indices)

    # Compute the maximum polynomial degree needed.
    # The spectral function v_ℓ^{(p,q)}(x) after clearing the common denominator
    # D(x) = x^{q_max} (x + 4m²)^{p_max + 1} becomes a polynomial of degree
    # (q_max - q) + (p_max - p).
    #
    # We need to track all (p,q) that contribute to each z_n,
    # including substituted variables.

    # Collect all original indices that contribute
    all_original_indices = set()
    for z_idx in z_indices:
        all_original_indices.add(z_idx)
    for elim_idx, sub_dict in substitution.items():
        all_original_indices.add(elim_idx)
        for free_idx in sub_dict:
            all_original_indices.add(free_idx)

    # Only keep indices within max_order
    all_original_indices = {
        (p, q) for p, q in all_original_indices if p + q <= max_order
    }

    if not all_original_indices:
        return None

    p_max = max(p for p, _ in all_original_indices)
    q_max = max(q for _, q in all_original_indices)

    # Maximum polynomial degree after clearing denominator
    poly_deg = q_max + p_max + 1

    # Compute spectral functions for all original indices
    # v_ℓ^{(p,q)}(x) × D(x) = coeff × x^{q_max-q} × (x+4m²)^{p_max-p}
    # where coeff = [Taylor_q of C_ℓ at 1] × 2^q

    geg_coeffs = gegenbauer_coefficients(ell, d)
    spectral_polys = {}  # (p,q) -> polynomial coefficients in x

    for p, q in all_original_indices:
        if q >= len(geg_coeffs) + 1:
            # Gegenbauer doesn't have enough terms; coefficient is 0
            spectral_polys[(p, q)] = [Fraction(0)] * (poly_deg + 1)
            continue

        taylor_q = _taylor_coeff_at_one(geg_coeffs, q)
        coeff = taylor_q * Fraction(2) ** q

        if coeff == 0:
            spectral_polys[(p, q)] = [Fraction(0)] * (poly_deg + 1)
            continue

        # Build x^{q_max - q}
        x_power = q_max - q
        poly_x = [Fraction(0)] * (x_power + 1)
        poly_x[x_power] = Fraction(1)

        # Build (x + 4m²)^{p_max - p}
        shift_power = p_max - p
        poly_shift = _expand_binomial_power(mu_threshold, shift_power)

        # Multiply
        product = _poly_multiply(poly_x, poly_shift)

        # Scale
        result = [c * coeff for c in product]

        # Pad to poly_deg + 1
        while len(result) < poly_deg + 1:
            result.append(Fraction(0))
        result = result[: poly_deg + 1]

        spectral_polys[(p, q)] = result

    # Now build the polynomial vector for each z-component.
    # For free variable z_k corresponding to index z_indices[k]:
    #   The effective spectral function is v_ℓ^{(z_k)}(x) plus
    #   contributions from eliminated variables that depend on z_k.
    #
    #   v_eff_ℓ^{(k)}(x) = v_ℓ^{(z_k)}(x) + Σ_{elim} sub[elim][z_k] × v_ℓ^{(elim)}(x)

    poly_vectors = []
    for k, z_idx in enumerate(z_indices):
        # Start with direct contribution
        if z_idx in spectral_polys:
            eff_poly = list(spectral_polys[z_idx])
        else:
            eff_poly = [Fraction(0)] * (poly_deg + 1)

        # Add contributions from eliminated variables
        for elim_idx, sub_dict in substitution.items():
            if z_idx in sub_dict and elim_idx in spectral_polys:
                factor = sub_dict[z_idx]
                for i in range(min(len(eff_poly), len(spectral_polys[elim_idx]))):
                    eff_poly[i] += factor * spectral_polys[elim_idx][i]

        poly_vectors.append(eff_poly)

    # Check if all polynomials are identically zero (skip this block)
    all_zero = all(
        all(c == 0 for c in pv)
        for pv in poly_vectors
    )
    if all_zero:
        return None

    # Build the block JSON
    # DampedRational prefactor: we use 1/D(x) as the prefactor
    # D(x) = x^{q_max} (x + 4m²)^{p_max + 1}
    # This has poles at x = 0 (mult q_max) and x = -4m² (mult p_max+1)
    # and the exponential damping e^{-x} for convergence.

    poles = []
    for _ in range(q_max):
        poles.append(fraction_to_str(Fraction(0), precision))
    for _ in range(p_max + 1):
        poles.append(fraction_to_str(-mu_threshold, precision))

    # Base for damped rational: e^{-1} ≈ 0.36788...
    e_inv = Fraction(1, 3)  # Approximate; for high precision use mpfr
    # Better: compute e^{-1} to high precision
    e_inv_str = _compute_e_inverse(precision)

    prefactor = {
        "base": e_inv_str,
        "constant": "1",
        "poles": poles,
    }

    # Format polynomials for JSON
    # Structure: polynomials[row][col][n] = polynomial coefficients
    # For 1×1 matrix: polynomials[0][0][n] = coefficients for z_n

    poly_json = [[
        [fraction_to_str(c, precision) for c in pv]
        for pv in poly_vectors
    ]]

    block = {
        "prefactor": prefactor,
        "polynomials": [poly_json],
    }

    return block


def _compute_e_inverse(precision: int) -> str:
    """Compute 1/e = e^{-1} to the given number of decimal digits."""
    from decimal import Decimal, getcontext
    getcontext().prec = precision + 20

    # e^{-1} via Taylor series: e^{-1} = Σ (-1)^n / n!
    result = Decimal(0)
    term = Decimal(1)
    for n in range(1, precision + 100):
        term = term / n
        if n % 2 == 0:
            result += term
        else:
            result -= term
        if abs(term) < Decimal(10) ** (-(precision + 10)):
            break
    result += Decimal(1)  # n=0 term

    # Actually compute properly: e = Σ 1/n!, then 1/e
    e_val = Decimal(1)
    factorial = Decimal(1)
    for n in range(1, precision + 100):
        factorial *= n
        term = Decimal(1) / factorial
        e_val += term
        if term < Decimal(10) ** (-(precision + 10)):
            break

    e_inv = Decimal(1) / e_val
    return format(e_inv, f'.{precision}f').rstrip('0').rstrip('.')


def write_pmp_json(pmp: Dict[str, Any], filepath: str) -> None:
    """Write PMP JSON to file."""
    with open(filepath, 'w') as f:
        json.dump(pmp, f, indent=2)


def generate_pmp_for_ratio_bound(
    numerator_index: Tuple[int, int],
    denominator_index: Tuple[int, int],
    max_spin: int = 10,
    max_order: int = 6,
    use_null_constraints: bool = True,
    crossing_type: str = "su",
    m_sq: Fraction = Fraction(1),
    d: int = 4,
    precision: int = 200,
    bound_direction: str = "lower",
) -> Dict[str, Any]:
    """
    Generate PMP for bounding the ratio W_{num} / W_{den}.

    This normalizes W_{den} = 1 and then minimizes W_{num}.

    Parameters
    ----------
    numerator_index : (int, int)
        Index of the Wilson coefficient in the numerator.
    denominator_index : (int, int)
        Index of the Wilson coefficient in the denominator (normalized to 1).
    ... (other parameters same as generate_pmp_json)

    Returns
    -------
    dict
        PMP JSON structure.
    """
    return generate_pmp_for_linear_functional_bound(
        objective_terms={numerator_index: Fraction(1)},
        normalization_terms={denominator_index: Fraction(1)},
        max_spin=max_spin,
        max_order=max_order,
        use_null_constraints=use_null_constraints,
        crossing_type=crossing_type,
        m_sq=m_sq,
        d=d,
        precision=precision,
        bound_direction=bound_direction,
    )
