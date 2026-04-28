"""
pmp_generator.py — Generate PMP JSON input files for SDPB.

This module translates the EFT dual optimization problem — formulated using
the CSDR (Crossing-Symmetric Dispersion Relations) scheme — into SDPB's
Polynomial Matrix Program (PMP) JSON format.

Physics setup (following user notes eq.(2))
-------------------------------------------
The heavy average of a kernel F(s1, ell) is:

  <F> = sum_{ell even} n^{(d)}_ell int_{delta0}^inf ds1/s1 * s1^{4-d}/pi
        * rho_ell(s1) * F(s1,ell)

Wilson coefficients W_{p,q} (first index p = n-m, second q = m):

  Objectives  (p >= 1, q >= 0): W_{p,q} = < C^alpha_ell(1)*(2ell+d-3)/s1^{2p+3q} >
  Null constr (p = n-m < 0):    W_{p,q} = < D^{(n,m)}_ell*C^alpha_ell(1)*(2ell+d-3)/s1^{2n+m} > = 0

  where alpha = (d-3)/2.

SDPB PMP structure
------------------
Decision variables z = (z_{p,q} for obj, c_{n,m} for null constraints).

For each even spin ell in {0,2,...,ell_max}, one 1×1 polynomial block:

  P^{(ell)}(x) = sum_{obj}  z_{p,q}  * obj_kernel(p,q,ell) * (1+x)^{K-2p-3q}
               + sum_{null} c_{n,m}  * null_kernel(n,m,ell) * (1+x)^{K-2n-m}

  where K = max spectral power, and the DampedRational prefactor is e^{-x}/(1+x)^K.

The SDPB problem maximizes b·z subject to c·z = 1 and P^{(ell)}(x) >= 0 for all x >= 0.

Key properties:
  * obj_kernel >= 0 for all ell, d => objectives always have positive kernels.
  * null_kernel can be negative (D can be negative) => null constraints are
    non-trivially enforced by the Lagrange multipliers c_{n,m}.
  * Different from Extremal EFT (fixed-t dispersion) spectral functions.
    The two subtraction schemes give different bounds.

IMPORTANT: d is a free parameter. Do NOT hardcode d=4.
"""

import json
from fractions import Fraction
from typing import Any, Dict, List, Optional, Tuple

from .csdr import (
    alpha_from_d,
    enumerate_null_pairs,
    enumerate_obj_pairs_nm,
    null_kernel_coeff,
    obj_kernel_coeff,
    poly_expand_1px,
    poly_pad,
    poly_scale,
    spectral_power_nm,
)
from .physics import fraction_to_str


# ---------------------------------------------------------------------------
# e^{-1} helper
# ---------------------------------------------------------------------------

def _compute_e_inverse(precision: int) -> str:
    """Compute 1/e to precision decimal digits."""
    from decimal import Decimal, getcontext
    getcontext().prec = precision + 20
    e_val = Decimal(1)
    fac = Decimal(1)
    for k in range(1, precision + 100):
        fac *= k
        term = Decimal(1) / fac
        e_val += term
        if term < Decimal(10) ** (-(precision + 10)):
            break
    return format(Decimal(1) / e_val, f".{precision}f").rstrip("0").rstrip(".")


# ---------------------------------------------------------------------------
# Single spin block
# ---------------------------------------------------------------------------

def _csdr_spin_block(
    ell: int,
    obj_pairs: List[Tuple[int, int]],
    null_pairs: List[Tuple[int, int]],
    K: int,
    d: int,
    precision: int,
) -> Optional[Dict[str, Any]]:
    """
    Build one SDPB PositiveMatrixWithPrefactor block for spin ell.

    Simplification record for this block
    -------------------------------------
    1. Compute alpha = (d-3)/2.
    2. For each objective (n,m) with n>m>=0:
         kappa = D^{(n,m)}_alpha * C^alpha_ell(1) * (2ell+d-3)
                 [obj_kernel_coeff, user's objective D formula, ell-independent]
         poly  = kappa * (1+x)^{K-(2n+m)}
    3. For each null constraint (n,m) with m>n>=1:
         kappa = D^{(n,m)}_{ell,alpha} * C^alpha_ell(1) * (2ell+d-3)
                 [null_kernel_coeff, CSDR eq.(11), ell-dependent]
         poly  = kappa * (1+x)^{K-(2n+m)}
    4. Prefactor: e^{-x} / (1+x)^K  (DampedRational, pole at x=-1, multiplicity K).
    5. Assemble 1×1 polynomial block.
    """
    poly_vectors: List[List[Fraction]] = []

    for (n, m) in obj_pairs:
        kappa = obj_kernel_coeff(n, m, ell, d)
        exp = K - spectral_power_nm(n, m)
        poly = poly_pad(poly_scale(poly_expand_1px(exp), kappa), K + 1)
        poly_vectors.append(poly)

    for (n, m) in null_pairs:
        kappa = null_kernel_coeff(n, m, ell, d)
        exp = K - (2 * n + m)
        poly = poly_pad(poly_scale(poly_expand_1px(exp), kappa), K + 1)
        poly_vectors.append(poly)

    if all(all(c == 0 for c in pv) for pv in poly_vectors):
        return None

    e_inv = _compute_e_inverse(precision)
    # SDPB PMP format uses "DampedRational" key (not "prefactor").
    # See test/data/end-to-end_tests/1d/input/pmp.json for the canonical example.
    damped_rational = {
        "base": e_inv,
        "constant": "1",
        "poles": ["-1"] * K,
    }

    poly_json = [[
        [fraction_to_str(c, precision) for c in pv]
        for pv in poly_vectors
    ]]

    return {
        "DampedRational": damped_rational,
        "polynomials": [poly_json],
    }


# ---------------------------------------------------------------------------
# Main PMP generators
# ---------------------------------------------------------------------------

def generate_csdr_pmp(
    obj_index: Tuple[int, int],
    norm_index: Tuple[int, int],
    d: int,
    K: int = 8,
    max_spin: int = 10,
    precision: int = 200,
    bound_direction: str = "upper",
) -> Dict[str, Any]:
    """
    Generate SDPB PMP JSON for bounding W_{obj_index} / W_{norm_index}
    using CSDR dispersion relations (user notes eq.(2)).

    Uses unified (n,m) notation for both objectives and null constraints:
      - Objectives: n > m >= 0, spectral power 2n+m, D from D_coeff_obj
      - Null:       m > n >= 1, spectral power 2n+m, D from D_coeff (CSDR eq.11)

    Simplification record
    ---------------------
    1. Enumerate objective pairs (n>m>=0, 2n+m<=K) via csdr.enumerate_obj_pairs_nm.
    2. Enumerate null constraint pairs (m>n>=1, 2n+m<=K) via csdr.enumerate_null_pairs.
    3. Build decision variable vector z = [z_{obj pairs...}, c_{null pairs...}].
       - Positions 0..N_obj-1: objective Wilson coefficients (indexed by n,m).
       - Positions N_obj..N_obj+N_null-1: null constraint Lagrange multipliers.
    4. Objective vector b: b[i] = +1 (upper) or -1 (lower) at obj_index position.
    5. Normalization vector c: c[i] = 1 at norm_index position.
    6. For each even spin ell in {0,2,...,max_spin}: build spin block via _csdr_spin_block.
    7. Output PMP JSON with metadata.

    Parameters
    ----------
    obj_index : (n, m)      Wilson coefficient W_{n-m,m} to bound (n>m>=0).
    norm_index : (n, m)     Wilson coefficient to normalize to 1 (n>m>=0).
    d : int                 Spacetime dimension (required, no default).
    K : int                 Maximum spectral power 2n+m (default 8).
    max_spin : int          Maximum even spin (default 10).
    precision : int         Output decimal precision (default 200).
    bound_direction : str   "upper" (maximize) or "lower" (minimize).

    Returns
    -------
    dict  SDPB PMP JSON.
    """
    obj_pairs = enumerate_obj_pairs_nm(K)
    null_pairs = enumerate_null_pairs(K)

    if obj_index not in obj_pairs:
        raise ValueError(
            f"obj_index {obj_index} not in objective pairs for K={K} "
            f"(need p>=1, q>=0, 2p+3q<={K})."
        )
    if norm_index not in obj_pairs:
        raise ValueError(
            f"norm_index {norm_index} not in objective pairs for K={K} "
            f"(need p>=1, q>=0, 2p+3q<={K})."
        )

    N_obj = len(obj_pairs)
    N_null = len(null_pairs)
    N = N_obj + N_null

    # Objective vector
    b: List[Fraction] = [Fraction(0)] * N
    obj_pos = obj_pairs.index(obj_index)
    b[obj_pos] = Fraction(1) if bound_direction == "upper" else Fraction(-1)

    # Normalization vector
    c: List[Fraction] = [Fraction(0)] * N
    norm_pos = obj_pairs.index(norm_index)
    c[norm_pos] = Fraction(1)

    # Spin blocks
    pmp_array = []
    for ell in range(0, max_spin + 1, 2):
        block = _csdr_spin_block(
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
            "description": "CSDR-based SDPB PMP (user notes eq.(2)); NOT Extremal EFT.",
            "d": d,
            "alpha": fraction_to_str(alpha_from_d(d), precision),
            "K_max_order": K,
            "max_spin": max_spin,
            "n_obj_pairs": N_obj,
            "n_null_pairs": N_null,
            "n_decision_variables": N,
            "n_spin_blocks": len(pmp_array),
            "bound_direction": bound_direction,
            "obj_index": list(obj_index),
            "norm_index": list(norm_index),
            "note": (
                "S1,S2,S3 are Mandelstam minus mu/3 so S1+S2+S3=0. "
                "The W_{p,q} basis uses x=-(S1*S2+S2*S3+S3*S1), y=-S1*S2*S3. "
                "d is a free parameter, not fixed to 4."
            ),
        },
    }


def write_pmp_json(pmp: Dict[str, Any], filepath: str) -> None:
    """Write PMP JSON to file."""
    with open(filepath, "w") as fh:
        json.dump(pmp, fh, indent=2)


def generate_pmp_for_ratio_bound(
    numerator_index: Tuple[int, int],
    denominator_index: Tuple[int, int],
    max_spin: int = 10,
    max_order: int = 6,
    d: int = 4,
    precision: int = 200,
    bound_direction: str = "lower",
    use_null_constraints: bool = True,
    **_kwargs: Any,
) -> Dict[str, Any]:
    """
    Generate PMP for bounding W_{numerator_index} / W_{denominator_index}.

    Both indices use (n,m) notation where W_{n-m,m} with n>m>=0.
    Uses CSDR kernels (user notes eq.(2)) with objective D from user's formula.

    Parameters
    ----------
    numerator_index : (n, m)    Wilson coefficient in numerator (n>m>=0).
    denominator_index : (n, m)  Wilson coefficient in denominator (n>m>=0).
    max_spin : int               Maximum even spin.
    max_order : int              Maximum spectral power K = 2n+m.
    d : int                      Spacetime dimension.
    precision : int              Output precision.
    bound_direction : str        "upper" or "lower".
    use_null_constraints : bool  If False, skips null pairs (testing only).
    **_kwargs : ignored         For backward compatibility.

    Returns
    -------
    dict  PMP JSON.
    """
    return generate_csdr_pmp(
        obj_index=numerator_index,
        norm_index=denominator_index,
        d=d,
        K=max_order,
        max_spin=max_spin,
        precision=precision,
        bound_direction=bound_direction,
    )


def generate_pmp_json(
    objective_index: Tuple[int, int],
    max_spin: int = 10,
    max_order: int = 6,
    d: int = 4,
    precision: int = 200,
    **_kwargs: Any,
) -> Dict[str, Any]:
    """
    Generate a PMP JSON for bounding a specific objective Wilson coefficient.

    objective_index uses (n,m) notation: W_{n-m,m} with n>m>=0.
    Normalizes the first available pair != objective_index.

    Parameters
    ----------
    objective_index : (n, m)  Wilson coefficient to minimize (n>m>=0).
    max_spin, max_order, d, precision: as in generate_csdr_pmp.
    **_kwargs : ignored for backward compatibility.

    Returns
    -------
    dict  PMP JSON.
    """
    obj_pairs = enumerate_obj_pairs_nm(max_order)
    if objective_index not in obj_pairs:
        raise ValueError(
            f"objective_index {objective_index} not found in pairs for K={max_order}."
        )
    norm_index = next(
        (p for p in obj_pairs if p != objective_index),
        None,
    )
    if norm_index is None:
        raise ValueError("Not enough objective pairs to form normalization.")
    return generate_csdr_pmp(
        obj_index=objective_index,
        norm_index=norm_index,
        d=d,
        K=max_order,
        max_spin=max_spin,
        precision=precision,
        bound_direction="lower",
    )


# ---------------------------------------------------------------------------
# Backward-compatibility shim for crossing_pipeline.py
# ---------------------------------------------------------------------------
# The crossing_pipeline uses fixed-t Extremal EFT spectral functions,
# which are a different calculation from the CSDR approach.  The shim
# below restores the original function so crossing_pipeline tests pass.

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
    Backward-compatibility wrapper used by crossing_pipeline.py.

    This function implements the OLD fixed-t Extremal EFT spectral function
    approach.  It is NOT the CSDR approach.  For new code use generate_csdr_pmp.

    For crossing_pipeline.py use with crossing_type="crossing_basis", this
    function builds SDPB blocks using fixed-t spectral functions and the
    variable-elimination approach for null constraints.
    """
    from collections import defaultdict
    from .physics import (
        fraction_to_str as _fts,
        gegenbauer_coefficients,
        _taylor_coeff_at_one,
        _expand_binomial_power,
        _poly_multiply,
    )
    from .null_constraints import (
        eliminate_variables,
        get_null_constraints,
        get_crossing_symmetric_null_constraints,
    )

    mu_threshold = 4 * m_sq

    # Enumerate all (p,q) indices
    all_indices = []
    for total in range(max_order + 1):
        for p in range(total + 1):
            all_indices.append((p, total - p))

    # Apply null constraints
    if use_null_constraints:
        if crossing_type == "full":
            constraints = get_crossing_symmetric_null_constraints(max_order, m_sq)
        elif crossing_type == "crossing_basis":
            from .crossing_pipeline import (
                derive_null_constraints as _derive_crossing_basis_nc,
            )
            constraints = _derive_crossing_basis_nc(max_order, m_sq)
        else:
            constraints = get_null_constraints(max_order, m_sq)
        free_indices, substitution = eliminate_variables(constraints, all_indices)
    else:
        free_indices = list(all_indices)
        substitution = {}

    z_indices = list(free_indices)

    def _reduce(terms):
        reduced = defaultdict(Fraction)
        for idx, coeff in terms.items():
            if idx in substitution:
                for fi, sc in substitution[idx].items():
                    reduced[fi] += coeff * sc
            else:
                reduced[idx] += coeff
        return [reduced.get(idx, Fraction(0)) for idx in z_indices]

    objective = _reduce(objective_terms)
    if bound_direction == "lower":
        objective = [-c for c in objective]
    normalization = _reduce(normalization_terms)
    if all(c == 0 for c in normalization):
        raise ValueError("Normalization vanishes after applying null constraints.")

    # Build spin blocks (old fixed-t approach)
    def _build_old_block(ell):
        all_orig = set(z_indices) | set(substitution.keys())
        for sub in substitution.values():
            all_orig |= set(sub.keys())
        all_orig = {(p, q) for p, q in all_orig if p + q <= max_order}
        if not all_orig:
            return None
        p_max = max(p for p, _ in all_orig)
        q_max = max(q for _, q in all_orig)
        poly_deg = q_max + p_max + 1
        geg = gegenbauer_coefficients(ell, d)
        spec = {}
        for p, q in all_orig:
            tq = _taylor_coeff_at_one(geg, q)
            coeff = tq * Fraction(2) ** q
            if coeff == 0:
                spec[(p, q)] = [Fraction(0)] * (poly_deg + 1)
                continue
            px = [Fraction(0)] * (q_max - q + 1)
            px[q_max - q] = Fraction(1)
            ps = _expand_binomial_power(mu_threshold, p_max - p)
            prod = _poly_multiply(px, ps)
            res = [c * coeff for c in prod]
            while len(res) < poly_deg + 1:
                res.append(Fraction(0))
            spec[(p, q)] = res[: poly_deg + 1]
        pvecs = []
        for k, zidx in enumerate(z_indices):
            ep = list(spec.get(zidx, [Fraction(0)] * (poly_deg + 1)))
            for eidx, sd in substitution.items():
                if zidx in sd and eidx in spec:
                    fac = sd[zidx]
                    for i in range(min(len(ep), len(spec[eidx]))):
                        ep[i] += fac * spec[eidx][i]
            pvecs.append(ep)
        if all(all(c == 0 for c in pv) for pv in pvecs):
            return None
        poles = (
            [fraction_to_str(Fraction(0), precision)] * q_max
            + [fraction_to_str(-mu_threshold, precision)] * (p_max + 1)
        )
        e_inv = _compute_e_inverse(precision)
        pf = {"base": e_inv, "constant": "1", "poles": poles}
        pj = [[
            [fraction_to_str(c, precision) for c in pv]
            for pv in pvecs
        ]]
        return {"DampedRational": pf, "polynomials": [pj]}

    pmp_array = []
    for ell in range(0, max_spin + 1, 2):
        blk = _build_old_block(ell)
        if blk is not None:
            pmp_array.append(blk)

    return {
        "objective": [fraction_to_str(v, precision) for v in objective],
        "normalization": [fraction_to_str(v, precision) for v in normalization],
        "PositiveMatrixWithPrefactorArray": pmp_array,
    }
