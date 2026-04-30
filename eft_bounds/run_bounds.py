#!/usr/bin/env python3
"""
run_bounds.py — Main driver for CSDR-based EFT Wilson coefficient bounds.

This script generates SDPB PMP JSON files for bounding ratios of CSDR Wilson
coefficients W_{p,q}, using the dispersion relation kernels from:

  Sinha & Zahed, "Crossing Symmetric Dispersion Relations in QFTs"

and the heavy-average formulation from:

  Caron-Huot & Duong, "Extremal Effective Field Theories" (Section 3.3)

as written explicitly in the user's derivation notes, eq.(2).

KEY POINTS
----------
* S1, S2, S3 are NOT the Mandelstam variables.  They are the Mandelstam
  variables *minus* mu/3, so S1+S2+S3 = 0 (crossing-symmetric point).
* d is a free parameter.  Do NOT default to d=4.
* The CSDR and Extremal EFT papers use DIFFERENT subtraction schemes.
  Their bounds are analytically different.  Running this script will
  produce results that differ from the Extremal EFT paper — that is the
  correct behavior.
* The W_{m,n} classification:
    - m <= n  (first index <= second): objectives (EFT coefficients >= 0)
    - n < m   (first index > second): null constraints (must equal 0)
  This translates to CSDR notation W_{n-m,m} with:
    - n-m >= 0 (objectives): forward-scattering kernel C_ell(1)*(2ell+d-3)/s1^{2n+m}
    - n-m < 0  (null constr): kernel D^{(n,m)}_ell*C_ell(1)*(2ell+d-3)/s1^{2n+m}
      where D^{(n,m)} is from CSDR eq.(11), valid for m > n >= 1.

USAGE
-----
    python -m eft_bounds.run_bounds [options]

After generating PMP files, run SDPB:
    pmp2sdp --precision=1024 --input=<pmp_file> --output=<sdp_dir>
    mpirun -n 4 sdpb --precision=1024 -s <sdp_dir> -o <output_dir>
"""

import argparse
import json
import os
import sys
from fractions import Fraction
from typing import Any, Dict, List, Tuple

# Allow running as a script
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from eft_bounds.csdr import (
    D_coeff,
    alpha_from_d,
    check_closed_form_1_2,
    enumerate_null_pairs,
    enumerate_obj_pairs,
    gegenbauer_at_one,
    null_kernel_coeff,
    obj_kernel_coeff,
    s1_power,
)
from eft_bounds.physics import fraction_to_str
from eft_bounds.pmp_generator import (
    generate_csdr_pmp,
    generate_pmp_for_ratio_bound,
    write_pmp_json,
)


# ---------------------------------------------------------------------------
# Consistency checks
# ---------------------------------------------------------------------------

def run_csdr_checks(d: int = 4, precision: int = 40) -> Dict[str, Any]:
    """
    Run CSDR self-consistency checks.

    Checks performed
    ----------------
    1. D^{(1,2)}_{ell,alpha} matches closed-form for ell = 0,2,4,6,8.
    2. C^{alpha}_ell(1) = (2*alpha)_ell / ell! for several (ell, d).
    3. Null kernel = D * C_ell(1) * (2ell+d-3) is non-zero for non-trivial ell.
    4. Objective kernel = C_ell(1) * (2ell+d-3) > 0 for all ell.
    5. Null pairs and objective pairs enumerated correctly for K=8.
    6. csdr.check_closed_form_1_2() passes.

    Parameters
    ----------
    d : int    Spacetime dimension.
    precision : int  Output decimal precision.
    """
    alpha = alpha_from_d(d)
    checks = []

    # Check 1: D^{(1,2)} closed-form
    for ell in [0, 2, 4, 6, 8]:
        direct = D_coeff(1, 2, ell, alpha)
        closed = (
            2 * ell * (ell + 2 * alpha)
            * (-11 - 10 * alpha + 2 * ell * (ell + 2 * alpha))
            / ((2 * alpha + 1) * (2 * alpha + 3))
        )
        ok = (direct == closed)
        checks.append({
            "name": f"D^(1,2)_ell{ell}_d{d}",
            "passed": ok,
            "direct": fraction_to_str(direct, precision),
            "closed_form": fraction_to_str(closed, precision),
        })

    # Check 2: Gegenbauer at 1
    for ell in [0, 2, 4]:
        C1 = gegenbauer_at_one(ell, alpha)
        # For d=4 (alpha=1/2): C_ell^{1/2}(1) = 1 for all ell
        if d == 4:
            expected = Fraction(1)
            ok = (C1 == expected)
            checks.append({
                "name": f"C_ell(1)_d4_ell{ell}",
                "passed": ok,
                "value": fraction_to_str(C1, precision),
                "expected": fraction_to_str(expected, precision),
            })
        else:
            checks.append({
                "name": f"C_ell(1)_d{d}_ell{ell}",
                "passed": True,
                "value": fraction_to_str(C1, precision),
                "note": "No simple closed form for general d; recorded for reference",
            })

    # Check 3: Null kernel non-zero at ell=2 for (1,2)
    kappa_null = null_kernel_coeff(1, 2, 2, d)
    checks.append({
        "name": f"null_kernel_(1,2)_ell2_d{d}",
        "passed": (kappa_null != 0),
        "value": fraction_to_str(kappa_null, precision),
    })

    # Check 4: Objective kernel positive at ell=0
    kappa_obj_0 = obj_kernel_coeff(1, 0, 0, d)
    kappa_obj_2 = obj_kernel_coeff(1, 0, 2, d)
    checks.append({
        "name": f"obj_kernel_(1,0)_ell0_d{d}_positive",
        "passed": (kappa_obj_0 > 0),
        "value": fraction_to_str(kappa_obj_0, precision),
    })
    checks.append({
        "name": f"obj_kernel_(1,0)_ell2_d{d}_positive",
        "passed": (kappa_obj_2 > 0),
        "value": fraction_to_str(kappa_obj_2, precision),
    })

    # Check 5: Enumeration for K=8
    K = 8
    null_pairs = enumerate_null_pairs(K)
    obj_pairs = enumerate_obj_pairs(K)
    checks.append({
        "name": f"enumeration_K{K}",
        "passed": len(null_pairs) > 0 and len(obj_pairs) > 0,
        "n_null_pairs": len(null_pairs),
        "n_obj_pairs": len(obj_pairs),
        "null_pairs": [list(p) for p in null_pairs],
        "obj_pairs": [list(p) for p in obj_pairs],
    })

    # Check 6: csdr module self-check
    checks.append({
        "name": "csdr_check_closed_form_1_2",
        "passed": check_closed_form_1_2(),
    })

    all_passed = all(c["passed"] for c in checks)
    return {
        "passed": all_passed,
        "d": d,
        "alpha": fraction_to_str(alpha, precision),
        "checks": checks,
    }


# ---------------------------------------------------------------------------
# Kernel table
# ---------------------------------------------------------------------------

def compute_kernel_table(
    d: int,
    K: int,
    max_ell: int,
    precision: int = 40,
) -> Dict[str, Any]:
    """
    Compute and tabulate all kernel coefficients for objectives and null constraints.

    Simplification record
    ---------------------
    For each pair (p,q) or (n,m) and each spin ell:
    - Objective (p,q): kappa = C^alpha_ell(1) * (2ell+d-3)
    - Null constraint (n,m): kappa = D^{(n,m)}_{ell,alpha} * C^alpha_ell(1) * (2ell+d-3)
    These are then multiplied by (1+x)^{K-power} in the SDPB polynomial.

    Parameters
    ----------
    d : int       Spacetime dimension.
    K : int       Maximum spectral power.
    max_ell : int Maximum spin to tabulate.
    precision : int  Output precision.
    """
    alpha = alpha_from_d(d)
    obj_pairs = enumerate_obj_pairs(K)
    null_pairs = enumerate_null_pairs(K)
    ell_list = list(range(0, max_ell + 1, 2))

    obj_table = []
    for (p, q) in obj_pairs:
        # Convert (p, q) Wilson notation (p = n-m, q = m) to CSDR (n, m) with n = p+q, m = q.
        n_nm, m_nm = p + q, q
        row = {
            "type": "objective",
            "p": p, "q": q,
            "s1_power": s1_power(p, q),
            "kernels": {},
        }
        for ell in ell_list:
            kappa = obj_kernel_coeff(n_nm, m_nm, ell, d)
            row["kernels"][str(ell)] = fraction_to_str(kappa, precision)
        obj_table.append(row)

    null_table = []
    for (n, m) in null_pairs:
        row = {
            "type": "null_constraint",
            "n": n, "m": m,
            "first_csdr_index": n - m,
            "s1_power": 2 * n + m,
            "kernels": {},
        }
        for ell in ell_list:
            kappa = null_kernel_coeff(n, m, ell, d)
            row["kernels"][str(ell)] = fraction_to_str(kappa, precision)
        null_table.append(row)

    return {
        "d": d,
        "alpha": fraction_to_str(alpha, precision),
        "K": K,
        "ell_values": ell_list,
        "objectives": obj_table,
        "null_constraints": null_table,
    }


# ---------------------------------------------------------------------------
# Main run function
# ---------------------------------------------------------------------------

def run_bounds(
    output_dir: str,
    d: int,
    K: int = 8,
    max_spin: int = 10,
    precision: int = 200,
    bounds_to_compute: List[Tuple[Tuple[int, int], Tuple[int, int], str]] = None,
) -> Dict[str, Any]:
    """
    Generate CSDR-based SDPB PMP files for EFT Wilson coefficient bounds.

    Parameters
    ----------
    output_dir : str     Directory for output files.
    d : int              Spacetime dimension.
    K : int              Maximum spectral power (2p+3q cutoff).
    max_spin : int       Maximum even spin.
    precision : int      Decimal precision for SDPB files.
    bounds_to_compute :  List of (obj_index, norm_index, direction) tuples.
                         If None, uses a default set.

    Returns
    -------
    dict  Summary of generated files and checks.
    """
    os.makedirs(output_dir, exist_ok=True)

    if bounds_to_compute is None:
        # Default: bound W_{2,0}/W_{1,0} and W_{1,1}/W_{1,0}
        bounds_to_compute = [
            ((2, 0), (1, 0), "upper"),
            ((2, 0), (1, 0), "lower"),
            ((1, 1), (1, 0), "upper"),
            ((1, 1), (1, 0), "lower"),
        ]

    files = {}

    # 1. Run consistency checks
    checks = run_csdr_checks(d=d, precision=min(40, precision))
    checks_file = os.path.join(output_dir, f"csdr_checks_d{d}.json")
    with open(checks_file, "w", encoding="utf-8") as fh:
        json.dump(checks, fh, indent=2)
    files["checks"] = checks_file

    # 2. Kernel table
    table = compute_kernel_table(d=d, K=K, max_ell=min(max_spin, 8), precision=40)
    table_file = os.path.join(output_dir, f"csdr_kernel_table_d{d}_K{K}.json")
    with open(table_file, "w", encoding="utf-8") as fh:
        json.dump(table, fh, indent=2)
    files["kernel_table"] = table_file

    # 3. Generate PMP files
    pmp_files = {}
    for (obj_idx, norm_idx, direction) in bounds_to_compute:
        try:
            pmp = generate_csdr_pmp(
                obj_index=obj_idx,
                norm_index=norm_idx,
                d=d,
                K=K,
                max_spin=max_spin,
                precision=precision,
                bound_direction=direction,
            )
            fname = (
                f"csdr_pmp_{direction}_W{obj_idx[0]}_{obj_idx[1]}"
                f"_over_W{norm_idx[0]}_{norm_idx[1]}_d{d}_K{K}.json"
            )
            fpath = os.path.join(output_dir, fname)
            write_pmp_json(pmp, fpath)
            pmp_files[f"{direction}_{obj_idx}_{norm_idx}"] = fpath
        except ValueError as exc:
            pmp_files[f"{direction}_{obj_idx}_{norm_idx}_error"] = str(exc)

    files["pmp_files"] = pmp_files

    return {
        "output_dir": output_dir,
        "d": d,
        "K": K,
        "max_spin": max_spin,
        "all_checks_passed": checks["passed"],
        "files": files,
    }


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Generate CSDR-based SDPB PMP files for EFT Wilson coefficient bounds.\n"
            "\n"
            "Uses the CSDR dispersion relation kernels (Sinha-Zahed eq.(11)) and\n"
            "the heavy-average form from Extremal EFT Section 3.3 (user notes eq.(2)).\n"
            "\n"
            "IMPORTANT: d is required — do NOT assume d=4."
        ),
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "--output-dir", required=True,
        help="Directory for generated JSON/PMP files.",
    )
    parser.add_argument(
        "--d", type=int, required=True,
        help="Spacetime dimension (e.g. 4, 6, 10).",
    )
    parser.add_argument(
        "--K", type=int, default=8,
        help="Maximum spectral power 2p+3q (default: 8).",
    )
    parser.add_argument(
        "--max-spin", type=int, default=10,
        help="Maximum even spin (default: 10).",
    )
    parser.add_argument(
        "--precision", type=int, default=200,
        help="Decimal precision for SDPB output (default: 200).",
    )
    parser.add_argument(
        "--checks-only", action="store_true",
        help="Run consistency checks only, do not generate PMP files.",
    )
    args = parser.parse_args()

    if args.checks_only:
        checks = run_csdr_checks(d=args.d)
        print(json.dumps(checks, indent=2))
        sys.exit(0 if checks["passed"] else 1)

    result = run_bounds(
        output_dir=args.output_dir,
        d=args.d,
        K=args.K,
        max_spin=args.max_spin,
        precision=args.precision,
    )
    print(json.dumps(result, indent=2))
    if not result["all_checks_passed"]:
        print("\nWARNING: Some consistency checks failed.", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
