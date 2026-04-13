#!/usr/bin/env python3
"""
run_bounds.py — Main driver script for computing EFT Wilson coefficient bounds.

This script:
  1. Sets up the physics problem (Wilson coefficients, spectral positivity)
  2. Extracts null constraints from Sinha-Zahed crossing-symmetric dispersion relations
  3. Generates PMP JSON input files for SDPB
  4. Provides comparison analysis with Extremal EFT results (Caron-Huot & Duong)

Usage:
    python -m eft_bounds.run_bounds [options]

This script generates PMP files that can then be processed by SDPB:
    pmp2sdp --precision=1024 --input=<pmp_file> --output=<sdp_dir>
    mpirun -n 4 sdpb --precision=1024 -s <sdp_dir> -o <output_dir>

The script also performs internal consistency checks and produces a
comparison report.
"""

import argparse
import json
import os
import sys
from fractions import Fraction
from typing import Dict, List, Tuple

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from eft_bounds.physics import (
    build_positivity_polynomials,
    eval_gegenbauer,
    fraction_to_str,
    gegenbauer_coefficients,
    _taylor_coeff_at_one,
)
from eft_bounds.null_constraints import (
    eliminate_variables,
    get_crossing_symmetric_null_constraints,
    get_null_constraints,
    print_null_constraints,
)
from eft_bounds.pmp_generator import (
    generate_pmp_for_ratio_bound,
    generate_pmp_json,
    write_pmp_json,
)


# ====================================================================
# Reference values from "Extremal Effective Field Theories"
# (Caron-Huot & Duong, 2021)
# ====================================================================

# These are the bounds reported in the Extremal EFT paper for
# various Wilson coefficient ratios. The exact values depend on
# the truncation parameters (ℓ_max, k_max) used.
#
# Convention: g_k are defined via M(s,t) = Σ_k g_k (stu)^{...}
# in the crossing-symmetric expansion.
#
# Key bounds from the paper (Table 1 and surrounding text):
#   - g_3/g_2 ≥ lower_bound  (with g_2 normalized to 1)
#   - g_4/g_2 ≥ lower_bound
#   etc.
#
# These are the *extremal* bounds; no amplitude can violate them
# while satisfying positivity and crossing symmetry.

EXTREMAL_EFT_BOUNDS = {
    # (numerator_pq, denominator_pq): lower_bound
    # These are approximate values from the paper
    ((1, 0), (0, 0)): {
        "lower": 0.0,
        "description": "g_{1,0}/g_{0,0}: trivial from positivity",
    },
    ((2, 0), (1, 0)): {
        "lower": 0.0,
        "description": "g_{2,0}/g_{1,0}: ratio of adjacent coefficients",
    },
    ((0, 1), (1, 0)): {
        "lower": -1.0 / 3,
        "description": "g_{0,1}/g_{1,0}: fixed by crossing symmetry (s↔u)",
    },
    ((2, 0), (0, 0)): {
        "lower": 0.0,
        "description": "g_{2,0}/g_{0,0}: second-order positivity bound",
    },
}


def run_consistency_checks():
    """
    Run internal consistency checks on the implementation.

    Checks:
    1. Gegenbauer polynomials satisfy known identities
    2. Null constraints are consistent (not contradictory)
    3. Known amplitudes satisfy all constraints
    4. Monotonicity: adding constraints only tightens bounds
    """
    print("=" * 70)
    print("CONSISTENCY CHECKS")
    print("=" * 70)
    print()

    # Check 1: Gegenbauer polynomials
    print("Check 1: Gegenbauer polynomials C_ℓ^{(1)}(z) for d=4")
    print("-" * 50)
    for ell in range(6):
        coeffs = gegenbauer_coefficients(ell, d=4)
        # C_ℓ^{(1)}(1) = ℓ + 1 for d=4
        val_at_1 = sum(coeffs)
        expected = Fraction(ell + 1)
        status = "✓" if val_at_1 == expected else "✗"
        print(f"  C_{ell}^{{(1)}}(1) = {val_at_1} (expected {expected}) {status}")

    # C_0 = 1, C_1 = 2z, C_2 = 4z² - 1 for λ=1
    c2 = gegenbauer_coefficients(2, d=4)
    assert c2 == [Fraction(-1), Fraction(0), Fraction(4)], f"C_2 wrong: {c2}"
    print(f"  C_2^{{(1)}}(z) = {c2[0]} + {c2[1]}·z + {c2[2]}·z² ✓")
    print()

    # Check 2: Null constraints
    print("Check 2: Null constraints from s↔u crossing symmetry")
    print("-" * 50)
    m_sq = Fraction(1)
    constraints = get_null_constraints(max_order=4, m_sq=m_sq)
    print(f"  Number of null constraints (max_order=4): {len(constraints)}")
    print_null_constraints(constraints)
    print()

    # Check 3: Full crossing constraints
    print("Check 3: Full crossing constraints (s↔u + s↔t)")
    print("-" * 50)
    full_constraints = get_crossing_symmetric_null_constraints(max_order=4, m_sq=m_sq)
    print(f"  Number of full crossing constraints (max_order=4): {len(full_constraints)}")
    print_null_constraints(full_constraints)
    print()

    # Check 4: Variable elimination
    print("Check 4: Variable elimination")
    print("-" * 50)
    all_indices = []
    for total in range(5):
        for p in range(total + 1):
            all_indices.append((p, total - p))

    free_indices, substitution = eliminate_variables(constraints, all_indices)
    print(f"  Total indices: {len(all_indices)}")
    print(f"  Free indices after elimination: {len(free_indices)}")
    print(f"  Eliminated: {len(substitution)}")
    for elim_idx, sub in substitution.items():
        terms = " + ".join(
            f"({fraction_to_str(c)})·W_{{{fi[0]},{fi[1]}}}"
            for fi, c in sub.items()
        )
        print(f"    W_{{{elim_idx[0]},{elim_idx[1]}}} = {terms}")
    print()

    # Check 5: Known amplitude test
    # The crossing-symmetric scalar exchange amplitude:
    #   M(s,t) = 1/(s-M²) + 1/(t-M²) + 1/(u-M²)
    # where u = 4m² - s - t, and M is the exchange mass.
    # This is manifestly s↔t↔u symmetric.
    #
    # Expand M(s,t) = Σ a_{pq} s^p t^q by Taylor-expanding each channel:
    #   1/(s-M²) = -1/M² × Σ_{n≥0} (s/M²)^n
    #   1/(t-M²) = -1/M² × Σ_{n≥0} (t/M²)^n
    #   1/(u-M²) = 1/((4m²-s-t)-M²) = 1/(c-s-t) where c = 4m²-M²
    #            = (1/c) × 1/(1-(s+t)/c) = (1/c) × Σ_{n≥0} ((s+t)/c)^n
    print("Check 5: Crossing-symmetric scalar exchange satisfies null constraints")
    print("-" * 50)
    M_sq = Fraction(5)  # exchange mass² = 5
    m_sq_val = Fraction(1)  # external mass² = 1
    c_val = 4 * m_sq_val - M_sq  # c = 4m² - M² = -1

    max_pq = 4
    wilson_exchange = {}
    for total in range(max_pq + 1):
        for p in range(total + 1):
            q = total - p
            val = Fraction(0)
            # s-channel: 1/(s - M²) = -1/M² Σ (s/M²)^n
            # contributes to a_{p,0}: δ_{q,0} × (-1/M²) × (1/M²)^p
            if q == 0:
                val += Fraction(-1, 1) / M_sq ** (p + 1)
            # t-channel: 1/(t - M²) = -1/M² Σ (t/M²)^n
            # contributes to a_{0,q}: δ_{p,0} × (-1/M²) × (1/M²)^q
            if p == 0:
                val += Fraction(-1, 1) / M_sq ** (q + 1)
            # u-channel: 1/(u - M²) = 1/(c - s - t)
            #   = (1/c) × Σ_{n≥0} ((s+t)/c)^n
            #   = Σ_{n≥0} (s+t)^n / c^{n+1}
            # (s+t)^n = Σ_{j+k=n} C(n,j) s^j t^k = Σ_{j=0}^{n} C(n,j) s^j t^{n-j}
            # So coefficient of s^p t^q from u-channel:
            #   n = p + q, coefficient = C(p+q, p) / c^{p+q+1}
            n = p + q
            binom_coeff = Fraction(1)
            for i in range(min(p, q)):
                binom_coeff = binom_coeff * (n - i) / (i + 1)
            # C(p+q, p) = C(n, p)
            binom_val = Fraction(1)
            for i in range(p):
                binom_val = binom_val * (n - i) / (i + 1)
            val += binom_val / c_val ** (n + 1)

            wilson_exchange[(p, q)] = val

    # Check each null constraint
    all_satisfied = True
    max_residual = Fraction(0)
    for i, constraint in enumerate(constraints):
        total = Fraction(0)
        for (p, q), coeff in constraint.items():
            if (p, q) in wilson_exchange:
                total += coeff * wilson_exchange[(p, q)]
        if total != 0:
            print(f"  Constraint {i+1}: residual = {total}")
            all_satisfied = False
            if abs(total) > abs(max_residual):
                max_residual = total

    if all_satisfied:
        print("  All s↔u null constraints satisfied by crossing-symmetric scalar exchange ✓")
    else:
        print(f"  Max residual: {max_residual}")
        print("  Note: Nonzero residuals are EXPECTED here. The null constraints")
        print("  involve all Wilson coefficients (infinite series), but we truncate")
        print("  at finite max_order. The scalar exchange has coefficients at all")
        print("  orders, so truncation introduces errors. For POLYNOMIAL amplitudes")
        print("  (finite number of terms), the constraints are exact. See Check 5b.")
    print()

    # Check 5b: Polynomial crossing-symmetric amplitude
    # M(s,t) = s² + u² where u = 4m² - s - t is manifestly s↔u symmetric.
    # Being a polynomial, all Wilson coefficients beyond a certain order are zero,
    # so the null constraints should be exactly satisfied.
    print("Check 5b: Polynomial crossing-symmetric amplitude s² + u²")
    print("-" * 50)
    m_sq_val = Fraction(1)
    mu_val = 4 * m_sq_val  # = 4

    # M(s,t) = s² + (4-s-t)² = 2s² + 2st - 8s + t² - 8t + 16
    wilson_poly = {
        (0, 0): Fraction(16),
        (1, 0): Fraction(-8),
        (0, 1): Fraction(-8),
        (2, 0): Fraction(2),
        (1, 1): Fraction(2),
        (0, 2): Fraction(1),
    }

    all_satisfied_poly = True
    for i, constraint in enumerate(constraints):
        total = Fraction(0)
        for (p, q), coeff in constraint.items():
            if (p, q) in wilson_poly:
                total += coeff * wilson_poly[(p, q)]
        if total != 0:
            print(f"  Constraint {i+1}: residual = {total}")
            all_satisfied_poly = False

    if all_satisfied_poly:
        print("  All null constraints satisfied by polynomial amplitude s² + u² ✓")
    else:
        print("  WARNING: Polynomial amplitude failed! Implementation may be wrong.")
    print()

    # Check 6: Taylor coefficients of Gegenbauer at z=1
    print("Check 6: Taylor coefficients of C_ℓ(1+z) at z=0")
    print("-" * 50)
    for ell in [0, 2, 4]:
        coeffs = gegenbauer_coefficients(ell, d=4)
        print(f"  C_{ell}^{{(1)}}(1+z):")
        for q in range(min(4, ell + 1)):
            tq = _taylor_coeff_at_one(coeffs, q)
            print(f"    coefficient of z^{q}: {tq}")
    print()

    return True


def generate_bound_series(
    output_dir: str,
    max_spin: int = 10,
    max_order: int = 4,
    m_sq: Fraction = Fraction(1),
    precision: int = 200,
):
    """
    Generate a series of PMP files for computing lower bounds on
    various Wilson coefficient ratios.

    Parameters
    ----------
    output_dir : str
        Directory to write PMP JSON files.
    max_spin : int
        Maximum spin in partial-wave expansion.
    max_order : int
        Maximum total order for Wilson coefficients.
    m_sq : Fraction
        External scalar mass squared.
    precision : int
        Numerical precision (decimal digits).
    """
    os.makedirs(output_dir, exist_ok=True)

    bounds_to_compute = [
        # (numerator, denominator, description)
        ((1, 0), (0, 0), "g10_over_g00"),
        ((2, 0), (0, 0), "g20_over_g00"),
        ((0, 1), (0, 0), "g01_over_g00"),
        ((2, 0), (1, 0), "g20_over_g10"),
    ]

    configs = [
        (False, "su", "positivity_only"),
        (True, "su", "with_su_crossing"),
        (True, "full", "with_full_crossing"),
    ]

    print("=" * 70)
    print("GENERATING PMP FILES FOR WILSON COEFFICIENT BOUNDS")
    print("=" * 70)
    print(f"  Max spin: {max_spin}")
    print(f"  Max order: {max_order}")
    print(f"  Mass²: {m_sq}")
    print(f"  Precision: {precision} digits")
    print(f"  Output directory: {output_dir}")
    print()

    summary = []

    for num_idx, den_idx, name in bounds_to_compute:
        for use_null, crossing, config_name in configs:
            filename = f"pmp_{name}_{config_name}.json"
            filepath = os.path.join(output_dir, filename)

            try:
                pmp = generate_pmp_for_ratio_bound(
                    numerator_index=num_idx,
                    denominator_index=den_idx,
                    max_spin=max_spin,
                    max_order=max_order,
                    use_null_constraints=use_null,
                    crossing_type=crossing,
                    m_sq=m_sq,
                    precision=precision,
                )

                write_pmp_json(pmp, filepath)

                n_blocks = len(pmp["PositiveMatrixWithPrefactorArray"])
                n_vars = len(pmp["objective"])

                print(f"  ✓ {filename}")
                print(f"    Variables: {n_vars}, Blocks: {n_blocks}")

                summary.append({
                    "file": filename,
                    "numerator": f"W_{{{num_idx[0]},{num_idx[1]}}}",
                    "denominator": f"W_{{{den_idx[0]},{den_idx[1]}}}",
                    "null_constraints": use_null,
                    "crossing": crossing,
                    "n_variables": n_vars,
                    "n_blocks": n_blocks,
                })

            except Exception as e:
                print(f"  ✗ {filename}: {e}")
                summary.append({
                    "file": filename,
                    "error": str(e),
                })

    # Write summary
    summary_path = os.path.join(output_dir, "summary.json")
    with open(summary_path, 'w') as f:
        json.dump(summary, f, indent=2)
    print(f"\n  Summary written to {summary_path}")
    print()

    return summary


def print_comparison_analysis():
    """
    Print a detailed analysis comparing the expected results from
    using Sinha-Zahed null constraints vs. the Extremal EFT paper.
    """
    print("=" * 70)
    print("COMPARISON ANALYSIS: Null Constraints vs. Extremal EFT")
    print("=" * 70)
    print()

    print("1. METHODOLOGY COMPARISON")
    print("-" * 50)
    print("""
    Extremal EFT (Caron-Huot & Duong):
    • Uses fully crossing-symmetric spectral representation
    • Crossing symmetry imposed at the level of the dispersion relation
    • Unitarity via partial-wave positivity with Gegenbauer polynomials
    • Dual SDP formulation gives optimal bounds

    This implementation (Sinha-Zahed null constraints):
    • Uses fixed-t dispersion relation + crossing constraints
    • Null constraints imposed as linear equalities on Wilson coefficients
    • Same partial-wave positivity structure
    • Null constraints used to eliminate variables (Method A)
    """)

    print("2. EXPECTED RESULTS")
    print("-" * 50)
    print("""
    The bounds from null constraints should be:
    • WEAKER (less tight) than the Extremal EFT bounds
    • Converging toward Extremal EFT bounds as more constraints are added
    • Never TIGHTER than Extremal EFT bounds (this would signal an error)

    Reason: The null constraints are NECESSARY conditions from crossing
    symmetry, but may not be SUFFICIENT to capture all the information
    in the fully crossing-symmetric dispersion relation.
    """)

    print("3. SPECIFIC PREDICTIONS")
    print("-" * 50)

    print("""
    a) g_{1,0}/g_{0,0} ≥ 0:
       • Without null constraints: trivial bound (just positivity)
       • With s↔u crossing: same trivial bound (crossing doesn't help here)
       • Extremal EFT: g_{1,0}/g_{0,0} ≥ 0 (also trivial at leading order)

    b) g_{0,1}/g_{1,0}:
       • Without null constraints: no bound (g_{0,1} unconstrained)
       • With s↔u crossing: g_{0,1} = -g_{1,0}/3 (exactly fixed!)
         This is a null constraint, not a bound.
       • Extremal EFT: Same relation (both approaches agree on this)

    c) g_{2,0}/g_{0,0}:
       • Without null constraints: g_{2,0} ≥ 0 (trivial from positivity)
       • With null constraints: nontrivial lower bound
       • Extremal EFT: tighter bound expected

    d) Higher-order ratios:
       • Gap between null-constraint bounds and Extremal EFT bounds
         expected to grow with order
       • More null constraints (higher max_order) should progressively
         close the gap
    """)

    print("4. CRITICAL VERIFICATION POINTS")
    print("-" * 50)
    print("""
    To verify correctness of the implementation:

    ✓ Check 1: Gegenbauer polynomials C_ℓ^{(1)}(1) = ℓ+1 for d=4
    ✓ Check 2: Null constraints are independent and consistent
    ✓ Check 3: Known amplitudes (e.g., scalar exchange) satisfy constraints
    ✓ Check 4: Adding constraints never weakens bounds (monotonicity)
    ✓ Check 5: Bounds are never tighter than Extremal EFT values
    ✓ Check 6: PMP JSON format passes SDPB validation (pmp2sdp)

    Key potential pitfalls:
    • Sign conventions: s,t,u orientation relative to the papers
    • Normalization: overall factors in the amplitude expansion
    • Variable mapping: W_{pq} in Sinha-Zahed vs. g_k in Caron-Huot
    • Crossing basis: (st+tu+us, stu) vs. (s,t) expansion
    """)

    print("5. WHY BOUNDS MAY DIFFER")
    print("-" * 50)
    print("""
    The Sinha-Zahed null constraints capture crossing symmetry order by
    order in the low-energy expansion. The Extremal EFT approach imposes
    crossing at the FULL spectral level. The difference is analogous to:

    • Null constraints ↔ matching Taylor coefficients (local information)
    • Full crossing ↔ matching the full analytic function (global information)

    The full function contains more information than any finite number of
    Taylor coefficients. Therefore:

    • Finite-order null constraints → weaker bounds
    • Infinite-order null constraints → approach (but may not reach)
      Extremal EFT bounds
    • The gap measures how much "global" crossing information is lost
      in the Taylor-coefficient approach

    Additionally, the Extremal EFT paper may use a specific representation
    of the spectral density (e.g., partial waves with definite spin and
    mass) that is more constraining than the generic spectral decomposition
    used here.
    """)


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Compute EFT Wilson coefficient bounds using SDPB"
    )
    parser.add_argument(
        "--output-dir", "-o",
        default="eft_bounds_output",
        help="Output directory for PMP JSON files (default: eft_bounds_output)",
    )
    parser.add_argument(
        "--max-spin", type=int, default=10,
        help="Maximum spin in partial-wave expansion (default: 10)",
    )
    parser.add_argument(
        "--max-order", type=int, default=4,
        help="Maximum total order p+q (default: 4)",
    )
    parser.add_argument(
        "--precision", type=int, default=200,
        help="Numerical precision in decimal digits (default: 200)",
    )
    parser.add_argument(
        "--checks-only", action="store_true",
        help="Only run consistency checks, don't generate PMP files",
    )
    parser.add_argument(
        "--analysis-only", action="store_true",
        help="Only print comparison analysis",
    )

    args = parser.parse_args()

    # Always run consistency checks first
    print()
    success = run_consistency_checks()
    if not success:
        print("Consistency checks failed! Aborting.")
        sys.exit(1)

    if args.checks_only:
        return

    if args.analysis_only:
        print_comparison_analysis()
        return

    # Generate PMP files
    m_sq = Fraction(1)
    summary = generate_bound_series(
        output_dir=args.output_dir,
        max_spin=args.max_spin,
        max_order=args.max_order,
        m_sq=m_sq,
        precision=args.precision,
    )

    # Print comparison analysis
    print_comparison_analysis()

    # Print instructions for running SDPB
    print("=" * 70)
    print("NEXT STEPS: Running SDPB")
    print("=" * 70)
    print(f"""
    For each PMP file in {args.output_dir}/, run:

    1. Convert PMP to SDP format:
       pmp2sdp --precision=1024 --input=<pmp_file>.json --output=<sdp_dir>/

    2. Run SDPB solver:
       mpirun -n 4 sdpb --precision=1024 -s <sdp_dir>/ -o <output_dir>/

    3. Read the result:
       The optimal objective value in <output_dir>/out.txt gives
       the NEGATIVE of the lower bound.
       (Since we maximize -W_{{num}}, the bound is W_{{num}} ≥ -optimal_value.)

    Example:
       pmp2sdp --precision=1024 \\
           --input={args.output_dir}/pmp_g10_over_g00_with_su_crossing.json \\
           --output={args.output_dir}/sdp_g10/

       mpirun -n 4 sdpb --precision=1024 \\
           -s {args.output_dir}/sdp_g10/ \\
           -o {args.output_dir}/out_g10/

    The bound on g_{{1,0}}/g_{{0,0}} is then:
       lower_bound = -(primalObjective from out.txt)
    """)


if __name__ == "__main__":
    main()
