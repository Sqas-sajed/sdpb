#!/usr/bin/env python3
"""
generate_examples.py — Generate all example visualization outputs.

Produces 800×800 JPEG and BMP images illustrating:
  - Gegenbauer polynomials
  - Spectral functions v_ℓ^{(p,q)}(x)
  - Null constraint formulas (text rendering)
  - Null constraint coefficient heatmaps (s↔u and full crossing)
  - Variable elimination diagrams (s↔u and full crossing)
  - Bound structure comparison charts
  - 2D allowed-region projections
  - Summary dashboard

All outputs are written to the `examples/outputs/` directory alongside
this script.

Usage:
    python -m eft_bounds.examples.generate_examples
    # or directly:
    python eft_bounds/examples/generate_examples.py
"""

import os
import sys

# Ensure parent package is importable regardless of how the script is invoked
_here = os.path.dirname(os.path.abspath(__file__))
_pkg_root = os.path.dirname(os.path.dirname(_here))
if _pkg_root not in sys.path:
    sys.path.insert(0, _pkg_root)

from fractions import Fraction

from eft_bounds.visualize import (
    plot_gegenbauer_polynomials,
    plot_spectral_functions,
    render_null_constraint_formulas,
    plot_null_constraints_heatmap,
    plot_variable_elimination,
    plot_bound_comparison,
    plot_allowed_region_2d,
    plot_summary_dashboard,
)

# ---------------------------------------------------------------------------
# Output directory — all images land here
# ---------------------------------------------------------------------------

OUTPUT_DIR = os.path.join(_here, "outputs")
os.makedirs(OUTPUT_DIR, exist_ok=True)


def p(name: str) -> str:
    """Return full output path for a given filename stem."""
    return os.path.join(OUTPUT_DIR, name)


# ---------------------------------------------------------------------------
# Main generation routine
# ---------------------------------------------------------------------------

def main() -> None:
    print("=" * 65)
    print("EFT BOUNDS — EXAMPLE VISUALIZATION GENERATOR")
    print("Output directory:", OUTPUT_DIR)
    print("=" * 65)
    print()

    m_sq = Fraction(1)

    # ------------------------------------------------------------------
    # 1. Gegenbauer polynomials
    # ------------------------------------------------------------------
    print("1. Gegenbauer polynomials …")
    plot_gegenbauer_polynomials(
        output_path=p("gegenbauer_polynomials.jpg"),
        spins=[0, 1, 2, 3, 4, 5],
        d=4,
    )

    # ------------------------------------------------------------------
    # 2. Spectral functions v_ℓ^{(p,q)}(x)
    # ------------------------------------------------------------------
    print("2. Spectral functions …")
    plot_spectral_functions(
        output_path=p("spectral_functions.jpg"),
        spins=[0, 2, 4],
        pq_pairs=[(1, 0), (2, 0), (0, 1), (3, 0)],
        m_sq=m_sq,
        x_max=12.0,
    )

    # ------------------------------------------------------------------
    # 3. Null constraint formulas — s↔u crossing
    # ------------------------------------------------------------------
    print("3a. Null constraint formulas (s↔u crossing) …")
    render_null_constraint_formulas(
        output_path=p("null_constraint_formulas_su.jpg"),
        max_order=4,
        m_sq=m_sq,
        crossing_type="su",
    )

    # ------------------------------------------------------------------
    # 3b. Null constraint formulas — full S₃ crossing
    # ------------------------------------------------------------------
    print("3b. Null constraint formulas (full S₃ crossing) …")
    render_null_constraint_formulas(
        output_path=p("null_constraint_formulas_full.jpg"),
        max_order=4,
        m_sq=m_sq,
        crossing_type="full",
    )

    # ------------------------------------------------------------------
    # 4. Null constraint heatmaps
    # ------------------------------------------------------------------
    print("4a. Null constraint heatmap (s↔u crossing) …")
    plot_null_constraints_heatmap(
        output_path=p("null_constraints_heatmap_su.jpg"),
        max_order=4,
        m_sq=m_sq,
        crossing_type="su",
    )

    print("4b. Null constraint heatmap (full S₃ crossing) …")
    plot_null_constraints_heatmap(
        output_path=p("null_constraints_heatmap_full.jpg"),
        max_order=4,
        m_sq=m_sq,
        crossing_type="full",
    )

    # ------------------------------------------------------------------
    # 5. Variable elimination diagrams
    # ------------------------------------------------------------------
    print("5a. Variable elimination (s↔u crossing) …")
    plot_variable_elimination(
        output_path=p("variable_elimination_su.jpg"),
        max_order=4,
        m_sq=m_sq,
        crossing_type="su",
    )

    print("5b. Variable elimination (full S₃ crossing) …")
    plot_variable_elimination(
        output_path=p("variable_elimination_full.jpg"),
        max_order=4,
        m_sq=m_sq,
        crossing_type="full",
    )

    # ------------------------------------------------------------------
    # 6. Bound structure comparison (degrees of freedom chart)
    # ------------------------------------------------------------------
    print("6. Bound comparison chart …")
    plot_bound_comparison(
        output_path=p("bound_comparison.jpg"),
        max_order=3,
        m_sq=m_sq,
        max_spin=6,
    )

    # ------------------------------------------------------------------
    # 7. 2D allowed-region projections
    # ------------------------------------------------------------------
    print("7a. 2D allowed region — W_{1,0} vs W_{2,0} …")
    plot_allowed_region_2d(
        output_path=p("allowed_region_W10_vs_W20.jpg"),
        max_order=3,
        m_sq=m_sq,
        x_index=(1, 0),
        y_index=(2, 0),
    )

    print("7b. 2D allowed region — W_{1,0} vs W_{0,1} …")
    plot_allowed_region_2d(
        output_path=p("allowed_region_W10_vs_W01.jpg"),
        max_order=3,
        m_sq=m_sq,
        x_index=(1, 0),
        y_index=(0, 1),
    )

    # ------------------------------------------------------------------
    # 8. Summary dashboard
    # ------------------------------------------------------------------
    print("8. Summary dashboard …")
    plot_summary_dashboard(
        output_path=p("summary_dashboard.jpg"),
        max_order=3,
        m_sq=m_sq,
    )

    # ------------------------------------------------------------------
    # Done
    # ------------------------------------------------------------------
    print()
    print("=" * 65)
    print("All example images generated successfully.")
    print()

    # List produced files
    files = sorted(os.listdir(OUTPUT_DIR))
    print(f"Files in {OUTPUT_DIR}:")
    for f in files:
        full = os.path.join(OUTPUT_DIR, f)
        size_kb = os.path.getsize(full) / 1024
        print(f"  {f:<55s}  {size_kb:7.1f} kB")


if __name__ == "__main__":
    main()
