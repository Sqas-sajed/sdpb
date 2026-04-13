"""
visualize.py — Visualization for EFT Wilson coefficient bounds and constraints.

This module provides functions to generate 800×800 JPEG and BMP images for:
  - Gegenbauer polynomial plots
  - Spectral function (v_ℓ^{(p,q)}(x)) plots
  - Null constraint formula renderings (text images)
  - Null constraint coefficient heatmaps
  - Variable elimination diagrams
  - Wilson coefficient bound comparison charts
  - 2D allowed-region plots in Wilson coefficient space

All outputs are 800×800 pixels saved as JPEG (.jpg) and/or BMP (.bmp).

Usage example:
    from eft_bounds.visualize import (
        plot_gegenbauer_polynomials,
        plot_spectral_functions,
        render_null_constraint_formulas,
        plot_null_constraints_heatmap,
        plot_variable_elimination,
        plot_bound_comparison,
        plot_allowed_region_2d,
    )
    plot_gegenbauer_polynomials("outputs/gegenbauer.jpg")
"""

from __future__ import annotations

import os
from fractions import Fraction
from typing import Dict, List, Optional, Tuple

import matplotlib
matplotlib.use("Agg")  # non-interactive backend
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
from matplotlib.colors import LinearSegmentedColormap
from matplotlib import rcParams
import numpy as np

# ---------------------------------------------------------------------------
# Global style settings
# ---------------------------------------------------------------------------

# Physics-paper style: use serif fonts and larger text
rcParams.update({
    "font.family": "serif",
    "font.size": 12,
    "axes.titlesize": 14,
    "axes.labelsize": 13,
    "xtick.labelsize": 11,
    "ytick.labelsize": 11,
    "legend.fontsize": 11,
    "figure.dpi": 100,         # so 800 pixels → 8 inch figure
    "savefig.dpi": 100,
    "axes.spines.top": False,
    "axes.spines.right": False,
})

_FIG_SIZE = (8.0, 8.0)        # inches → 800×800 at 100 dpi
_COLORS = plt.rcParams["axes.prop_cycle"].by_key()["color"]

# Custom blue-white-red diverging colormap for constraint heatmaps
_BWR = LinearSegmentedColormap.from_list(
    "bwr_soft",
    [(0.1, 0.3, 0.8), (1.0, 1.0, 1.0), (0.8, 0.1, 0.1)],
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _save(fig: plt.Figure, path: str, fmt: Optional[str] = None) -> None:
    """
    Save figure as JPEG (.jpg) and BMP (.bmp), each 800×800 pixels.

    Matplotlib does not support BMP natively, so we:
      1. Save the figure as JPEG via matplotlib (forced to 8×8 in at 100 dpi).
      2. Open the JPEG with Pillow, resize to exactly 800×800, and save as BMP.
    """
    from PIL import Image as PilImage

    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    base = os.path.splitext(path)[0]
    jpg_path = base + ".jpg"
    bmp_path = base + ".bmp"

    # Force figure to 8×8 inches before saving so output is exactly 800×800
    fig.set_size_inches(8.0, 8.0, forward=True)
    fig.savefig(jpg_path, format="jpeg", dpi=100, bbox_inches=None)
    plt.close(fig)

    # Convert JPEG → BMP via Pillow (ensure exactly 800×800)
    try:
        img = PilImage.open(jpg_path).resize((800, 800), PilImage.LANCZOS)
        img.save(bmp_path, format="BMP")
    except Exception as exc:  # pragma: no cover
        print(f"  Warning: BMP conversion failed: {exc}")

    print(f"  Saved: {jpg_path}")
    print(f"  Saved: {bmp_path}")


def _fraction_to_float(f: Fraction) -> float:
    return float(f.numerator) / float(f.denominator)


# ---------------------------------------------------------------------------
# 1. Gegenbauer polynomial plots
# ---------------------------------------------------------------------------

def plot_gegenbauer_polynomials(
    output_path: str = "outputs/gegenbauer.jpg",
    spins: List[int] = None,
    d: int = 4,
) -> None:
    """
    Plot Gegenbauer polynomials C_ℓ^{(d/2-1)}(z) for several spins.

    Produces an 800×800 JPEG (and BMP) showing the polynomials on [-1, 1]
    together with their values at z=1 (= ℓ+1 for d=4).

    Parameters
    ----------
    output_path : str
        Output file path (extension .jpg or .bmp).
    spins : list of int
        Spins to plot (default [0, 1, 2, 3, 4, 5]).
    d : int
        Spacetime dimension (default 4).
    """
    from .physics import gegenbauer_coefficients

    if spins is None:
        spins = [0, 1, 2, 3, 4, 5]

    z = np.linspace(-1.0, 1.0, 600)

    fig, axes = plt.subplots(1, 2, figsize=_FIG_SIZE)
    fig.suptitle(
        rf"Gegenbauer Polynomials $C_{{\ell}}^{{(\lambda)}}(z)$, "
        rf"$\lambda = (d-2)/2 = {(d-2)//2}$",
        fontsize=15,
    )

    ax_main = axes[0]
    ax_vals = axes[1]

    # Left panel: C_ℓ(z) on [-1, 1]
    for i, ell in enumerate(spins):
        coeffs = gegenbauer_coefficients(ell, d=d)
        c_float = [_fraction_to_float(c) for c in coeffs]
        y = np.polyval(c_float[::-1], z)  # polyval expects highest-degree first
        ax_main.plot(z, y, color=_COLORS[i % len(_COLORS)],
                     linewidth=2, label=rf"$\ell={ell}$")

    ax_main.axhline(0, color="black", linewidth=0.7, linestyle="--")
    ax_main.axvline(0, color="black", linewidth=0.7, linestyle=":")
    ax_main.set_xlabel(r"$z = \cos\theta$")
    ax_main.set_ylabel(r"$C_{\ell}^{(\lambda)}(z)$")
    ax_main.set_title(r"Polynomial values on $[-1,\,1]$")
    ax_main.legend(loc="upper left", frameon=False)
    ax_main.set_xlim(-1, 1)

    # Right panel: C_ℓ(1) = ℓ+1 for d=4 and coefficient table
    ell_list = list(spins)
    vals_at_1 = []
    for ell in ell_list:
        coeffs = gegenbauer_coefficients(ell, d=d)
        v = _fraction_to_float(sum(coeffs))
        vals_at_1.append(v)
        # Verify: C_ℓ^{(1)}(1) = ℓ+1
    expected = [ell + 1 for ell in ell_list]

    x_pos = np.arange(len(ell_list))
    ax_vals.bar(x_pos - 0.18, vals_at_1, width=0.35, label=r"$C_\ell^{(\lambda)}(1)$",
                color="#4477AA", alpha=0.85)
    ax_vals.bar(x_pos + 0.18, expected, width=0.35, label=r"Expected $\ell+1$",
                color="#EE6677", alpha=0.65, hatch="//")
    ax_vals.set_xticks(x_pos)
    ax_vals.set_xticklabels([rf"$\ell={e}$" for e in ell_list])
    ax_vals.set_ylabel(r"Value at $z=1$")
    ax_vals.set_title(r"$C_\ell^{(\lambda)}(1) = \ell+1$ identity check")
    ax_vals.legend(frameon=False)

    # Add annotation
    for xi, (computed, exp) in enumerate(zip(vals_at_1, expected)):
        match = "[OK]" if abs(computed - exp) < 1e-10 else "[FAIL]"
        ax_vals.text(xi, max(computed, exp) + 0.2, match,
                     ha="center", fontsize=9, color="green" if "[OK]" in match else "red")

    fig.tight_layout(rect=[0, 0, 1, 0.95])
    _save(fig, output_path)


# ---------------------------------------------------------------------------
# 2. Spectral function plots
# ---------------------------------------------------------------------------

def plot_spectral_functions(
    output_path: str = "outputs/spectral_functions.jpg",
    spins: List[int] = None,
    pq_pairs: List[Tuple[int, int]] = None,
    m_sq: Fraction = Fraction(1),
    x_max: float = 10.0,
) -> None:
    """
    Plot spectral functions v_ℓ^{(p,q)}(x) vs x = μ - 4m².

    Parameters
    ----------
    output_path : str
        Output file path.
    spins : list of int
        Spins to show (default [0, 2, 4]).
    pq_pairs : list of (int, int)
        (p, q) Wilson coefficient indices (default [(1,0), (2,0), (0,1)]).
    m_sq : Fraction
        Mass squared (default 1).
    x_max : float
        Upper limit for x axis.
    """
    from .physics import gegenbauer_coefficients, _taylor_coeff_at_one

    if spins is None:
        spins = [0, 2, 4]
    if pq_pairs is None:
        pq_pairs = [(1, 0), (2, 0), (0, 1), (3, 0)]

    mu_th = float(4 * m_sq)
    x = np.linspace(0.01, x_max, 400)

    n_pairs = len(pq_pairs)
    ncols = 2
    nrows = (n_pairs + 1) // 2

    fig, axes = plt.subplots(nrows, ncols, figsize=_FIG_SIZE)
    fig.suptitle(
        r"Spectral functions $v_\ell^{(p,q)}(x)$ vs $x = \mu - 4m^2$",
        fontsize=15,
    )
    axes = np.array(axes).flatten()

    for ax_idx, (p, q) in enumerate(pq_pairs):
        ax = axes[ax_idx]
        for i, ell in enumerate(spins):
            # v_ℓ^{(p,q)}(x) = c_q × 2^q / (x^q × (x + 4m²)^{p+1})
            # where c_q = [coeff of z^q in C_ℓ(1+z)]
            coeffs = gegenbauer_coefficients(ell, d=4)
            c_q = _fraction_to_float(_taylor_coeff_at_one(coeffs, q))
            factor = c_q * (2.0 ** q)
            denom = (x ** q) * ((x + mu_th) ** (p + 1))
            # Avoid zero denominator
            safe_denom = np.where(np.abs(denom) < 1e-12, np.nan, denom)
            v = factor / safe_denom
            ax.plot(x, v, color=_COLORS[i % len(_COLORS)],
                    linewidth=1.8, label=rf"$\ell={ell}$")

        ax.axhline(0, color="black", linewidth=0.7, linestyle="--")
        ax.set_xlabel(r"$x = \mu - 4m^2$")
        ax.set_ylabel(rf"$v_{{\ell}}^{{({p},{q})}}(x)$")
        ax.set_title(rf"$W_{{({p},{q})}}$ spectral function")
        ax.legend(frameon=False, loc="upper right")
        # Clip extreme values for readability
        finite_vals = v[np.isfinite(v)]
        if len(finite_vals) > 0:
            ymax = np.percentile(np.abs(finite_vals), 95) * 1.3
            ax.set_ylim(-0.05 * ymax, ymax)

    # Hide unused axes
    for ax_idx in range(len(pq_pairs), len(axes)):
        axes[ax_idx].set_visible(False)

    fig.tight_layout(rect=[0, 0, 1, 0.95])
    _save(fig, output_path)


# ---------------------------------------------------------------------------
# 3. Null constraint formula rendering
# ---------------------------------------------------------------------------

def render_null_constraint_formulas(
    output_path: str = "outputs/null_constraint_formulas.jpg",
    max_order: int = 4,
    m_sq: Fraction = Fraction(1),
    crossing_type: str = "su",
) -> None:
    """
    Render the null constraints as a typeset text image (800×800 JPEG/BMP).

    Each null constraint Σ C_{pq}^{(k)} W_{pq} = 0 is displayed as a
    formatted equation, with colour-coded positive (blue) and negative (red)
    terms.

    Parameters
    ----------
    output_path : str
        Output file path.
    max_order : int
        Maximum order for null constraints.
    m_sq : Fraction
        Mass squared.
    crossing_type : str
        "su" for s↔u only, "full" for full S₃ crossing.
    """
    from .null_constraints import (
        get_null_constraints,
        get_crossing_symmetric_null_constraints,
    )

    if crossing_type == "full":
        constraints = get_crossing_symmetric_null_constraints(max_order, m_sq)
        title = rf"Full S₃ Crossing Null Constraints (max order = {max_order})"
    else:
        constraints = get_null_constraints(max_order, m_sq)
        title = rf"Null Constraints from $s\leftrightarrow u$ Crossing (max order = {max_order})"

    fig, ax = plt.subplots(figsize=_FIG_SIZE)
    ax.set_axis_off()
    fig.patch.set_facecolor("white")

    ax.text(0.5, 0.98, title,
            transform=ax.transAxes, ha="center", va="top",
            fontsize=14, fontweight="bold", color="#222222")

    # Header line
    ax.text(0.5, 0.92,
            r"Each line is one constraint: $\sum_{p,q}\,C_{pq}^{(k)}\,W_{pq}=0$",
            transform=ax.transAxes, ha="center", va="top",
            fontsize=12, color="#444444", style="italic")

    ax.plot([0.02, 0.98], [0.89, 0.89], transform=ax.transAxes,
            color="#888888", linewidth=0.8, clip_on=False)

    n = len(constraints)
    # Vertical positions for each constraint
    y_start = 0.875
    y_step = min(0.055, 0.80 / max(n, 1))

    for k, constraint in enumerate(constraints):
        y = y_start - k * y_step
        if y < 0.02:
            ax.text(0.5, y + y_step * 0.5,
                    f"  … {n - k} more constraints not shown …",
                    transform=ax.transAxes, ha="center", va="center",
                    fontsize=10, color="#888888")
            break

        # Build the constraint string with colour segments
        # We'll render each term separately with colour
        terms = []
        for (p, q), coeff in sorted(constraint.items(), key=lambda item: item[0]):
            c_float = float(coeff)
            if coeff == Fraction(1):
                terms.append((f"$W_{{{p},{q}}}$", c_float))
            elif coeff == Fraction(-1):
                terms.append((f"$-W_{{{p},{q}}}$", c_float))
            else:
                # Format fraction nicely
                if coeff.denominator == 1:
                    c_str = str(int(coeff))
                else:
                    c_str = rf"\tfrac{{{coeff.numerator}}}{{{coeff.denominator}}}"
                sign = "+" if c_float >= 0 else ""
                terms.append((rf"${sign}{c_str}\,W_{{{p},{q}}}$", c_float))

        # Place constraint number label
        ax.text(0.02, y, rf"$k={k+1}$:", transform=ax.transAxes,
                ha="left", va="center", fontsize=10, color="#555555")

        # Place terms with colour
        x_cursor = 0.10
        for t_str, t_val in terms:
            color = "#1155AA" if t_val >= 0 else "#AA1111"
            ax.text(x_cursor, y, t_str + " ",
                    transform=ax.transAxes, ha="left", va="center",
                    fontsize=10.5, color=color)
            # Estimate width: rough heuristic based on string length
            x_cursor += max(0.055, len(t_str) * 0.011)
            if x_cursor > 0.88:
                # wrap to next visual line (just show "...")
                ax.text(x_cursor, y, "…", transform=ax.transAxes,
                        ha="left", va="center", fontsize=10, color="#888888")
                break

        ax.text(x_cursor, y, "$= 0$",
                transform=ax.transAxes, ha="left", va="center",
                fontsize=10.5, color="#222222")

    # Footer
    ax.text(0.5, 0.02,
            rf"Total constraints: {n} | "
            rf"$m^2={m_sq}$ | crossing: {'full S₃' if crossing_type=='full' else 's↔u'}",
            transform=ax.transAxes, ha="center", va="bottom",
            fontsize=10, color="#777777")

    fig.tight_layout()
    _save(fig, output_path)


# ---------------------------------------------------------------------------
# 4. Null constraint coefficient heatmap
# ---------------------------------------------------------------------------

def plot_null_constraints_heatmap(
    output_path: str = "outputs/null_constraints_heatmap.jpg",
    max_order: int = 4,
    m_sq: Fraction = Fraction(1),
    crossing_type: str = "su",
) -> None:
    """
    Plot a heatmap of the null constraint coefficient matrix C_{pq}^{(k)}.

    Rows = constraint index k, Columns = Wilson coefficient (p,q).
    Colour encodes the coefficient value (blue = positive, red = negative).

    Parameters
    ----------
    output_path : str
        Output file path.
    max_order : int
        Maximum order.
    m_sq : Fraction
        Mass squared.
    crossing_type : str
        "su" or "full".
    """
    from .null_constraints import (
        get_null_constraints,
        get_crossing_symmetric_null_constraints,
    )

    if crossing_type == "full":
        constraints = get_crossing_symmetric_null_constraints(max_order, m_sq)
        title = rf"Null Constraint Matrix $C_{{pq}}^{{(k)}}$ — Full S₃ Crossing"
    else:
        constraints = get_null_constraints(max_order, m_sq)
        title = rf"Null Constraint Matrix $C_{{pq}}^{{(k)}}$ — $s\leftrightarrow u$ Crossing"

    # Collect all (p,q) indices that appear
    all_pq = set()
    for c in constraints:
        all_pq.update(c.keys())
    # Sort: by total order, then q
    sorted_pq = sorted(all_pq, key=lambda x: (x[0] + x[1], x[1]))

    n_k = len(constraints)
    n_pq = len(sorted_pq)
    pq_to_col = {pq: i for i, pq in enumerate(sorted_pq)}

    # Build matrix
    matrix = np.zeros((n_k, n_pq))
    for k, constraint in enumerate(constraints):
        for pq, coeff in constraint.items():
            matrix[k, pq_to_col[pq]] = float(coeff)

    fig, ax = plt.subplots(figsize=_FIG_SIZE)

    vmax = np.abs(matrix).max() or 1.0
    im = ax.imshow(matrix, aspect="auto", cmap=_BWR,
                   vmin=-vmax, vmax=vmax, interpolation="nearest")

    # Axis labels
    col_labels = [rf"$W_{{{p},{q}}}$" for p, q in sorted_pq]
    row_labels = [rf"$k={k+1}$" for k in range(n_k)]

    ax.set_xticks(range(n_pq))
    ax.set_xticklabels(col_labels, rotation=45, ha="right", fontsize=9)
    ax.set_yticks(range(n_k))
    ax.set_yticklabels(row_labels, fontsize=9)
    ax.set_xlabel(r"Wilson coefficient $W_{pq}$")
    ax.set_ylabel(r"Constraint index $k$")
    ax.set_title(title, fontsize=13, pad=12)

    # Annotate cells with non-zero values
    for k in range(n_k):
        for j in range(n_pq):
            val = matrix[k, j]
            if abs(val) > 1e-12:
                frac_val = constraints[k].get(sorted_pq[j], Fraction(0))
                text = _format_fraction_short(frac_val)
                txt_color = "white" if abs(val) > 0.4 * vmax else "black"
                ax.text(j, k, text, ha="center", va="center",
                        fontsize=8.5, color=txt_color, fontweight="bold")

    cbar = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.04)
    cbar.set_label(r"Coefficient $C_{pq}^{(k)}$", fontsize=11)

    fig.tight_layout()
    _save(fig, output_path)


def _format_fraction_short(f: Fraction) -> str:
    """Format a Fraction as a short string for heatmap annotations."""
    if f == 0:
        return ""
    if f.denominator == 1:
        return str(int(f))
    # Shorten: if denominator ≤ 9, write as fraction
    if abs(f.denominator) <= 9:
        return f"{f.numerator}/{f.denominator}"
    # Otherwise scientific-ish
    v = float(f)
    if abs(v) >= 10:
        return f"{v:.0f}"
    return f"{v:.2g}"


# ---------------------------------------------------------------------------
# 5. Variable elimination diagram
# ---------------------------------------------------------------------------

def plot_variable_elimination(
    output_path: str = "outputs/variable_elimination.jpg",
    max_order: int = 4,
    m_sq: Fraction = Fraction(1),
    crossing_type: str = "su",
) -> None:
    """
    Visualize the variable elimination tree from null constraints.

    Shows which Wilson coefficients are FREE (independent) and which are
    ELIMINATED (expressed as linear combinations of free variables).
    Also displays the explicit substitution for each eliminated variable.

    Parameters
    ----------
    output_path : str
        Output file path.
    max_order : int
        Maximum order.
    m_sq : Fraction
        Mass squared.
    crossing_type : str
        "su" or "full".
    """
    from .null_constraints import (
        get_null_constraints,
        get_crossing_symmetric_null_constraints,
        eliminate_variables,
    )

    if crossing_type == "full":
        constraints = get_crossing_symmetric_null_constraints(max_order, m_sq)
        title_suffix = "Full S₃"
    else:
        constraints = get_null_constraints(max_order, m_sq)
        title_suffix = r"$s\leftrightarrow u$"

    all_indices = []
    for total in range(max_order + 1):
        for p in range(total + 1):
            q = total - p
            all_indices.append((p, q))

    free_indices, substitution = eliminate_variables(constraints, all_indices)

    fig, axes = plt.subplots(1, 2, figsize=_FIG_SIZE)
    fig.suptitle(
        rf"Variable Elimination via Null Constraints ({title_suffix} crossing, order ≤ {max_order})",
        fontsize=13,
    )

    # Left panel: grid showing free vs. eliminated
    ax_grid = axes[0]
    max_tot = max_order
    grid_size = max_tot + 1

    for p in range(grid_size):
        for q in range(grid_size):
            if p + q > max_tot:
                continue
            pq = (p, q)
            if pq in substitution:
                color = "#EE4444"  # eliminated — red
                label = "elim."
            else:
                color = "#44AA44"  # free — green
                label = "free"

            rect = mpatches.FancyBboxPatch(
                (q - 0.42, p - 0.42), 0.84, 0.84,
                boxstyle="round,pad=0.05",
                facecolor=color, alpha=0.75, edgecolor="white", linewidth=1.5,
            )
            ax_grid.add_patch(rect)
            ax_grid.text(q, p, rf"$W_{{{p},{q}}}$",
                         ha="center", va="center", fontsize=9.5,
                         color="white", fontweight="bold")

    ax_grid.set_xlim(-0.6, max_tot + 0.6)
    ax_grid.set_ylim(-0.6, max_tot + 0.6)
    ax_grid.set_xlabel(r"$q$ (power of $t$)")
    ax_grid.set_ylabel(r"$p$ (power of $s$)")
    ax_grid.set_title(
        rf"Free ({len(free_indices)}) vs. Eliminated ({len(substitution)})",
        fontsize=12,
    )
    ax_grid.set_xticks(range(grid_size))
    ax_grid.set_yticks(range(grid_size))
    ax_grid.spines[:].set_visible(False)
    ax_grid.tick_params(length=0)

    # Legend
    legend_patches = [
        mpatches.Patch(facecolor="#44AA44", alpha=0.75, label="Free (independent)"),
        mpatches.Patch(facecolor="#EE4444", alpha=0.75, label="Eliminated (dependent)"),
    ]
    ax_grid.legend(handles=legend_patches, loc="upper right", frameon=True, fontsize=10)

    # Right panel: substitution formulas as text
    ax_sub = axes[1]
    ax_sub.set_axis_off()
    ax_sub.set_title("Substitution Rules\n(eliminated = linear combination of free)", fontsize=11)

    n_sub = len(substitution)
    if n_sub == 0:
        ax_sub.text(0.5, 0.5, "No eliminations needed.",
                    transform=ax_sub.transAxes, ha="center", va="center",
                    fontsize=12, color="#888888")
    else:
        y_step = min(0.08, 0.88 / n_sub)
        y = 0.94

        for elim_idx, sub in sorted(substitution.items(),
                                     key=lambda x: (x[0][0] + x[0][1], x[0][1])):
            p_e, q_e = elim_idx
            lhs = rf"$W_{{{p_e},{q_e}}} =$"

            rhs_parts = []
            for (p_f, q_f), coeff in sorted(sub.items(),
                                             key=lambda x: (x[0][0] + x[0][1], x[0][1])):
                if coeff.denominator == 1:
                    c_str = str(int(coeff))
                else:
                    c_str = rf"\frac{{{coeff.numerator}}}{{{coeff.denominator}}}"
                rhs_parts.append(rf"${c_str}\,W_{{{p_f},{q_f}}}$")

            # Wrap long lines
            rhs = ("  " + " + ".join(rhs_parts)) if rhs_parts else "  0"

            if y < 0.04:
                ax_sub.text(0.05, y + y_step, f"… {n_sub} total …",
                            transform=ax_sub.transAxes, ha="left", va="center",
                            fontsize=9, color="#888888")
                break

            ax_sub.text(0.05, y, lhs, transform=ax_sub.transAxes,
                        ha="left", va="top", fontsize=10, color="#CC2222")
            ax_sub.text(0.30, y, rhs if len(rhs) <= 80 else rhs[:77] + "…",
                        transform=ax_sub.transAxes,
                        ha="left", va="top", fontsize=9.5, color="#1155AA")

            y -= y_step

    fig.tight_layout(rect=[0, 0, 1, 0.95])
    _save(fig, output_path)


# ---------------------------------------------------------------------------
# 6. Bound comparison bar chart
# ---------------------------------------------------------------------------

def plot_bound_comparison(
    output_path: str = "outputs/bound_comparison.jpg",
    max_order: int = 3,
    m_sq: Fraction = Fraction(1),
    max_spin: int = 6,
) -> None:
    """
    Plot a comparison of the *objective function coefficients* for each
    Wilson coefficient ratio across three constraint configurations:
      - Positivity only (no null constraints)
      - With s↔u crossing null constraints
      - With full S₃ crossing null constraints

    Since running SDPB requires the compiled binary, this function instead
    visualises the *structure* of the bounds problem: how many free variables
    remain, and what the objective direction (minimise/maximise) is.

    Additionally, it shows the reduction in degrees of freedom graphically
    as a proxy for how much the null constraints tighten the allowed region.

    Parameters
    ----------
    output_path : str
        Output file path.
    max_order : int
        Maximum expansion order.
    m_sq : Fraction
        Mass squared.
    max_spin : int
        Maximum spin for partial-wave expansion.
    """
    from .null_constraints import (
        get_null_constraints,
        get_crossing_symmetric_null_constraints,
        eliminate_variables,
    )

    all_indices = []
    for total in range(max_order + 1):
        for p in range(total + 1):
            q = total - p
            all_indices.append((p, q))

    total_vars = len(all_indices)

    # Count free variables for each configuration
    configs = [
        ("Positivity\nonly", None, "#4477CC"),
        (r"$s\leftrightarrow u$" + "\ncrossing", "su", "#EE7722"),
        ("Full S₃\ncrossing", "full", "#33AA55"),
    ]

    free_counts = []
    elim_counts = []
    for _, ctype, _ in configs:
        if ctype is None:
            free_counts.append(total_vars)
            elim_counts.append(0)
        elif ctype == "su":
            c = get_null_constraints(max_order, m_sq)
            free, sub = eliminate_variables(c, all_indices)
            free_counts.append(len(free))
            elim_counts.append(len(sub))
        else:
            c = get_crossing_symmetric_null_constraints(max_order, m_sq)
            free, sub = eliminate_variables(c, all_indices)
            free_counts.append(len(free))
            elim_counts.append(len(sub))

    # Objective ratios to compare
    ratio_labels = [
        (r"$W_{1,0}/W_{0,0}$", (1, 0), (0, 0)),
        (r"$W_{2,0}/W_{0,0}$", (2, 0), (0, 0)),
        (r"$W_{0,1}/W_{0,0}$", (0, 1), (0, 0)),
        (r"$W_{2,0}/W_{1,0}$", (2, 0), (1, 0)),
    ]

    fig = plt.figure(figsize=_FIG_SIZE)
    gs = gridspec.GridSpec(2, 2, figure=fig, hspace=0.45, wspace=0.38)

    fig.suptitle(
        "Impact of Null Constraints on Wilson Coefficient Bounds\n"
        "(Reduction in Degrees of Freedom)",
        fontsize=14,
    )

    # Panel A: Free variables remaining
    ax_dof = fig.add_subplot(gs[0, :])
    x_pos = np.arange(len(configs))
    bars_free = ax_dof.bar(x_pos, free_counts, color=[c[2] for c in configs],
                           alpha=0.85, edgecolor="white", linewidth=1.5, width=0.55)
    bars_elim = ax_dof.bar(x_pos, elim_counts, bottom=free_counts,
                           color="lightgrey", alpha=0.7, edgecolor="white",
                           linewidth=1, width=0.55, label="Eliminated variables")

    for xi, (fc, ec, cfg) in enumerate(zip(free_counts, elim_counts, configs)):
        ax_dof.text(xi, fc - 0.4, str(fc), ha="center", va="top",
                    color="white", fontsize=12, fontweight="bold")
        if ec > 0:
            ax_dof.text(xi, fc + ec / 2.0, f"(-{ec})", ha="center", va="center",
                        color="#444444", fontsize=10)

    ax_dof.set_xticks(x_pos)
    ax_dof.set_xticklabels([c[0] for c in configs], fontsize=12)
    ax_dof.set_ylabel("Number of Wilson coefficient variables")
    ax_dof.set_title(
        rf"Degrees of freedom: {total_vars} total $W_{{p,q}}$ indices (order ≤ {max_order})",
        fontsize=11,
    )
    ax_dof.legend(frameon=False)
    ax_dof.set_ylim(0, total_vars + 2)

    # Panel B–D: For each ratio, show which variables are in the objective
    # under each configuration (free vs. captured by the objective via substitution)
    from .null_constraints import get_null_constraints, get_crossing_symmetric_null_constraints

    for ratio_idx, (rlabel, num_pq, den_pq) in enumerate(ratio_labels):
        row, col = divmod(ratio_idx, 2)
        ax = fig.add_subplot(gs[1, col])

        obj_var_counts = []
        for _, ctype, clr in configs:
            if ctype is None:
                c_list = []
            elif ctype == "su":
                c_list = get_null_constraints(max_order, m_sq)
            else:
                c_list = get_crossing_symmetric_null_constraints(max_order, m_sq)

            free, sub = eliminate_variables(c_list, all_indices)

            # How many free variables does the objective touch?
            if num_pq in sub:
                n_touched = len(sub[num_pq])
            elif num_pq in free:
                n_touched = 1
            else:
                n_touched = 0
            obj_var_counts.append(n_touched)

        ax.bar(x_pos, obj_var_counts, color=[c[2] for c in configs],
               alpha=0.85, edgecolor="white", linewidth=1, width=0.55)
        for xi, cnt in enumerate(obj_var_counts):
            ax.text(xi, cnt + 0.05, str(cnt), ha="center", va="bottom",
                    fontsize=11, fontweight="bold", color="#333333")

        ax.set_xticks(x_pos)
        ax.set_xticklabels([c[0] for c in configs], fontsize=9)
        ax.set_ylabel("Free vars in objective")
        ax.set_title(rlabel, fontsize=12)
        ax.set_ylim(0, max(obj_var_counts) + 1.5)

    fig.tight_layout(rect=[0, 0, 1, 0.93])
    _save(fig, output_path)


# ---------------------------------------------------------------------------
# 7. 2D allowed-region plot
# ---------------------------------------------------------------------------

def plot_allowed_region_2d(
    output_path: str = "outputs/allowed_region_2d.jpg",
    max_order: int = 3,
    m_sq: Fraction = Fraction(1),
    x_index: Tuple[int, int] = (1, 0),
    y_index: Tuple[int, int] = (2, 0),
    n_points: int = 400,
) -> None:
    """
    Visualise the 2D projected allowed region for two Wilson coefficient ratios.

    Since solving SDPB is needed for exact bounds, this function instead
    draws:
      (a) The trivial positivity-only allowed region (quadrant W_{p,q} ≥ 0)
      (b) Sample crossing-symmetric amplitudes (polynomial toy models)
          that satisfy the null constraints
      (c) Analytic constraints derived by setting the substitution rules
          (linear hyperplanes in the 2D space)

    Parameters
    ----------
    output_path : str
        Output file path.
    max_order : int
        Maximum order for null constraints.
    m_sq : Fraction
        Mass squared.
    x_index, y_index : (int, int)
        The two Wilson coefficients to plot (normalised by W_{0,0}).
    n_points : int
        Number of sample points for boundary curves.
    """
    from .null_constraints import (
        get_null_constraints,
        get_crossing_symmetric_null_constraints,
        eliminate_variables,
    )

    all_indices = []
    for total in range(max_order + 1):
        for p in range(total + 1):
            q = total - p
            all_indices.append((p, q))

    fig, axes = plt.subplots(1, 2, figsize=_FIG_SIZE)
    fig.suptitle(
        rf"Wilson Coefficient Space: $W_{{{x_index[0]},{x_index[1]}}}$ vs "
        rf"$W_{{{y_index[0]},{y_index[1]}}}$  (normalised by $W_{{0,0}}=1$)",
        fontsize=13,
    )

    xlims = (-2, 6)
    ylims = (-2, 8)

    for ax_idx, (ax, crossing_type) in enumerate(zip(axes, ["su", "full"])):
        if crossing_type == "su":
            constraints = get_null_constraints(max_order, m_sq)
            ctype_label = r"$s\leftrightarrow u$ crossing"
        else:
            constraints = get_crossing_symmetric_null_constraints(max_order, m_sq)
            ctype_label = "Full S₃ crossing"

        free, sub = eliminate_variables(constraints, all_indices)

        # Shade positivity-only region (trivial: all W_{p,q} ≥ 0 means first quadrant)
        xx = np.linspace(xlims[0], xlims[1], 300)
        yy = np.linspace(ylims[0], ylims[1], 300)

        ax.fill_between(xx, np.zeros_like(xx), ylims[1],
                        color="#AADDFF", alpha=0.3, label="Positivity only ($W_{pq}≥0$)")

        # Draw null-constraint hyperplanes that constrain the 2D slice
        # For each constraint, if it involves exactly x_index or y_index,
        # draw the implied boundary line.
        for k, constraint in enumerate(constraints):
            has_x = x_index in constraint
            has_y = y_index in constraint

            if not (has_x or has_y):
                continue

            # Get constraints involving both x and y
            cx = float(constraint.get(x_index, Fraction(0)))
            cy = float(constraint.get(y_index, Fraction(0)))
            # Sum of all other terms (treated as constants on the plot)
            other = sum(
                float(c)
                for pq, c in constraint.items()
                if pq not in (x_index, y_index)
            )
            # c_x * x + c_y * y + other = 0
            # If both cx and cy are non-zero:
            if abs(cx) > 1e-10 and abs(cy) > 1e-10:
                # y = -(cx * x + other) / cy
                x_line = np.linspace(xlims[0], xlims[1], 200)
                y_line = -(cx * x_line + other) / cy
                ax.plot(x_line, y_line, "--", linewidth=1.2,
                        color=f"C{k % 10}", alpha=0.7,
                        label=rf"Constraint $k={k+1}$")

        # Generate sample crossing-symmetric points (toy amplitudes)
        np.random.seed(42)
        sample_x = []
        sample_y = []
        n_norm_free = 50  # sample free-parameter values
        for _ in range(n_norm_free * 20):
            # Random free-parameter vector
            free_vals = {idx: np.random.uniform(-1, 3) for idx in free}
            free_vals[(0, 0)] = 1.0  # normalise

            # Compute x_index and y_index via substitution if eliminated
            def _get_w(pq):
                if pq in sub:
                    return sum(free_vals.get(f, 0.0) * float(c)
                               for f, c in sub[pq].items())
                return free_vals.get(pq, 0.0)

            w_x = _get_w(x_index)
            w_y = _get_w(y_index)

            # Apply positivity filter: all free and derived coefficients ≥ 0
            positivity_ok = True
            for idx in all_indices:
                w_val = _get_w(idx)
                if w_val < -0.5:  # loose filter
                    positivity_ok = False
                    break

            if positivity_ok and xlims[0] <= w_x <= xlims[1] and ylims[0] <= w_y <= ylims[1]:
                sample_x.append(w_x)
                sample_y.append(w_y)

            if len(sample_x) >= n_norm_free:
                break

        if sample_x:
            ax.scatter(sample_x, sample_y, s=20, c="#FF6633", alpha=0.6,
                       zorder=5, label="Sample crossing-sym. points")

        ax.axhline(0, color="black", linewidth=0.8)
        ax.axvline(0, color="black", linewidth=0.8)
        ax.set_xlim(xlims)
        ax.set_ylim(ylims)
        ax.set_xlabel(rf"$W_{{{x_index[0]},{x_index[1]}}}$")
        ax.set_ylabel(rf"$W_{{{y_index[0]},{y_index[1]}}}$")
        ax.set_title(f"Constraints: {ctype_label}", fontsize=11)
        ax.legend(loc="upper left", frameon=True, fontsize=8, ncol=1)

    fig.tight_layout(rect=[0, 0, 1, 0.93])
    _save(fig, output_path)


# ---------------------------------------------------------------------------
# 8. Comprehensive summary figure
# ---------------------------------------------------------------------------

def plot_summary_dashboard(
    output_path: str = "outputs/summary_dashboard.jpg",
    max_order: int = 3,
    m_sq: Fraction = Fraction(1),
) -> None:
    """
    Produce a single 800×800 overview figure with four panels:
      - Gegenbauer polynomials (top-left)
      - Null constraint heatmap (top-right)
      - Variable elimination grid (bottom-left)
      - Bound structure bar chart (bottom-right)

    Parameters
    ----------
    output_path : str
        Output path.
    max_order : int
        Maximum expansion order.
    m_sq : Fraction
        Mass squared.
    """
    from .physics import gegenbauer_coefficients, _taylor_coeff_at_one
    from .null_constraints import (
        get_null_constraints,
        get_crossing_symmetric_null_constraints,
        eliminate_variables,
    )

    fig = plt.figure(figsize=_FIG_SIZE)
    gs = gridspec.GridSpec(2, 2, figure=fig, hspace=0.42, wspace=0.38)
    fig.suptitle(
        rf"EFT Wilson Coefficient Bounds — Summary Dashboard (order ≤ {max_order})",
        fontsize=14, fontweight="bold",
    )

    # ---- Panel A: Gegenbauer polynomials ----
    ax_geg = fig.add_subplot(gs[0, 0])
    z = np.linspace(-1.0, 1.0, 300)
    for i, ell in enumerate([0, 2, 4]):
        coeffs = gegenbauer_coefficients(ell, d=4)
        c_f = [_fraction_to_float(c) for c in coeffs]
        y_geg = np.polyval(c_f[::-1], z)
        ax_geg.plot(z, y_geg, color=_COLORS[i], linewidth=2,
                    label=rf"$\ell={ell}$")
    ax_geg.axhline(0, color="black", linewidth=0.7, linestyle="--")
    ax_geg.set_title(r"Gegenbauer $C_\ell^{(1)}(z)$", fontsize=11)
    ax_geg.set_xlabel(r"$z$", fontsize=10)
    ax_geg.legend(frameon=False, fontsize=9)

    # ---- Panel B: Null constraint heatmap ----
    ax_heat = fig.add_subplot(gs[0, 1])
    constraints = get_null_constraints(max_order, m_sq)
    all_pq = set()
    for c in constraints:
        all_pq.update(c.keys())
    sorted_pq = sorted(all_pq, key=lambda x: (x[0] + x[1], x[1]))
    n_k = len(constraints)
    n_pq = len(sorted_pq)
    pq_to_col = {pq: i for i, pq in enumerate(sorted_pq)}
    mat = np.zeros((n_k, n_pq))
    for k, constraint in enumerate(constraints):
        for pq, coeff in constraint.items():
            mat[k, pq_to_col[pq]] = float(coeff)
    vmax = np.abs(mat).max() or 1.0
    ax_heat.imshow(mat, aspect="auto", cmap=_BWR, vmin=-vmax, vmax=vmax)
    ax_heat.set_xticks(range(n_pq))
    ax_heat.set_xticklabels(
        [rf"$W_{{{p},{q}}}$" for p, q in sorted_pq],
        rotation=45, ha="right", fontsize=7,
    )
    ax_heat.set_yticks(range(n_k))
    ax_heat.set_yticklabels([rf"$k={k+1}$" for k in range(n_k)], fontsize=7)
    ax_heat.set_title(
        rf"Null constraint matrix ($s\leftrightarrow u$, n={n_k})", fontsize=10
    )

    # ---- Panel C: Free vs. eliminated ----
    ax_elim = fig.add_subplot(gs[1, 0])
    all_indices = []
    for total in range(max_order + 1):
        for p in range(total + 1):
            q = total - p
            all_indices.append((p, q))
    free, sub = eliminate_variables(constraints, all_indices)

    configs = ["Positivity\nonly", r"$s\leftrightarrow u$", "Full S₃"]
    free_ns, elim_ns = [], []
    for cfg_type in [None, "su", "full"]:
        if cfg_type is None:
            free_ns.append(len(all_indices)); elim_ns.append(0)
        elif cfg_type == "su":
            c = get_null_constraints(max_order, m_sq)
            fr, sb = eliminate_variables(c, all_indices)
            free_ns.append(len(fr)); elim_ns.append(len(sb))
        else:
            c = get_crossing_symmetric_null_constraints(max_order, m_sq)
            fr, sb = eliminate_variables(c, all_indices)
            free_ns.append(len(fr)); elim_ns.append(len(sb))

    x_pos = np.arange(len(configs))
    ax_elim.bar(x_pos, free_ns, color=["#4477CC", "#EE7722", "#33AA55"],
                alpha=0.85, edgecolor="white", width=0.55)
    ax_elim.bar(x_pos, elim_ns, bottom=free_ns, color="lightgrey",
                alpha=0.7, edgecolor="white", width=0.55)
    for xi, (fn, en) in enumerate(zip(free_ns, elim_ns)):
        ax_elim.text(xi, fn / 2.0, str(fn), ha="center", va="center",
                     color="white", fontsize=10, fontweight="bold")
    ax_elim.set_xticks(x_pos)
    ax_elim.set_xticklabels(configs, fontsize=9)
    ax_elim.set_ylabel("# variables")
    ax_elim.set_title("Degrees of freedom per config.", fontsize=10)

    # ---- Panel D: Spectral function for ℓ=0,2 ----
    ax_spec = fig.add_subplot(gs[1, 1])
    x_spec = np.linspace(0.05, 8.0, 300)
    mu_th = float(4 * m_sq)
    for i, ell in enumerate([0, 2, 4]):
        coeffs = gegenbauer_coefficients(ell, d=4)
        # v_ℓ^{(1,0)}(x) = C_ℓ(1) / (x + 4m²)^2   [q=0, p=1]
        c0 = _fraction_to_float(_taylor_coeff_at_one(coeffs, 0))
        v = c0 / ((x_spec + mu_th) ** 2)
        ax_spec.plot(x_spec, v, color=_COLORS[i], linewidth=2,
                     label=rf"$\ell={ell}$")
    ax_spec.set_xlabel(r"$x = \mu - 4m^2$", fontsize=10)
    ax_spec.set_ylabel(r"$v_\ell^{(1,0)}(x)$", fontsize=10)
    ax_spec.set_title(r"Spectral fn. $v_\ell^{(p=1,q=0)}$", fontsize=10)
    ax_spec.legend(frameon=False, fontsize=9)

    fig.tight_layout(rect=[0, 0, 1, 0.93])
    _save(fig, output_path)
