# EFT Wilson Coefficient Bounds via SDPB

## Overview

This package computes **lower bounds on Wilson coefficients** in Effective Field
Theories (EFTs) for 2→2 identical scalar scattering, using:

- **Spectral positivity** from unitarity (partial-wave expansion)
- **Null constraints** from crossing-symmetric dispersion relations
  (Sinha & Zahed, "Crossing Symmetric Dispersion Relations in QFTs", 2021)
- **SDPB** semidefinite program solver

The approach follows the dual problem formulation of **"Extremal Effective Field
Theories"** (Caron-Huot & Duong, 2021).

## Quick Start

### 1. Generate PMP input files

```bash
# Run from the repository root
python -m eft_bounds.run_bounds \
    -o eft_bounds_output \
    --max-spin 10 \
    --max-order 4 \
    --precision 200
```

This generates PMP JSON files for various Wilson coefficient ratios,
with three configurations each:
- **positivity_only**: No crossing constraints (baseline)
- **with_su_crossing**: s↔u null constraints from Sinha-Zahed
- **with_full_crossing**: Full S₃ crossing (s↔u + s↔t) constraints

### 1b. Generate the exact crossing-basis automation pipeline

To derive the full-crossing null constraints from the exact
crossing-symmetric basis

\[
M(s,t) = \sum_{p,q} W_{p,q}\,x^p y^q,\qquad
x=-(st+tu+us),\quad y=-stu,
\]

expand them into ordinary `(s,t)` polynomials, export JSON artifacts, and
generate paired-inequality PMP files plus lower/upper bound PMP examples, run:

```bash
python -m eft_bounds.crossing_pipeline \
    --output-dir /tmp/eft_crossing_pipeline \
    --max-degree 6 \
    --mass-squared 1 \
    --precision 80 \
    --max-spin 10
```

This writes:

- `crossing_to_ordinary.json` — exact basis expansion of each `x^p y^q`
- `null_constraints.json` — the full list of exact null constraints on ordinary coefficients
- `null_constraints_pmp.json` — the same constraints encoded as paired SDPB PMP blocks
- `eq3_to_eq23_translation.json` — exact map from the Sinha-Zahed eq.(3) basis to the
  Caron-Huot/Duong eq.(2.3) basis
- `lower_*.json`, `upper_*.json` — example two-sided SDPB PMP files
- `reliability_checks.json` — verification results on at least five examples

### 2. Convert and solve with SDPB

```bash
# For each PMP file:
pmp2sdp --precision=1024 \
    --input=eft_bounds_output/pmp_g10_over_g00_with_su_crossing.json \
    --output=eft_bounds_output/sdp_g10/

mpirun -n 4 sdpb --precision=1024 \
    -s eft_bounds_output/sdp_g10/ \
    -o eft_bounds_output/out_g10/
```

### 3. Read the bound

The optimal objective value gives the **negative** of the lower bound:
```
lower_bound = -(primalObjective from out.txt)
```

## Physics Background

### Amplitude Expansion

For 2→2 scattering of identical scalars with mass $m$, the amplitude is
expanded as:

$$\mathcal{M}(s,t) = \sum_{p,q \geq 0} W_{pq} \, s^p \, t^q$$

where $s, t$ are Mandelstam variables and $u = 4m^2 - s - t$.

### Spectral Decomposition

Unitarity gives the Wilson coefficients as positive sums over the spectral
density:

$$W_{pq} = \sum_\ell \int_{4m^2}^{\infty} d\mu \, \rho_\ell(\mu) \, v_\ell^{(p,q)}(\mu)$$

where $\rho_\ell(\mu) \geq 0$ and $v_\ell^{(p,q)}(\mu)$ are known
kinematic functions involving Gegenbauer polynomials.

### Null Constraints

Crossing symmetry $\mathcal{M}(s,t) = \mathcal{M}(u,t)$ imposes linear
relations among the Wilson coefficients:

$$\sum_{p,q} C_{pq}^{(k)} W_{pq} = 0, \quad k = 1, 2, \ldots$$

These **null constraints** reduce the number of independent Wilson
coefficients.

### SDPB Formulation

The bound on $W_{\text{num}} / W_{\text{den}}$ is computed by solving:

$$\text{maximize} \quad -W_{\text{num}}$$
$$\text{subject to} \quad \sum_n z_n \, v_\ell^{(n)}(x) \geq 0 \quad \forall x \geq 0, \;\forall \ell$$
$$\text{and} \quad W_{\text{den}} = 1$$

## Method for Incorporating Null Constraints

### Method A: Variable Elimination (Used)

Null constraints are used to express dependent Wilson coefficients in
terms of independent ones. This:
- Reduces the problem dimensionality
- Maintains mathematical exactness
- Avoids numerical issues from enforcing equalities as paired inequalities

The spectral functions $v_\ell^{(n)}(x)$ are combined accordingly:

$$v_{\ell,\text{eff}}^{(k)}(x) = v_\ell^{(z_k)}(x) + \sum_{\text{eliminated}} \text{sub}[\text{elim}][z_k] \cdot v_\ell^{(\text{elim})}(x)$$

### Method B: Paired Inequalities (Alternative)

Each null constraint can also be encoded as two opposing inequality
constraints:
- $+\sum C_{pq}^{(k)} W_{pq} \geq 0$
- $-\sum C_{pq}^{(k)} W_{pq} \geq 0$

These appear as 1×1 constant polynomial matrix blocks in the PMP.
This is less efficient but useful for cross-checking.

### Exact full-crossing derivation

The new `eft_bounds.crossing_pipeline` module derives the null constraints
without guessing relations directly in the ordinary `(s,t)` basis:

1. enumerate the symmetric basis monomials `x^p y^q` up to the chosen EFT cutoff,
2. expand each monomial exactly into ordinary `s^a t^b` coefficients using
   `u = 4m² - s - t`,
3. compute the exact left-nullspace of that basis-expansion matrix,
4. export the resulting null constraints and optional PMP encodings.

This is the recommended route when you want the exact full-crossing relations
associated with the Sinha-Zahed eq.(3) basis.

## Comparison with Extremal EFT Results

### Expected Behavior

| Configuration | Expected bounds |
|---------------|----------------|
| Positivity only | Trivial ($W_{pq} \geq 0$) |
| With finite-order null constraints | Nontrivial but weaker than Extremal EFT |
| With all null constraints (infinite) | Should approach Extremal EFT bounds |
| Extremal EFT (full crossing) | Optimal bounds |

### Why Bounds Differ

The null constraints capture crossing symmetry **order by order** in the
low-energy expansion. The Extremal EFT approach imposes crossing at the
**full spectral level**. The difference is analogous to:

- Null constraints ↔ matching Taylor coefficients (local information)
- Full crossing ↔ matching the full analytic function (global information)

### Key Verification Points

1. **Monotonicity**: Adding more null constraints should only tighten bounds
2. **Consistency**: Bounds should never be tighter than Extremal EFT values
3. **Known amplitudes**: Scalar exchange should satisfy all constraints
4. **Polynomial test**: Finite-order crossing-symmetric amplitudes (like
   $s^2 + u^2$) should exactly satisfy all null constraints

### Truncation Effects

The null constraints derived from $\mathcal{M}(s,t) = \mathcal{M}(u,t)$
involve all Wilson coefficients in principle. When truncated to a finite
order $p+q \leq K$, the constraints are:

- **Exact** for polynomial amplitudes (finitely many nonzero $W_{pq}$)
- **Approximate** for amplitudes with infinitely many terms (e.g., scalar
  exchange), with truncation errors from higher-order contributions

For the SDPB optimization, where Wilson coefficients are the decision
variables, the truncated constraints are valid constraints on the
finite-dimensional parameter space being optimized over.

## File Structure

```
eft_bounds/
├── __init__.py           # Package initialization
├── physics.py            # Gegenbauer polynomials, spectral functions
├── null_constraints.py   # Crossing-symmetric null constraints
├── pmp_generator.py      # PMP JSON file generation for SDPB
├── run_bounds.py         # Main driver script
└── README.md             # This file
```

## Module Documentation

### `physics.py`
- `gegenbauer_coefficients(ell, d)` — Gegenbauer polynomial $C_\ell^{(d/2-1)}(z)$
- `build_positivity_polynomials(ell, indices, ...)` — Spectral functions for SDPB
- `fraction_to_str(f, precision)` — High-precision decimal string output

### `null_constraints.py`
- `get_null_constraints(max_order, m_sq)` — s↔u crossing null constraints
- `get_crossing_symmetric_null_constraints(max_order, m_sq)` — Full S₃ crossing
- `eliminate_variables(constraints, indices)` — Variable elimination (Method A)

### `crossing_pipeline.py`
- `derive_null_constraints(max_degree, m_sq)` — exact full-crossing null constraints from the `x,y` basis
- `build_crossing_to_ordinary_map(...)` — exact expansion of `x^p y^q` into ordinary polynomials
- `build_translation_metadata(...)` — exact map from Sinha-Zahed eq.(3) coefficients to Extremal-EFT eq.(2.3) coefficients
- `generate_two_sided_bound_pmps(...)` — lower/upper SDPB PMP files for one symmetric-basis ratio
- `run_reliability_checks(...)` — exact checks on at least five amplitudes

### `pmp_generator.py`
- `generate_pmp_json(objective_index, ...)` — Generate PMP for single coefficient
- `generate_pmp_for_ratio_bound(num, den, ...)` — Generate PMP for ratio bound
- `write_pmp_json(pmp, filepath)` — Write PMP to JSON file

### `run_bounds.py`
- `run_consistency_checks()` — Validate implementation correctness
- `generate_bound_series(output_dir, ...)` — Generate all PMP files
- `print_comparison_analysis()` — Print comparison with Extremal EFT

## References

1. S. Caron-Huot & V.-D. Duong, "Extremal Effective Field Theories",
   JHEP 05 (2021) 280, [arXiv:2011.02957]

2. A. Sinha & A. Zahed, "Crossing Symmetric Dispersion Relations in QFTs",
   Phys. Rev. Lett. 126 (2021) 181601, [arXiv:2012.04877]

3. D. Simmons-Duffin, "A Semidefinite Program Solver for the Conformal
   Bootstrap", JHEP 06 (2015) 174, [arXiv:1502.02033]
