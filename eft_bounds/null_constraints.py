"""
null_constraints.py — Null constraints from crossing-symmetric dispersion relations.

Implements the null constraints derived in:
  Sinha & Zahed, "Crossing Symmetric Dispersion Relations in QFTs" (2021)

These are exact linear relations among Wilson coefficients W_{pq} that
follow from crossing symmetry of the 2→2 scattering amplitude of
identical scalars. They take the form:
    Σ_{p,q} C_{pq}^{(k)} W_{pq} = 0

for each null constraint index k.

Convention:
  The amplitude is expanded as M(s,t,u) = Σ_{p,q≥0} W_{pq} x^p y^q
  where (using the crossing-symmetric variables from Sinha-Zahed):
    x = -(st + tu + us)
    y = -stu
  or equivalently, the amplitude is expanded around the crossing-symmetric
  point s = t = u = 4m²/3 in suitable variables.

  For the simpler fixed-t dispersion relation approach, the Wilson coefficients
  are defined via:
    M(s,t) = Σ_{p,q≥0} g_{p,q} s^p t^q

  The crossing symmetry s ↔ u (where u = 4m² - s - t) imposes:
    M(s,t) = M(u,t) = M(4m² - s - t, t)

  This generates null constraints order by order in the low-energy expansion.

Reference equations from Sinha & Zahed (2021):
  - Eq. (3.10)-(3.20): Null constraints at various orders
  - Table 1: Summary of null constraints

The null constraints used here are derived from the crossing relation
  M(s,t) = M(u,t) with u = 4m² - s - t
for identical scalar scattering in d=4.
"""

from fractions import Fraction
from typing import Dict, List, Tuple


def get_null_constraints(
    max_order: int = 6,
    m_sq: Fraction = Fraction(1),
) -> List[Dict[Tuple[int, int], Fraction]]:
    """
    Compute null constraints from crossing symmetry.

    For identical scalar scattering with M(s,t) = M(u,t), u = 4m² - s - t,
    the amplitude must be symmetric under s ↔ u.

    Expanding M(s,t) = Σ_{p,q} a_{pq} s^p t^q, crossing symmetry s ↔ u gives:
      Σ_{p,q} a_{pq} s^p t^q = Σ_{p,q} a_{pq} (4m² - s - t)^p t^q

    Matching coefficients of s^p t^q on both sides generates linear relations
    among the a_{pq}.

    The crossing condition is that M(s,t) - M(4m²-s-t, t) = 0 identically.

    Expanding (4m² - s - t)^p in powers of s and t using the binomial theorem:
      (4m² - s - t)^p = Σ_{j+k≤p} (-1)^{j+k} C(p,j,k) (4m²)^{p-j-k} s^j t^k

    where C(p,j,k) = p! / (j! k! (p-j-k)!) is the multinomial coefficient.

    The null constraints come from requiring that the coefficient of each
    s^a t^b monomial vanishes in M(s,t) - M(u,t) = 0.

    For the coefficient of s^a t^b:
      a_{a,b} - Σ_{p≥a, q≤b} a_{p,q} × [coeff of s^a t^{b-q} in (4m²-s-t)^p]
    must vanish.

    The terms where p is ODD in (the expansion of (-1)^p) give nontrivial
    null constraints.

    Parameters
    ----------
    max_order : int
        Maximum total order p+q to include (default 6).
    m_sq : Fraction
        Mass squared (default 1).

    Returns
    -------
    list of dict
        Each dict maps (p, q) -> coefficient, representing
        Σ_{p,q} coeff × a_{pq} = 0.
    """
    mu = 4 * m_sq  # threshold = 4m²
    constraints = []

    # Crossing symmetry: M(s,t) = M(4m²-s-t, t)
    # We expand (4m²-s-t)^p and collect monomials s^a t^b.
    #
    # The crossing relation gives, for each monomial s^a t^b:
    #   a_{a,b} = Σ_p Σ_q [coeff of s^a t^{b-q} in (4m²-s-t)^p] × a_{p,q}
    #
    # This is an identity if a_{a,b} appears on both sides consistently.
    # The null constraints arise when terms with ODD powers of s have
    # a sign flip under s → -s (in the shifted variable).
    #
    # More precisely, let's define σ = s - 4m²/3 (shift to crossing-symmetric
    # point for stu-symmetric case), but for s↔u symmetry it's simpler:
    #
    # Under s ↔ u = 4m² - s - t:
    #   s → 4m² - s - t
    # So s changes sign (up to shifts).
    #
    # The simplest approach: expand M(s,t) in powers of s and t,
    # and impose M(s,t) = M(4m²-s-t, t) order by order.

    # Build the crossing matrix: for each target monomial s^a t^b,
    # what is the coefficient in terms of the original a_{p,q}?
    #
    # M(u,t) = Σ_{p,q} a_{p,q} (4m²-s-t)^p t^q
    # The coefficient of s^a t^b in (4m²-s-t)^p t^q is:
    #   [coeff of s^a t^{b-q} in (4m²-s-t)^p]
    # = [coeff of s^a t^{b-q} in Σ_{j+k≤p} C(p;j,k) (4m²)^{p-j-k} (-s)^j (-t)^k]
    # = (-1)^{a} (-1)^{b-q} × C(p; a, b-q) × (4m²)^{p-a-(b-q)}
    #   if a + (b-q) ≤ p, else 0.
    #
    # where C(p; j, k) = p! / (j! k! (p-j-k)!)

    for total_order in range(1, max_order + 1):
        for a in range(total_order + 1):
            b = total_order - a
            if b < 0:
                continue

            # Coefficient of s^a t^b in M(s,t) is simply a_{a,b} (coefficient 1).
            # Coefficient of s^a t^b in M(u,t) = Σ_{p,q} a_{p,q} × T(p,q,a,b)
            # where T(p,q,a,b) = coefficient of s^a t^{b-q} in (4m²-s-t)^p
            #                  = (-1)^{a+b-q} × multinomial(p; a, b-q) × μ^{p-a-b+q}

            constraint = {}

            for p in range(max_order + 1):
                for q in range(max_order + 1):
                    if p + q > max_order:
                        continue

                    # Direct coefficient: δ_{(a,b),(p,q)} from M(s,t)
                    direct = Fraction(1) if (p == a and q == b) else Fraction(0)

                    # Crossed coefficient from M(u,t):
                    r = b - q  # power of t from (4m²-s-t)^p needed
                    if r < 0 or a + r > p:
                        crossed = Fraction(0)
                    else:
                        sign = (-1) ** (a + r)
                        multi = _multinomial(p, a, r)
                        mu_power = mu ** (p - a - r)
                        crossed = Fraction(sign) * multi * mu_power

                    net = direct - crossed
                    if net != 0:
                        constraint[(p, q)] = net

            if constraint and len(constraint) > 1:
                # Only keep nontrivial constraints (more than one term)
                constraints.append(constraint)

    # Remove duplicate/proportional constraints
    constraints = _remove_proportional(constraints)

    return constraints


def _multinomial(n: int, j: int, k: int) -> Fraction:
    """
    Multinomial coefficient n! / (j! k! (n-j-k)!)
    """
    if j < 0 or k < 0 or j + k > n:
        return Fraction(0)
    from math import factorial
    return Fraction(factorial(n), factorial(j) * factorial(k) * factorial(n - j - k))


def _remove_proportional(
    constraints: List[Dict[Tuple[int, int], Fraction]]
) -> List[Dict[Tuple[int, int], Fraction]]:
    """Remove constraints that are proportional to each other."""
    unique = []
    for c in constraints:
        if not c:
            continue
        is_proportional = False
        for u in unique:
            if _are_proportional(c, u):
                is_proportional = True
                break
        if not is_proportional:
            unique.append(c)
    return unique


def _are_proportional(
    c1: Dict[Tuple[int, int], Fraction],
    c2: Dict[Tuple[int, int], Fraction],
) -> bool:
    """Check if two constraints are proportional."""
    if set(c1.keys()) != set(c2.keys()):
        return False
    if not c1:
        return True

    # Find ratio from first nonzero pair
    key0 = next(iter(c1))
    if c2[key0] == 0:
        return c1[key0] == 0
    ratio = c1[key0] / c2[key0]

    for key in c1:
        if c2[key] == 0:
            if c1[key] != 0:
                return False
        elif c1[key] / c2[key] != ratio:
            return False
    return True


def eliminate_variables(
    constraints: List[Dict[Tuple[int, int], Fraction]],
    all_indices: List[Tuple[int, int]],
) -> Tuple[
    List[Tuple[int, int]],
    Dict[Tuple[int, int], Dict[Tuple[int, int], Fraction]],
]:
    """
    Use null constraints to eliminate dependent Wilson coefficients.

    Given null constraints of the form Σ C_{pq}^{(k)} W_{pq} = 0,
    solve for some W_{pq} in terms of others.

    This implements Method A from the plan: variable elimination.

    The algorithm maintains the invariant that all substitutions only
    reference truly free (non-eliminated) variables.

    Parameters
    ----------
    constraints : list of dict
        Null constraints, each mapping (p,q) -> coefficient.
    all_indices : list of (int, int)
        All Wilson coefficient indices (p,q) being considered.

    Returns
    -------
    free_indices : list of (int, int)
        The independent (free) Wilson coefficient indices after elimination.
    substitution : dict
        Maps eliminated index → dict of {free_index: coefficient}.
        So W_{eliminated} = Σ coeff × W_{free}.
        All values in the dict reference only free (non-eliminated) variables.
    """
    remaining_constraints = [dict(c) for c in constraints if c]
    eliminated = set()
    substitution = {}
    available_set = set(all_indices)

    for constraint in remaining_constraints:
        # Step 1: Substitute all previously eliminated variables.
        # Since we maintain the invariant that all substitutions
        # only reference free variables, this produces a reduced
        # constraint purely in terms of free variables.
        reduced = {}
        for idx, coeff in constraint.items():
            if idx in substitution:
                for free_idx, sub_coeff in substitution[idx].items():
                    reduced[free_idx] = reduced.get(free_idx, Fraction(0)) + coeff * sub_coeff
            else:
                reduced[idx] = reduced.get(idx, Fraction(0)) + coeff

        # Remove zero entries
        reduced = {k: v for k, v in reduced.items() if v != 0}

        if not reduced:
            continue

        # Step 2: Pick a free variable to eliminate.
        candidates = [
            idx for idx in reduced
            if idx in available_set and idx not in eliminated
            and reduced[idx] != 0
        ]
        if not candidates:
            continue

        elim_idx = max(candidates, key=lambda idx: (idx[0] + idx[1], idx[1]))
        elim_coeff = reduced[elim_idx]

        # Build the new substitution: W_{elim} = Σ (-coeff/elim_coeff) W_{free}
        new_sub = {}
        for idx, coeff in reduced.items():
            if idx == elim_idx:
                continue
            new_sub[idx] = -coeff / elim_coeff

        # Step 3: Update ALL existing substitutions to remove references
        # to elim_idx (maintaining the invariant).
        for prev_elim in list(substitution.keys()):
            prev_sub = substitution[prev_elim]
            if elim_idx in prev_sub:
                factor = prev_sub.pop(elim_idx)
                for free_idx, free_coeff in new_sub.items():
                    prev_sub[free_idx] = prev_sub.get(free_idx, Fraction(0)) + factor * free_coeff
                # Clean up zeros
                substitution[prev_elim] = {k: v for k, v in prev_sub.items() if v != 0}

        substitution[elim_idx] = new_sub
        eliminated.add(elim_idx)

    free_indices = [idx for idx in all_indices if idx not in eliminated]
    return free_indices, substitution


def print_null_constraints(
    constraints: List[Dict[Tuple[int, int], Fraction]],
) -> None:
    """Print null constraints in human-readable form."""
    for i, c in enumerate(constraints):
        terms = []
        for (p, q), coeff in sorted(c.items()):
            if coeff == 1:
                terms.append(f"W_{{{p},{q}}}")
            elif coeff == -1:
                terms.append(f"-W_{{{p},{q}}}")
            else:
                terms.append(f"({coeff})·W_{{{p},{q}}}")
        print(f"  Null constraint {i+1}: {' + '.join(terms)} = 0")


def get_crossing_symmetric_null_constraints(
    max_order: int = 6,
    m_sq: Fraction = Fraction(1),
) -> List[Dict[Tuple[int, int], Fraction]]:
    """
    Compute a simple full-crossing approximation on ordinary s,t coefficients.

    This helper combines:
      1. the direct s↔u relations derived in `get_null_constraints`, and
      2. the additional s↔t symmetry conditions a_{a,b} = a_{b,a}.

    It is useful as a lightweight approximation in the ordinary polynomial
    basis, but it is not the exact crossing-basis derivation from the
    Sinha-Zahed eq.(3) variables. For the exact derivation from the
    crossing-symmetric basis x = -(st+tu+us), y = -stu, use:

        eft_bounds.crossing_pipeline.derive_null_constraints(...)

    Parameters
    ----------
    max_order : int
        Maximum total order 2p+3q in the expansion.
    m_sq : Fraction
        Mass squared (default 1).

    Returns
    -------
    list of dict
        Each dict maps (p, q) → coefficient for the null constraint.
    """
    # For the crossing-symmetric variable expansion, the null constraints
    # arise from the requirement that the dispersive representation
    # (which uses fixed-t or crossing-symmetric dispersion relations)
    # reproduces the correct analytic structure.
    #
    # From Sinha-Zahed (2021), the key null constraints are:
    #
    # At leading orders in the σ₂, σ₃ expansion:
    #
    # 1) W_{01} relates to W_{10}:
    #    The coefficient of σ₃ in the amplitude is completely fixed by
    #    lower-order coefficients via the dispersion relation.
    #
    # 2) At higher orders, each new order in σ₃ introduces constraints.
    #
    # We implement the explicit null constraints from the paper.
    # These are obtained by expanding the crossing-symmetric dispersion
    # relation and matching coefficients.

    # For the simpler s↔u crossing, we use the function above.
    # Here we add the additional t↔u constraints.
    #
    # The full set of crossing constraints for M(s,t) = M(t,s) (s↔t)
    # gives: for each monomial s^a t^b with a ≠ b,
    #   a_{a,b} = a_{b,a}
    # These are simple symmetry constraints.
    #
    # Combined with s↔u (which we already have), this gives the full S₃.

    constraints_su = get_null_constraints(max_order, m_sq)

    # Add s↔t symmetry constraints: a_{a,b} = a_{b,a} for a ≠ b
    constraints_st = []
    for a in range(max_order + 1):
        for b in range(a + 1, max_order + 1):
            if a + b > max_order:
                continue
            # a_{a,b} - a_{b,a} = 0
            constraints_st.append({(a, b): Fraction(1), (b, a): Fraction(-1)})

    all_constraints = constraints_su + constraints_st
    return _remove_proportional(all_constraints)
