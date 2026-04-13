"""
physics.py — Spectral functions and partial-wave kinematics for EFT bounds.

This module implements the core physics needed to set up the dual problem
for bounding Wilson coefficients in 2→2 identical scalar scattering.

Key objects:
  - Gegenbauer polynomials C_ℓ(cos θ) for partial-wave expansion
  - Spectral functions v_ℓ^{(p,q)}(μ) that give the Wilson coefficient
    contributions from a spin-ℓ state at mass μ
  - Polynomial-in-x representations for SDPB input

References:
  - Caron-Huot & Duong, "Extremal Effective Field Theories" (2021)
  - The spectral decomposition is:
      g_{p,q} = Σ_ℓ ∫ dμ ρ_ℓ(μ) v_ℓ^{(p,q)}(μ)
    where ρ_ℓ(μ) ≥ 0 is the spectral density.
"""

from fractions import Fraction
from functools import lru_cache
from typing import List, Tuple

# We use exact rational arithmetic where possible to avoid precision loss.
# For the final JSON output, we convert to high-precision decimal strings.


# ---------------------------------------------------------------------------
# Gegenbauer polynomials
# ---------------------------------------------------------------------------

@lru_cache(maxsize=256)
def gegenbauer_coefficients(ell: int, d: int = 4) -> List[Fraction]:
    """
    Compute coefficients of the Gegenbauer polynomial C_ℓ^{(d/2-1)}(z)
    as a polynomial in z, returned as [c_0, c_1, ..., c_ℓ].

    For d=4, this is C_ℓ^{(1)}(z) = U_ℓ(z) (Chebyshev of second kind
    normalized so C_ℓ^{(1)}(1) = ℓ+1).

    Uses the recurrence:
      C_0^λ(z) = 1
      C_1^λ(z) = 2λ z
      (n+1) C_{n+1}^λ(z) = 2(n+λ) z C_n^λ(z) - (n+2λ-1) C_{n-1}^λ(z)

    Parameters
    ----------
    ell : int
        Spin (angular momentum quantum number), must be non-negative.
    d : int
        Spacetime dimension (default 4).

    Returns
    -------
    list of Fraction
        Coefficients [c_0, c_1, ..., c_ℓ] such that
        C_ℓ(z) = c_0 + c_1 z + c_2 z² + ... + c_ℓ z^ℓ.
    """
    lam = Fraction(d - 2, 2)  # λ = (d-2)/2

    if ell == 0:
        return [Fraction(1)]
    if ell == 1:
        return [Fraction(0), 2 * lam]

    # Recurrence: (n+1) C_{n+1} = 2(n+λ) z C_n - (n+2λ-1) C_{n-1}
    c_prev = [Fraction(1)]  # C_0
    c_curr = [Fraction(0), 2 * lam]  # C_1

    for n in range(1, ell):
        # Multiply c_curr by z: shift coefficients right
        z_times_curr = [Fraction(0)] + list(c_curr)

        # Pad to same length
        max_len = len(z_times_curr)
        c_prev_padded = list(c_prev) + [Fraction(0)] * (max_len - len(c_prev))

        a_coeff = 2 * (n + lam)
        b_coeff = n + 2 * lam - 1
        denom = n + 1

        c_next = []
        for i in range(max_len):
            val = (a_coeff * z_times_curr[i] - b_coeff * c_prev_padded[i]) / denom
            c_next.append(val)

        c_prev = c_curr
        c_curr = c_next

    return c_curr


def eval_gegenbauer(ell: int, z: Fraction, d: int = 4) -> Fraction:
    """Evaluate Gegenbauer polynomial C_ℓ^{(d/2-1)}(z) at a point z."""
    coeffs = gegenbauer_coefficients(ell, d)
    result = Fraction(0)
    z_power = Fraction(1)
    for c in coeffs:
        result += c * z_power
        z_power *= z
    return result


# ---------------------------------------------------------------------------
# Scattering kinematics for 2→2 identical scalars
# ---------------------------------------------------------------------------

def mandelstam_to_cos_theta(s: Fraction, t: Fraction, m_sq: Fraction) -> Fraction:
    """
    Compute cos(θ) from Mandelstam variables for 2→2 scattering of
    identical scalars with mass² = m_sq.

    cos(θ) = 1 + 2t / (s - 4m²)

    Parameters
    ----------
    s : Fraction
        Center-of-mass energy squared.
    t : Fraction
        Momentum transfer squared.
    m_sq : Fraction
        Mass squared of the external scalars.
    """
    return 1 + 2 * t / (s - 4 * m_sq)


# ---------------------------------------------------------------------------
# Spectral functions v_ℓ^{(p,q)}(μ)
# ---------------------------------------------------------------------------

def spectral_function_coefficients(
    ell: int,
    p: int,
    q: int,
    m_sq: Fraction = Fraction(1),
    d: int = 4,
    x_shift: Fraction = Fraction(0),
) -> List[Fraction]:
    """
    Compute the polynomial coefficients of v_ℓ^{(p,q)}(x) where x = μ - μ_threshold.

    The Wilson coefficient g_{p,q} in the low-energy expansion
        M(s,t) = Σ_{p,q} g_{p,q} s^p t^q
    receives contributions from each partial wave:
        g_{p,q} = Σ_ℓ ∫ dμ ρ_ℓ(μ) v_ℓ^{(p,q)}(μ)

    For fixed-t dispersion relations, the absorptive part of the amplitude
    at mass μ and spin ℓ contributes:
        v_ℓ^{(p,q)}(μ) = (coefficient from expanding μ^{-n} C_ℓ(1 + 2t/(μ-4m²))
                           in powers of s^p t^q)

    We work in units where the threshold is μ_th = 4m², and define x = μ - 4m² ≥ 0.

    The function v_ℓ^{(p,q)}(x) for the spectral decomposition in the
    forward-limit dispersion relation approach is:

    For the s-channel partial wave at mass² = μ = x + 4m²:
      The partial wave contributes to the t-expansion of the amplitude as:
        A_ℓ(s) C_ℓ(cos θ)
      where cos θ = 1 + 2t/(s - 4m²).

    The contribution to g_{p,q} from a delta-function spectral density
    at mass μ is obtained by expanding:
      1/(μ^{p+1}) × [coefficient of t^q in C_ℓ(1 + 2t/(μ - 4m²))]

    In terms of x = μ - 4m²:
      v_ℓ^{(p,q)}(x) = 1/(x + 4m²)^{p+1} × [coeff of t^q in C_ℓ(1 + 2t/x)]

    For the Gegenbauer expansion in t:
      C_ℓ(1 + 2t/x) = Σ_k c_k (2t/x)^k → coefficient of t^q = c_q × 2^q / x^q

    So:
      v_ℓ^{(p,q)}(x) = c_q × 2^q / (x^q × (x + 4m²)^{p+1})

    where c_q is the q-th coefficient when expanding C_ℓ(1+z) around z=0.

    Since C_ℓ(1+z) involves expanding the Gegenbauer about z=0, we need
    the Taylor coefficients of C_ℓ(1+z).

    We return polynomial coefficients of the *numerator* of v_ℓ^{(p,q)}(x),
    after clearing the denominator x^q (x + 4m²)^{p+1}.

    For SDPB, we need a polynomial in x. Since v_ℓ^{(p,q)}(x) has a
    rational dependence on x, we can multiply by a suitable denominator
    (absorbed into the DampedRational prefactor or handled via the
    polynomial structure).

    In practice, for the simplest setup we work at fixed (p,q) and
    the positivity condition for each spin ℓ becomes:
      Σ_n y_n × v_ℓ^{(n)}(x) ≥ 0 for all x ≥ 0

    where the index n runs over the Wilson coefficients being optimized.

    Parameters
    ----------
    ell : int
        Spin.
    p : int
        Power of s in the Wilson coefficient.
    q : int
        Power of t in the Wilson coefficient.
    m_sq : Fraction
        Mass squared (default 1).
    d : int
        Spacetime dimension (default 4).
    x_shift : Fraction
        Additional shift (default 0), so μ = x + 4m² + x_shift.

    Returns
    -------
    list of Fraction
        Polynomial coefficients [a_0, a_1, ...] of v_ℓ^{(p,q)}(x).
    """
    # Get Gegenbauer coefficients
    geg_coeffs = gegenbauer_coefficients(ell, d)

    # Taylor coefficients of C_ℓ(1 + z) around z = 0
    # If C_ℓ(w) = Σ_k c_k w^k, then C_ℓ(1+z) = Σ_k c_k (1+z)^k
    # The coefficient of z^q in C_ℓ(1+z):
    taylor_q = _taylor_coeff_at_one(geg_coeffs, q)

    # v_ℓ^{(p,q)}(x) = taylor_q × 2^q / (x^q × (x + 4m²)^{p+1})
    #
    # For SDPB, we multiply through by the denominator and incorporate
    # the 1/x^q part into the prefactor (as poles).
    #
    # The polynomial part (after clearing denominators) is:
    #   numerator(x) = taylor_q × 2^q    (a constant)
    # The denominator is x^q × (x + 4m²)^{p+1}
    #
    # We encode this as a degree-0 polynomial times a DampedRational
    # with poles at x=0 (multiplicity q) and x=-4m² (multiplicity p+1).

    # For the simplified approach (Method A: polynomial in x after
    # multiplying by the common denominator), we return the overall
    # numerical coefficient.
    numerator = taylor_q * Fraction(2) ** q

    return numerator


def _taylor_coeff_at_one(poly_coeffs: List[Fraction], q: int) -> Fraction:
    """
    Given polynomial P(w) = Σ_k c_k w^k, compute the coefficient of z^q
    in the Taylor expansion of P(1+z) around z=0.

    P(1+z) = Σ_k c_k (1+z)^k
    Coefficient of z^q = Σ_{k≥q} c_k × C(k, q)

    where C(k,q) is the binomial coefficient.
    """
    result = Fraction(0)
    for k in range(q, len(poly_coeffs)):
        result += poly_coeffs[k] * _binomial(k, q)
    return result


def _binomial(n: int, k: int) -> Fraction:
    """Compute binomial coefficient C(n, k) as a Fraction."""
    if k < 0 or k > n:
        return Fraction(0)
    if k == 0 or k == n:
        return Fraction(1)
    result = Fraction(1)
    for i in range(min(k, n - k)):
        result = result * (n - i) / (i + 1)
    return result


# ---------------------------------------------------------------------------
# Building the positivity polynomial matrix for SDPB
# ---------------------------------------------------------------------------

def build_positivity_polynomials(
    ell: int,
    wilson_indices: List[Tuple[int, int]],
    m_sq: Fraction = Fraction(1),
    d: int = 4,
    poly_degree: int = 10,
) -> List[List[Fraction]]:
    """
    Build the 1×1 polynomial matrix entries for spin ℓ.

    For each Wilson coefficient index (p_n, q_n) in wilson_indices,
    we compute v_ℓ^{(p_n, q_n)}(x) as a polynomial in x.

    The positivity constraint is:
      M_ℓ(x) = Σ_n z_n × v_ℓ^{(n)}(x) ≥ 0  for all x ≥ 0

    After clearing the common denominator D(x) = x^{q_max} × (x+4m²)^{p_max+1},
    we get a polynomial positivity constraint.

    Parameters
    ----------
    ell : int
        Spin.
    wilson_indices : list of (int, int)
        List of (p, q) pairs for the Wilson coefficients to include.
    m_sq : Fraction
        Mass squared.
    d : int
        Spacetime dimension.
    poly_degree : int
        Maximum polynomial degree for the x-polynomial representation.

    Returns
    -------
    list of list of Fraction
        For each Wilson coefficient, the polynomial coefficients of
        v_ℓ^{(n)}(x) in the common-denominator form.
    """
    mu_threshold = 4 * m_sq
    q_max = max(q for _, q in wilson_indices) if wilson_indices else 0
    p_max = max(p for p, _ in wilson_indices) if wilson_indices else 0

    polynomials = []

    for p_n, q_n in wilson_indices:
        # The spectral function is:
        #   v_ℓ^{(p,q)}(x) = [taylor_q of C_ℓ at 1] × 2^q / (x^q (x+4m²)^{p+1})
        #
        # After multiplying by the common denominator D(x) = x^{q_max} (x+4m²)^{p_max+1}:
        #   v_ℓ^{(p,q)}(x) × D(x) = taylor_q × 2^q × x^{q_max - q} × (x+4m²)^{p_max - p}
        #
        # This is a polynomial in x of degree (q_max - q) + (p_max - p).

        geg_coeffs = gegenbauer_coefficients(ell, d)
        taylor_q = _taylor_coeff_at_one(geg_coeffs, q_n)
        coeff = taylor_q * Fraction(2) ** q_n

        # Build x^{q_max - q_n} as polynomial
        x_power = q_max - q_n
        poly_x = [Fraction(0)] * (x_power + 1)
        poly_x[x_power] = Fraction(1)

        # Build (x + 4m²)^{p_max - p_n} as polynomial
        shift_power = p_max - p_n
        poly_shift = _expand_binomial_power(mu_threshold, shift_power)

        # Multiply the two polynomials
        product = _poly_multiply(poly_x, poly_shift)

        # Scale by the coefficient
        result = [c * coeff for c in product]

        # Pad or truncate to poly_degree + 1 coefficients
        while len(result) < poly_degree + 1:
            result.append(Fraction(0))
        result = result[: poly_degree + 1]

        polynomials.append(result)

    return polynomials


def _expand_binomial_power(a: Fraction, n: int) -> List[Fraction]:
    """
    Expand (x + a)^n as a polynomial in x.
    Returns coefficients [c_0, c_1, ..., c_n] where (x+a)^n = Σ c_k x^k.
    """
    if n == 0:
        return [Fraction(1)]

    coeffs = [Fraction(0)] * (n + 1)
    for k in range(n + 1):
        coeffs[k] = _binomial(n, k) * (a ** (n - k))
    return coeffs


def _poly_multiply(p: List[Fraction], q: List[Fraction]) -> List[Fraction]:
    """Multiply two polynomials represented as coefficient lists."""
    if not p or not q:
        return [Fraction(0)]
    result = [Fraction(0)] * (len(p) + len(q) - 1)
    for i, pi in enumerate(p):
        for j, qj in enumerate(q):
            result[i + j] += pi * qj
    return result


# ---------------------------------------------------------------------------
# Utility: high-precision string output
# ---------------------------------------------------------------------------

def fraction_to_str(f: Fraction, precision: int = 200) -> str:
    """
    Convert a Fraction to a decimal string with the given number of
    significant digits.
    """
    if f == 0:
        return "0"

    from decimal import Decimal, getcontext
    getcontext().prec = precision + 10

    d = Decimal(f.numerator) / Decimal(f.denominator)
    # Format with enough digits
    s = format(d, f'.{precision}f')

    # Strip trailing zeros but keep at least one decimal place for non-integers
    if '.' in s:
        s = s.rstrip('0').rstrip('.')

    # If the result is just "-0", return "0"
    if s == '-0':
        return '0'

    return s
