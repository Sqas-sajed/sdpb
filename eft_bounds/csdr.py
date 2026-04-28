"""
csdr.py — CSDR kernel formulas for general spacetime dimension d.

Implements the kernel functions needed to express all W_{n-m,m} as
"heavy averages" following user notes eq.(2) and:
  - Sinha & Zahed, "Crossing Symmetric Dispersion Relations in QFTs" eq.(11)
  - Caron-Huot & Duong, "Extremal Effective Field Theories" Section 3.3

Convention
----------
  alpha = (d-3)/2  (Gegenbauer parameter, chosen so 2*alpha = d-3 and
                   (2*ell + 2*alpha) = (2*ell + d - 3) matches eq.(1))

  S1, S2, S3 = Mandelstam variables shifted so S1+S2+S3 = 0.
    Concretely: S_i = (Mandelstam variable)_i - mu/3,  sum S_i = 0.
  x = -(S1*S2 + S2*S3 + S3*S1),  y = -S1*S2*S3.
  Amplitude: M = sum_{p,q>=0} W_{p,q} x^p y^q

  W_{n-m,m}: first index = n-m, second index = m.
  - n > m  (first index > 0): EFT objectives.
  - n = m  (first index = 0): EFT objective W_{0,m}.
  - n < m  (first index < 0): CSDR null constraints (W = 0). Kernel from eq.(11).

User notes eq.(2) [derived for null constraints m > n >= 1]:

  W_{n-m,m} = < D^{(n,m)}_{ell,alpha} / n^{(d)}_ell
                 * C^{(alpha)}_ell(1) * (2*ell+d-3)
                 / s1^{2n+m} >

  where the heavy average is:
    <F> = sum_{ell even} n^{(d)}_ell * int_{delta0}^inf ds1/s1
          * s1^{4-d}/pi * rho_ell(s1) * F(s1,ell)

  and   n^{(d)}_J = (4*pi)^{d/2} * (d+2J-3) * Gamma(d+J-3)
                    / (pi * Gamma((d-2)/2) * Gamma(J+1))

WHY n^{(d)}_ell CANCELS IN SDPB BLOCKS
----------------------------------------
Expanding the heavy average explicitly:

  <F> = sum_{ell even} n^{(d)}_ell * int ρ_ell * F(s1,ell) * (measure)
      = sum_{ell even} n^{(d)}_ell * int ρ_ell
          * [D^{(n,m)} / n^{(d)}_ell * C_ell(1) * (2ell+d-3) / s1^{2n+m}]
          * (measure)
      = sum_{ell even} int ρ_ell
          * D^{(n,m)} * C_ell(1) * (2ell+d-3) / s1^{2n+m}
          * (measure)

So n^{(d)}_ell cancels between the prefactor and the denominator of F.

For SDPB, we need: P^{(ell)}(x) = Σ_k y_k * (spectral contribution of z_k at spin ell) >= 0.
Since n^{(d)}_ell > 0, the condition ⟨F⟩ ≥ 0 for all ρ_ell ≥ 0 is equivalent to:
  D^{(n,m)} * C_ell(1) * (2ell+d-3) * (polynomial in s1) ≥ 0
for all even ell and all s1 ≥ delta0.
This is precisely what SDPB enforces.  The factor n^{(d)}_ell is absorbed into
the measure weights and is positive, so it does NOT change the feasibility problem.
Therefore, SDPB polynomial blocks do NOT need to include n^{(d)}_ell explicitly.

OBJECTIVE D FORMULA (TODO — formula not yet provided by user)
--------------------------------------------------------------
For objectives (n >= m, first CSDR index n-m >= 0), the CSDR paper (Sinha-Zahed)
uses a different kernel formula from eq.(11), which applies only to null constraints
(m > n).  The user has indicated the correct objective formula in their notes but
the explicit formula was not yet transmitted.

Current placeholder for objectives:
  kappa^{obj}_{p,q,ell} = C^{(alpha)}_ell(1) * (2*ell+d-3)

This is the leading forward-scattering contribution (positive-definite).
Once the user provides the correct formula, replace obj_kernel_coeff below.

Note on CSDR vs Extremal EFT:
The forms of null constraints and objectives in CSDR may or may not differ from
the fixed-t dispersion relation result.  Do not assume they differ; verify numerically.
"""

from fractions import Fraction
from math import factorial
from typing import List, Tuple


# ---------------------------------------------------------------------------
# Basic special functions
# ---------------------------------------------------------------------------

def pochhammer(a: Fraction, n: int) -> Fraction:
    """
    Pochhammer symbol (a)_n = a*(a+1)*...*(a+n-1).

    (a)_0 = 1 by convention.
    """
    if n < 0:
        raise ValueError(f"pochhammer requires n >= 0, got n={n}")
    result = Fraction(1)
    for i in range(n):
        result *= a + i
    return result


def alpha_from_d(d: int) -> Fraction:
    """
    Return the Gegenbauer parameter alpha = (d-3)/2.

    This ensures (2*ell + 2*alpha) = (2*ell + d - 3), which matches
    the spin measure in the CSDR paper and user notes eq.(1).
    """
    if d < 3:
        raise ValueError(f"Spacetime dimension d must be >= 3, got d={d}")
    return Fraction(d - 3, 2)


def gegenbauer_at_one(ell: int, alpha: Fraction) -> Fraction:
    """
    C^{(alpha)}_ell(1) = (2*alpha)_ell / ell!

    For d=4, alpha=1/2: C^{(1/2)}_ell(1) = 1 for all ell.
    For d=6, alpha=3/2: C^{(3/2)}_ell(1) = (ell+1)*(ell+2)/2.
    """
    return pochhammer(2 * alpha, ell) / Fraction(factorial(ell))


# ---------------------------------------------------------------------------
# D coefficient — CSDR eq.(11)
# ---------------------------------------------------------------------------

def D_coeff(n: int, m: int, ell: int, alpha: Fraction) -> Fraction:
    """
    Compute D^{(n,m)}_{ell,alpha} from CSDR eq.(11).

    Valid for m > n >= 1.  Returns exact rational result.

    Formula:
      D^{(n,m)}_{ell,alpha} = Gamma(m-n) * sum_{j=n}^{m}
          (-4)^j * (-ell/2)_j * (alpha + ell/2)_j * (3j - m - 2n)
          / (j! * (m-j)! * (j-n)! * (alpha + 1/2)_j)

    The factor Gamma(m-n) = (m-n-1)! appears outside the sum.

    Parameters
    ----------
    n : int   First CSDR index, n >= 1.
    m : int   Second CSDR index, m > n.
    ell : int  Spin (angular momentum, non-negative integer).
    alpha : Fraction  Gegenbauer parameter = (d-3)/2.

    Returns
    -------
    Fraction  Exact rational value of D^{(n,m)}_{ell,alpha}.

    Notes
    -----
    For odd ell, the Pochhammer (-ell/2)_j contains a zero factor when
    j > ell/2, so the sum naturally terminates at j = ell/2 for odd ell.
    Even so, D is well-defined; for identical scalars only even ell matter.
    """
    if n < 1:
        raise ValueError(
            f"D_coeff requires n >= 1 (null constraint index), got n={n}. "
            "For objectives (n >= m >= 0), use obj_kernel_coeff instead."
        )
    if m <= n:
        raise ValueError(
            f"D_coeff requires m > n, got n={n}, m={m}. "
            "This is a null constraint (first CSDR index < 0); for objectives use obj_kernel_coeff."
        )

    # Gamma(m - n) = (m - n - 1)!  since m > n >= 1 implies m - n >= 1
    gamma_factor = Fraction(factorial(m - n - 1))

    total = Fraction(0)
    for j in range(n, m + 1):
        poch_neg_ell_half = pochhammer(Fraction(-ell, 2), j)
        poch_alpha_ell_half = pochhammer(alpha + Fraction(ell, 2), j)
        three_j_term = Fraction(3 * j - m - 2 * n)
        poch_alpha_half = pochhammer(alpha + Fraction(1, 2), j)

        if poch_alpha_half == 0:
            # Denominator zero: term is zero (or skip)
            continue

        numerator = (
            Fraction((-4) ** j)
            * poch_neg_ell_half
            * poch_alpha_ell_half
            * three_j_term
            * gamma_factor
        )
        denominator = (
            Fraction(factorial(j))
            * Fraction(factorial(m - j))
            * Fraction(factorial(j - n))
            * poch_alpha_half
        )
        total += numerator / denominator

    return total


def D_coeff_batch(n: int, m: int, ell_values: List[int], alpha: Fraction) -> List[Fraction]:
    """
    Compute D^{(n,m)}_{ell,alpha} for a list of ell values.

    Parameters
    ----------
    n, m : int  CSDR indices (m > n >= 1).
    ell_values : list of int  Spin values to evaluate.
    alpha : Fraction  Gegenbauer parameter.

    Returns
    -------
    list of Fraction
    """
    return [D_coeff(n, m, ell, alpha) for ell in ell_values]


# ---------------------------------------------------------------------------
# SDPB kernel coefficients
# ---------------------------------------------------------------------------

def null_kernel_coeff(n: int, m: int, ell: int, d: int) -> Fraction:
    """
    Kernel coefficient for null constraint W_{n-m,m} at spin ell.

    From user notes eq.(2), the heavy-average kernel (after dividing out
    n^{(d)}_ell which cancels between measure and integrand) is:

      kappa_{n,m,ell} = D^{(n,m)}_{ell,alpha} * C^{(alpha)}_ell(1) * (2*ell + d - 3)

    where alpha = (d-3)/2.

    Parameters
    ----------
    n : int   CSDR first argument (n >= 1).
    m : int   CSDR second argument (m > n).
    ell : int  Spin.
    d : int   Spacetime dimension.

    Returns
    -------
    Fraction  Exact rational kernel coefficient.
    """
    alpha = alpha_from_d(d)
    D = D_coeff(n, m, ell, alpha)
    C1 = gegenbauer_at_one(ell, alpha)
    spin_measure = Fraction(2 * ell + d - 3)
    return D * C1 * spin_measure


def obj_kernel_coeff(p: int, q: int, ell: int, d: int) -> Fraction:
    """
    Forward-scattering kernel coefficient for objective W_{p,q} at spin ell.

    For objectives (p = n-m >= 0, q = m >= 0), the CSDR large-s1 kernel
    from eq.(11) gives D = 0 (empty sum). We instead use the forward-limit
    spectral function, which is the leading contribution from each partial
    wave to W_{p,q} in the CSDR dispersion relation:

      kappa_{p,q,ell}^obj = C^{(alpha)}_ell(1) * (2*ell + d - 3)

    This kernel is strictly positive for all even ell >= 0 and d >= 4,
    consistent with W_{p,q} >= 0 from unitarity in the CSDR scheme.

    Parameters
    ----------
    p : int   First Wilson-coefficient index (p = n-m >= 0).
    q : int   Second Wilson-coefficient index (q = m >= 0).
    ell : int  Spin.
    d : int   Spacetime dimension.

    Returns
    -------
    Fraction  Exact rational kernel coefficient.
    """
    alpha = alpha_from_d(d)
    C1 = gegenbauer_at_one(ell, alpha)
    spin_measure = Fraction(2 * ell + d - 3)
    return C1 * spin_measure


# ---------------------------------------------------------------------------
# Spectral power
# ---------------------------------------------------------------------------

def s1_power(p: int, q: int) -> int:
    """
    Denominator power of s1 for W_{p,q}: returns 2p+3q.

    Equivalently, for CSDR parameters n = p+q, m = q: returns 2n+m = 2p+3q.

    Parameters
    ----------
    p : int  First Wilson-coefficient index (can be negative for null constraints).
    q : int  Second Wilson-coefficient index (>= 0).

    Returns
    -------
    int  Power of s1 in the kernel denominator.
    """
    return 2 * p + 3 * q


# ---------------------------------------------------------------------------
# Pair enumeration
# ---------------------------------------------------------------------------

def enumerate_null_pairs(K: int) -> List[Tuple[int, int]]:
    """
    Enumerate CSDR null constraint pairs (n, m) with m > n >= 1 and 2n+m <= K.

    These correspond to W_{n-m,m} with first index n-m < 0.

    Parameters
    ----------
    K : int  Maximum allowed 2n+m.

    Returns
    -------
    list of (n, m)  Sorted by 2n+m, then m.
    """
    pairs = []
    for n in range(1, K + 1):
        for m in range(n + 1, K + 1):
            if 2 * n + m <= K:
                pairs.append((n, m))
    pairs.sort(key=lambda x: (2 * x[0] + x[1], x[1]))
    return pairs


def enumerate_obj_pairs(K: int) -> List[Tuple[int, int]]:
    """
    Enumerate objective Wilson-coefficient pairs (p, q) with p >= 1, q >= 0
    and 2p+3q <= K.  (These have positive first CSDR index p = n-m >= 1.)

    The case p=0 (first CSDR index = 0) is excluded because D^{(n,n)} from
    eq.(11) is degenerate (n=m gives an empty sum), requiring a separate
    treatment.  For the SDPB problem, p=0 objectives can be included via
    W_{0,q} = < C^{(alpha)}_ell(1)*(2ell+d-3) / s1^{3q} > if needed.

    Parameters
    ----------
    K : int  Maximum allowed 2p+3q.

    Returns
    -------
    list of (p, q)  Sorted by 2p+3q, then q.
    """
    pairs = []
    for p in range(1, K + 1):
        for q in range(0, K + 1):
            if 2 * p + 3 * q <= K and 2 * p + 3 * q > 0:
                pairs.append((p, q))
    pairs.sort(key=lambda x: (2 * x[0] + 3 * x[1], x[1]))
    return pairs


def enumerate_obj_pairs_including_zero(K: int) -> List[Tuple[int, int]]:
    """
    Enumerate objective Wilson-coefficient pairs (p, q) with p >= 0, q >= 0
    and 2p+3q <= K and 2p+3q > 0.

    Includes p=0 cases (W_{0,q} for q >= 1).

    Parameters
    ----------
    K : int  Maximum allowed 2p+3q.

    Returns
    -------
    list of (p, q)  Sorted by 2p+3q, then q.
    """
    pairs = []
    for p in range(0, K + 1):
        for q in range(0, K + 1):
            val = 2 * p + 3 * q
            if 0 < val <= K:
                pairs.append((p, q))
    pairs.sort(key=lambda x: (2 * x[0] + 3 * x[1], x[1]))
    return pairs


# ---------------------------------------------------------------------------
# Polynomial utilities
# ---------------------------------------------------------------------------

def poly_expand_1px(k: int) -> List[Fraction]:
    """
    Return coefficients of (1+x)^k as a polynomial in x: [C(k,0), C(k,1), ..., C(k,k)].

    Parameters
    ----------
    k : int  Non-negative integer exponent.

    Returns
    -------
    list of Fraction  Polynomial coefficients [a_0, a_1, ..., a_k].
    """
    if k < 0:
        raise ValueError(f"poly_expand_1px requires k >= 0, got k={k}")
    from math import comb
    return [Fraction(comb(k, j)) for j in range(k + 1)]


def poly_scale(poly: List[Fraction], scale: Fraction) -> List[Fraction]:
    """Multiply a polynomial by a scalar."""
    return [c * scale for c in poly]


def poly_pad(poly: List[Fraction], length: int) -> List[Fraction]:
    """Pad polynomial with trailing zeros to the given length."""
    result = list(poly)
    while len(result) < length:
        result.append(Fraction(0))
    return result


# ---------------------------------------------------------------------------
# Self-consistency check
# ---------------------------------------------------------------------------

def check_closed_form_1_2(precision: int = 20) -> bool:
    """
    Verify D^{(1,2)}_{ell,1/2} against the closed form quoted in CSDR paper:

        D^{(1,2)}_{ell,alpha} = 2*ell*(ell+2*alpha)*(-11-10*alpha+2*ell*(ell+2*alpha))
                                 / ((2*alpha+1)*(2*alpha+3))

    Returns True if all test cases pass.
    """
    alpha = Fraction(1, 2)
    n, m = 1, 2
    passed = True
    for ell in [0, 2, 4, 6, 8]:
        direct = D_coeff(n, m, ell, alpha)
        closed = (
            2 * ell * (ell + 2 * alpha)
            * (-11 - 10 * alpha + 2 * ell * (ell + 2 * alpha))
            / ((2 * alpha + 1) * (2 * alpha + 3))
        )
        if direct != closed:
            passed = False
    return passed
