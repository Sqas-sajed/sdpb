"""
csdr.py — CSDR kernel formulas for general spacetime dimension d.

Implements the kernel functions needed to express all W_{n-m,m} as
"heavy averages" following user notes eq.(2) and:
  - Sinha & Zahed, "Crossing Symmetric Dispersion Relations in QFTs" eq.(11)
  - Caron-Huot & Duong, "Extremal Effective Field Theories" Section 3.3

UNIFIED (n,m) NOTATION
-----------------------
All Wilson coefficients W_{n-m,m} are labelled by a single pair (n, m):
  - n > m >= 0 : EFT objectives (first CSDR index n-m > 0).
  - n = m >= 1 : EFT objective W_{0,m} (first CSDR index = 0).
  - n < m, m > n >= 1 : CSDR null constraints W_{n-m,m} = 0.

Both objectives and null constraints use the same expression (user notes eq.(2)):

  W_{n-m,m} = < D^{(n,m)}_{ell,alpha} / n^{(d)}_ell
                 * C^{(alpha)}_ell(1) * (2*ell+d-3)
                 / s1^{2n+m} >

The spectral power is 2n+m for ALL cases.

GEGENBAUER PARAMETER
---------------------
  alpha = (d-3)/2

CSDR MANDELSTAM BASIS
---------------------
  S1, S2, S3 = (s - mu/3, t - mu/3, u - mu/3)  with S1+S2+S3 = 0, mu = s+t+u.
  x = -(S1*S2 + S2*S3 + S3*S1),  y = -S1*S2*S3.
  Amplitude: M = sum_{n',m'>=0} W_{n',m'} x^{n'} y^{m'}

HEAVY AVERAGE DEFINITION (Extremal EFT, Section 3.3)
-----------------------------------------------------
  <F> = sum_{ell even} n^{(d)}_ell * int_{delta0}^inf ds1/s1
        * s1^{4-d}/pi * rho_ell(s1) * F(s1,ell)

  n^{(d)}_J = (4*pi)^{d/2} * (d+2J-3) * Gamma(d+J-3)
              / (pi * Gamma((d-2)/2) * Gamma(J+1))

WHY n^{(d)}_ell CANCELS — EXPLICIT PROOF
-----------------------------------------
The integrand F = D/n_ell^(d) * C_ell(1) * (2ell+d-3) / s1^{2n+m}, so:

  <F> = sum_ell n_ell^(d) * integral rho_ell * [D / n_ell^(d) * C_ell(1) * ...] * (measure)
      = sum_ell integral rho_ell * D * C_ell(1) * (2ell+d-3) * (measure)

n_ell^(d) > 0 cancels exactly. The SDPB polynomial P^(ell)(x) encodes the
positivity condition on D * C_ell(1) * (2ell+d-3) as a polynomial in x = s1/delta0 - 1 >= 0.
Multiplying each spin block by n_ell^(d) > 0 is a positive rescaling and does not
change the feasibility region. SDPB blocks therefore do NOT include n_ell^(d).

D COEFFICIENT — TWO CASES
--------------------------
Case 1: NULL CONSTRAINTS (m > n >= 1)
  D^{(n,m)}_{ell,alpha} from CSDR eq.(11) [Sinha-Zahed]:
    D = Gamma(m-n) * sum_{j=n}^{m}
          (-4)^j * (-ell/2)_j * (alpha + ell/2)_j * (3j-m-2n)
          / (j! * (m-j)! * (j-n)! * (alpha + 1/2)_j)
  Implemented in D_coeff(n, m, ell, alpha).
  D depends on the spin ell.

Case 2: OBJECTIVES (n >= m >= 0)
  Two sub-cases:

  2a. n > m (first CSDR index n-m > 0, "standard objectives" like W_{1,0}, W_{2,0}):
    D = sum_{j=0}^{m}
          (-4)^j * (-ell/2)_j * (alpha+ell/2)_j * (3j-m-2n) * (n-j)!
          / ((alpha+1/2)_j * j! * (m-j)! * (n-m)!) * (-1)^{m+j+1}
    For m=0 the sum reduces to the single j=0 term: D = 2n (independent of ell).
    Implemented in D_coeff_obj(n, m, ell, alpha).

  2b. n = m (first CSDR index = 0, "on-diagonal objectives" like W_{0,1}):
    The n > m formula has a 0/0 indeterminate form at j = m = n (the last term
    has both (3n-n-2n)=0 and the implicit Gamma(m-n) diverge).  The correct value
    is the regularised limit of the null formula as m -> n from above:

      D^{(n,n)}_{ell,alpha} = (-1)^{n+1} * 4^n * (-ell/2)_n * (alpha+ell/2)_n
                               / (n! * (alpha+1/2)_n)

    Key property: (-ell/2)_n = 0 whenever ell = 0, 2, ..., 2(n-1), so the kernel
    vanishes for the lowest 2n/2 = n even spins.  In particular for n=1 (W_{0,1}):
      D^{(1,1)}_{ell,alpha} = -2*ell*(alpha+ell/2)/(alpha+1/2)
    which equals 0 at ell=0 and is negative for all ell >= 2.
    Implemented in D_coeff_obj(n, m, ell, alpha) with a special n==m branch.

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
    Compute D^{(n,m)}_{ell,alpha} from CSDR eq.(11).  NULL CONSTRAINTS only.

    Valid for m > n >= 1.  Returns exact rational result.

    Formula (Sinha-Zahed eq.(11)):
      D^{(n,m)}_{ell,alpha} = Gamma(m-n) * sum_{j=n}^{m}
          (-4)^j * (-ell/2)_j * (alpha + ell/2)_j * (3j - m - 2n)
          / (j! * (m-j)! * (j-n)! * (alpha + 1/2)_j)

    This coefficient DEPENDS on the spin ell.

    Parameters
    ----------
    n : int   First CSDR index, n >= 1.
    m : int   Second CSDR index, m > n.
    ell : int  Spin (angular momentum, non-negative integer).
    alpha : Fraction  Gegenbauer parameter = (d-3)/2.

    Returns
    -------
    Fraction  Exact rational value of D^{(n,m)}_{ell,alpha}.
    """
    if n < 1:
        raise ValueError(
            f"D_coeff (null) requires n >= 1, got n={n}. "
            "For objectives (n >= m >= 0), use D_coeff_obj."
        )
    if m <= n:
        raise ValueError(
            f"D_coeff (null) requires m > n, got n={n}, m={m}. "
            "For objectives (n >= m), use D_coeff_obj."
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


def D_coeff_obj(n: int, m: int, ell: int, alpha: Fraction) -> Fraction:
    """
    Compute D^{(n,m)}_{ell,alpha} for OBJECTIVES (n >= m >= 0).

    Two sub-cases depending on whether n == m or n > m.

    Case n > m (includes m == 0)
    ----------------------------
    Formula (derived from CSDR eq.(4), valid for n > m >= 0):

      D^{(n,m)}_{ell,alpha} = sum_{j=0}^{m}
          (-4)^j * (-ell/2)_j * (alpha+ell/2)_j * (3j-m-2n) * (n-j)!
          / ((alpha+1/2)_j * j! * (m-j)! * (n-m)!) * (-1)^{m+j+1}

    For m == 0 the sum collapses to D = 2n (independent of ell).

    Case n == m (on-diagonal objectives, first CSDR index = 0)
    -----------------------------------------------------------
    The n > m formula has a 0/0 indeterminate form at j = n = m (the last
    summand has (3n-n-2n) = 0 while the Gamma(m-n) prefactor of the parent
    null formula diverges).  The correct value is obtained as the regularised
    limit of the CSDR null formula as m -> n from above:

      D^{(n,n)}_{ell,alpha} = (-1)^{n+1} * 4^n
                               * (-ell/2)_n * (alpha+ell/2)_n
                               / (n! * (alpha+1/2)_n)

    This vanishes whenever (-ell/2)_n = 0, i.e. for ell = 0, 2, ..., 2(n-1).
    For n = 1 (W_{0,1}):
      D^{(1,1)}_{ell,alpha} = -2*ell*(alpha+ell/2)/(alpha+1/2)
    which is 0 at ell=0 and strictly negative for ell >= 2.

    Parameters
    ----------
    n : int    First CSDR index (n >= m).
    m : int    Second CSDR index (m >= 0).
    ell : int  Spin (ell-dependent coefficient).
    alpha : Fraction  Gegenbauer parameter = (d-3)/2.

    Returns
    -------
    Fraction  Exact rational value of D^{(n,m)}_{ell,alpha}.
    """
    if m < 0:
        raise ValueError(f"D_coeff_obj requires m >= 0, got m={m}.")
    if n < m:
        raise ValueError(
            f"D_coeff_obj requires n >= m, got n={n}, m={m}. "
            "For null constraints (m > n >= 1), use D_coeff."
        )

    # --- n == m: use the regularised limit formula ---
    if n == m and n >= 1:
        # D^{(n,n)} = (-1)^{n+1} * 4^n * (-ell/2)_n * (alpha+ell/2)_n / [n! * (alpha+1/2)_n]
        poch_neg = pochhammer(Fraction(-ell, 2), n)
        poch_pos = pochhammer(alpha + Fraction(ell, 2), n)
        poch_denom = pochhammer(alpha + Fraction(1, 2), n)
        if poch_denom == 0:
            return Fraction(0)
        sign = Fraction((-1) ** (n + 1))
        return (
            sign
            * Fraction(4 ** n)
            * poch_neg
            * poch_pos
            / (Fraction(factorial(n)) * poch_denom)
        )

    # --- n > m (includes m == 0): standard sum formula ---
    # For m = 0: sum has only j=0 term:
    #   (-1)^{0+0+1} * (-4)^0 * 1 * 1 * (0-0-2n) * n! / (1 * 1 * 1 * n!) = (-1)*(-2n) = 2n
    total = Fraction(0)
    for j in range(m + 1):  # j = 0 to m
        poch_neg_ell_half = pochhammer(Fraction(-ell, 2), j)
        poch_alpha_ell_half = pochhammer(alpha + Fraction(ell, 2), j)
        poch_alpha_half = pochhammer(alpha + Fraction(1, 2), j)

        if poch_alpha_half == 0:
            continue

        three_j_term = Fraction(3 * j - m - 2 * n)
        n_minus_j_fac = Fraction(factorial(n - j))   # n >= m >= j, so n-j >= 0
        sign = Fraction((-1) ** (m + j + 1))

        numerator = (
            Fraction((-4) ** j)
            * poch_neg_ell_half
            * poch_alpha_ell_half
            * three_j_term
            * n_minus_j_fac
            * sign
        )
        denominator = (
            poch_alpha_half
            * Fraction(factorial(j))
            * Fraction(factorial(m - j))
            * Fraction(factorial(n - m))
        )
        total += numerator / denominator

    return total


def D_coeff_batch(n: int, m: int, ell_values: List[int], alpha: Fraction) -> List[Fraction]:
    """
    Compute D^{(n,m)}_{ell,alpha} (null, eq.11) for a list of ell values.

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

    From user notes eq.(2), the SDPB block coefficient (after n_ell^(d) cancels) is:

      kappa_{n,m,ell} = D^{(n,m)}_{ell,alpha} * C^{(alpha)}_ell(1) * (2*ell + d - 3)

    where D^{(n,m)}_{ell,alpha} is from CSDR eq.(11) (depends on ell),
    and alpha = (d-3)/2.

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


def obj_kernel_coeff(n: int, m: int, ell: int, d: int) -> Fraction:
    """
    Kernel coefficient for objective W_{n-m,m} at spin ell.

    Uses the corrected objective D formula (D_coeff_obj), which is ell-dependent
    (same Pochhammer structure as the null formula: (-ell/2)_j (alpha+ell/2)_j / (alpha+1/2)_j).

    From user notes eq.(2):

      kappa_{n,m,ell}^obj = D^{(n,m)}_{ell,alpha} * C^{(alpha)}_ell(1) * (2*ell + d - 3)

    Parameters
    ----------
    n : int   CSDR first index (n >= m).
    m : int   CSDR second index (m >= 0).
    ell : int  Spin.
    d : int   Spacetime dimension.

    Returns
    -------
    Fraction  Exact rational kernel coefficient.
    """
    alpha = alpha_from_d(d)
    D = D_coeff_obj(n, m, ell, alpha)
    C1 = gegenbauer_at_one(ell, alpha)
    spin_measure = Fraction(2 * ell + d - 3)
    return D * C1 * spin_measure


# ---------------------------------------------------------------------------
# Spectral power
# ---------------------------------------------------------------------------

def spectral_power_nm(n: int, m: int) -> int:
    """
    Denominator power of s1 for W_{n-m,m}: returns 2n+m.

    Valid for both objectives (n >= m) and null constraints (m > n >= 1).

    Parameters
    ----------
    n : int  First CSDR index.
    m : int  Second CSDR index.

    Returns
    -------
    int  Power of s1 in the kernel denominator (2n+m).
    """
    return 2 * n + m


def s1_power(p: int, q: int) -> int:
    """
    Denominator power of s1 for W_{p,q}: returns 2p+3q.

    This is equivalent to spectral_power_nm(p+q, q) = 2(p+q)+q = 2p+3q.

    Kept for backward compatibility.  For new code use spectral_power_nm(n,m).

    Parameters
    ----------
    p : int  First Wilson-coefficient index (p = n-m).
    q : int  Second Wilson-coefficient index (q = m).

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


def enumerate_obj_pairs_nm(K: int) -> List[Tuple[int, int]]:
    """
    Enumerate CSDR objective pairs (n, m) with n >= m >= 0 and 2n+m <= K, 2n+m > 0.

    Includes both:
    - W_{n-m,m} with first index n-m > 0 (n > m): standard objectives.
    - W_{0,m} with first index 0, n = m >= 1: physically valid operators.

    The spectral power is 2n+m.

    Parameters
    ----------
    K : int  Maximum allowed 2n+m.

    Returns
    -------
    list of (n, m)  Sorted by 2n+m, then m.
    """
    pairs = []
    for m in range(0, K + 1):
        for n in range(m, K + 1):  # n >= m (includes n == m for W_{0,m})
            val = 2 * n + m
            if 0 < val <= K:
                pairs.append((n, m))
    pairs.sort(key=lambda x: (2 * x[0] + x[1], x[1]))
    return pairs
    return pairs


def enumerate_obj_pairs(K: int) -> List[Tuple[int, int]]:
    """
    Enumerate objective pairs in (p, q) notation: p = n-m >= 1, q = m >= 0, 2p+3q <= K.

    Kept for backward compatibility.  For new code use enumerate_obj_pairs_nm(K).

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
# Self-consistency checks
# ---------------------------------------------------------------------------

def check_closed_form_1_2(precision: int = 20) -> bool:
    """
    Verify D^{(1,2)}_{ell,1/2} (null, eq.11) against the closed form quoted in CSDR paper:

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


def check_obj_D_basic() -> bool:
    """
    Basic sanity checks for D_coeff_obj (objectives, corrected ell-dependent formula).

    Verifies:
      - m=0 cases: D(n, 0, ell=0) = 2n  (ell-independent forward dispersion).
      - n=m cases: D(n, n, ell=0) = 0   (the limit formula vanishes at ell=0
          because (-ell/2)_n = 0 when ell=0 and n>=1).
      - n=m, ell>=2: D is ell-dependent and negative.
        D(1,1,2,alpha=1/2) = -6  [= -2*(alpha=1/2)*(1+2)/(alpha+1/2=1) ]

    Returns True if all checks pass.
    """
    alpha = Fraction(1, 2)
    passed = True

    # m=0 forward-dispersion cases (ell-independent)
    # D(1,0,any_ell) = 2
    val = D_coeff_obj(1, 0, 0, alpha)
    if val != Fraction(2):
        passed = False

    # D(2,0,any_ell) = 4
    val = D_coeff_obj(2, 0, 0, alpha)
    if val != Fraction(4):
        passed = False

    # n=m limit formula: D(n,n,ell=0) = 0 because (-ell/2)_n|_{ell=0} = 0 for n>=1.
    # Old buggy sum formula returned -3 here; the corrected limit gives 0.
    val = D_coeff_obj(1, 1, 0, alpha)
    if val != Fraction(0):
        passed = False

    # n=m, ell=2: D(1,1,2,1/2) = -2*2*(1/2+1)/1 = -6
    val = D_coeff_obj(1, 1, 2, alpha)
    if val != Fraction(-6):
        passed = False

    # Verify D is a Fraction (no errors) for a few more cases.
    for ell in [4, 6]:
        v = D_coeff_obj(1, 1, ell, alpha)
        if not isinstance(v, Fraction):
            passed = False

    return passed
