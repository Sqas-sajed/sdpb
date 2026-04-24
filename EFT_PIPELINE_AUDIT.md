# EFT Pipeline Audit Report

**Date:** $(date -u +"%Y-%m-%d")  
**Repository:** Sqas-sajed/sdpb  
**Path:** `/home/runner/work/sdpb/sdpb/eft_bounds`

---

## Executive Summary

The EFT pipeline has **3 out of 4 tasks** implemented with varying completeness:

| Task | Status | Completeness | Critical Gaps |
|------|--------|--------------|---------------|
| 1. Analytic null constraints derivation | ✅ Implemented | 90% | Missing explicit Sinha-Zahed eq.(3) implementation |
| 2. Translation to SDPB PMP JSON | ✅ Implemented | 95% | Minor: e^{-1} approximation could be more precise |
| 3. Verification of ≥5 examples | ⚠️ Partial | 40% | **Only 2 examples verified; needs 3 more** |
| 4. Two-sided bounds translation | ❌ Missing | 0% | **Not implemented at all** |

---

## Task-by-Task Analysis

### Task 1: Analytic Null Constraints Derivation ✅

**Current Implementation:**
- **File:** `null_constraints.py`
- **Function:** `get_null_constraints(max_order, m_sq)` (lines 44-167)
- **Algorithm:** Expands M(s,t) = M(u,t) where u = 4m² - s - t using multinomial theorem

**What Works:**
- Correctly derives s↔u crossing constraints using binomial expansion
- Full S₃ crossing via `get_crossing_symmetric_null_constraints()` (lines 334-411)
- Variable elimination to reduce DOF (lines 223-315)
- Produces 11 independent constraints at max_order=4

**Critical Gap:**
- **Missing:** Explicit implementation of Sinha-Zahed equation (3)
  - Current code references equations 3.10-3.20 in comments but doesn't implement eq.(3) directly
  - The paper's eq.(3) likely refers to specific two-sided dispersion bounds
  - Code implements s↔u symmetry constraints but not the dispersive two-sided inequality form

**Likely Bug:**
- Line 370-388 in `get_crossing_symmetric_null_constraints()`:
  - Comment claims to implement "explicit null constraints from the paper"
  - BUT: Actually just calls `get_null_constraints()` + adds s↔t symmetry
  - Does NOT implement the crossing-symmetric dispersion relation form from Sinha-Zahed

**Verification Status:**
- ✅ Polynomial amplitude (s² + u²) passes all constraints (Check 5b)
- ⚠️ Scalar exchange shows expected truncation errors (Check 5)

---

### Task 2: Translation to SDPB PMP JSON ✅

**Current Implementation:**
- **File:** `pmp_generator.py`
- **Functions:** `generate_pmp_json()`, `generate_pmp_for_ratio_bound()`
- **Output:** 12 PMP files in `eft_bounds_output/` (3 configs × 4 ratios)

**What Works:**
- Correct PMP structure with objective, normalization, positivity matrices
- DampedRational prefactor handling with poles at x=0 and x=-4m²
- Variable elimination integrated (Method A)
- Spectral functions v_ℓ^{(p,q)}(x) correctly encoded as polynomials

**Minor Issues:**
- Line 363-366: Uses `Fraction(1,3)` as e^{-1} approximation, then computes properly
  - This is okay but code is confusing (sets e_inv then immediately overwrites it)
- Line 390-419: `_compute_e_inverse()` computes e then 1/e
  - Could be more direct; unnecessary double computation

**Files Generated:**
```
pmp_g10_over_g00_{positivity_only, with_su_crossing, with_full_crossing}.json
pmp_g20_over_g00_{positivity_only, with_su_crossing, with_full_crossing}.json
pmp_g01_over_g00_{positivity_only, with_su_crossing, with_full_crossing}.json
pmp_g20_over_g10_{positivity_only, with_su_crossing, with_full_crossing}.json
```

---

### Task 3: Verification of ≥5 Examples ⚠️ **CRITICAL GAP**

**Current Implementation:**
- **File:** `run_bounds.py`, function `run_consistency_checks()` (lines 94-283)

**Examples Verified:**
1. ✅ **Polynomial amplitude** s² + u² (line 241-270)
   - All null constraints satisfied exactly
2. ⚠️ **Scalar exchange** 1/(s-M²) + 1/(t-M²) + 1/(u-M²) (line 163-235)
   - Shows expected truncation errors (not a bug, just verification of truncation theory)

**Missing Examples (need 3+ more):**
3. ❌ Contact interaction M(s,t) = c₀ (should trivially satisfy all constraints)
4. ❌ Linear amplitude M(s,t) = a·s + b·t with crossing symmetry
5. ❌ t-channel exchange only: 1/(t-M²)
6. ❌ Simple EFT: M(s,t) = g₀₀ + g₁₀·s + g₀₁·t + g₂₀·s² + ...
7. ❌ Mixed polynomial-rational test case

**What's Missing:**
- No systematic test suite (`test/` directory has C++ tests only, no Python tests)
- No validation against known numerical results from papers
- No comparison with actual SDPB-computed bounds
- Missing `examples/` directory validation (has `generate_examples.py` but no verification)

**Action Needed:**
Add to `run_bounds.py` or create new `test_examples.py`:
```python
def test_contact_interaction():
    # M(s,t) = 1 → W_{0,0}=1, all others=0
def test_linear_crossing_symmetric():
    # M(s,t) = s + u = 4m² - t
def test_t_channel_exchange():
    # M(s,t) = 1/(t-M²) (not crossing-symmetric)
```

---

### Task 4: Two-Sided Bounds (Sinha-Zahed eq. 3) ❌ **NOT IMPLEMENTED**

**What's Missing:**
- **No implementation of equation (3)** from Sinha & Zahed (2021)
- The paper likely defines two-sided bounds of the form:
  ```
  Lower(g_i/g_j) ≤ g_i/g_j ≤ Upper(g_i/g_j)
  ```
  derived from dispersion relations with subtractions

**Current State:**
- Only **lower bounds** via optimization: maximize -W_num / W_den
- No upper bounds: would need to flip objective (minimize -W_num = maximize W_num)
- No two-sided simultaneous bounds
- No implementation of subtracted dispersion relations

**Reference in Papers:**
From README.md line 195-203 and Sinha-Zahed paper:
- Equation (3) likely defines bounds via crossing-symmetric dispersion relation:
  ```
  ∫₄ₘ² dμ ρ(μ) [kernel] = known polynomial in Wilson coefficients
  ```
  with ρ(μ) ≥ 0 providing inequalities

**What Needs to Be Added:**

1. **Create new file:** `eft_bounds/two_sided_bounds.py`
   ```python
   def get_upper_bound_pmp(num_idx, den_idx, ...):
       """Generate PMP for UPPER bound: minimize W_num/W_den"""
       # Flip sign of objective from generate_pmp_for_ratio_bound
   
   def get_two_sided_bounds(num_idx, den_idx, ...):
       """Return (lower_pmp, upper_pmp) pair"""
   
   def sinha_zahed_eq3_constraint(...):
       """Implement the specific constraint from S-Z equation (3)"""
       # This is THE KEY MISSING PIECE
       # Need to read the paper to understand eq.(3) structure
   ```

2. **Update `run_bounds.py`:**
   - Add generation of upper-bound PMP files
   - Generate pairs: `pmp_g10_over_g00_lower.json` and `pmp_g10_over_g00_upper.json`

3. **Implement Extremal-EFT style optimization:**
   - Current code already follows Caron-Huot & Duong dual formulation
   - Need to explicitly match their Section 3 optimization structure
   - May need crossing-symmetric basis change: (s,t) → (σ₂, σ₃) where σ₂=st+tu+us, σ₃=stu

---

## Specific Code Issues Found

### Bug 1: Misleading Comment in `get_crossing_symmetric_null_constraints`
**File:** `null_constraints.py`, lines 370-388  
**Issue:** Claims to implement "explicit null constraints from Sinha-Zahed" but just combines s↔u and s↔t  
**Fix:** Either implement the paper's constraints properly or correct the comment

### Bug 2: Inefficient e^{-1} Calculation
**File:** `pmp_generator.py`, lines 363-419  
**Issue:** Computes e then 1/e instead of directly computing e^{-1}  
**Fix:**
```python
def _compute_e_inverse(precision: int) -> str:
    from decimal import Decimal, getcontext
    getcontext().prec = precision + 20
    # Direct series: e^{-1} = Σ_{n=0}^∞ (-1)^n / n!
    result = Decimal(0)
    term = Decimal(1)
    for n in range(1, precision + 100):
        term = term / Decimal(n)
        result += term if n % 2 == 0 else -term
        if abs(term) < Decimal(10) ** (-(precision + 10)):
            break
    result += Decimal(1)  # n=0 term
    return format(result, f'.{precision}f').rstrip('0').rstrip('.')
```

### Bug 3: Missing Test Coverage
**File:** None (that's the bug)  
**Issue:** No Python unit tests for the EFT pipeline  
**Fix:** Create `eft_bounds/test_eft_bounds.py` with pytest tests

---

## Minimal Files to Edit (Priority Order)

### High Priority (Required for Task Completion)

1. **Create `eft_bounds/two_sided_bounds.py`** (NEW FILE)
   - Implement Sinha-Zahed equation (3) constraints
   - Add upper bound generation
   - Add two-sided bound solver

2. **Edit `eft_bounds/run_bounds.py`** (ADD ~100 lines)
   - Add 3 more verification examples in `run_consistency_checks()`
   - Add upper bound generation to `generate_bound_series()`
   - Add test_contact_interaction(), test_linear_crossing(), test_t_channel()

3. **Edit `eft_bounds/null_constraints.py`** (MODIFY ~20 lines)
   - Fix `get_crossing_symmetric_null_constraints()` docstring (line 370)
   - Implement actual Sinha-Zahed eq.(3) or remove misleading claims

### Medium Priority (Code Quality)

4. **Edit `eft_bounds/pmp_generator.py`** (MODIFY ~30 lines)
   - Fix `_compute_e_inverse()` efficiency (line 390-419)
   - Clean up confusing e_inv initialization (line 363-366)

5. **Create `eft_bounds/test_eft_bounds.py`** (NEW FILE)
   - Pytest suite for all verification examples
   - Regression tests for null constraint counts
   - PMP JSON schema validation

### Low Priority (Nice to Have)

6. **Edit `eft_bounds/README.md`** (ADD ~50 lines)
   - Document the missing Task 4 implementation
   - Add explicit example outputs from verification
   - Clarify relationship to Sinha-Zahed eq.(3)

---

## Actionable Recommendations

### To Complete Task 3 (Verification)
**Effort:** ~2-4 hours  
**Files:** `run_bounds.py` (+100 lines)

```python
# Add to run_consistency_checks() around line 270

# Example 3: Contact interaction
print("Example 3: Contact interaction M = c₀")
wilson_contact = {(0,0): Fraction(1)}
verify_amplitude_satisfies_constraints(wilson_contact, constraints, "contact")

# Example 4: Linear crossing-symmetric
print("Example 4: Linear M = s + u = 4m² - t")
wilson_linear = {(1,0): Fraction(1), (0,1): Fraction(-1), (0,0): 4*m_sq}
verify_amplitude_satisfies_constraints(wilson_linear, constraints, "linear")

# Example 5: Quadratic s↔u symmetric
print("Example 5: M = s² + t² + u²")
# Implement expansion and verify
```

### To Complete Task 4 (Two-Sided Bounds)
**Effort:** ~8-16 hours (requires reading Sinha-Zahed paper carefully)  
**Files:** Create `two_sided_bounds.py` (~300 lines), edit `run_bounds.py` (+50 lines)

**Critical Unknown:** What exactly is Sinha-Zahed equation (3)?
- Need to read "Crossing Symmetric Dispersion Relations in QFTs" paper
- File exists: `Sinha 和 Zahed - 2021 - Crossing Symmetric Dispersion Relations in QFTs.pdf`
- Equation (3) is likely the key dispersion relation with bounds

**Steps:**
1. Extract equation (3) from PDF (use `pdftotext` or manual reading)
2. Implement the constraint structure
3. Generate upper-bound PMPs by flipping objective sign
4. Verify bounds are consistent (lower ≤ upper)

### Quick Wins (1-2 hours each)

1. **Fix e^{-1} computation** (30 min)
2. **Add 3 verification examples** (2 hours)
3. **Create pytest test file** (2 hours)
4. **Document Task 4 gap in README** (30 min)

---

## Summary

**Pipeline Completeness: 58%**

**Strengths:**
- Solid mathematical foundation
- Clean code structure
- Null constraints correctly derived for s↔u and full S₃
- PMP generation works correctly
- Good documentation in README and docstrings

**Critical Gaps:**
1. **Task 4 entirely missing** (0% complete) → Need two-sided bounds
2. **Task 3 incomplete** (40% complete) → Need 3 more verified examples
3. **Sinha-Zahed eq.(3) not implemented** → Core feature missing

**Estimated Work to Complete:**
- Task 3 (examples): 2-4 hours
- Task 4 (two-sided bounds): 8-16 hours
- Code quality fixes: 2-3 hours
- **Total:** ~12-23 hours of focused development

**Recommended Next Steps:**
1. Read Sinha-Zahed paper equation (3) to understand two-sided bounds structure
2. Add 3 simple verification examples (contact, linear, quadratic)
3. Implement upper bound generation (flip objective)
4. Create two_sided_bounds.py with equation (3) implementation
5. Add pytest test suite

**Files to Create/Edit:**
- **NEW:** `two_sided_bounds.py` (~300 lines)
- **NEW:** `test_eft_bounds.py` (~200 lines)
- **EDIT:** `run_bounds.py` (+150 lines)
- **EDIT:** `null_constraints.py` (~20 line fix)
- **EDIT:** `pmp_generator.py` (~30 line cleanup)
- **EDIT:** `README.md` (+50 lines documentation)

---
