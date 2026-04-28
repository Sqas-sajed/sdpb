# SDPB Workflow Guide: EFT Bounds from CSDR Dispersion Relations

This guide explains how to use the tools in this repository to produce plots
of allowed regions for Wilson coefficients, using SDPB on Windows Docker.
No prior coding experience is assumed.

---

## What you are doing (physics summary)

You want to find the **allowed region** of normalized Wilson coefficients in
an Effective Field Theory (EFT). The idea is:

1. The EFT amplitude has coefficients W_{p,q} that are "heavy averages" of
   spectral data via CSDR (Crossing-Symmetric Dispersion Relations).
2. Unitarity and crossing symmetry constrain which (W_{p,q}) values are
   physically possible.
3. SDPB (Semidefinite Program Bootstrap) computes the **extremal boundary**
   of this allowed region.

The Wilson coefficients are labeled W_{n−m, m}, indexed by (n, m):

| (n, m) | Wilson coeff | Physics analog |
|--------|-------------|----------------|
| (1, 0) | W_{1,0}     | g₂ (sub-leading forward) |
| (1, 1) | W_{0,1}     | g₃ (crossing) |
| (2, 0) | W_{2,0}     | g₄ (sub-leading) |
| (1, 2) | W_{-1,2}    | **null constraint** (= 0 by crossing) |

A "plot of allowed (g̃₃, g̃₄)" means the 2D region:

    g̃₃ = W_{0,1} / W_{1,0},    g̃₄ = W_{2,0} / W_{1,0}

---

## Setup on Windows with Docker

### Step 1: Install Docker Desktop

Download from https://www.docker.com/products/docker-desktop and install.
After installation, open "Docker Desktop" and wait for it to start.

### Step 2: Pull the SDPB Docker image

Open a **PowerShell** or **Command Prompt** window and run:

```powershell
docker pull davidsimmons1979/sdpb:latest
```

Wait for the download to complete (it may take several minutes).

Test it works:
```powershell
docker run --rm davidsimmons1979/sdpb:latest sdpb --help
```
You should see the SDPB help text.

### Step 3: Get Python (if you don't have it)

Download Python 3.9+ from https://www.python.org/downloads/.
During installation, check the box **"Add Python to PATH"**.

---

## Extract the required Python files

You need **four files** from this repository. Copy them to a folder on your
computer, e.g. `C:\sdpb_eft\`.

```
eft_bounds/
    __init__.py
    csdr.py
    physics.py
    pmp_generator.py
    generate_pmp.py
```

Create a folder `C:\sdpb_eft\eft_bounds\` and put all five `.py` files there.
Also copy `generate_pmp.py` to `C:\sdpb_eft\` (one level above).

Your folder should look like:
```
C:\sdpb_eft\
    generate_pmp.py          ← the main script you run
    eft_bounds\
        __init__.py
        csdr.py
        physics.py
        pmp_generator.py
        generate_pmp.py      ← also here (for import)
```

---

## Quick test: generate one PMP file

Open **PowerShell** and navigate to your folder:
```powershell
cd C:\sdpb_eft
```

List available operators (to understand the notation):
```powershell
python generate_pmp.py --d 4 --K 4 --obj 1 0 --norm 1 0 --list-pairs
```

Generate an upper bound on g̃₃ = W_{0,1}/W_{1,0}:
```powershell
python generate_pmp.py --obj 1 1 --norm 1 0 --d 4 --K 4 --max-spin 10 --delta0 40 --direction upper --output upper_g3.json
```

This creates `upper_g3.json` — the PMP input file for SDPB.

---

## Run SDPB with Docker (step by step)

### Step 1: Convert PMP → SDP binary format

SDPB needs its own binary format. The tool `pmp2sdp` does this conversion.

```powershell
# From C:\sdpb_eft\, run:
docker run --rm -v "${PWD}:/work" davidsimmons1979/sdpb:latest `
    pmp2sdp --input /work/upper_g3.json --output /work/sdp_upper_g3/
```

This creates a folder `sdp_upper_g3/` with binary files.

### Step 2: Run SDPB

```powershell
docker run --rm -v "${PWD}:/work" davidsimmons1979/sdpb:latest `
    sdpb --sdpDir /work/sdp_upper_g3/ --outDir /work/out_upper_g3/
```

This may take a few minutes. When done, the result is in `out_upper_g3/`.

### Step 3: Read the result

```powershell
type out_upper_g3\out.txt
```

Look for the line:
```
primalObjective = 1.23456...
```

This number is the **upper bound on g̃₃**.

---

## Produce a 2D allowed-region plot of (g̃₃, g̃₄)

The allowed region is a convex set in the (g̃₃, g̃₄) plane. To plot it,
you need to find the boundary by running SDPB many times with different
"scanning" directions.

### Strategy: scan the boundary

The boundary can be traced by computing, for many angles θ:

    max (cos θ · g̃₃ + sin θ · g̃₄)

For each θ, this is one SDPB run where:
- Objective = cos(θ) * W_{0,1} + sin(θ) * W_{2,0}
- Normalization = W_{1,0}

Currently `generate_pmp.py` supports single-ratio bounds. For the 2D scan,
a simple approach is to bound one coordinate while fixing the other.

**Simple 4-point approach:**

Run these 4 bounds to get approximate corners of the allowed region:

```powershell
# Upper bound on g̃₃:
python generate_pmp.py --obj 1 1 --norm 1 0 --d 4 --K 4 --max-spin 10 --delta0 40 --direction upper -o ub_g3.json

# Lower bound on g̃₃:
python generate_pmp.py --obj 1 1 --norm 1 0 --d 4 --K 4 --max-spin 10 --delta0 40 --direction lower -o lb_g3.json

# Upper bound on g̃₄:
python generate_pmp.py --obj 2 0 --norm 1 0 --d 4 --K 4 --max-spin 10 --delta0 40 --direction upper -o ub_g4.json

# Lower bound on g̃₄:
python generate_pmp.py --obj 2 0 --norm 1 0 --d 4 --K 4 --max-spin 10 --delta0 40 --direction lower -o lb_g4.json
```

For each file, run `pmp2sdp` and `sdpb` as described above.

### Collect results

After running all 4:
- g̃₃ ∈ [lb_g3, ub_g3]
- g̃₄ ∈ [lb_g4, ub_g4]

The allowed region is **inside** these intervals.

### Plot in Python (on your computer)

Once you have the bound values (e.g., lb_g3 = -2.5, ub_g3 = 3.0, etc.):

```python
import matplotlib.pyplot as plt
import matplotlib.patches as patches

lb_g3, ub_g3 = -2.5, 3.0   # replace with your values
lb_g4, ub_g4 = -1.0, 4.0

fig, ax = plt.subplots()
rect = patches.Rectangle(
    (lb_g3, lb_g4),
    ub_g3 - lb_g3, ub_g4 - lb_g4,
    linewidth=2, edgecolor='blue', facecolor='lightblue', alpha=0.5
)
ax.add_patch(rect)
ax.set_xlim(lb_g3 - 1, ub_g3 + 1)
ax.set_ylim(lb_g4 - 1, ub_g4 + 1)
ax.set_xlabel(r'$\tilde{g}_3 = W_{0,1}/W_{1,0}$')
ax.set_ylabel(r'$\tilde{g}_4 = W_{2,0}/W_{1,0}$')
ax.set_title('Allowed EFT region (CSDR bounds, K=4, d=4)')
plt.tight_layout()
plt.savefig('allowed_region.png', dpi=150)
plt.show()
```

Install matplotlib if needed: `pip install matplotlib`

---

## Understanding the parameters

| Parameter | Meaning | Effect |
|-----------|---------|--------|
| `--K` | Max spectral power 2n+m | More operators → tighter bounds, slower |
| `--d` | Spacetime dimension | Physical: d=4 for 4D, d=3 for 3D CFT, etc. |
| `--delta0` | IR cutoff for s₁ | s₁ starts at delta0; larger = less constraining |
| `--max-spin` | Max spin ℓ summed | More spins → tighter bounds, slower |
| `--precision` | Decimal digits | 200 is enough; reduce to 50 for speed tests |

**Recommended values for a quick test:**
```
K=4, max_spin=6, delta0=40, d=4, precision=50
```
This runs in seconds on a modern laptop.

**Recommended values for publication-quality bounds:**
```
K=10, max_spin=20, delta0=40, d=4, precision=200
```
This may take hours. Run on a computing cluster if available.

---

## Relation to Extremal EFT (Caron-Huot & Duong 2021)

Figure 8 of "Extremal EFT" uses **fixed-t dispersion relations**, not CSDR.
The two approaches give different operator bases and different bounds.

In Extremal EFT notation:
- g₃ corresponds to the amplitude coefficient at order s·t
- g₄ corresponds to the amplitude coefficient at order s²

These map to CSDR Wilson coefficients W_{n−m,m} through the equations in the
user's LaTeX notes. The exact mapping depends on the normalization convention.

This code uses the CSDR formula (Sinha & Zahed 2021) with:
    W_{n−m,m} = ⟨ D^{(n,m)}_{ℓ,α} · C^α_ℓ(1) · (2ℓ+d−3) / s₁^{2n+m} ⟩

where ⟨·⟩ denotes the spectral average with positivity-preserving measure.

---

## Troubleshooting

**"pmp2sdp: command not found"**
→ The Docker image may use a different binary name. Try:
```powershell
docker run --rm davidsimmons1979/sdpb:latest ls /usr/local/bin/
```
And use whichever executable converts PMP files.

**"SDPB did not converge"**
→ The problem may be infeasible (no solution exists). Try:
- Increasing `--K` (more decision variables)
- Increasing `--max-spin`
- Using `--direction lower` instead of `upper`

**"The bound value seems wrong"**
→ Check that `--delta0` matches the physical setup. The cutoff δ₀ affects
the numerical values of W_{n,m}.

**File size is large**
→ Reduce `--precision` to 50 for tests. The SDPB solver uses arbitrary-precision
arithmetic internally, so 50 digits in the input is usually fine.

---

## Example files provided

The `eft_bounds/examples/` directory contains pre-generated PMP files:

| File | What it computes |
|------|-----------------|
| `upper_W01_over_W10_K4_d4_delta40.json` | Upper bound on W_{0,1}/W_{1,0} (∝ g̃₃) |
| `lower_W01_over_W10_K4_d4_delta40.json` | Lower bound on W_{0,1}/W_{1,0} |
| `upper_W20_over_W10_K4_d4_delta40.json` | Upper bound on W_{2,0}/W_{1,0} (∝ g̃₄) |
| `lower_W20_over_W10_K4_d4_delta40.json` | Lower bound on W_{2,0}/W_{1,0} |

Parameters: K=4, d=4, δ₀=40, max_spin=10, precision=50.

Run each through `pmp2sdp` + `sdpb` to get the bound values.
