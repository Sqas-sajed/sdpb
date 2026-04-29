# SDPB Workflow Guide: EFT Bounds from CSDR Dispersion Relations

This guide explains how to use the tools in this repository to produce plots
of allowed regions for Wilson coefficients, using SDPB on Windows with Docker.
No prior coding experience is assumed.

The official SDPB documentation is in [`sdpbUsage.md`](../sdpbUsage.md) and
[`docs/Docker.md`](../docs/Docker.md) — read those alongside this guide.

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
| (1, 2) | W_{−1,2}    | **null constraint** (= 0 by crossing) |

A "plot of allowed (g̃₃, g̃₄)" means the 2D region:

    g̃₃ = W_{0,1} / W_{1,0},    g̃₄ = W_{2,0} / W_{1,0}

---

## Setup on Windows with Docker

### Step 1: Install Docker Desktop

Download from https://www.docker.com/products/docker-desktop and install.
After installation, open "Docker Desktop" and wait for it to start (green
icon in the taskbar).

### Step 2: Pull the official SDPB Docker image

Open a **PowerShell** window and run:

```powershell
docker pull bootstrapcollaboration/sdpb:master
```

This downloads the official SDPB image (see [`docs/Docker.md`](../docs/Docker.md)).
You can also use a specific release tag, e.g. `bootstrapcollaboration/sdpb:3.1.0`.

Test it works:
```powershell
docker run --rm bootstrapcollaboration/sdpb:master sdpb --help
```
You should see the SDPB option list.

### Step 3: Get Python (if you don't have it)

Download Python 3.9+ from https://www.python.org/downloads/.
During installation, check the box **"Add Python to PATH"**.

---

## Extract the required Python files

You need **five files** from `eft_bounds/` in this repository. Copy them to a
folder on your computer, e.g. `C:\sdpb_eft\`.

Your folder structure must look like this:

```
C:\sdpb_eft\
    generate_pmp.py              ← the script you run from the command line
    eft_bounds\
        __init__.py
        csdr.py
        physics.py
        pmp_generator.py
        generate_pmp.py          ← keep a copy here too (for internal imports)
```

> **Where to find these files:** In the GitHub repository, go to the
> `eft_bounds/` directory, click each `.py` file, then click the download
> button (or copy the raw text into a new file with the same name).

---

## Quick test: generate one PMP file

Open **PowerShell** and navigate to your working folder:
```powershell
cd C:\sdpb_eft
```

List available operators (useful to understand the notation):
```powershell
python generate_pmp.py --d 4 --K 4 --obj 1 0 --norm 1 0 --list-pairs
```

Generate an upper-bound PMP file for g̃₃ = W_{0,1}/W_{1,0}:
```powershell
python generate_pmp.py --obj 1 1 --norm 1 0 --d 4 --K 4 --max-spin 10 --delta0 40 --direction upper --output upper_g3.json
```

This creates `upper_g3.json` — the PMP input file ready for SDPB.

> The `--precision` flag in `generate_pmp.py` is **decimal digits** of output
> coefficients (default 200).  This is different from SDPB's `--precision`,
> which is in **bits** (1024 bits ≈ 308 decimal digits).  They do not need
> to match exactly — SDPB's precision must be large enough to handle the
> coefficient values.

---

## Run SDPB with Docker (step by step)

The workflow has two stages (as described in [`sdpbUsage.md`](../sdpbUsage.md)):

1. **`pmp2sdp`**: convert the JSON PMP file into SDPB's internal binary SDP format.
2. **`sdpb`**: run the semidefinite program solver on the SDP.

Docker mounts your local folder into the container at the path
`/usr/local/share/sdpb/` (see [`docs/Docker.md`](../docs/Docker.md)).

Suppose your PMP files are in `C:\sdpb_eft\`.

### Step 1: Convert PMP → SDP with `pmp2sdp`

```powershell
docker run --rm `
    -v "C:/sdpb_eft/:/usr/local/share/sdpb/" `
    bootstrapcollaboration/sdpb:master `
    mpirun --allow-run-as-root -n 4 pmp2sdp --precision 1024 `
        -i /usr/local/share/sdpb/upper_g3.json `
        -o /usr/local/share/sdpb/sdp_upper_g3
```

- `-v "C:/sdpb_eft/:/usr/local/share/sdpb/"` — mounts your folder.
- `mpirun --allow-run-as-root -n 4` — runs with 4 CPU cores inside Docker.
  Change `4` to match the number of cores on your machine.
- `--precision 1024` — working precision in **bits** (1024 is a safe default).
- `-i …` — input JSON file (inside the container's mounted path).
- `-o …` — output directory for the SDP binary files.

This creates `C:\sdpb_eft\sdp_upper_g3\` containing `control.json`,
`objectives.json`, and binary block files.

### Step 2: Run SDPB

```powershell
docker run --rm `
    -v "C:/sdpb_eft/:/usr/local/share/sdpb/" `
    bootstrapcollaboration/sdpb:master `
    mpirun --allow-run-as-root -n 4 sdpb --precision=1024 `
        -s /usr/local/share/sdpb/sdp_upper_g3 `
        -o /usr/local/share/sdpb/out_upper_g3 `
        -c /usr/local/share/sdpb/out_upper_g3/ck
```

- `-s …` — path to the SDP directory produced by `pmp2sdp`.
- `-o …` — output directory for the solver result.
- `-c …` — checkpoint directory (used to resume interrupted runs).

This may take a few minutes. When complete, `C:\sdpb_eft\out_upper_g3\` is created.

### Step 3: Read the result

```powershell
type C:\sdpb_eft\out_upper_g3\out.txt
```

The output looks like (see the example in [`sdpbUsage.md`](../sdpbUsage.md)):
```
terminateReason = "found primal-dual optimal solution";
primalObjective = 1.23456...;
dualObjective   = 1.23456...;
dualityGap      = 3.5e-31;
primalError     = 2.8e-309;
dualError       = 7.7e-305;
Solver runtime  = 42;
```

The **upper bound** on g̃₃ is the value of `primalObjective`.
A `terminateReason` of `"found primal-dual optimal solution"` means the
solver converged. If it says `"dual infeasible"` the problem may be
infeasible (no allowed region exists for those parameters).

> **Note on newly created files:** Files written by Docker may be owned by
> `root`.  If you cannot delete them in Windows Explorer, run:
> `docker run --rm -v "C:/sdpb_eft/:/usr/local/share/sdpb/" bootstrapcollaboration/sdpb:master rm -rf /usr/local/share/sdpb/out_upper_g3`

---

## Produce a 2D allowed-region plot of (g̃₃, g̃₄)

The allowed region is a convex set in the (g̃₃, g̃₄) plane.

### Step 1: Generate 4 PMP files

Run these four commands to get the extreme bounds:

```powershell
# Upper bound on g̃₃ = W_{0,1}/W_{1,0}:
python generate_pmp.py --obj 1 1 --norm 1 0 --d 4 --K 4 --max-spin 10 --delta0 40 --direction upper -o ub_g3.json

# Lower bound on g̃₃:
python generate_pmp.py --obj 1 1 --norm 1 0 --d 4 --K 4 --max-spin 10 --delta0 40 --direction lower -o lb_g3.json

# Upper bound on g̃₄ = W_{2,0}/W_{1,0}:
python generate_pmp.py --obj 2 0 --norm 1 0 --d 4 --K 4 --max-spin 10 --delta0 40 --direction upper -o ub_g4.json

# Lower bound on g̃₄:
python generate_pmp.py --obj 2 0 --norm 1 0 --d 4 --K 4 --max-spin 10 --delta0 40 --direction lower -o lb_g4.json
```

### Step 2: Run pmp2sdp + sdpb for each file

Repeat the two Docker commands from the previous section for each of the 4
JSON files.  Use different directory names each time, e.g.
`sdp_ub_g3`, `sdp_lb_g3`, etc.  Example for `lb_g3.json`:

```powershell
docker run --rm -v "C:/sdpb_eft/:/usr/local/share/sdpb/" bootstrapcollaboration/sdpb:master mpirun --allow-run-as-root -n 4 pmp2sdp --precision 1024 -i /usr/local/share/sdpb/lb_g3.json -o /usr/local/share/sdpb/sdp_lb_g3

docker run --rm -v "C:/sdpb_eft/:/usr/local/share/sdpb/" bootstrapcollaboration/sdpb:master mpirun --allow-run-as-root -n 4 sdpb --precision=1024 -s /usr/local/share/sdpb/sdp_lb_g3 -o /usr/local/share/sdpb/out_lb_g3 -c /usr/local/share/sdpb/out_lb_g3/ck
```

### Step 3: Collect the results

Open each `out.txt` and record the `primalObjective` value:

| File | What it is | Value (example) |
|------|-----------|-----------------|
| `out_ub_g3/out.txt` | upper bound g̃₃ | `3.12` |
| `out_lb_g3/out.txt` | lower bound −g̃₃ | `2.45` → g̃₃ ≥ −2.45 |
| `out_ub_g4/out.txt` | upper bound g̃₄ | `4.50` |
| `out_lb_g4/out.txt` | lower bound −g̃₄ | `1.20` → g̃₄ ≥ −1.20 |

> **For lower bounds:** when `--direction lower` is used, `generate_pmp.py`
> negates the objective.  The bound on g̃ is therefore
> **negative** of `primalObjective`.  For example, if `primalObjective = 2.45`
> in the lower-bound run, the true lower bound is g̃₃ ≥ −2.45.

### Step 4: Plot in Python

Install matplotlib (open PowerShell):
```powershell
pip install matplotlib
```

Create a file `plot_region.py` in `C:\sdpb_eft\`:

```python
import matplotlib.pyplot as plt
import matplotlib.patches as patches

# ── Replace these values with what you read from out.txt ──
ub_g3  =  3.12   # primalObjective from out_ub_g3/out.txt
lb_g3  = -2.45   # NEGATIVE of primalObjective from out_lb_g3/out.txt
ub_g4  =  4.50   # primalObjective from out_ub_g4/out.txt
lb_g4  = -1.20   # NEGATIVE of primalObjective from out_lb_g4/out.txt

fig, ax = plt.subplots(figsize=(6, 5))
rect = patches.Rectangle(
    (lb_g3, lb_g4),
    ub_g3 - lb_g3,
    ub_g4 - lb_g4,
    linewidth=2, edgecolor='blue', facecolor='lightblue', alpha=0.5,
    label='Allowed region (bounding box)'
)
ax.add_patch(rect)
ax.set_xlim(lb_g3 - 0.5, ub_g3 + 0.5)
ax.set_ylim(lb_g4 - 0.5, ub_g4 + 0.5)
ax.set_xlabel(r'$\tilde{g}_3 = W_{0,1}/W_{1,0}$', fontsize=12)
ax.set_ylabel(r'$\tilde{g}_4 = W_{2,0}/W_{1,0}$', fontsize=12)
ax.set_title('Allowed EFT region (CSDR, K=4, d=4, $\\delta_0$=40)', fontsize=11)
ax.legend()
plt.tight_layout()
plt.savefig('allowed_region.png', dpi=150)
print('Saved: allowed_region.png')
plt.show()
```

Run it:
```powershell
python plot_region.py
```

This produces a rectangular **bounding box** of the allowed region.  For a
more accurate boundary (convex hull), you would run SDPB for many more
directions and connect the boundary points — the same approach as Figure 8
in "Extremal EFT".

---

## One-line combined commands (copy-paste ready)

These combine pmp2sdp and sdpb on a single line (no line continuation):

```powershell
docker run --rm -v "C:/sdpb_eft/:/usr/local/share/sdpb/" bootstrapcollaboration/sdpb:master mpirun --allow-run-as-root -n 4 pmp2sdp --precision 1024 -i /usr/local/share/sdpb/upper_g3.json -o /usr/local/share/sdpb/sdp_upper_g3

docker run --rm -v "C:/sdpb_eft/:/usr/local/share/sdpb/" bootstrapcollaboration/sdpb:master mpirun --allow-run-as-root -n 4 sdpb --precision=1024 -s /usr/local/share/sdpb/sdp_upper_g3 -o /usr/local/share/sdpb/out_upper_g3 -c /usr/local/share/sdpb/out_upper_g3/ck
```

---

## Understanding the parameters

### Parameters for `generate_pmp.py`

| Flag | Meaning | Notes |
|------|---------|-------|
| `--K` | Max spectral power 2n+m | More operators → tighter bounds, slower |
| `--d` | Spacetime dimension | Physical: d=4 for 4D |
| `--delta0` | IR cutoff δ₀ | s₁ starts at δ₀; matches the paper |
| `--max-spin` | Max even spin ℓ | More spins → tighter bounds, slower |
| `--precision` | Decimal digits in output | 200 is safe; use 50 for quick tests |
| `--obj n m` | Wilson coeff to bound | (n, m) pair with n ≥ m ≥ 0 |
| `--norm n m` | Normalization coeff | Set this W to 1 |
| `--direction` | `upper` or `lower` | Which extremum to compute |

### Parameters for `pmp2sdp` and `sdpb`

These take `--precision` in **bits** (not decimal digits):

| `--precision` bits | Decimal digits | Recommendation |
|--------------------|----------------|----------------|
| 512  | ~154 | Minimum safe for 50-digit PMP input |
| 1024 | ~308 | Safe default for 200-digit PMP input |
| 2048 | ~616 | For high-precision production runs |

Use `sdpb --help` (via Docker) to see all available options.

---

## Relation to Extremal EFT (Caron-Huot & Duong 2021)

Figure 8 of "Extremal EFT" uses **fixed-t dispersion relations**, not CSDR.
The two approaches give different operator bases and different bounds.

This code uses the CSDR formula (Sinha & Zahed 2021):

    W_{n−m,m} = ⟨ D^{(n,m)}_{ℓ,α} · C^α_ℓ(1) · (2ℓ+d−3) / s₁^{2n+m} ⟩

where D is spin-dependent and α = (d−3)/2.

---

## Troubleshooting

**`terminateReason = "dual infeasible"`**
→ The SDP may be infeasible.  Try increasing `--K` or `--max-spin`.

**The bound value looks wrong / very large**
→ Check `--delta0` matches your physical setup; changing δ₀ rescales W_{n,m}.
→ Also check `--direction`: lower-bound runs negate the objective.

**Docker command fails with `path not found`**
→ Use forward slashes in the `-v` flag path on Windows:
   `-v "C:/sdpb_eft/:/usr/local/share/sdpb/"` (not backslashes).

**Files created by Docker cannot be deleted**
→ They are owned by root inside the container. Delete them with:
   `docker run --rm -v "C:/sdpb_eft/:/usr/local/share/sdpb/" bootstrapcollaboration/sdpb:master rm -rf /usr/local/share/sdpb/sdp_upper_g3`

**`'A was not numerically HPD'` error from SDPB**
→ Increase `--precision` in `pmp2sdp` and `sdpb` (e.g. try 2048 instead of 1024),
   as described in [`sdpbUsage.md`](../sdpbUsage.md).

**Out-of-memory error**
→ Reduce parallel cores (`-n 2` instead of `-n 4`) or add
  `--maxSharedMemory=2G` to the `sdpb` command.

---

## Example files provided

The `eft_bounds/examples/` directory contains pre-generated PMP files
(K=4, d=4, δ₀=40, max_spin=10, precision=50 decimal digits):

| File | Computes |
|------|---------|
| `upper_W01_over_W10_K4_d4_delta40.json` | Upper bound on g̃₃ = W_{0,1}/W_{1,0} |
| `lower_W01_over_W10_K4_d4_delta40.json` | Lower bound on g̃₃ |
| `upper_W20_over_W10_K4_d4_delta40.json` | Upper bound on g̃₄ = W_{2,0}/W_{1,0} |
| `lower_W20_over_W10_K4_d4_delta40.json` | Lower bound on g̃₄ |

**To use them directly** (no Python needed), copy the `examples/` folder to
`C:\sdpb_eft\examples\` and run (example for upper g̃₃):

```powershell
docker run --rm -v "C:/sdpb_eft/:/usr/local/share/sdpb/" bootstrapcollaboration/sdpb:master mpirun --allow-run-as-root -n 4 pmp2sdp --precision 1024 -i /usr/local/share/sdpb/examples/upper_W01_over_W10_K4_d4_delta40.json -o /usr/local/share/sdpb/sdp_ub_g3

docker run --rm -v "C:/sdpb_eft/:/usr/local/share/sdpb/" bootstrapcollaboration/sdpb:master mpirun --allow-run-as-root -n 4 sdpb --precision=1024 -s /usr/local/share/sdpb/sdp_ub_g3 -o /usr/local/share/sdpb/out_ub_g3 -c /usr/local/share/sdpb/out_ub_g3/ck

type C:\sdpb_eft\out_ub_g3\out.txt
```

The `primalObjective` line gives the upper bound on g̃₃.
