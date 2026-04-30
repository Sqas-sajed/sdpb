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
icon in the taskbar).  Make sure the **"Use WSL 2 based engine"** setting is
enabled (recommended, default on most Windows installs).

### Step 2: Switch Docker Desktop to Linux containers mode

The SDPB image is a **Linux** container.  Docker Desktop on Windows can run
either Linux containers or Windows containers, but **not both at the same
time**.  It must be in **Linux containers mode** before any `docker run`
command will work; if it is in Windows containers mode every `docker run`
fails with `exec format error` regardless of the `--platform` flag.

To check and switch:

1. Find the Docker Desktop icon in the Windows taskbar notification area
   (bottom-right).
2. **Right-click** the icon.
3. If the menu shows **"Switch to Linux containers…"**, click it and wait for
   Docker to restart.  (If it shows "Switch to Windows containers…" instead,
   you are already in Linux mode — nothing to do.)

After switching, verify the setting is correct:
Open Docker Desktop → **Settings** → **General** → confirm
**"Use the WSL 2 based engine"** is ticked.

### Step 3: Pull the official SDPB Docker image

Open a **PowerShell** window and run:

```powershell
docker pull --platform linux/amd64 bootstrapcollaboration/sdpb:3.1.0
```

The `--platform linux/amd64` flag ensures the x86-64 Linux image is pulled
explicitly.  We use the stable `3.1.0` release tag; avoid the `master` tag
as it may have binary-format issues that cause `exec format error` on some
platforms.

Test it works:

```powershell
docker run --rm --platform linux/amd64 bootstrapcollaboration/sdpb:3.1.0 sdpb --help
```

You should see the SDPB option list.

### Step 4: Get Python (if you don't have it)

Download Python 3.9+ from https://www.python.org/downloads/.
During installation, check the box **"Add Python to PATH"**.

---

## Extract the required Python files

You need **six files** from `eft_bounds/` in this repository. Copy the entire
`eft_bounds/` folder to your working directory on your computer,
e.g. `C:\sdpb_eft\`.

Your folder structure must look exactly like this:

```
C:\sdpb_eft\
    eft_bounds\
        __init__.py
        csdr.py
        null_constraints.py
        physics.py
        pmp_generator.py
        generate_pmp.py          ← this is the script you run from the command line
```

> **Where to find these files:** In the GitHub repository, open the
> `eft_bounds/` directory and download all six `.py` files listed above.
> Click each file, then click the **Raw** button and save the page (or copy
> the text into a new file with the same name).  Keep them inside the
> `eft_bounds\` sub-folder — **do not** move `generate_pmp.py` to the parent
> folder.

> **Why `null_constraints.py` is required:** `pmp_generator.py` imports from
> it for the backward-compatibility shim used by `crossing_pipeline.py`.
> Without this file present, Python will raise an `ImportError` as soon as
> `pmp_generator` is loaded.

> **Note on the JSON output:** `pmp_generator.py` automatically strips any
> internal `_metadata` key from the file it writes to disk.  This is
> intentional — `pmp2sdp` rejects unknown top-level keys, so the metadata is
> kept in memory only and never written to the `.json` file that Docker reads.

---

## Quick test: generate one PMP file

Open **PowerShell** and navigate to your working folder:
```powershell
cd C:\sdpb_eft
```

List available operators (useful to understand the notation):
```powershell
python eft_bounds\generate_pmp.py --d 4 --K 4 --obj 1 0 --norm 1 0 --list-pairs
```

Generate an upper-bound PMP file for g̃₃ = W_{0,1}/W_{1,0}:
```powershell
python eft_bounds\generate_pmp.py --obj 1 1 --norm 1 0 --d 4 --K 4 --max-spin 10 --delta0 40 --direction upper --output upper_g3.json
```

This creates `upper_g3.json` and prints **Next steps** with the exact Docker
commands ready to copy and paste.

> The `--precision` flag in `generate_pmp.py` is **decimal digits** of output
> coefficients (default 200).  This is different from SDPB's `--precision`,
> which is in **bits** (1024 bits ≈ 308 decimal digits).  They do not need
> to match exactly.

---

## Run SDPB with Docker (step by step)

The workflow has two stages (as described in [`sdpbUsage.md`](../sdpbUsage.md)):

1. **`pmp2sdp`**: convert the JSON PMP file into SDPB's internal binary SDP format.
2. **`sdpb`**: run the semidefinite program solver on the SDP.

Docker mounts your local folder into the container at `/usr/local/share/sdpb/`
(see [`docs/Docker.md`](../docs/Docker.md)).

**Important:** always include `--platform linux/amd64` in every `docker run`
command.  Without it, Docker may try to run the wrong binary format and fail
with `exec format error`.

All commands below are **single-line** — copy each as one line.

### Step 1: Convert PMP → SDP with `pmp2sdp`

```powershell
docker run --rm --platform linux/amd64 -v "C:/sdpb_eft/:/usr/local/share/sdpb/" bootstrapcollaboration/sdpb:3.1.0 mpirun --allow-run-as-root -n 4 pmp2sdp --precision 1024 -i /usr/local/share/sdpb/upper_g3.json -o /usr/local/share/sdpb/sdp_upper_g3
```

Options explained:
- `--platform linux/amd64` — use the x86-64 image (required on ARM machines).
- `-v "C:/sdpb_eft/:/usr/local/share/sdpb/"` — mount your folder.
- `mpirun --allow-run-as-root -n 4` — run with 4 CPU cores inside Docker.
  Change `4` to match the number of cores on your machine.
- `--precision 1024` — working precision in **bits**.
- `-i …` — input JSON file path inside the container.
- `-o …` — output directory for the SDP binary files.

This creates `C:\sdpb_eft\sdp_upper_g3\` with `control.json`,
`objectives.json`, and binary block files.

### Step 2: Run SDPB

```powershell
docker run --rm --platform linux/amd64 -v "C:/sdpb_eft/:/usr/local/share/sdpb/" bootstrapcollaboration/sdpb:3.1.0 mpirun --allow-run-as-root -n 4 sdpb --precision=1024 -s /usr/local/share/sdpb/sdp_upper_g3 -o /usr/local/share/sdpb/out_upper_g3 -c /usr/local/share/sdpb/out_upper_g3/ck
```

Options explained:
- `-s …` — path to the SDP directory produced by `pmp2sdp`.
- `-o …` — output directory for the solver result.
- `-c …` — checkpoint directory (allows resuming interrupted runs).

### Step 3: Read the result

```powershell
type C:\sdpb_eft\out_upper_g3\out.txt
```

The output looks like (see [`sdpbUsage.md`](../sdpbUsage.md) for details):
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
A `terminateReason` of `"found primal-dual optimal solution"` means success.

> **Note on root-owned files:** Files written by Docker may be owned by
> root and cannot be deleted in Windows Explorer.  Delete them with:
> `docker run --rm --platform linux/amd64 -v "C:/sdpb_eft/:/usr/local/share/sdpb/" bootstrapcollaboration/sdpb:3.1.0 rm -rf /usr/local/share/sdpb/out_upper_g3`

---

## Produce a 2D allowed-region plot of (g̃₃, g̃₄) — recommended workflow

The recommended way to generate all four bounds at once is the **`run_bounds.py`
batch driver**.  It generates the four PMP JSON files in one call and has the
same `--K` convention as `generate_pmp.py` (maximum spectral power `2n+m`).

### Step 1: Generate all 4 PMP files with `run_bounds.py`

```powershell
cd C:\sdpb_eft
python -m eft_bounds.run_bounds --output-dir pmp_files --d 4 --K 8 --max-spin 50 --precision 1024
```

This creates four files in `C:\sdpb_eft\pmp_files\`:

| File | Computes |
|------|---------|
| `csdr_pmp_lower_W2_0_over_W1_0_d4_K8.json` | Lower bound on g̃₄ |
| `csdr_pmp_upper_W2_0_over_W1_0_d4_K8.json` | Upper bound on g̃₄ |
| `csdr_pmp_lower_W1_1_over_W1_0_d4_K8.json` | Lower bound on g̃₃ |
| `csdr_pmp_upper_W1_1_over_W1_0_d4_K8.json` | Upper bound on g̃₃ |

> **Parameters**: `--K 8` sets the spectral cutoff; `--max-spin 50` includes
> even spins ℓ = 0, 2, …, 50 per block; `--precision 1024` writes 1024-digit
> coefficients (matching `--precision=1024` in the Docker step).

### Step 2: Run pmp2sdp for each PMP file

```powershell
docker run --rm --platform linux/amd64 -v "C:/sdpb_eft/:/usr/local/share/sdpb/" bootstrapcollaboration/sdpb:3.1.0 mpirun --allow-run-as-root -n 4 pmp2sdp --precision 1024 -i /usr/local/share/sdpb/pmp_files/csdr_pmp_lower_W2_0_over_W1_0_d4_K8.json -o /usr/local/share/sdpb/sdp_lb_g4

docker run --rm --platform linux/amd64 -v "C:/sdpb_eft/:/usr/local/share/sdpb/" bootstrapcollaboration/sdpb:3.1.0 mpirun --allow-run-as-root -n 4 pmp2sdp --precision 1024 -i /usr/local/share/sdpb/pmp_files/csdr_pmp_upper_W2_0_over_W1_0_d4_K8.json -o /usr/local/share/sdpb/sdp_ub_g4

docker run --rm --platform linux/amd64 -v "C:/sdpb_eft/:/usr/local/share/sdpb/" bootstrapcollaboration/sdpb:3.1.0 mpirun --allow-run-as-root -n 4 pmp2sdp --precision 1024 -i /usr/local/share/sdpb/pmp_files/csdr_pmp_lower_W1_1_over_W1_0_d4_K8.json -o /usr/local/share/sdpb/sdp_lb_g3

docker run --rm --platform linux/amd64 -v "C:/sdpb_eft/:/usr/local/share/sdpb/" bootstrapcollaboration/sdpb:3.1.0 mpirun --allow-run-as-root -n 4 pmp2sdp --precision 1024 -i /usr/local/share/sdpb/pmp_files/csdr_pmp_upper_W1_1_over_W1_0_d4_K8.json -o /usr/local/share/sdpb/sdp_ub_g3
```

### Step 3: Run SDPB for each SDP

```powershell
docker run --rm --platform linux/amd64 -v "C:/sdpb_eft/:/usr/local/share/sdpb/" bootstrapcollaboration/sdpb:3.1.0 mpirun --allow-run-as-root -n 4 sdpb --precision=1024 -s /usr/local/share/sdpb/sdp_lb_g4 -o /usr/local/share/sdpb/out_lb_g4 --checkpointDir /usr/local/share/sdpb/ck_lb_g4

docker run --rm --platform linux/amd64 -v "C:/sdpb_eft/:/usr/local/share/sdpb/" bootstrapcollaboration/sdpb:3.1.0 mpirun --allow-run-as-root -n 4 sdpb --precision=1024 -s /usr/local/share/sdpb/sdp_ub_g4 -o /usr/local/share/sdpb/out_ub_g4 --checkpointDir /usr/local/share/sdpb/ck_ub_g4

docker run --rm --platform linux/amd64 -v "C:/sdpb_eft/:/usr/local/share/sdpb/" bootstrapcollaboration/sdpb:3.1.0 mpirun --allow-run-as-root -n 4 sdpb --precision=1024 -s /usr/local/share/sdpb/sdp_lb_g3 -o /usr/local/share/sdpb/out_lb_g3 --checkpointDir /usr/local/share/sdpb/ck_lb_g3

docker run --rm --platform linux/amd64 -v "C:/sdpb_eft/:/usr/local/share/sdpb/" bootstrapcollaboration/sdpb:3.1.0 mpirun --allow-run-as-root -n 4 sdpb --precision=1024 -s /usr/local/share/sdpb/sdp_ub_g3 -o /usr/local/share/sdpb/out_ub_g3 --checkpointDir /usr/local/share/sdpb/ck_ub_g3
```

> **Do not** pass `-c` (checkpoint load) for a fresh run.  Only `--checkpointDir`
> is needed — it tells SDPB where to save checkpoints if the run is interrupted.

### Step 4: Collect results

```powershell
type C:\sdpb_eft\out_lb_g4\out.txt
type C:\sdpb_eft\out_ub_g4\out.txt
type C:\sdpb_eft\out_lb_g3\out.txt
type C:\sdpb_eft\out_ub_g3\out.txt
```

#### Which runs converge and which diverge — and why

**Only two of the four SDPB runs will converge.  The other two are expected to
diverge.**  This is not a bug.

| Run | `terminateReason` | What it means |
|-----|-------------------|---------------|
| `out_lb_g4` | `"found primal-dual optimal solution"` ✓ | **g̃₄ has a finite lower bound** |
| `out_ub_g4` | `"maxComplementarity exceeded"` | g̃₄ has **no** finite upper bound from CSDR alone |
| `out_ub_g3` | `"found primal-dual optimal solution"` ✓ | **g̃₃ has a finite upper bound** |
| `out_lb_g3` | `"maxComplementarity exceeded"` | g̃₃ has **no** finite lower bound from CSDR alone |

**Physically**, crossing-symmetric dispersion relations (CSDR) impose:
- A *lower* bound on g̃₄ (the forward-scattering integral is non-negative).
- An *upper* bound on g̃₃ (a positivity condition at finite spin).

They do not, by themselves, give an upper bound on g̃₄ or a lower bound on g̃₃.
The `"maxComplementarity exceeded"` message means the SDP found that the
feasible region is unbounded in that direction — the solver's interior-point
parameter µ diverged, exactly as expected for an unbounded problem.

> **`"maxComplementarity exceeded"` ≠ error.**  It means the bound is ±∞
> in that direction.  Do not try to read a physical value from `primalObjective`
> in those runs — the number printed is meaningless (the SDP never converged).

#### Reading the physical bounds from converged runs

**Important:** SDPB *always* maximizes `b·y`.  The sign of `b` depends on which
direction you are bounding:

| `bound_direction` | `b[obj]` | SDPB maximizes | `primalObjective` | Physical bound |
|-------------------|----------|----------------|-------------------|----------------|
| `"upper"` | **+1** | `+z_obj` | max of `z_obj` | bound = **+primalObjective** |
| `"lower"` | **−1** | `−z_obj` | max of `−z_obj` | bound = **−primalObjective** |

Because the lower-bound run negates the objective (`b = −1`), the solver
maximizes `−g̃₄`.  You must **negate `primalObjective`** to recover the actual
lower bound:

    lower bound on g̃₄  =  −primalObjective(lb_g4 run)
    upper bound on g̃₃  =  +primalObjective(ub_g3 run)

**Verified results from `sdpb_data2/` (K=8, d=4, max-spin=50):**

| Run | `primalObjective` | How to convert | Physical bound |
|-----|-------------------|----------------|----------------|
| `out_lb_g4/out.txt` | **0.5000** | −primalObjective | **g̃₄ ≥ −0.5000** |
| `out_ub_g3/out.txt` | **0.3017** | +primalObjective | **g̃₃ ≤ +0.3017** |

> **Why g̃₄ ≥ −1/2 is exact and analytical:**  The ℓ=0 spin block (the only
> block where the g̃₃ kernel is zero and g̃₄ has a non-zero kernel) reads
>
>     P⁰(0) = 2·z_{1,0} + 4·z_{2,0} ≥ 0
>
> With normalization `z_{1,0} = 1` this gives `z_{2,0} ≥ −1/2` exactly.
> SDPB correctly finds that the ℓ=0 block is the binding constraint and
> returns `primalObjective = max(−z_{2,0}) = 1/2` (the ℓ=0 block saturated
> at x=0).  Increasing K or max-spin does **not** change this bound because
> no higher-spin block adds a tighter constraint on the ratio W_{2,0}/W_{1,0}.

> `terminateReason = "found primal-dual optimal solution"` means the solver
> converged.  The `dualityGap` should be ≲ 10⁻²⁰ for a well-converged result.

---

## Running at higher cutoff K=12 for tighter bounds

Larger K includes more spectral operators and can produce tighter bounds.
The commands below use K=12; the workflow is identical.

### Generate PMP files at K=12

```powershell
python -m eft_bounds.run_bounds --output-dir pmp_files_K12 --d 4 --K 12 --max-spin 50 --precision 1024
```

### Run pmp2sdp for K=12

```powershell
docker run --rm --platform linux/amd64 -v "C:/sdpb_eft/:/usr/local/share/sdpb/" bootstrapcollaboration/sdpb:3.1.0 mpirun --allow-run-as-root -n 4 pmp2sdp --precision 1024 -i /usr/local/share/sdpb/pmp_files_K12/csdr_pmp_lower_W2_0_over_W1_0_d4_K12.json -o /usr/local/share/sdpb/sdp_K12_lb_g4

docker run --rm --platform linux/amd64 -v "C:/sdpb_eft/:/usr/local/share/sdpb/" bootstrapcollaboration/sdpb:3.1.0 mpirun --allow-run-as-root -n 4 pmp2sdp --precision 1024 -i /usr/local/share/sdpb/pmp_files_K12/csdr_pmp_upper_W2_0_over_W1_0_d4_K12.json -o /usr/local/share/sdpb/sdp_K12_ub_g4

docker run --rm --platform linux/amd64 -v "C:/sdpb_eft/:/usr/local/share/sdpb/" bootstrapcollaboration/sdpb:3.1.0 mpirun --allow-run-as-root -n 4 pmp2sdp --precision 1024 -i /usr/local/share/sdpb/pmp_files_K12/csdr_pmp_lower_W1_1_over_W1_0_d4_K12.json -o /usr/local/share/sdpb/sdp_K12_lb_g3

docker run --rm --platform linux/amd64 -v "C:/sdpb_eft/:/usr/local/share/sdpb/" bootstrapcollaboration/sdpb:3.1.0 mpirun --allow-run-as-root -n 4 pmp2sdp --precision 1024 -i /usr/local/share/sdpb/pmp_files_K12/csdr_pmp_upper_W1_1_over_W1_0_d4_K12.json -o /usr/local/share/sdpb/sdp_K12_ub_g3
```

### Run SDPB for K=12

```powershell
docker run --rm --platform linux/amd64 -v "C:/sdpb_eft/:/usr/local/share/sdpb/" bootstrapcollaboration/sdpb:3.1.0 mpirun --allow-run-as-root -n 4 sdpb --precision=1024 -s /usr/local/share/sdpb/sdp_K12_lb_g4 -o /usr/local/share/sdpb/out_K12_lb_g4 --checkpointDir /usr/local/share/sdpb/ck_K12_lb_g4

docker run --rm --platform linux/amd64 -v "C:/sdpb_eft/:/usr/local/share/sdpb/" bootstrapcollaboration/sdpb:3.1.0 mpirun --allow-run-as-root -n 4 sdpb --precision=1024 -s /usr/local/share/sdpb/sdp_K12_ub_g4 -o /usr/local/share/sdpb/out_K12_ub_g4 --checkpointDir /usr/local/share/sdpb/ck_K12_ub_g4

docker run --rm --platform linux/amd64 -v "C:/sdpb_eft/:/usr/local/share/sdpb/" bootstrapcollaboration/sdpb:3.1.0 mpirun --allow-run-as-root -n 4 sdpb --precision=1024 -s /usr/local/share/sdpb/sdp_K12_lb_g3 -o /usr/local/share/sdpb/out_K12_lb_g3 --checkpointDir /usr/local/share/sdpb/ck_K12_lb_g3

docker run --rm --platform linux/amd64 -v "C:/sdpb_eft/:/usr/local/share/sdpb/" bootstrapcollaboration/sdpb:3.1.0 mpirun --allow-run-as-root -n 4 sdpb --precision=1024 -s /usr/local/share/sdpb/sdp_K12_ub_g3 -o /usr/local/share/sdpb/out_K12_ub_g3 --checkpointDir /usr/local/share/sdpb/ck_K12_ub_g3
```

**Verified results from `sdpb_data2_K12/` (K=12, d=4, max-spin=50):**

| Run | `terminateReason` | `primalObjective` | Physical bound |
|-----|-------------------|-------------------|----------------|
| `out_K12_lb_g4` | optimal ✓ | **0.5000** | g̃₄ ≥ **−primalObjective** = **−0.5000** |
| `out_K12_ub_g4` | maxComplementarity exceeded | (meaningless) | no finite upper bound (expected) |
| `out_K12_ub_g3` | optimal ✓ | **0.3018** | g̃₃ ≤ **+primalObjective** = **+0.3018** |
| `out_K12_lb_g3` | maxComplementarity exceeded | (meaningless) | no finite lower bound (expected) |

The K=12 bounds are the same as K=8 for these operators at max-spin=50,
confirming that the ℓ=0 block gives the tight lower bound on g̃₄ and increasing
K does not change it.

---

## Plotting the allowed region

### What the CSDR bounds give

The CSDR method provides **two half-plane constraints** in the (g̃₃, g̃₄) plane.
To read the physical bounds from `primalObjective`:

- **lower-bound run** (`lb_g4`): physical lower bound = **−primalObjective**
- **upper-bound run** (`ub_g3`): physical upper bound = **+primalObjective**

From the verified K=8 results:

    g̃₄  ≥  −0.5000    (= −primalObjective from lb_g4 run; K=8 and K=12 give the same)
    g̃₃  ≤  +0.3017    (= +primalObjective from ub_g3 run; K=8)

The allowed region is the **intersection of these two half-planes**, which is
an infinite wedge.

### Relation to Figure 8 of "Extremal EFT" (Caron-Huot & Duong 2021)

Figure 8 of the paper shows a **bounded** polygon.  That figure uses
**fixed-t dispersion relations** (Sections 3.3–3.4 of the paper), which give
additional constraints (e.g. a finite upper bound on g̃₄ from a different
subtraction scheme) that are *not* present in CSDR.

The CSDR approach implemented here correctly recovers the two constraints
above.  To reproduce the full bounded polygon of figure 8, additional
fixed-t constraints would need to be added as separate SDP blocks —
this is a different calculation that is not yet implemented in this codebase.

### Step 5: Plot in Python

Install matplotlib if needed:

```powershell
pip install matplotlib
```

Create `plot_region.py` in `C:\sdpb_eft\`.

#### Single-cutoff plot (K=8)

```python
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

# ---------------------------------------------------------------
# Fill in primalObjective from the two CONVERGED runs only.
# (The two "maxComplementarity exceeded" runs give no finite bound.)
# ---------------------------------------------------------------
lb_g4 = 0.5000   # primalObjective from out_lb_g4/out.txt
ub_g3 = 0.3017   # primalObjective from out_ub_g3/out.txt

# Plot window: extend well beyond the bounds to show the half-planes
g3_min, g3_max = -0.5, 0.7
g4_min, g4_max =  0.0, 2.0

fig, ax = plt.subplots(figsize=(7, 5))

# Shade the ALLOWED region (g3 <= ub_g3  AND  g4 >= lb_g4)
g3_vals = np.linspace(g3_min, ub_g3, 500)
ax.fill_betweenx(
    [lb_g4, g4_max],
    g3_min, ub_g3,
    color='lightblue', alpha=0.6, label='CSDR allowed region (K=8, d=4)'
)

# Draw bound lines
ax.axvline(ub_g3, color='blue',   lw=2, linestyle='--', label=rf'$\tilde{{g}}_3 \leq {ub_g3:.4f}$')
ax.axhline(lb_g4, color='navy',   lw=2, linestyle='-',  label=rf'$\tilde{{g}}_4 \geq {lb_g4:.4f}$')

ax.set_xlim(g3_min, g3_max)
ax.set_ylim(g4_min, g4_max)
ax.set_xlabel(r'$\tilde{g}_3 = W_{0,1}/W_{1,0}$', fontsize=13)
ax.set_ylabel(r'$\tilde{g}_4 = W_{2,0}/W_{1,0}$', fontsize=13)
ax.set_title('CSDR allowed region (K=8, d=4, max-spin=50)', fontsize=12)
ax.legend(fontsize=11)
plt.tight_layout()
plt.savefig('allowed_region_K8.png', dpi=150)
print('Saved: allowed_region_K8.png')
plt.show()
```

#### K=8 and K=12 overlay (to see convergence)

```python
import matplotlib.pyplot as plt
import numpy as np

# K=8 bounds (from sdpb_data2/)
lb_g4_K8  = 0.5000
ub_g3_K8  = 0.3017

# K=12 bounds (from sdpb_data2_K12/)
lb_g4_K12 = 0.5000
ub_g3_K12 = 0.3018

g3_min, g3_max = -0.5, 0.7
g4_min, g4_max =  0.0, 2.0

fig, ax = plt.subplots(figsize=(7, 5))

# K=8 allowed region
ax.fill_betweenx(
    [lb_g4_K8, g4_max], g3_min, ub_g3_K8,
    color='lightblue', alpha=0.5, label='K=8 allowed'
)
# K=12 allowed region (slightly different shading to show overlap)
ax.fill_betweenx(
    [lb_g4_K12, g4_max], g3_min, ub_g3_K12,
    color='lightyellow', alpha=0.5, label='K=12 allowed'
)

# Bound lines
ax.axvline(ub_g3_K8,  color='blue',   lw=2, linestyle='--',
           label=rf'$\tilde{{g}}_3 \leq {ub_g3_K8:.4f}$ (K=8)')
ax.axvline(ub_g3_K12, color='green',  lw=2, linestyle=':',
           label=rf'$\tilde{{g}}_3 \leq {ub_g3_K12:.4f}$ (K=12)')
ax.axhline(lb_g4_K8,  color='navy',   lw=2, linestyle='-',
           label=rf'$\tilde{{g}}_4 \geq {lb_g4_K8:.4f}$ (K=8, K=12)')

ax.set_xlim(g3_min, g3_max)
ax.set_ylim(g4_min, g4_max)
ax.set_xlabel(r'$\tilde{g}_3 = W_{0,1}/W_{1,0}$', fontsize=13)
ax.set_ylabel(r'$\tilde{g}_4 = W_{2,0}/W_{1,0}$', fontsize=13)
ax.set_title('CSDR allowed region, K=8 vs K=12 (d=4, max-spin=50)', fontsize=11)
ax.legend(fontsize=10)
plt.tight_layout()
plt.savefig('allowed_region_K8_K12.png', dpi=150)
print('Saved: allowed_region_K8_K12.png')
plt.show()
```

Run either script:

```powershell
python plot_region.py
```

The plot shows the **allowed half-plane wedge** in (g̃₃, g̃₄) space.
Points to the left of the vertical dashed line and above the horizontal line
are allowed by CSDR unitarity constraints.

---

## Understanding the parameters

### Parameters for `run_bounds.py` (batch driver — recommended)

| Flag | Meaning | Recommended value |
|------|---------|-------------------|
| `--output-dir` | Output directory for PMP JSON files | e.g. `pmp_files` |
| `--d` | Spacetime dimension | 4 (for 4D QFT) |
| `--K` | Max spectral power 2n+m | 8 (default); 12 for tighter bounds |
| `--max-spin` | Max even spin ℓ | **50** (use at least 20) |
| `--precision` | Decimal digits in PMP coefficients | 1024 (matches SDPB `--precision=1024` bits ≈ 308 digits; use 1024 for safety) |

### Parameters for `generate_pmp.py` (single-bound driver)

| Flag | Meaning | Notes |
|------|---------|-------|
| `--K` | Max spectral power 2n+m | Same convention as `run_bounds.py` |
| `--d` | Spacetime dimension | Physical: d=4 for 4D |
| `--delta0` | IR cutoff δ₀ | s₁ starts at δ₀; matches the paper |
| `--max-spin` | Max even spin ℓ | **50** recommended; more spins → tighter bounds, slower |
| `--precision` | Decimal digits in output | 200 is safe; use 50 for quick tests |
| `--obj n m` | Wilson coeff to bound | (n, m) pair with n ≥ m ≥ 0 |
| `--norm n m` | Normalization coeff | Set this W to 1 |
| `--direction` | `upper` or `lower` | Which extremum to compute |

### Parameters for `pmp2sdp` and `sdpb`

These take `--precision` in **bits**:

| `--precision` bits | Decimal digits | Recommendation |
|--------------------|----------------|----------------|
| 512  | ~154 | Quick tests |
| 1024 | ~308 | Safe default (matches 200-digit PMP input) |
| 2048 | ~616 | High-precision production runs |

Run `docker run --rm --platform linux/amd64 bootstrapcollaboration/sdpb:3.1.0 sdpb --help` to see all available `sdpb` options.

---

## Relation to Extremal EFT (Caron-Huot & Duong 2021) and Figure 8

Figure 8 of "Extremal EFT" shows a **bounded** polygon in the (g̃₃, g̃₄) plane.
That figure uses **fixed-t dispersion relations** (Sections 3.3–3.4 of the
paper), which provide *more* constraints than CSDR alone:

| Method | Constraints obtained | Region shape |
|--------|---------------------|-------------|
| **CSDR** (this code) | g̃₄ ≥ 0.5, g̃₃ ≤ 0.3017 | Open wedge (2 half-planes) |
| **Fixed-t** (paper) | upper+lower bounds on both g̃₃ and g̃₄ | Bounded polygon |

The CSDR approach in this code implements the Sinha-Zahed (2021) formula:

    W_{n−m,m} = ⟨ D^{(n,m)}_{ℓ,α} · C^α_ℓ(1) · (2ℓ+d−3) / s₁^{2n+m} ⟩

where D is spin-dependent and α = (d−3)/2.

The CSDR method gives physically correct bounds for g̃₄ ≥ 0 (lower bound) and
g̃₃ ≤ 0.3017 (upper bound) at K=8, d=4, max-spin=50.  These are not the same
as the bounds in the paper because the two formalisms use different subtraction
schemes and measure different combinations of Wilson coefficients.

---

## Troubleshooting

**`exec /usr/bin/mpirun: exec format error`**
→ Docker is trying to run a Linux binary with the wrong container mode.
  Work through these checks in order:

  **Check 1 — Docker is in Windows containers mode (most common cause):**
  Right-click the Docker Desktop tray icon.  If it offers
  **"Switch to Linux containers…"**, click it and wait for Docker to restart,
  then retry.  The SDPB image is a Linux container and will not run while
  Docker is in Windows containers mode.

  **Check 2 — wrong image cached:** even in Linux mode, an old cached image
  for the wrong architecture can cause this.  Re-pull explicitly:
  ```powershell
  docker pull --platform linux/amd64 bootstrapcollaboration/sdpb:3.1.0
  ```
  Then retry your `docker run --platform linux/amd64 …` command.

**`terminateReason = "maxComplementarity exceeded"` (mu diverges)**
→ This is **expected** for the `ub_g4` and `lb_g3` runs.  It means the SDP
  is unbounded in that direction — no finite upper bound on g̃₄ or finite lower
  bound on g̃₃ exists from CSDR alone.  The primalObjective printed is not
  meaningful for these runs; ignore it.
→ If this happens for `lb_g4` or `ub_g3` (which should converge), it indicates
  a bug in the PMP generation.  Re-generate the PMP files making sure you are
  using the latest `pmp_generator.py` (the version in this repo).



**Docker command fails with "path not found"**
→ Use forward slashes in the `-v` flag on Windows:
  `-v "C:/sdpb_eft/:/usr/local/share/sdpb/"` (not backslashes).

**`'A was not numerically HPD'` error**
→ Increase `--precision` in `pmp2sdp` and `sdpb` (e.g. try `2048` instead
  of `1024`), as described in [`sdpbUsage.md`](../sdpbUsage.md).

**Out-of-memory error (`std::bad_alloc`)**
→ Reduce the number of MPI cores (use `-n 2` instead of `-n 4`) or add
  `--maxSharedMemory=2G` to the `sdpb` command, as described in
  [`sdpbUsage.md`](../sdpbUsage.md).

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

**To use them directly** (copy `examples/` to `C:\sdpb_eft\examples\`):

```powershell
docker run --rm --platform linux/amd64 -v "C:/sdpb_eft/:/usr/local/share/sdpb/" bootstrapcollaboration/sdpb:3.1.0 mpirun --allow-run-as-root -n 4 pmp2sdp --precision 1024 -i /usr/local/share/sdpb/examples/upper_W01_over_W10_K4_d4_delta40.json -o /usr/local/share/sdpb/sdp_ub_g3

docker run --rm --platform linux/amd64 -v "C:/sdpb_eft/:/usr/local/share/sdpb/" bootstrapcollaboration/sdpb:3.1.0 mpirun --allow-run-as-root -n 4 sdpb --precision=1024 -s /usr/local/share/sdpb/sdp_ub_g3 -o /usr/local/share/sdpb/out_ub_g3 -c /usr/local/share/sdpb/out_ub_g3/ck

type C:\sdpb_eft\out_ub_g3\out.txt
```

The `primalObjective` line gives the upper bound on g̃₃.
