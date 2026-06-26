# Sensitivity Training (SuperNEMO)

This repository contains training material to compute background expectations and Bayesian sensitivity estimates in the SuperNEMO context.

## Repository contents

- `data/`: input ROOT files and `data_info.csv` index.
- `examples/Instantiate.ipynb`: environment bootstrap notebook (activate, instantiate, precompile).
- `examples/ex1_background.ipynb`: expected background workflow and basic outputs.
- `examples/ex2_bayes.ipynb`: Bayesian fit and sensitivity workflow.
- `examples/ex2_background.jl`, `examples/ex2_bayes.jl`: script versions of the notebook examples.
- `examples/ex3_freq1D.jl`, `examples/ex3_freq1D.jl`: Frequentist 1D example
- `examples/ex4_freqND.jl`, `examples/ex4_freqND.jl`: Frequentist ND example


# First steps:

## 1. Clone repository to your sps scratch

To clone the repository use

```bash
cd /sps/nemo/scratch/MY_NAME/PATH_TO_WHERE_YOU_WANT_ME
git clone https://github.com/Shoram444/SN_Sensitivity_training.git
```

## 2. Copy example data files to your directory

Before you start: copy large assets

Large files are not tracked in git (`data/` ROOT files).
Before running the training, copy them from:

`/sps/nemo/scratch/mpetro/Sensitivity_training`

From your local clone root, for example:

```bash
# copy ROOT inputs and metadata index
cp -r /sps/nemo/scratch/mpetro/Sensitivity_training/data ./

```

If you prefer, use `rsync` instead of `cp`.


## 3. Install julia on CC or your own computer

### 3.1 Julia setup on CC / SLURM

Load Julia from modules:

```bash
module load julia
julia --version
```

On CC clusters with heterogeneous CPUs, set a generic CPU target for precompilation by adding this line to `~/.profile`:

```bash
export JULIA_CPU_TARGET="generic"
```

Then reload your profile in each new shell (or log out/in):

```bash
source ~/.profile
```

### 3.2 Julia on own PC
Go to https://github.com/JuliaLang/juliaup

**windows** 
```bash
winget install --name Julia --id 9NJNWW8PVKMN -e -s msstore
```

**linux/mac**
```bash
curl -fsSL https://install.julialang.org | sh
```

**homebrew**
```bash
brew install juliaup
```


## 4. First-time environment initialization

To install all the dependencies and start the project (this may take some 5-10 minutes)

```bash
cd PATH_TO_THIS_DIRECTORY
julia --project=. examples/Instantiate.jl
```

Notes:

- `examples/Instantiate.ipynb` is intended to be run first.
- `Manifest.toml` is pinned to Julia `1.10.3`.

## 5. IJulia setup (for notebook.cc.in2p3.fr)

For learning I find it best to use the jupyter notebooks. Although everything can be done via simple `.jl` scripts in equivalent way.

Install kernel for this repository use:

Open julia REPL by running

```bash
julia
```

Then inside REPL:

```julia
import Pkg
Pkg.activate(".")
Pkg.add("IJulia")
using IJulia
IJulia.installkernel(
	"Julia 1.10.3 (Sensitivity_training)",
	"--project=$(abspath("."))";
	env = Dict("JULIA_CPU_TARGET" => "generic")
)
```

After installing, restart your Jupyter session on `https://notebook.cc.in2p3.fr/` and select the kernel named `Julia 1.10.3 (Sensitivity_training)`.

**Now we are ready to run the examples**

-------

## Running modes (choose one)

Users do not have to use `https://notebook.cc.in2p3.fr/`.

### A) CC notebook server (Jupyter)

- Use the IJulia kernel steps above.
- Open notebooks in `examples/` and run in order.

### B) VS Code (local or CC remote session)

1. Install the VS Code extensions `Julia` and `Jupyter`.
2. Open this repository folder in VS Code.
3. In a terminal inside the repo:

```bash
module load julia
source ~/.profile
julia --project=. -e 'import Pkg; Pkg.instantiate(); Pkg.precompile()'
```

4. Open an `.ipynb` notebook and select the Julia kernel you installed.
5. If using scripts (`.jl`), run them in the integrated terminal with:

```bash
julia --project=. examples/ex1_background.jl
```

6. In notebook UI, select the Julia kernel from the kernel picker.

### C) Script-only mode (`.jl` files)

You can skip notebooks and run:

```bash
julia --project=. examples/ex1_background.jl
julia --project=. examples/ex2_bayes.jl
```

----
## Notebook guide

### 1) `examples/Instantiate.ipynb` or `examples/Instantiate.jl` 

Purpose:

- Activate this project environment.
- Instantiate all dependencies from `Project.toml` / `Manifest.toml`.
- Precompile packages.

Recommended use:

- Run this notebook before running any other example notebook.

### 2) `examples/ex1_background.ipynb` or `examples/ex1_background.jl`

Purpose:

- Load process metadata from `data/data_info.csv`.
- Open ROOT files with `UnROOT` and convert to `LazyTree` tables.
- Build `SNSensitivityEstimate.DataProcess` objects.
- Compute expected background histograms via `get_bkg_counts_1D`.
- Produce stacked background plots and counts tables.

Key outputs produced by this workflow:

- `data/out/ex1/stacked_background_histogram.png`
- `data/out/ex1/counts_per_process.csv`
- `data/out/ex1/counts_per_process_roi.csv`

### 3) `examples/ex2_bayes.ipynb`

Purpose:

- Start from the same data-loading stage as `ex1`.
- Apply analysis cuts (`sumE`, `dy`, `dz`, `Pint`, `Pext`).
- Build signal/background templates and pseudo-data.
- Define prior + likelihood + posterior with `BAT.jl`.
- Run MCMC (`BAT.bat_sample`) and inspect posterior distributions.
- Convert posterior signal limit to half-life sensitivity.
- Optionally repeat pseudo-experiments to get a median sensitivity.

Key outputs produced by this workflow:

- `data/out/ex2/posterior_distributions.png`
- `data/out/ex2/fit_result.png`
- `data/out/ex2/sensitivity_distribution.png`


Runtime note:

- MCMC parts are compute-heavy and can take significant time.

### 3) `examples/ex3_freq1D.ipynb`

Purpose:

- Start from the same data-loading stage as `ex1`.
- Apply analysis cuts (`sumE`, `dy`, `dz`, `Pint`, `Pext`).
- Calculate sensitivity with `my_roi` (your own proposal for ROI)
- Build map of the sensitivity per ROI in the selected single channel
- Calculate the best obtainable sensitivity

Key outputs produced by this workflow:

- `data/out/ex3/my_roi_plot.png`
- `data/out/ex3/frequentist_sensitivity_map.png`
- `data/out/ex3/best_roi_plot.png`


### 4) `examples/ex4_freqND.ipynb`

Purpose:

- Start from the same data-loading stage as `ex1`.
- 

## Recommended execution order

1. `examples/Instantiate.ipynb`
2. `examples/ex1_background.ipynb`
3. `examples/ex2_bayes.ipynb`

## Path note

Some notebook cells currently use absolute paths under `/sps/nemo/scratch/mpetro/Sensitivity_training/`.
For portability, prefer your own relative paths, for example:

- use `Pkg.activate(".")` instead of an absolute project path;
- use `CSV.read("data/data_info.csv", DataFrame)`;
- use `joinpath("data", file_name)` for ROOT files;
- write outputs to `data/out/...` relative to repository root.

## Julia quick intro (for first-time users)

### Minimal Julia basics

- Comments start with `#`.
- Variables are created by assignment: `x = 3`.
- Arrays use square brackets: `a = [1, 2, 3]`.
- Functions are defined with `function ... end` (or short form `f(x) = x^2`).
- Indexing starts at `1` (not `0`).

Example:

```julia
x = [10, 20, 30]
println(x[1])    # 10
println(length(x))
```

### Package/environment basics

- This repo is a Julia project (`Project.toml` + `Manifest.toml`).
- Always use `--project=.` or `Pkg.activate(".")` in this folder.
- Install/resolve packages once with `Pkg.instantiate()`.

### Workflow summary for this training

1. Load Julia module and environment (`module load julia`, `source ~/.profile`).
2. Activate project and instantiate.
3. Install the IJulia kernel first.
4. Run `examples/Instantiate.ipynb` once.
5. Continue with background example (`ex1`) then Bayesian example (`ex2`).

