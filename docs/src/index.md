```@meta
CurrentModule = qPCRModels
```

# qPCRModels.jl

`qPCRModels.jl` is a Julia package that fits qPCR fluorescence amplification curves with non-linear least squares. It implements:

| Model | Type | Fitted range | Non-linear parameters |
|---|---|---|---|
| `BirchModel` | mechanistic | **whole curve** | ``b``, ``t_m`` |
| `MAK2Model` | mechanistic | exponential phase only (up to the SDM) | ``x_0`` |
| `LogisticModel` | empirical (4PL) | whole curve | ``a``, ``t_m`` |
| `RichardsModel` | empirical (5PL) | whole curve | ``a``, ``t_m``, ``\nu`` |

Every model has the form ``F(c) = F_\mathrm{bg} + F_\mathrm{max} x(c)``, where ``c`` is the cycle number. The baseline ``F_\mathrm{bg}`` and scale ``F_\mathrm{max}`` are linear parameters, eliminated by Variable Projection, so the Levenberg–Marquardt optimiser only searches over the non-linear parameters.

See [The Birch model](@ref) for the mechanistic derivation of the default model.

## Installation

```julia
using Pkg
Pkg.add(url = "https://github.com/gschivre/qPCRModels.jl")
```

For development:

```julia
Pkg.develop(path = "path/to/qPCRModels.jl")
```

## Quick start

```julia
using qPCRModels

# Load a built-in dataset
ds = qpcRDataset("batsch1")

# Fit the Birch model to the first curve
m = BirchModel(ds[:, 2])
fit!(m)

hasconverged(m)   # true
get_params(m)     # (b, tₘ)
get_x0(m)         # initial quantity n₀/K
get_cy0(m)        # Cy₀ cycle
AICc(m)           # model selection criterion
get_stderrors(m)  # standard errors of b, tₘ, F_bg, F_max and ln x₀
```

A fitted model prints a summary in the REPL, with standard errors:

```
BirchModel · 45 cycles · converged
  b          1.133 ± 0.061
  tₘ         25.389 ± 0.046
  F_bg       5477 ± 12
  F_max      4946 ± 22
  log₁₀ x₀  -7.302 ± 0.021
  Cy₀        21.01
  RSE = 51.34 · AICc = -206.2
```

`RSE` is the residual standard error, in fluorescence units.

### Absolute quantification

``x_0`` is a relative quantity (``n_0/K``). To get absolute quantities (e.g. copies), fit a standard curve ``\log_{10} x_0 = \alpha + \beta \log_{10} q`` on wells of known quantity ``q``, then apply it to new samples:

```julia
standards = [BirchModel(y) for y in standard_curves] # wells with 10⁵, 10⁴, 10³ and 10² copies
foreach(fit!, standards)
cal = calibration([1e5, 1e4, 1e3, 1e2], standards)
quantify(cal, sample) # copies in a fitted sample well
```

The line is fitted robustly (Tukey bisquare weights), so outlier wells are down-weighted. Models that did not converge are left out. With a mechanistic model, ``x_0`` is proportional to the template, so the slope should be close to 1. `fixed_slope = true` fixes it to 1 and fits only the intercept.

## Command-line interface

```bash
# Fit the default Birch model to a built-in dataset
julia --project=. bin/qpcr_fit.jl batsch1

# Fit a custom dataset with MAK2 and export a CSV
julia --project=. bin/qpcr_fit.jl path/to/data.csv --model mak2 -o results.csv

# Absolute quantification with standard wells of the same run
julia --project=. bin/qpcr_fit.jl path/to/plate.csv --standards standards.csv -o results.csv
```

Options:

- `-m, --model <name>`: `birch` (default), `logistic`, `richards` or `mak2`.
- `-o, --output <path>`: output CSV file.
- `-s, --standards <path>`: standards file, for absolute quantification (see below).
- `--fixed-slope`: fix the slope of the standard curve to 1.
- `--no-robust`: fit the standard curve by ordinary least squares, with no outlier rejection.
- `--fit-noise`: also fit curves that are flagged as non-amplifying.
- `-h, --help`: show help.


The input file has one column per curve and one row per cycle, with a header line. The first column (the cycle number) is skipped unless it is the only column. Columns can be separated by tabs, commas, semicolons or spaces.

The standards file lists the standard curves of the input data and their known quantities, one per line. The header line is optional:

```text
curve,copies
F1.1,100000
F1.2,100000
F2.1,10000
F2.2,10000
F3.1,1000
```

With `--standards`, the tool fits the standard curve on these wells (as [`calibration`](@ref) does: robust fit unless `--no-robust` is given, non-converged wells left out). It prints the slope, the intercept, R² and the rejected outliers. Then it reports the quantity of every curve, standards included, in the units of the file. The CSV output gets two more columns: `Known_quantity` (standards only) and `Quantity`. You need at least 3 converged standards with at least two different quantities, or 2 with `--fixed-slope`.

### Standalone executable

A compiled version of the command-line tool runs without Julia. Download the archive for your system (Linux x64 or macOS ARM) from the [Releases page](https://github.com/gschivre/qPCRModels.jl/releases), extract it, and run `bin/qpcr_fit`:

```bash
tar xzf qpcr_fit-v0.2.0-linux-x64.tar.gz
cd qpcr_fit-v0.2.0-linux-x64
bin/qpcr_fit batsch1
bin/qpcr_fit path/to/data.csv --model mak2 -o results.csv
```

It takes the same arguments and options as `bin/qpcr_fit.jl`. Keep the extracted folder together: the executable needs the `lib/` folder next to it, and reads the built-in datasets from `share/qpcR_datasets/`. To call it from anywhere, add the `bin/` folder to your `PATH` instead of moving the executable.

The executable starts in a few milliseconds: fitting all 375 curves of `dil4reps94` takes about 0.2 s, compared with about 8 s for `julia bin/qpcr_fit.jl`, which must compile the code first.

#### Building it yourself

The executable is built with [JuliaC](https://github.com/JuliaLang/JuliaC.jl), which needs Julia 1.12 or later. JuliaC is not a dependency of the package: it has its own environment in `juliac/`.

```bash
julia --project=juliac -e 'using Pkg; Pkg.instantiate()'   # first time only
julia --project=juliac juliac/build.jl
```

The bundle is written to `juliac/qpcr_fit/` (about 115 MB, 39 MB compressed). It is built with `--trim=safe`, which keeps only the code that the program can reach. This only works if every call has a concrete type. That is why `bin/qpcr_fit.jl` avoids `Dict{Symbol,Any}`, chooses the model with a branch instead of a `Dict` of types, and writes to `Core.stdout`. The `Binary` workflow builds the executable on every push, so a change that breaks trimming fails in CI.
