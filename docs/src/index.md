```@meta
CurrentModule = qPCRModels
```

# qPCRModels.jl

`qPCRModels.jl` is a Julia package that fits qPCR fluorescence amplification curves with non-linear least squares. It implements:

| Model | Type | Fitted range | Non-linear parameters |
|---|---|---|---|
| `BirchModel` | mechanistic | **whole curve** | `b`, `tₘ` |
| `MAK2Model` | mechanistic | exponential phase only (up to the SDM) | `x₀` |
| `LogisticModel` | empirical (4PL) | whole curve | `a`, `tₘ` |
| `RichardsModel` | empirical (5PL) | whole curve | `a`, `tₘ`, `ν` |

Every model has the form `F(c) = F_bg + F_max · x(c)`. The baseline `F_bg` and scale `F_max` are linear parameters, eliminated by Variable Projection, so the Levenberg–Marquardt optimiser only searches over the non-linear parameters.

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
```

## Command-line interface

```bash
# Fit the default Birch model to a built-in dataset
julia --project=. bin/qpcr_fit.jl batsch1

# Fit a custom dataset with MAK2 and export a CSV
julia --project=. bin/qpcr_fit.jl path/to/data.csv --model mak2 -o results.csv
```

Options:

- `-m, --model <name>`: `birch` (default), `logistic`, `richards` or `mak2`.
- `-o, --output <path>`: output CSV file.
- `--fit-noise`: also fit curves that are flagged as non-amplifying.
- `-h, --help`: show help.

