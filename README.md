# qPCRModels.jl

`qPCRModels.jl` is a Julia package for non-linear least squares model fitting of qPCR fluorescence amplification curves. It implements the Birch model, 4-parameter Logistic, 5-parameter Richards, and MAK2 models.

## Package Architecture

```
qPCRModels/
├── Project.toml              # Package definition & dependencies
├── src/
│   ├── qPCRModels.jl         # Primary module entrypoint (exports API)
│   ├── models.jl             # NonLinearLeastSquare structs & methods
│   ├── solver.jl             # Levenberg-Marquardt & rootfinding solvers
│   ├── utils.jl              # Thresholding & amplification metrics
│   └── datasets.jl           # Built-in qPCR dataset management
├── bin/
│   └── qpcr_fit.jl           # Command-line interface executable
├── test/
│   └── runtests.jl           # Package test suite
└── qpcR_datasets/            # Built-in qPCR datasets
```

## Quick Start (As a Package)

```julia
using Pkg
Pkg.activate(".")
using qPCRModels

# Load a built-in dataset
ds = qpcRDataset("batsch1")

# Fit Birch model to the first curve
m = BirchModel(ds[:, 2])
fit!(m)

# Inspect fit results
println("Converged: ", hasconverged(m))
println("Initial Concentration (x0): ", get_x0(m))
println("AICc: ", AICc(m))
```

## Running Tests

Execute the package test suite:
```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

## CLI Usage

Fit the default Birch model to a built-in dataset:
```bash
julia --project=. qpcr_fit.jl batsch1
```

Fit a custom dataset with the MAK2 model and export CSV output:
```bash
julia --project=. qpcr_fit.jl path/to/data.csv --model mak2 -o results.csv
```

### CLI Options
- `-m, --model <name>`: Model to fit (`birch`, `logistic`, `richards`, or `mak2`). Default: `birch`.
- `-o, --output <path>`: Output CSV file path.
- `--fit-noise`: Force fitting of noisy curves that show no amplification.
- `-h, --help`: Show help menu.
