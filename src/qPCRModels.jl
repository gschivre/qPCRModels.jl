"""
    qPCRModels

`qPCRModels` is a Julia package for non-linear least squares modeling and analysis of real-time quantitative PCR (qPCR) fluorescence amplification curves.

# Key Features
- **Model Implementations**:
  - `BirchModel`: Mechanistic logistic model accounting for reaction kinetics.
  - `LogisticModel`: Standard 4-parameter logistic curve model.
  - `RichardsModel`: 5-parameter generalized logistic model with asymmetry factor `ν`.
  - `MAK2Model`: Mechanistic mass-action kinetic qPCR model.
- **Solvers & Root-finding**:
  - `LevenbergMarquardt!`: Dampened non-linear least squares optimization.
  - `newtsafe_rootfinding`: Safe Newton-Raphson rootfinder with bisection fallback.
- **Diagnostics & Metrics**:
  - `cycletreshold`, `hac`, `sdm`: Threshold, half-amplification cycle, and second-derivative maximum estimators.
  - `isnoise`: Noise and non-amplification detection heuristic based on Median Absolute Deviation.
  - `AICc`, `chisq`, `reducedchisq`, `get_x0`, `get_logx0`: Model quality criteria and initial concentration estimators.
- **Built-in Datasets**:
  - `qpcRDataset`, `load_qpcR_dataset`: Access to benchmark qPCR datasets (e.g. `batsch1`, `boggy`, `testdat`).
"""
module qPCRModels

using LinearAlgebra
using LogExpFunctions
using SpecialFunctions
using Statistics
using Printf

# Include submodules and component files (models before solver for NonLinearLeastSquare type)
include("utils.jl")
include("models.jl")
include("solver.jl")
include("datasets.jl")

# Export public interface
export NonLinearLeastSquare, BirchModel, LogisticModel, RichardsModel, MAK2Model
export LevenbergMarquardt!, newtsafe_rootfinding
export filter_kwargs, cycletreshold, hac, sdm, mad, isnoise
export fit!, refit!, update_n!, hasamplification, hasconverged, isfitted
export residualsJacobian, residualsJacobian!, residuals!, residuals
export chisq, reducedchisq, get_predicted, get_params, get_linearparams, get_linearparams!
export get_x0, get_logx0, AICc
export qpcRDataset, load_qpcR_dataset, DatasetDetails

end # module qPCRModels
