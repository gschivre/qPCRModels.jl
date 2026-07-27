# Non-linear least squares qPCR model definitions and evaluation methods

"""
    NonLinearLeastSquare

Abstract base type for all non-linear least squares qPCR amplification models in `qPCRModels`.

Subtypes must implement parameter transformations (`get_params`), linear variable projection parameter computation (`get_linearparams!`), and residual/Jacobian calculations (`residualsJacobian!`).
"""
abstract type NonLinearLeastSquare end

### Birch logistic model

"""
    BirchModel{T} <: NonLinearLeastSquare

Mechanistic logistic qPCR model based on the Birch model kinetics for fluorescence growth curves.

Uses Variable Projection (VarPro) to decouple non-linear rate/position parameters from linear scale and offset parameters.

# Fields
- `rawdata::Vector{T}`: Original unscaled input fluorescence vector.
- `data::Vector{Float64}`: Scaled data normalized by standard deviation.
- `σ::Float64`: Standard deviation of the raw dataset.
- `params::Vector{Float64}`: Non-linear model parameters `[params[1], tₘ]`.
- `linearparams::Vector{Float64}`: Linear scale and offset parameters.
- `chisq::Float64`: Chi-squared merit function result after fitting.
- `fitted::Bool`: `true` if model fitting has been executed.
- `converged::Bool`: `true` if model optimization converged.
- `noamp::Bool`: `true` if input curve is detected as noise.
- `xmat::Matrix{Float64}`: Work buffer matrix for Variable Projection.

# Constructor
```julia
m = BirchModel(data; detection_threshold = 5.0)
```
"""
mutable struct BirchModel{T} <: NonLinearLeastSquare
    rawdata::Vector{T}
    data::Vector{Float64}
    σ::Float64 # data standard deviation
    params::Vector{Float64}
    linearparams::Vector{Float64}
    chisq::Float64
    fitted::Bool
    converged::Bool
    noamp::Bool
    xmat::Matrix{Float64} # matrix buffer for Variable Projection

    function BirchModel(data::Vector{T}; detection_threshold::Real = 5.0) where {T}
        init_params = [0.0; hac(data; linear_interp = true)] # logistic start
        n = length(data)
        σ = std(data)
        noamp = isnoise(data; detection_threshold)
        return new{T}(data, data ./ σ, σ, init_params, NaN * Vector{Float64}(undef, 2), NaN, false, false, noamp, ones(n, 2))
    end
end

"""
    get_params(nlm::BirchModel)

Extract transformed physical parameters `(b, tₘ)` from a `BirchModel`.
- `b`: Reaction efficiency parameter (`b > -1`).
- `tₘ`: Cycle number at half-amplification inflection point.
"""
function get_params(nlm::BirchModel)
    b = max(eps() - 1, expm1(nlm.params[1])) # b > -1
    tₘ = nlm.params[2]
    return (b, tₘ)
end

# Heuristics to find approriate initial value for x₀
function _get_x0_start(b, tₘ)
    if b > 50
        return 0.5
    elseif b < -0.5
        if tₘ > 0
            return max(1e-9, exp(-log(2) * (tₘ + 1)))
        else
            return 1.0 - 1e-9
        end
    else
        return max(1e-9, min(1.0 - 1e-9, logistic(-log(2) * tₘ)))
    end
end
function _get_logx0_start(b, tₘ)
    if b > 50
        return -log(2)
    elseif b < -0.5
        if tₘ > 0
            return -log(2) * (tₘ + 1)
        else
            return -1e-9
        end
    else
        return min(-1e-9, loglogistic(-log(2) * tₘ))
    end
end

"""
    get_linearparams!(nlm::BirchModel; RFtol = 10 * eps(), RFmaxiter = 100)

Compute linear scale and baseline offset parameters for `BirchModel` via linear least squares on the feature matrix `xmat`.
"""
function get_linearparams!(nlm::BirchModel; RFtol = 10 * eps(), RFmaxiter = 100)
    b, tₘ = get_params(nlm)
    f = let b = b
        x -> log(x) - (1 + b) * log1p(-x)
    end
    df = let b = b
        x -> 1 / x + (1 + b) / (1 - x)
    end
    x, _ = newtsafe_rootfinding(f, df, 0.0, 1.0; x₀ = _get_x0_start(b, tₘ), RFtol = RFtol, RFmaxiter = RFmaxiter, offset = log(2) * (tₘ - b - 1))
    nlm.xmat[1, 2] = x
    i = 2
    while (i ≤ length(nlm.data)) && (RFtol < x < (1 - RFtol))
        x, _ = newtsafe_rootfinding(f, df, x, 1.0; x₀ = x, RFtol = RFtol, RFmaxiter = RFmaxiter, offset = log(2) * (tₘ - b - i))
        nlm.xmat[i, 2] = x
        i += 1
    end
    for j in i:length(nlm.data)
        nlm.xmat[j, 2] = x
    end
    nlm.linearparams .= nlm.xmat \ nlm.data
    return nlm.linearparams
end
get_linearparams_filtkwargs!(nlm::BirchModel; kwargs...) = get_linearparams!(nlm; filter_kwargs(kwargs, (:RFtol, :RFmaxiter))...)

"""
    get_x0(nlm::BirchModel; kwargs...)

Estimate initial target concentration `x₀` (at cycle 0) from fitted `BirchModel`.
"""
function get_x0(nlm::BirchModel; RFtol = 10 * eps(), RFmaxiter = 100)
    b, tₘ = get_params(nlm)
    f = let b = b
        x -> log(x) - (1 + b) * log1p(-x)
    end
    df = let b = b
        x -> 1 / x + (1 + b) / (1 - x)
    end
    x₀, _ = newtsafe_rootfinding(f, df, 0.0, 1.0; x₀ = _get_x0_start(b, tₘ), RFtol = RFtol, RFmaxiter = RFmaxiter, offset = log(2) * (tₘ - b))
    return x₀
end

"""
    get_logx0(nlm::BirchModel; kwargs...)

Estimate log-transformed initial target concentration `ln(x₀)` from fitted `BirchModel`.
"""
function get_logx0(nlm::BirchModel; RFtol = 10 * eps(), RFmaxiter = 100)
    b, tₘ = get_params(nlm)
    f = let b = b
        z -> z - (1 + b) * log1mexp(z)
    end
    df = let b = b
        z -> 1 + (1 + b) * exp(z) / (1 - exp(z)) # z < 0
    end
    z₀, _ = newtsafe_rootfinding(f, df, -log(2) * (max(0.0, tₘ) + 1), 0.0; x₀ = _get_logx0_start(b, tₘ), RFtol = RFtol, RFmaxiter = RFmaxiter, offset = log(2) * (tₘ - b))
    return z₀
end

# Compute residuals and Jacobian
function residualsJacobian!(r, J, nlm::BirchModel; RFtol = 10 * eps(), RFmaxiter = 100)
    b, tₘ = get_params(nlm)
    f = let b = b
        x -> log(x) - (1 + b) * log1p(-x)
    end
    df = let b = b
        x -> 1 / x + (1 + b) / (1 - x)
    end
    x, _ = newtsafe_rootfinding(f, df, 0.0, 1.0; x₀ = _get_x0_start(b, tₘ), RFtol = RFtol, RFmaxiter = RFmaxiter, offset = log(2) * (tₘ - b - 1))
    nlm.xmat[1, 2] = x
    dF = df(x)
    J[1, 1] = -(1 + b) * (log1p(-x) + log(2)) / dF
    J[1, 2] = log(2) / dF
    i = 2
    while (i ≤ length(nlm.data)) && (RFtol < x < (1 - RFtol))
        x, _ = newtsafe_rootfinding(f, df, x, 1.0; x₀ = x, RFtol = RFtol, RFmaxiter = RFmaxiter, offset = log(2) * (tₘ - b - i))
        nlm.xmat[i, 2] = x
        dF = df(x)
        J[i, 1] = -(1 + b) * (log1p(-x) + log(2)) / dF
        J[i, 2] = log(2) / dF
        i += 1
    end
    for j in i:length(nlm.data)
        nlm.xmat[j, 2] = x
        J[j, 1] = -(1 + b) * (log1p(-x) + log(2)) / dF
        J[j, 2] = log(2) / dF
    end
    copyto!(r, nlm.data)
    F = svd(nlm.xmat)
    Finv = pinv(F)
    nlm.linearparams .= Finv * nlm.data
    mul!(r, nlm.xmat, nlm.linearparams, -1.0, 1.0)
    rc = 0.0 # residuals coupling term
    @views for j in 1:2
        rc = dot(J[:, j], r)
        J[:, j] .*= nlm.linearparams[2]
        J[:, j] .-= F.U * (F.U' * J[:, j])
        J[:, j] .+= Finv[2, :] * rc
    end
    return nothing
end

### 4-parameters logistic model

"""
    LogisticModel{T} <: NonLinearLeastSquare

Standard 4-parameter logistic model for qPCR curves with symmetric sigmoidal growth kinetics.

# Constructor
```julia
m = LogisticModel(data; detection_threshold = 5.0)
```
"""
mutable struct LogisticModel{T} <: NonLinearLeastSquare
    rawdata::Vector{T}
    data::Vector{Float64}
    σ::Float64 # data standard deviation
    params::Vector{Float64}
    linearparams::Vector{Float64}
    chisq::Float64
    fitted::Bool
    converged::Bool
    noamp::Bool
    xmat::Matrix{Float64} # matrix buffer for Variable Projection

    function LogisticModel(data::Vector{T}; detection_threshold::Real = 5.0) where {T}
        init_params = [log(log(2)); hac(data; linear_interp = true)] # starting rate ~ln(2)
        n = length(data)
        σ = std(data)
        noamp = isnoise(data; detection_threshold)
        return new{T}(data, data ./ σ, σ, init_params, NaN * Vector{Float64}(undef, 2), NaN, false, false, noamp, ones(n, 2))
    end
end

"""
    get_params(nlm::LogisticModel)

Return non-linear parameter tuple `(a, tₘ)` where `a` is amplification rate and `tₘ` is mid-point cycle.
"""
function get_params(nlm::LogisticModel)
    # we could have used log(1.5 + 0.5 * p / (1 + abs(p))) to enforce 0 < a < ln(2)
    # but the ln(2) rate would only be reached at the infinite limit
    a = exp(nlm.params[1]) # a > 0 so rate could be higher than the perfect doubling
    tₘ = nlm.params[2]
    return (a, tₘ)
end


"""
    get_linearparams!(nlm::LogisticModel)

Compute linear scale and baseline offset parameters for `LogisticModel` via linear least squares on the feature matrix `xmat`.
"""
function get_linearparams!(nlm::LogisticModel)
    a, tₘ = get_params(nlm)
    for i in eachindex(nlm.data)
        nlm.xmat[i, 2] = logistic(a * (i - tₘ))
    end
    nlm.linearparams .= nlm.xmat \ nlm.data
    return nlm.linearparams
end
get_linearparams_filtkwargs!(nlm::LogisticModel; kwargs...) = get_linearparams!(nlm; filter_kwargs(kwargs, ())...)

"""
    get_x0(nlm::LogisticModel)

Estimate initial concentration `x₀` from fitted 4-parameter `LogisticModel`.
"""
function get_x0(nlm::LogisticModel)
    a, tₘ = get_params(nlm)
    return logistic(-a * tₘ)
end

"""
    get_logx0(nlm::LogisticModel)

Estimate log initial concentration `ln(x₀)` from fitted 4-parameter `LogisticModel`.
"""
function get_logx0(nlm::LogisticModel)
    a, tₘ = get_params(nlm)
    return loglogistic(-a * tₘ)
end

# Compute residuals and Jacobian
function residualsJacobian!(r, J, nlm::LogisticModel)
    a, tₘ = get_params(nlm)
    for i in eachindex(nlm.data)
        x = logistic(a * (i - tₘ))
        nlm.xmat[i, 2] = x
        dFinv = x * (1 - x) # use the implicit form of the logistic ODE
        J[i, 1] = (tₘ - i) * a * dFinv
        J[i, 2] = a * dFinv
    end
    copyto!(r, nlm.data)
    F = svd(nlm.xmat)
    Finv = pinv(F)
    nlm.linearparams .= Finv * nlm.data
    mul!(r, nlm.xmat, nlm.linearparams, -1.0, 1.0)
    rc = 0.0 # residuals coupling term
    @views for j in 1:2
        rc = dot(J[:, j], r)
        J[:, j] .*= nlm.linearparams[2]
        J[:, j] .-= F.U * (F.U' * J[:, j])
        J[:, j] .+= Finv[2, :] * rc
    end
    return nothing
end

### Richards's logistic model (5-parameters logistic)

"""
    RichardsModel{T} <: NonLinearLeastSquare

5-parameter generalized Richards logistic model introducing an asymmetry parameter `ν`.

# Constructor
```julia
m = RichardsModel(data; detection_threshold = 5.0)
```
"""
mutable struct RichardsModel{T} <: NonLinearLeastSquare
    rawdata::Vector{T}
    data::Vector{Float64}
    σ::Float64 # data standard deviation
    params::Vector{Float64}
    linearparams::Vector{Float64}
    chisq::Float64
    fitted::Bool
    converged::Bool
    noamp::Bool
    xmat::Matrix{Float64} # matrix buffer for Variable Projection

    function RichardsModel(data::Vector{T}; detection_threshold::Real = 5.0) where {T}
        init_params = [log(log(2)); hac(data; linear_interp = true); 0.0] # logistic start
        n = length(data)
        σ = std(data)
        noamp = isnoise(data; detection_threshold)
        return new{T}(data, data ./ σ, σ, init_params, NaN * Vector{Float64}(undef, 2), NaN, false, false, noamp, ones(n, 2))
    end
end

"""
    RichardsLogistic(t, ν)

Evaluate Richards asymmetric logistic function at normalized time/cycle `t` with exponent `ν`.
"""
function RichardsLogistic(t, ν)
    return exp(logRichardsLogistic(t, ν))
end

"""
    logRichardsLogistic(t, ν)

Numerically stable log-evaluation of Richards asymmetric logistic curve.
"""
function logRichardsLogistic(t, ν)
    if t < 0
        return (t - log1pexp(t)) / ν
    else
        return -log1pexp(-t) / ν
    end
end

"""
    get_params(nlm::RichardsModel)

Return non-linear parameter tuple `(a, tₘ, ν)` for `RichardsModel`.
"""
function get_params(nlm::RichardsModel)
    # we could have used log(1.5 + 0.5 * p / (1 + abs(p))) to enforce 0 < a < ln(2)
    # but the ln(2) rate would only be reached at the infinite limit
    a = exp(nlm.params[1]) # a > 0 so rate could be higher than the perfect doubling
    tₘ = nlm.params[2]
    ν = exp(nlm.params[3]) # ν > 0
    return (a, tₘ, ν)
end

"""
    get_linearparams!(nlm::RichardsModel)

Compute linear scale and baseline offset parameters for `RichardsModel` via linear least squares on the feature matrix `xmat`.
"""
function get_linearparams!(nlm::RichardsModel)
    a, tₘ, ν = get_params(nlm)
    for i in eachindex(nlm.data)
        nlm.xmat[i, 2] = RichardsLogistic(a * (i - tₘ), ν)
    end
    nlm.linearparams .= nlm.xmat \ nlm.data
    return nlm.linearparams
end
get_linearparams_filtkwargs!(nlm::RichardsModel; kwargs...) = get_linearparams!(nlm; filter_kwargs(kwargs, ())...)

"""
    get_x0(nlm::RichardsModel)

Estimate initial target concentration `x₀` from fitted `RichardsModel`.
"""
function get_x0(nlm::RichardsModel)
    a, tₘ, ν = get_params(nlm)
    return RichardsLogistic(-a * tₘ, ν)
end

"""
    get_logx0(nlm::RichardsModel)

Estimate log initial target concentration `ln(x₀)` from fitted `RichardsModel`.
"""
function get_logx0(nlm::RichardsModel)
    a, tₘ, ν = get_params(nlm)
    return logRichardsLogistic(-a * tₘ, ν)
end

# Compute residuals and Jacobian
function residualsJacobian!(r, J, nlm::RichardsModel)
    a, tₘ, ν = get_params(nlm)
    for i in eachindex(nlm.data)
        logx = logRichardsLogistic(a * (i - tₘ), ν)
        x = exp(logx)
        nlm.xmat[i, 2] = x
        dFinv = x * (1 - x) # use the implicit form of the logistic ODE
        J[i, 1] = (tₘ - i) * a * dFinv
        J[i, 2] = a * dFinv
        J[i, 3] = ν * logx * exp(-log1mexp(ν * logx)) * dFinv
        if isnan(J[i, 3])
            J[i, 3] = 0.0
        end
    end
    copyto!(r, nlm.data)
    F = svd(nlm.xmat)
    Finv = pinv(F)
    nlm.linearparams .= Finv * nlm.data
    mul!(r, nlm.xmat, nlm.linearparams, -1.0, 1.0)
    rc = 0.0 # residuals coupling term
    @views for j in 1:3
        rc = dot(J[:, j], r)
        J[:, j] .*= nlm.linearparams[2]
        J[:, j] .-= F.U * (F.U' * J[:, j])
        J[:, j] .+= Finv[2, :] * rc
    end
    return nothing
end

### MAK2 model

"""
    MAK2Model{T} <: NonLinearLeastSquare

Mechanistic mass-action kinetic (MAK2) model for qPCR amplification curves.

Models amplification based on reactant accumulation up to cut-off cycle `n`.

# Constructor
```julia
m = MAK2Model(data; detection_threshold = 5.0)
```
"""
mutable struct MAK2Model{T} <: NonLinearLeastSquare
    rawdata::Vector{T}
    data::Vector{Float64}
    σ::Float64 # data standard deviation
    n::Int # last cycle included in the fit
    params::Vector{Float64}
    linearparams::Vector{Float64}
    chisq::Float64
    fitted::Bool
    converged::Bool
    noamp::Bool
    xmat::Matrix{Float64} # matrix buffer for Variable Projection

    function MAK2Model(data::Vector{T}; detection_threshold::Real = 5.0) where {T}
        ct = 4 + sdm(view(data, 4:length(data)); parabolic_interp = true) # minimum five point for model fitting
        init_params = [-log(2) * (ct + 1)]
        n::Int = floor(ct)
        σ = std(data)
        noamp = isnoise(data; detection_threshold)
        return new{T}(data, data ./ σ, σ, n, init_params, NaN * Vector{Float64}(undef, 2), NaN, false, false, noamp, ones(n, 2))
    end
end

"""
    get_x0(nlm::MAK2Model)

Get initial target concentration `x₀` from `MAK2Model`.
"""
function get_x0(nlm::MAK2Model)
    x₀ = exp(nlm.params[1]) # x₀ > 0
    return x₀
end

"""
    get_logx0(nlm::MAK2Model)

Get log initial concentration `ln(x₀)` from `MAK2Model`.
"""
get_logx0(nlm::MAK2Model) = nlm.params[1]

"""
    get_params(nlm::MAK2Model)

Return non-linear parameter `z₀ = ln(1 + x₀)` for `MAK2Model`.
"""
function get_params(nlm::MAK2Model)
    return log1pexp(nlm.params[1]) # ln(1 + x₀)
end

"""
    get_linearparams!(nlm::MAK2Model; RFtol = 10 * eps(), RFmaxiter = 100)

Compute linear scale and baseline offset parameters for `MAK2Model` via linear least squares on the feature matrix `xmat`.
"""
function get_linearparams!(nlm::MAK2Model; RFtol = 10 * eps(), RFmaxiter = 100)
    z₀ = get_params(nlm)
    f = let z₀ = z₀
        z -> expinti(z) - expinti(z₀)
    end
    df = z -> exp(z) / z
    fz₀ = f(z₀)
    z, _ = newtsafe_rootfinding(f, df, z₀, log1pexp(1 + logexpm1(z₀)); x₀ = z₀, RFtol = RFtol, RFmaxiter = RFmaxiter, offset = -(1 + fz₀))
    nlm.xmat[1, 2] = expm1(z)
    for i in 2:nlm.n
        z, _ = newtsafe_rootfinding(f, df, z, log1pexp(i + logexpm1(z₀)); x₀ = z, RFtol = RFtol, RFmaxiter = RFmaxiter, offset = -(i + fz₀))
        nlm.xmat[i, 2] = expm1(z)
    end
    nlm.linearparams .= nlm.xmat \ view(nlm.data, 1:nlm.n)
    return nlm.linearparams
end
get_linearparams_filtkwargs!(nlm::MAK2Model; kwargs...) = get_linearparams!(nlm; filter_kwargs(kwargs, (:RFtol, :RFmaxiter))...)

# Compute residuals and Jacobian
function residualsJacobian!(r, J, nlm::MAK2Model; RFtol = 10 * eps(), RFmaxiter = 100)
    # we use the implicit ODE, given it is autonomous this equivalent to fit the recurence relation
    z₀ = get_params(nlm)
    f = let z₀ = z₀
        z -> expinti(z) - expinti(z₀)
    end
    df = z -> exp(z) / z
    fz₀ = f(z₀)
    z, _ = newtsafe_rootfinding(f, df, z₀, log1pexp(1 + logexpm1(z₀)); x₀ = z₀, RFtol = RFtol, RFmaxiter = RFmaxiter, offset = -(1 + fz₀))
    nlm.xmat[1, 2] = expm1(z)
    J[1] = -expm1(z₀) * z / z₀
    for i in 2:nlm.n
        z, _ = newtsafe_rootfinding(f, df, z, log1pexp(i + logexpm1(z₀)); x₀ = z, RFtol = RFtol, RFmaxiter = RFmaxiter, offset = -(i + fz₀))
        nlm.xmat[i, 2] = expm1(z)
        J[i] = -expm1(z₀) * z / z₀
    end
    copyto!(r, view(nlm.data, 1:nlm.n))
    F = svd(nlm.xmat)
    Finv = pinv(F)
    nlm.linearparams .= Finv * view(nlm.data, 1:nlm.n)
    mul!(r, nlm.xmat, nlm.linearparams, -1.0, 1.0)
    rc = dot(J, r) # residuals coupling term
    rmul!(J, nlm.linearparams[2])
    J .-= F.U * (F.U' * J)
    J .+= view(Finv, 2, :) * rc
    return nothing
end

"""
    update_n!(nlm::MAK2Model, n::Int)

Update the maximum cycle count `n` considered during fitting of a `MAK2Model`.
"""
function update_n!(nlm::MAK2Model, n::Int)
    @assert (n > 4) "`MAK2Model` has 3 parameters (k) and AICc involves n-k-1 so n must be ≥ 5!"
    nlm.n = n
    nlm.xmat = ones(n, 2)
    nlm.params .= -log(2) * (n + 1)
    nlm.linearparams .= NaN
    nlm.chisq = NaN
    nlm.fitted = false
    nlm.converged = false
    return nothing
end

"""
    refit!(nlm::MAK2Model; kwargs...)

Iteratively fit a `MAK2Model` over candidate cycle bounds `n` and select the cycle length minimizing AICc.
"""
function refit!(nlm::MAK2Model; kwargs...)
    out = zeros(length(nlm.data) - 4)
    for (i, n) in Iterators.enumerate(5:length(nlm.data)) # MAK2Model has 3 parameters (k) and AICc involves n-k-1
        update_n!(nlm, n)
        fit!(nlm; kwargs...)
        out[i] = AICc(nlm; filter_kwargs(kwargs, (:RFtol, :RFmaxiter))...)
    end
    _, n = findmin(out)
    update_n!(nlm, n + 4)
    fit!(nlm; kwargs...)
    return out
end

### Generic evaluation functions

"""
    fit!(nlm::NonLinearLeastSquare; fit_noise::Bool = false, kwargs...)

Fit non-linear model `nlm` to its fluorescence data vector using `LevenbergMarquardt!`.

# Keyword Arguments
- `fit_noise::Bool`: If `false` (default), skips optimization if curve was detected as noise (`nlm.noamp == true`).
- `kwargs...`: Additional options passed to optimizer (e.g. `LMtol`, `LMmaxiter`).
"""
function fit!(nlm::NonLinearLeastSquare; fit_noise::Bool = false, kwargs...)
    (nlm.noamp && !fit_noise) && return nothing
    chisq, converged = LevenbergMarquardt!(nlm; kwargs...)
    get_linearparams_filtkwargs!(nlm; kwargs...)
    nlm.chisq = chisq
    nlm.fitted = true
    nlm.converged = converged
    return nothing
end

"""
    hasamplification(nlm::NonLinearLeastSquare)

Return `true` if curve contains real amplification signal (not noise).
"""
hasamplification(nlm::NonLinearLeastSquare) = !nlm.noamp

"""
    hasconverged(nlm::NonLinearLeastSquare)

Return `true` if the non-linear optimization solver converged cleanly.
"""
hasconverged(nlm::NonLinearLeastSquare) = nlm.converged

"""
    isfitted(nlm::NonLinearLeastSquare)

Return `true` if `fit!` has been called on model `nlm`.
"""
isfitted(nlm::NonLinearLeastSquare) = nlm.fitted

"""
    residualsJacobian(nlm::NonLinearLeastSquare; kwargs...)

Return a tuple `(r, J)` of the residual vector `r` and Jacobian matrix `J`.

# Keyword Arguments
- `scaled::Bool`: If `true`, returns residuals on standardized scale. If `false`, returns residuals on raw fluorescence scale.
"""
function residualsJacobian(nlm::NonLinearLeastSquare; kwargs...)
    n = size(nlm.xmat, 1)
    m = length(nlm.params)
    r = zeros(n)
    J = zeros(n, m)
    residualsJacobian!(r, J, nlm; kwargs...)
    return (r, J)
end

"""
    residuals(nlm::NonLinearLeastSquare, kwargs...)

Compute vector of model residuals.

# Keyword Arguments
- `scaled::Bool`: If `true`, returns residuals on standardized scale. If `false`, returns residuals on raw fluorescence scale.
"""
residuals(nlm::NonLinearLeastSquare, kwargs...) = residuals!(zeros(size(nlm.xmat, 1)), nlm, kwargs...)

"""
    residuals!(r, nlm::NonLinearLeastSquare; scaled::Bool = true, kwargs...)

Fill pre-allocated vector `r` with model residual values `y_obs - y_pred`.

# Keyword Arguments
- `scaled::Bool`: If `true`, returns residuals on standardized scale. If `false`, returns residuals on raw fluorescence scale.
"""
function residuals!(r, nlm::NonLinearLeastSquare; scaled::Bool = true, kwargs...)
    nlm.fitted || get_linearparams!(nlm; kwargs...)
    for i in eachindex(r)
        r[i] = nlm.data[i]
    end
    if scaled # on the standardized scale
        return mul!(r, nlm.xmat, nlm.linearparams, -1.0, 1.0)
    else # same scale as original data
        return mul!(r, nlm.xmat, nlm.linearparams, -nlm.σ, nlm.σ)
    end
end

"""
    chisq(nlm::NonLinearLeastSquare; kwargs...)

Compute Chi-squared goodness-of-fit merit sum `∑ rᵢ²`.
"""
function chisq(nlm::NonLinearLeastSquare; kwargs...)
    if nlm.fitted
        return nlm.chisq
    else
        return sum(abs2, residuals(nlm; kwargs...))
    end
end

"""
    reducedchisq(nlm::NonLinearLeastSquare; kwargs...)

Compute reduced Chi-squared statistic `χ² / (n - k)` accounting for degrees of freedom.
"""
function reducedchisq(nlm::NonLinearLeastSquare; kwargs...)
    n = size(nlm.xmat, 1)
    k = length(nlm.params) + length(nlm.linearparams)
    return chisq(nlm; kwargs...) / (n - k)
end

"""
    get_predicted(nlm::NonLinearLeastSquare; scaled::Bool = true, kwargs...)

Compute fitted curve predictions `y_pred` across all cycles.

# Keyword Arguments
- `scaled::Bool`: Return values on standardized scale if `true`, or raw scale if `false`.
"""
function get_predicted(nlm::NonLinearLeastSquare; scaled::Bool = true, kwargs...)
    nlm.fitted || get_linearparams!(nlm; kwargs...)
    if scaled
        return nlm.xmat * nlm.linearparams
    else
        return nlm.xmat * nlm.linearparams * nlm.σ
    end
end

"""
    get_linearparams(nlm::NonLinearLeastSquare; kwargs...)

Return fitted linear scale and offset parameter vector.
"""
function get_linearparams(nlm::NonLinearLeastSquare; kwargs...)
    if nlm.fitted
        return nlm.linearparams
    else
        return get_linearparams!(nlm; kwargs...)
    end
end

"""
    AICc(nlm::NonLinearLeastSquare; kwargs...)

Compute Corrected Akaike Information Criterion (AICc) for model selection and quality assessment.
"""
function AICc(nlm::NonLinearLeastSquare; kwargs...)
    chi2 = chisq(nlm; kwargs...)
    n = size(nlm.xmat, 1)
    k = length(nlm.params) + length(nlm.linearparams)
    return n * log(chi2 / n) + 2 * k + n * (1 + log(2π)) + 2 * (k ^ 2 + k) / (n - k - 1)
end
