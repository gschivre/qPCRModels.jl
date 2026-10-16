# Non-linear least squares qPCR model definitions and evaluation methods

"""
    NonLinearLeastSquare

Abstract base type for all non-linear least squares qPCR amplification models in `qPCRModels`.

Subtypes must implement parameter transformations (`get_params`), linear variable projection parameter computation (`get_linearparams!`), and the curve with its sensitivities to the non-linear parameters (`_curve_jacobian!`, used by the generic `residualsJacobian!`).
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
Trailing missing cycles (`NaN`) are dropped; missing values before the last measured cycle throw an `ArgumentError`.
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
        data = _trim_trailing_nonfinite(data)
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

# Heuristics to find approriate initial value for ln(x₀)
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
    if _birch_curve!(nlm; RFtol = RFtol, RFmaxiter = RFmaxiter)
        nlm.linearparams .= nlm.xmat \ nlm.data
    else
        nlm.linearparams .= NaN
    end
    return nlm.linearparams
end
get_linearparams_filtkwargs!(nlm::BirchModel; kwargs...) = get_linearparams!(nlm; filter_kwargs(kwargs, (:RFtol, :RFmaxiter))...)

"""
    get_x0(nlm::BirchModel; kwargs...)

Estimate initial target concentration `x₀` (at cycle 0) from fitted `BirchModel`.
"""
function get_x0(nlm::BirchModel; RFtol = 10 * eps(), RFmaxiter = 100)
    # solved in log space: the root finder's absolute tolerance would make x₀ inaccurate for late curves
    return exp(get_logx0(nlm; RFtol = RFtol, RFmaxiter = RFmaxiter))
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

"""
    get_cy0(nlm::BirchModel)

Compute the Cy₀ cycle for the `BirchModel` (see Guescini, M. et al. BMC Bioinformatics 9, 326 (2008)).
The Cy₀ cycle is defined here as the intersection of the tangent line at the inflexion point and the baseline.
"""
function get_cy0(nlm::BirchModel)
    b, tₘ = get_params(nlm)
    xᵢ = 1 / (1 + sqrt(1 + b)) # value at the inflexion point
    cᵢ = tₘ - b + (log(xᵢ) - (1 + b) * log1p(-xᵢ)) / log(2) # inflexion cycle
    return cᵢ - 1 / (log(2) * xᵢ)
end

# Fill `nlm.xmat[:, 2]` with the Birch curve, solving z - (1 + b) ln(1 - eᶻ) = ln(2) (c - tₘ + b) for z = ln x cycle by cycle,
# and, when `J` is given, the sensitivities of x to the non-linear parameters `(params[1], tₘ)`.
# Solving for ln x keeps the early cycles accurate in relative terms: in x, the root finder's absolute tolerance
# swamps values below ~1e-14 (late curves), which broke the monotone bracketing between cycles.
# Once x is within rounding of 1 the curve is frozen at its last value (plateau).
# Returns `false` when the parameters are not representable in Float64 (e.g. b = expm1(params[1]) overflows).
function _birch_curve!(nlm::BirchModel, J = nothing; RFtol = 10 * eps(), RFmaxiter = 100)
    b, tₘ = get_params(nlm)
    (isfinite(b) && isfinite(tₘ)) || return false
    g = let b = b
        z -> z - (1 + b) * log1mexp(z) # z < 0
    end
    dg = let b = b
        z -> 1 + (1 + b) / expm1(-z) # = 1 + (1 + b) eᶻ / (1 - eᶻ)
    end
    sens = let b = b # from implicit differentiation, in terms of z = ln x; vanishes at x = 1
        function (z)
            z < 0 || return (0.0, 0.0)
            dF = exp(-z) - (1 + b) / expm1(z) # 1/x + (1 + b)/(1 - x)
            return (-(1 + b) * (log1mexp(z) + log(2)) / dF, log(2) / dF)
        end
    end
    n = length(nlm.data)
    target = log(2) * (1 - tₘ + b)
    lo = min(target, 0.0) - 1.0 # g(z) ≈ z for z → -∞, step down until the cycle-1 root is bracketed
    while g(lo) > target
        lo = 2 * lo
        isfinite(lo) || return false
    end
    z = lo
    i = 1
    while i ≤ n
        target = log(2) * (i - tₘ + b)
        root_ok = g(lo) ≤ target # g(0) = +∞ so the upper bracket always holds
        if root_ok
            znew, _ = newtsafe_rootfinding(g, dg, lo, 0.0; x₀ = clamp(target, lo, -RFtol), RFtol = RFtol, RFmaxiter = RFmaxiter, offset = -target)
            root_ok = isfinite(znew)
        end
        if !root_ok # next root not representable above the previous one: plateau reached in Float64
            i == 1 && return false
            break
        end
        z = znew
        nlm.xmat[i, 2] = exp(z)
        isnothing(J) || ((J[i, 1], J[i, 2]) = sens(z))
        i += 1
        (z < -RFtol) || break # freeze on the plateau (x within rounding of 1)
        lo = z
    end
    for j in i:n
        nlm.xmat[j, 2] = exp(z)
        isnothing(J) || ((J[j, 1], J[j, 2]) = sens(z))
    end
    return true
end

_curve_jacobian!(J, nlm::BirchModel; RFtol = 10 * eps(), RFmaxiter = 100) = _birch_curve!(nlm, J; RFtol = RFtol, RFmaxiter = RFmaxiter)

### 4-parameters logistic model

"""
    LogisticModel{T} <: NonLinearLeastSquare

Standard 4-parameter logistic model for qPCR curves with symmetric sigmoidal growth kinetics.

# Constructor
```julia
m = LogisticModel(data; detection_threshold = 5.0)
```
Trailing missing cycles (`NaN`) are dropped; missing values before the last measured cycle throw an `ArgumentError`.
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
        data = _trim_trailing_nonfinite(data)
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

"""
    get_cy0(nlm::LogisticModel)

Compute the Cy₀ cycle for the `LogisticModel` (see Guescini, M. et al. BMC Bioinformatics 9, 326 (2008)).
The Cy₀ cycle is defined here as the intersection of the tangent line at the inflexion point and the baseline.
"""
function get_cy0(nlm::LogisticModel)
    a, tₘ = get_params(nlm)
    return tₘ - 2 / a
end

function _curve_jacobian!(J, nlm::LogisticModel)
    a, tₘ = get_params(nlm)
    for i in eachindex(nlm.data)
        x = logistic(a * (i - tₘ))
        nlm.xmat[i, 2] = x
        dFinv = x * (1 - x) # use the implicit form of the logistic ODE
        J[i, 1] = (tₘ - i) * a * dFinv
        J[i, 2] = a * dFinv
    end
    return true
end

### Richards's logistic model (5-parameters logistic)

"""
    RichardsModel{T} <: NonLinearLeastSquare

5-parameter generalized Richards logistic model introducing an asymmetry parameter `ν`.

# Constructor
```julia
m = RichardsModel(data; detection_threshold = 5.0)
```
Trailing missing cycles (`NaN`) are dropped; missing values before the last measured cycle throw an `ArgumentError`.
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
        data = _trim_trailing_nonfinite(data)
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

"""
    get_cy0(nlm::RichardsModel)

Compute the Cy₀ cycle for the `RichardsModel` (see Guescini, M. et al. BMC Bioinformatics 9, 326 (2008)).
The Cy₀ cycle is defined here as the intersection of the tangent line at the inflexion point and the baseline.
"""
function get_cy0(nlm::RichardsModel)
    a, tₘ, ν = get_params(nlm)
    return tₘ - (1 + ν + log(ν)) / a
end

function _curve_jacobian!(J, nlm::RichardsModel)
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
    return true
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
Trailing missing cycles (`NaN`) are dropped; missing values before the last measured cycle throw an `ArgumentError`.
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
        data = _trim_trailing_nonfinite(data)
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
    if _mak2_curve!(nlm, get_params(nlm); RFtol = RFtol, RFmaxiter = RFmaxiter)
        nlm.linearparams .= nlm.xmat \ view(nlm.data, 1:nlm.n)
    else
        nlm.linearparams .= NaN
    end
    return nlm.linearparams
end
get_linearparams_filtkwargs!(nlm::MAK2Model; kwargs...) = get_linearparams!(nlm; filter_kwargs(kwargs, (:RFtol, :RFmaxiter))...)

# Fill `nlm.xmat[:, 2]` with the MAK2 curve: we use the implicit solution of the ODE dx/dc = ln(2) ln(1 + x),
# i.e. Ei(z) - Ei(z₀) = ln(2) c with z = ln(1 + x).
# The ln(2) time scale makes the exponential phase double each cycle, as the MAK2 recursion x + ln(1 + x) does,
# so x₀ only differs from the recursion estimate by a constant factor (fold changes are preserved).
# Returns `false` when the curve is not representable in Float64 for this z₀ (Ei overflow or catastrophic
# cancellation for large z₀, underflow of z₀ to 0), which happens for absurd x₀ on non-amplifying windows.
function _mak2_curve!(nlm::MAK2Model, z₀; RFtol = 10 * eps(), RFmaxiter = 100)
    ((z₀ > 0) && isfinite(expinti(z₀))) || return false
    f = let z₀ = z₀
        z -> expinti(z) - expinti(z₀)
    end
    df = z -> exp(z) / z
    lxm₀ = logexpm1(z₀)
    z = z₀
    for i in 1:nlm.n
        offset = -log(2) * i
        hi = log1pexp(log(2) * (i + 1) + lxm₀) # x(c) ≤ x₀ 2ᶜ, one extra doubling keeps the root strictly inside
        flo, fhi = f(z) + offset, f(hi) + offset
        (isfinite(fhi) && (flo * fhi ≤ 0)) || return false # also false when NaN
        z, _ = newtsafe_rootfinding(f, df, z, hi; x₀ = z, RFtol = RFtol, RFmaxiter = RFmaxiter, offset = offset)
        nlm.xmat[i, 2] = expm1(z)
    end
    return true
end

function _curve_jacobian!(J, nlm::MAK2Model; RFtol = 10 * eps(), RFmaxiter = 100)
    z₀ = get_params(nlm)
    _mak2_curve!(nlm, z₀; RFtol = RFtol, RFmaxiter = RFmaxiter) || return false
    for i in 1:nlm.n
        J[i, 1] = -expm1(z₀) * log1p(nlm.xmat[i, 2]) / z₀
    end
    return true
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

# Data the model is fitted to (MAK2 only uses its window)
_fitdata(nlm::NonLinearLeastSquare) = nlm.data
_fitdata(nlm::MAK2Model) = view(nlm.data, 1:nlm.n)

# Compute residuals and Jacobian by Variable Projection.
# Each model's `_curve_jacobian!(J, nlm)` fills `nlm.xmat[:, 2]` with the curve x and `J` with -∂x/∂params,
# and returns `false` when the curve is not representable in Float64.
function residualsJacobian!(r, J, nlm::NonLinearLeastSquare; kwargs...)
    if !_curve_jacobian!(J, nlm; kwargs...)
        # trial step outside floating-point range: infinite χ² makes LevenbergMarquardt! reject it
        fill!(r, Inf)
        fill!(J, 0.0)
        return nothing
    end
    y = _fitdata(nlm)
    copyto!(r, y)
    F = svd(nlm.xmat)
    Finv = _pinv(F)
    nlm.linearparams .= Finv * y
    mul!(r, nlm.xmat, nlm.linearparams, -1.0, 1.0)
    @views for j in axes(J, 2)
        rc = dot(J[:, j], r) # residuals coupling term
        J[:, j] .*= nlm.linearparams[2]
        J[:, j] .-= F.U * (F.U' * J[:, j])
        J[:, j] .+= Finv[2, :] * rc
    end
    return nothing
end

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
    residuals(nlm::NonLinearLeastSquare; kwargs...)

Compute vector of model residuals.

# Keyword Arguments
- `scaled::Bool`: If `true`, returns residuals on standardized scale. If `false`, returns residuals on raw fluorescence scale.
"""
residuals(nlm::NonLinearLeastSquare; kwargs...) = residuals!(zeros(size(nlm.xmat, 1)), nlm; kwargs...)

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

### Standard errors

# Names of the physical non-linear parameters returned by `get_params` (MAK2 has none besides x₀)
_param_names(::BirchModel) = (:b, :tₘ)
_param_names(::LogisticModel) = (:a, :tₘ)
_param_names(::RichardsModel) = (:a, :tₘ, :ν)
_param_names(::MAK2Model) = ()

# Derivative of each physical parameter with respect to its internal parameter `nlm.params[k]`
function _dphysical(nlm::BirchModel)
    b, _ = get_params(nlm)
    return (1 + b, 1.0) # b = expm1(params[1])
end
_dphysical(nlm::LogisticModel) = (get_params(nlm)[1], 1.0) # a = exp(params[1])
function _dphysical(nlm::RichardsModel)
    a, _, ν = get_params(nlm)
    return (a, 1.0, ν) # a = exp(params[1]), ν = exp(params[3])
end
_dphysical(::MAK2Model) = ()

# Gradient of ln(x₀) with respect to `nlm.params`
function _logx0_gradient(nlm::BirchModel; kwargs...)
    b, _ = get_params(nlm)
    z₀ = get_logx0(nlm; kwargs...)
    dg = 1 + (1 + b) / expm1(-z₀) # implicit differentiation of z - (1 + b) ln(1 - eᶻ) = ln(2) (b - tₘ)
    return [(1 + b) * (log(2) + log1mexp(z₀)) / dg, -log(2) / dg]
end
function _logx0_gradient(nlm::LogisticModel)
    a, tₘ = get_params(nlm)
    d = logistic(a * tₘ) # d/du loglogistic(u) at u = -a tₘ
    return [-a * tₘ * d, -a * d]
end
function _logx0_gradient(nlm::RichardsModel)
    a, tₘ, ν = get_params(nlm)
    d = logistic(a * tₘ) / ν # ln x₀ = loglogistic(-a tₘ) / ν
    return [-a * tₘ * d, -a * d, -get_logx0(nlm)]
end
_logx0_gradient(::MAK2Model; kwargs...) = [1.0] # params[1] = ln(x₀)

# Gauss–Newton covariance s² (JᵀJ)⁻¹ of (params..., linearparams...) on the scaled data, where J is the Jacobian
# of the model c₁ + c₂ x in all parameters. `NaN` when JᵀJ is singular or the curve cannot be computed.
function _covariance(nlm::NonLinearLeastSquare; kwargs...)
    n, m = size(nlm.xmat, 1), length(nlm.params)
    J = zeros(n, m + 2)
    @views if _curve_jacobian!(J[:, 1:m], nlm; kwargs...)
        J[:, 1:m] .*= -nlm.linearparams[2] # _curve_jacobian! gives -∂x/∂params
        J[:, m + 1] .= 1.0
        J[:, m + 2] .= nlm.xmat[:, 2]
        H = cholesky!(Symmetric(J' * J); check = false)
        issuccess(H) && return reducedchisq(nlm; kwargs...) .* inv(H)
    end
    return fill(NaN, m + 2, m + 2)
end

"""
    get_stderrors(nlm::NonLinearLeastSquare; kwargs...)

Standard errors of a fitted model, as a `NamedTuple`: the non-linear parameters returned by `get_params`
(none for `MAK2Model`), the baseline `F_bg` and scale `F_max` (raw fluorescence units) and `logx0` (of `get_logx0`).

They come from the Gauss–Newton covariance `s² (JᵀJ)⁻¹` at the optimum, where `J` is the Jacobian of the model
in all its parameters and `s² = χ²/(n - k)`, and are propagated to the reported quantities by the delta method.
All values are `NaN` when the model is not fitted.
"""
function get_stderrors(nlm::NonLinearLeastSquare; kwargs...)
    names = (_param_names(nlm)..., :F_bg, :F_max, :logx0)
    isfitted(nlm) || return NamedTuple{names}(ntuple(_ -> NaN, length(names)))
    m = length(nlm.params)
    C = _covariance(nlm; kwargs...)
    se = sqrt.(max.(diag(C), 0.0))
    g = _logx0_gradient(nlm; kwargs...)
    se_logx0 = sqrt(max(dot(g, view(C, 1:m, 1:m), g), 0.0))
    se_nonlinear = _dphysical(nlm) .* Tuple(se[1:length(_param_names(nlm))])
    return NamedTuple{names}((se_nonlinear..., nlm.σ * se[m + 1], nlm.σ * se[m + 2], se_logx0))
end
