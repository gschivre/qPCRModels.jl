# Calibration (standard) curves for absolute quantification

"""
    CalibrationCurve

Standard curve `log₁₀ x₀ = intercept + slope · log₁₀ q` relating the model estimate `x₀` to a known template
quantity `q` (e.g. copies), fitted by [`calibration`](@ref). Use [`quantify`](@ref) to convert the `x₀` of a new
sample into an absolute quantity.

# Fields
- `slope`, `intercept`: fitted coefficients (`slope = 1` when it is fixed).
- `slope_se`, `intercept_se`: approximate standard errors (`slope_se = 0` when the slope is fixed).
- `r2`: coefficient of determination of the used points.
- `weights::Vector{Float64}`: final weight of each input point (bisquare weights in `[0, 1]` if robust, `0` for rejected or unused points).
- `used::BitVector`: inputs with a finite `x₀` (and, for models, a converged fit).
- `robust::Bool`, `fixed_slope::Bool`: fitting options.
- `converged::Bool`: `true` if the robust iterations converged.
"""
struct CalibrationCurve
    slope::Float64
    intercept::Float64
    slope_se::Float64
    intercept_se::Float64
    r2::Float64
    weights::Vector{Float64}
    used::BitVector
    robust::Bool
    fixed_slope::Bool
    converged::Bool
end

"""
    calibration(quantity, models; robust = true, fixed_slope = false)
    calibration(quantity, x0; robust = true, fixed_slope = false)

Fit the standard curve `log₁₀ x₀ = intercept + slope · log₁₀ q` from standards of known quantity `q`.

The second argument is either a vector of fitted models (their `get_logx0` is used, and models that did not
converge are left out) or a vector of `x₀` values (non-positive or non-finite values are left out).

For a mechanistic model, `x₀` is proportional to the template, so the slope should be close to 1: a slope far
from 1 points to a biased `x₀`. With `fixed_slope = true` the slope is set to 1 and only the intercept is fitted.

With `robust = true` (default) the line is fitted with Tukey's bisquare weights, so that outlier wells are
down-weighted (see the `weights` field); otherwise by ordinary least squares.
Standard errors come from the final weighted least-squares step and are approximate for the robust fit.

# Example
```julia
standards = [BirchModel(y) for y in standard_curves] # wells with 10⁵, 10⁴, 10³ and 10² copies
foreach(fit!, standards)
cal = calibration([1e5, 1e4, 1e3, 1e2], standards)
sample = BirchModel(sample_curve)
fit!(sample)
quantify(cal, sample) # copies in the sample well
```
"""
function calibration(quantity::AbstractVector{<:Real}, models::AbstractVector{<:NonLinearLeastSquare}; kwargs...)
    log10x0 = [hasconverged(m) ? get_logx0(m) / log(10) : NaN for m in models]
    return _calibration(quantity, log10x0; kwargs...)
end
function calibration(quantity::AbstractVector{<:Real}, x0::AbstractVector{<:Real}; kwargs...)
    return _calibration(quantity, [x > 0 ? log10(x) : NaN for x in x0]; kwargs...)
end

function _calibration(quantity, y; robust::Bool = true, fixed_slope::Bool = false)
    length(quantity) == length(y) || throw(DimensionMismatch("$(length(quantity)) quantities for $(length(y)) x₀ values"))
    all(>(0), quantity) || throw(ArgumentError("quantities must be positive"))
    x = log10.(quantity)
    used = BitVector(isfinite.(y) .& isfinite.(x))
    p = fixed_slope ? 1 : 2
    count(used) > p || throw(ArgumentError("calibration needs at least $(p + 1) usable points, got $(count(used))"))
    fixed_slope || length(unique(x[used])) > 1 || throw(ArgumentError("calibration needs at least two different quantities"))
    xu, yu = x[used], y[used]
    X = fixed_slope ? ones(length(xu), 1) : [ones(length(xu)) xu]
    yfit = fixed_slope ? yu .- xu : yu
    β, w, converged = robust ? _bisquare_irls(X, yfit) : (X \ yfit, ones(length(yu)), true)
    intercept, slope = β[1], fixed_slope ? 1.0 : β[2]
    r = yu .- (intercept .+ slope .* xu)
    dof = count(>(0), w) - p
    s² = dof > 0 ? sum(w .* abs2.(r)) / dof : NaN
    H = cholesky!(Symmetric(X' * (w .* X)); check = false) # singular if the outliers leave a single quantity
    se = issuccess(H) ? sqrt.(max.(s² .* diag(inv(H)), 0.0)) : fill(NaN, p)
    r2 = 1 - sum(abs2, r) / sum(abs2, yu .- mean(yu))
    weights = zeros(length(y))
    weights[used] .= w
    return CalibrationCurve(slope, intercept, fixed_slope ? 0.0 : se[2], se[1], r2, weights, used, robust, fixed_slope, converged)
end

"""
    quantify(cal::CalibrationCurve, m::NonLinearLeastSquare; kwargs...)
    quantify(cal::CalibrationCurve, x0::Real)

Absolute template quantity of a sample, in the units of the quantities used to build `cal`, from its fitted
model (via `get_logx0`) or its `x₀`. The model should be the same type as the one used for the calibration.
"""
quantify(cal::CalibrationCurve, m::NonLinearLeastSquare; kwargs...) = _quantify(cal, get_logx0(m; kwargs...) / log(10))
quantify(cal::CalibrationCurve, x0::Real) = _quantify(cal, x0 > 0 ? log10(x0) : NaN)
_quantify(cal::CalibrationCurve, log10x0) = exp10((log10x0 - cal.intercept) / cal.slope)
