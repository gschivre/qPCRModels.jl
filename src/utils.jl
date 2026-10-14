# Utility functions for qPCR models

# Filter keyword arguments
filter_kwargs(kwargs, allowed) = (; (k => v for (k, v) in kwargs if k in allowed)...)

"""
    cycletreshold(x::AbstractVector, val::Real; linear_interp::Bool = false)

Estimate the fractional cycle number at which fluorescence signal `x` reaches target threshold `val`.

# Arguments
- `x::AbstractVector`: Fluorescence values across cycles.
- `val::Real`: Target threshold value.
- `linear_interp::Bool`: If `true`, applies linear interpolation between neighboring cycles for sub-cycle resolution.

# Returns
- Fractional cycle number (`Float64`), or `NaN` if `val` is never reached.
"""
function cycletreshold(x::AbstractVector, val::Real; linear_interp::Bool = false)
    (maximum(x) < val) && return Nan
    ct = findfirst(x .≥ val)
    if linear_interp && (ct > 1)
        δ = (val - x[ct]) / (x[ct] - x[ct - 1])
    else
        δ = 0.0
    end
    return ct + δ # type stable
end

"""
    hac(x::AbstractVector; kwargs...)

Compute the Half Amplification Cycle (HAC), defined as the fractional cycle where fluorescence reaches 50% of the total dynamic range `[min(x), max(x)]`.

# Arguments
- `x::AbstractVector`: Vector of fluorescence readings per cycle.
- `kwargs...`: Additional options passed to `cycletreshold` (e.g. `linear_interp = true`).

# Returns
- Fractional cycle number corresponding to 50% max fluorescence.
"""
function hac(x::AbstractVector; kwargs...)
    xmin = minimum(x)
    xmax = maximum(x)
    h = xmin + 0.5 * (xmax - xmin)
    return cycletreshold(x, h; kwargs...)
end

"""
    sdm(x::AbstractVector; parabolic_interp::Bool = false)

Calculate the Second Derivative Maximum (SDM) cycle number using Savitzky-Golay smoothing (order 4, window size 7).

# Arguments
- `x::AbstractVector`: Fluorescence curve data.
- `parabolic_interp::Bool`: If `true`, applies parabolic interpolation around the peak derivative for sub-cycle accuracy.

# Returns
- Cycle number corresponding to maximum acceleration (inflection point derivative maximum).
"""
function sdm(x::AbstractVector; parabolic_interp::Bool = false)
    #=
    The Savitzky–Golay coefficient could be obtained via:
    hw = window_length >> 1
    V = [i ^ j for i in -hw:hw, j in 0:order]
    VVinv = inv(V' * V)
    sg_coeff = factorial(derivative_order) .* (sum([VVinv[derivative_order + 1, j + 1] * (i ^ j) for j in 0:order]) for i in -hw:hw)
    =#
    sg_coeff = (-13, 67, -19, -70, -19, 67, -13) ./ 132 # order 4, window length of 7
    left_padding = (x[1] + x[2] + x[3]) / 3 # first cycle should ~baseline so constant padding is the best
    b, a = ((2 * x[end - 1] + 5 * x[end] - x[end - 2]) / 6, (x[end] - x[end - 2]) / 2) # linear padding as the last cycle not always reached the plateau
    right_padding = b .+ a .* (1, 2, 3)
    n = length(x)
    d2x = zeros(n)
    for i in 1:n
        for (jj, j) in Iterators.enumerate(-3:3)
            if (i + j) < 1
                val = left_padding
            elseif (i + j) > n
                val = right_padding[j - (n - i)]
            else
                val = x[i + j]
            end
            d2x[i] += val * sg_coeff[jj]
        end
    end
    aₘ, ct = findmax(d2x)
    if parabolic_interp # use parabolic interpolation to get sub-cycle resolution
        aₘ₋₁ = (ct > 1 ? d2x[ct - 1] : 0.0) # use zero padding for the second order derivative
        aₘ₊₁ = (ct < (n - 1) ? d2x[ct + 1] : 0.0)
        δ = 0.5 * (aₘ₋₁ - aₘ₊₁) / (aₘ₋₁ - 2 * aₘ + aₘ₊₁)
    else
        δ = 0.0
    end
    return ct + δ # type stable
end

# Pseudo-inverse (as a matrix) from an SVD factorization.
# Kept private: extending `LinearAlgebra.pinv(::SVD)` clashes with the method LinearAlgebra itself defines in recent Julia versions.
function _count_svdvals(S, atol::Real, rtol::Real)
    isempty(S) && return 0
    tol = max(rtol * S[1], atol)
    return iszero(S[1]) ? 0 : searchsortedlast(S, tol, rev = true)
end
function _pinv(F::SVD{T}; atol::Real = 0, rtol::Real = (eps(real(float(oneunit(T))))*min(size(F)...)) * iszero(atol)) where {T}
    k = _count_svdvals(F.S, atol, rtol)
    return @views F.Vt[k:-1:1, :]' * Diagonal(inv.(F.S[k:-1:1])) * F.U[:,k:-1:1]'
end

"""
    mad(x)

Compute the Median Absolute Deviation (MAD) scaled by ~`1.4826` for normal consistency.

# Arguments
- `x`: Input vector or array.

# Returns
- Scaled median absolute deviation.
"""
function mad(x)
    m = median(x)
    xc = abs.(x .- m)
    return median!(xc) * 1.4826022185056018
end

"""
    isnoise(x; detection_threshold::Real = 5.0)

Determine whether fluorescence curve `x` represents background noise without true amplification.

Compares the total observed range `max(x) - min(x)` against expected noise level from Gaussian extreme value theory.

# Arguments
- `x`: Vector of cycle fluorescence data.
- `detection_threshold::Real`: Multiplier for noise threshold (default: `5.0`).

# Returns
- `true` if curve is classified as noise (no amplification); `false` if real amplification is present.
"""
function isnoise(x; detection_threshold::Real = 5.0)
    obs = maximum(x) - minimum(x) # observed range
    # expected range:
    # for N iid Gaussian variable the expected range max - min = 2 * σ * sqrt(2 * log(N))
    # for a pure noise amplification curve var(xᵢ - xᵢ₋₁) = 2σ² ⟹ σ = std(xᵢ - xᵢ₋₁) / √2
    # replace σ with the robust estimator ̂σ = mad(xᵢ - xᵢ₋₁)
    expect = 2 * mad(diff(x)) * sqrt(log(length(x)))
    return obs ≤ (detection_threshold * expect)
end
