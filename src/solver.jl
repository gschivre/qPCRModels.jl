# Non-linear optimization and root-finding solvers

"""
    LevenbergMarquardt!(nlm::NonLinearLeastSquare; kwargs...)

Fit parameters of a non-linear least squares model `nlm` using the Levenberg-Marquardt optimization algorithm with adaptive damping.

# Arguments
- `nlm::NonLinearLeastSquare`: Target qPCR non-linear model.

# Keyword Arguments
- `LMtol::Real = 1e-8`: Convergence tolerance on gradient norm.
- `LMtol_chisq::Real = 10 * eps()`: Tolerance for relative change in χ² merit function.
- `LMtol_step::Real = 10 * eps()`: Tolerance on step size Δ.
- `Ndone::Int = 10`: Number of consecutive small steps required for termination.
- `LMmaxiter::Int = 1_000`: Maximum allowed iterations.
- `λinit::Real = 1e-2`: Initial damping factor multiplier.
- `λdec::Real = 2`: Factor by which damping λ is divided upon successful step.
- `λinc::Real = 3`: Factor by which damping λ is multiplied upon failed step.
- `min_damping::Real = 10 * eps()`: Lower limit for damping parameter.
- `max_dampening::Real = 0.1 * inv(eps())`: Upper limit for damping parameter.

# Returns
- `(chisq::Float64, converged::Bool)`: Final chi-squared sum of squared residuals and convergence status flag.
"""
function LevenbergMarquardt!(nlm::NonLinearLeastSquare; LMtol = 1e-8, LMtol_chisq = 10 * eps(), LMtol_step = 10 * eps(), Ndone = 10, LMmaxiter = 1_000, λinit = 1e-2, λdec = 2, λinc = 3, min_damping = 10 * eps(), max_dampening = 0.1 * inv(eps()), kwargs...)
    p = nlm.params
    p_tmp = similar(p)
    copyto!(p_tmp, p)
    m = length(p)
    r, J = residualsJacobian(nlm; kwargs...) # residuals & Jacobian
    r_new = similar(r)
    J_new = similar(J)
    λ = max(λinit * maximum(x -> sum(abs2, x), eachcol(J)), min_damping)
    chisq = sum(abs2, r) # χ² merit function value
    g = similar(p) # gradient
    Δ = similar(p)
    δ = Inf
    done = 0
    H = zeros(eltype(p), m, m) # dampened approximate Hessian
    for _ in 1:LMmaxiter
        mul!(H, J', J)
        for i in 1:m
            H[i, i] += min(max(λ * H[i, i], min_damping), max_dampening)
        end
        mul!(g, J', r)
        Δ .= H \ g
        p .-= Δ
        residualsJacobian!(r_new, J_new, nlm; kwargs...) # new residuals & Jacobian
        chisq_new = sum(abs2, r_new)
        δ = chisq_new - chisq
        done += (maximum(abs, Δ) < LMtol_step) || (abs(δ) < max(LMtol_chisq, LMtol_chisq * chisq_new))
        if δ < 0 # success: decrease λ
            λ /= λdec # delayed gratification
            copyto!(p_tmp, p)
            copyto!(r, r_new)
            copyto!(J, J_new)
            chisq = chisq_new
        else # failure: increase λ
            λ *= λinc
            p .= p_tmp
        end
        (maximum(abs, g) < LMtol) && return (chisq, true) # terminate when gradient is small
        (done ≥ Ndone) && return (chisq, true) # terminate when step or χ² improvement are negligible
    end
    @warn "Reached maximum number of iterations" func = "LevenbergMarquardt!"
    return (chisq, false)
end

"""
    newtsafe_rootfinding(f, df, a, b; x₀ = (a + b) / 2, RFtol = 10 * eps(), RFmaxiter = 100, offset = 0.0)

Globally convergent root-finding solver combining Newton-Raphson iterations with bisection bounds.
Code adapted from Press, W. H. (2007), Numerical recipes 3rd edition: The art of scientific computing, p460-461.

Solves `f(x) + offset = 0` within bracket `[a, b]`. Automatically reverts to bisection step if a Newton iteration jumps outside bounds or fails to converge rapidly.

# Arguments
- `f`: Function returning scalar evaluation `f(x)`.
- `df`: Derivative function `f'(x)`.
- `a::Real`, `b::Real`: Search bracket bounds satisfying `(f(a)+offset)*(f(b)+offset) <= 0`.

# Keyword Arguments
- `x₀::Real`: Initial starting guess (defaults to interval midpoint `(a + b)/2`).
- `RFtol::Real`: Desired absolute tolerance on root position.
- `RFmaxiter::Int`: Maximum allowed iterations (default: `100`).
- `offset::Real`: Target offset value (default: `0.0`).

# Returns
- `(x::Float64, converged::Bool)`: Estimated root location and boolean success flag.
"""
function newtsafe_rootfinding(f, df, a, b; x₀ = (a + b) / 2, RFtol = 10 * eps(), RFmaxiter = 100, offset = 0.0)
    fa = f(a) + offset
    fb = f(b) + offset
    (abs(fa) < RFtol) && return (a, true)
    (abs(fb) < RFtol) && return (b, true)
    @assert (fa * fb) ≤ 0 "Root not bracketed"
    # orient the search so that f(xₗ) < 0
    xₗ, xₕ = ifelse(fa > 0, (b, a), (a, b))
    x = x₀
    dx_old = abs(xₗ - xₕ)
    dx = dx_old
    fx = f(x) + offset # use the offset to solve for f(x) = offset
    (abs(fx) < RFtol) && return (x, true)
    dfx = df(x)
    for _ in 1:RFmaxiter
        # fallback to bisection if Newton would go out of bounds or not converge fast enough
        reject_newton = ((((x - xₕ) * dfx - fx) * ((x - xₗ) * dfx - fx)) > 0) || ((abs(2.0 * fx) > abs(dx_old * dfx)))
        dx_old = dx
        dx = (reject_newton ? (0.5 * (xₕ - xₗ)) : (fx / dfx))
        x = ifelse(reject_newton, (xₗ + dx), (x - dx))
        (abs(dx) < RFtol) && return (x, true)
        fx = f(x) + offset
        dfx = df(x)
        # update brackets
        lower_bracket_update = fx < 0
        xₗ = ifelse(lower_bracket_update, x, xₗ)
        xₕ = ifelse(lower_bracket_update, xₕ, x)
    end
    @warn "Reached maximum number of iterations" func = "newtsafe_rootfinding"
    return (x, false)
end
