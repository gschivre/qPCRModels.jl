# Pretty printing of fitted models and calibration curves

# Decimals shown: two significant digits of the standard error, or four significant digits of the value without one
function _decimals(v, se)
    (isfinite(se) && se > 0) && return clamp(1 - floor(Int, log10(se)), 0, 8)
    (isfinite(v) && !iszero(v)) && return clamp(3 - floor(Int, log10(abs(v))), 0, 8)
    return 0
end
_fmt(x, d) = isfinite(x) ? Printf.format(Printf.Format("%.$(d)f"), x) : string(x)

# Padding by display width, so that labels with subscripts such as "log₁₀ x₀" line up
_rpad(s, w) = s * " "^max(w - textwidth(s), 0)

# Print rows (label, value, standard error or `nothing`) as a table: values start in the same column,
# with room for a minus sign so that their first digits line up
function _print_table(io::IO, rows)
    vals = [_fmt(v, _decimals(v, something(se, NaN))) for (_, v, se) in rows]
    sign = any(startswith('-'), vals) ? " " : ""
    wl = maximum(row -> textwidth(row[1]), rows)
    for ((label, v, se), val) in zip(rows, vals)
        print(io, "\n  ", _rpad(label, wl), "  ", startswith(val, '-') ? "" : sign, val)
        isnothing(se) || print(io, " ± ", _fmt(se, _decimals(v, se)))
    end
    return nothing
end

_sig(x) = @sprintf("%.4g", x)

function _status(nlm::NonLinearLeastSquare)
    isfitted(nlm) || return (nlm.noamp ? "no amplification" : "not fitted", :light_black)
    return hasconverged(nlm) ? ("converged", :green) : ("not converged", :yellow)
end

_cycles(nlm::NonLinearLeastSquare) = "$(length(nlm.data)) cycles"
_cycles(nlm::MAK2Model) = "cycles 1–$(nlm.n) of $(length(nlm.data))"

function Base.show(io::IO, ::MIME"text/plain", nlm::NonLinearLeastSquare)
    status, color = _status(nlm)
    printstyled(io, nameof(typeof(nlm)); bold = true)
    print(io, " · ", _cycles(nlm), " · ")
    printstyled(io, status; color = color)
    isfitted(nlm) || return nothing
    se = get_stderrors(nlm)
    p = get_params(nlm)
    rows = Tuple{String, Float64, Union{Nothing, Float64}}[(String(k), p[i], se[k]) for (i, k) in enumerate(_param_names(nlm))]
    push!(rows, ("F_bg", nlm.σ * nlm.linearparams[1], se.F_bg), ("F_max", nlm.σ * nlm.linearparams[2], se.F_max),
          ("log₁₀ x₀", get_logx0(nlm) / log(10), se.logx0 / log(10)))
    nlm isa MAK2Model || push!(rows, ("Cy₀", get_cy0(nlm), nothing))
    _print_table(io, rows)
    # residual standard error in fluorescence units, more telling than χ²ν on the standardised data
    print(io, "\n  RSE = ", _sig(nlm.σ * sqrt(reducedchisq(nlm))), " · AICc = ", @sprintf("%.1f", AICc(nlm)))
    return nothing
end

function Base.show(io::IO, nlm::NonLinearLeastSquare)
    print(io, nameof(typeof(nlm)), "(")
    if isfitted(nlm)
        p = get_params(nlm)
        for (i, k) in enumerate(_param_names(nlm))
            print(io, k, " = ", _sig(p[i]), ", ")
        end
        nlm isa MAK2Model && print(io, "log₁₀ x₀ = ", _sig(get_logx0(nlm) / log(10)), ", n = ", nlm.n, ", ")
    end
    print(io, first(_status(nlm)), ")")
    return nothing
end

function Base.show(io::IO, ::MIME"text/plain", cal::CalibrationCurve)
    nused = count(cal.used)
    nout = count(iszero, cal.weights[cal.used])
    printstyled(io, "CalibrationCurve"; bold = true)
    print(io, " · ", cal.robust ? "robust" : "least squares", cal.fixed_slope ? ", slope fixed to 1" : "", " · ")
    print(io, nused, nused == length(cal.used) ? "" : " of $(length(cal.used))", " points")
    cal.robust && print(io, " (", nout, nout == 1 ? " outlier)" : " outliers)")
    cal.converged || printstyled(io, " · not converged"; color = :yellow)
    print(io, "\n  log₁₀ x₀ = intercept + slope · log₁₀ q")
    rows = Tuple{String, Float64, Union{Nothing, Float64}}[]
    cal.fixed_slope || push!(rows, ("slope", cal.slope, cal.slope_se))
    push!(rows, ("intercept", cal.intercept, cal.intercept_se), ("R²", cal.r2, nothing))
    _print_table(io, rows)
    return nothing
end

Base.show(io::IO, cal::CalibrationCurve) = print(io, "CalibrationCurve(slope = ", _sig(cal.slope), ", intercept = ", _sig(cal.intercept), ")")
