#!/usr/bin/env julia
# bin/qpcr_fit.jl - CLI tool for qPCR model fitting.
#
# The script is also the entry point of the standalone executable built by juliac/build.jl,
# so it must stay trimmable: every value has a concrete type, output goes through Core.stdout
# and Core.stderr, and caught exceptions are never printed as such.

using Printf
using qPCRModels

# Constants
const MODEL_NAMES = ("birch", "logistic", "richards", "mak2")
const LINE = "========================================================================================"
const DASHES = "----------------------------------------------------------------------------------------"

const HELP = """
qPCR Model Fitting CLI Tool
===========================
Usage:
  qpcr_fit <dataset_or_filepath> [options]                         (standalone executable)
  julia --project=. bin/qpcr_fit.jl <dataset_or_filepath> [options]

Arguments:
  <dataset_or_filepath>   Path to a TSV/CSV file containing qPCR data, OR
                          the name of a built-in dataset (e.g., batsch1, boggy, testdat).

Options:
  -m, --model <name>      Model to fit: "birch" (default), "logistic", "richards", or "mak2"
  -o, --output <path>     Path to save results as a CSV file
  -s, --standards <path>  File of standard curves with known quantities (e.g. copies), one per line:
                          <curve name> <quantity>, same delimiters as the data, optional header.
                          A standard curve is fitted on them and the quantity of every curve is reported.
  --fixed-slope           Fix the slope of the standard curve to 1 (x0 proportional to the quantity)
  --no-robust             Fit the standard curve by ordinary least squares (no outlier rejection)
  --fit-noise             Force fitting even on curves detected as noise (no amplification)
  -h, --help              Show this help message and exit

Examples:
  qpcr_fit batsch1
  qpcr_fit qpcR_datasets/boggy.txt --model mak2 -o boggy_results.csv
  qpcr_fit plate.csv --standards standards.csv -o plate_results.csv
"""

struct Options
    input::String # empty if not given
    model::String
    output::String # empty if not given
    standards::String # empty if not given
    fixed_slope::Bool
    robust::Bool
    fit_noise::Bool
    help::Bool
end

struct CurveResult
    curve::String
    status::String
    x0::Float64 # NaN unless the fit converged
    log_x0::Float64
    aicc::Float64
    red_chisq::Float64
end

# Bad arguments or data. Its message is a `String`: the `AbstractString` message of
# `ErrorException` or `ArgumentError` cannot be printed in a trimmed executable.
struct CLIError <: Exception
    msg::String
end

out(s::String) = (print(Core.stdout, s); nothing)
err(s::String) = (print(Core.stderr, s); nothing)

# Message of a caught exception; anything else than bad input or a failed file operation is a bug and is rethrown
errmsg(e) = e isa CLIError ? e.msg : e isa SystemError ? e.prefix : rethrow()

function parse_args(args::Vector{String})
    input, model, output, standards, fixed_slope, robust, fit_noise, help = "", "birch", "", "", false, true, false, false
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "-h" || arg == "--help"
            help = true
        elseif arg == "-m" || arg == "--model"
            i + 1 <= length(args) || throw(CLIError("Missing value for option --model"))
            model = lowercase(args[i + 1])
            i += 1
        elseif arg == "-o" || arg == "--output"
            i + 1 <= length(args) || throw(CLIError("Missing value for option --output"))
            output = args[i + 1]
            i += 1
        elseif arg == "-s" || arg == "--standards"
            i + 1 <= length(args) || throw(CLIError("Missing value for option --standards"))
            standards = args[i + 1]
            i += 1
        elseif arg == "--fixed-slope"
            fixed_slope = true
        elseif arg == "--no-robust"
            robust = false
        elseif arg == "--fit-noise"
            fit_noise = true
        elseif startswith(arg, "-")
            throw(CLIError("Unknown option: " * arg))
        else
            isempty(input) || throw(CLIError("Multiple input arguments specified"))
            input = arg
        end
        i += 1
    end
    return Options(input, model, output, standards, fixed_slope, robust, fit_noise, help)
end

# Non-empty lines of a text file
function read_lines(filepath::String)
    isfile(filepath) || throw(CLIError("File not found: " * abspath(filepath)))
    lines = filter(!isempty, strip.(split(read(filepath, String), '\n')))
    isempty(lines) && throw(CLIError("File is empty: " * filepath))
    return lines
end

# Delimiter guessed from the first line: tab, semicolon, space or comma
function detect_delim(first_line::AbstractString)
    occursin('\t', first_line) && return '\t'
    occursin(';', first_line) && return ';'
    (occursin(' ', first_line) && !occursin(',', first_line)) && return ' '
    return ','
end

function load_custom_file(filepath::String)
    lines = read_lines(filepath)
    first_line = lines[1]
    delim = detect_delim(first_line)

    colnames = String[String(strip(h)) for h in split(first_line, delim; keepempty = false)]
    data = zeros(Float64, length(lines) - 1, length(colnames))
    for i in 2:length(lines)
        tokens = split(lines[i], delim; keepempty = false)
        if length(tokens) < length(colnames)
            throw(CLIError("Row " * string(i) * " has fewer columns than the header (" * string(length(tokens)) * " vs " * string(length(colnames)) * ")"))
        end
        for j in eachindex(colnames)
            v = tryparse(Float64, tokens[j])
            isnothing(v) && throw(CLIError("Row " * string(i) * ", column " * string(j) * ": '" * String(tokens[j]) * "' is not a number"))
            data[i - 1, j] = v
        end
    end

    return colnames, data
end

# Standards file: one `<curve name> <quantity>` line per standard, with an optional header line
function load_standards(filepath::String)
    lines = read_lines(filepath)
    delim = detect_delim(lines[1])
    names, quantities = String[], Float64[]
    for (i, line) in enumerate(lines)
        tokens = split(line, delim; keepempty = false)
        length(tokens) == 2 || throw(CLIError("Standards row " * string(i) * ": expected a curve name and a quantity, got " * string(length(tokens)) * " fields"))
        name = String(strip(tokens[1]))
        q = tryparse(Float64, tokens[2])
        if isnothing(q)
            i == 1 && continue # header
            throw(CLIError("Standards row " * string(i) * ": '" * String(tokens[2]) * "' is not a number"))
        end
        (isfinite(q) && q > 0) || throw(CLIError("Standards row " * string(i) * ": the quantity of '" * name * "' must be positive"))
        name in names && throw(CLIError("Standards row " * string(i) * ": '" * name * "' is listed twice"))
        push!(names, name)
        push!(quantities, q)
    end
    isempty(names) && throw(CLIError("No standard found in " * filepath))
    return names, quantities
end

# Standard curve from the fitted standards; the checks of `calibration` are done first, so that the error
# is a `CLIError` (the message of an `ArgumentError` cannot be printed in a trimmed executable)
function fit_calibration(results::Vector{CurveResult}, curves::Vector{String}, names::Vector{String}, quantities::Vector{Float64}, fixed_slope::Bool, robust::Bool)
    index = Dict{String, Int}(c => i for (i, c) in enumerate(curves))
    x0 = Float64[]
    for name in names
        haskey(index, name) || throw(CLIError("Standard '" * name * "' is not a curve of the input data"))
        push!(x0, results[index[name]].x0)
    end
    usable = [isfinite(x) && x > 0 for x in x0]
    needed = fixed_slope ? 2 : 3
    count(usable) >= needed || throw(CLIError("The standard curve needs at least " * string(needed) * " converged standards, got " * string(count(usable))))
    (fixed_slope || length(unique(quantities[usable])) > 1) || throw(CLIError("The standard curve needs converged standards with at least two different quantities"))
    return calibration(quantities, x0; fixed_slope = fixed_slope, robust = robust)
end

function calibration_summary(cal::CalibrationCurve, names::Vector{String})
    s = "Standard curve: log10(x0) = intercept + slope * log10(quantity), " * (cal.robust ? "robust" : "least-squares") * " fit on " * string(count(cal.used)) * " of " * string(length(names)) * " standards\n"
    s *= cal.fixed_slope ? "  slope     = 1 (fixed)\n" : @sprintf("  slope     = %.4f ± %.4f\n", cal.slope, cal.slope_se)
    s *= @sprintf("  intercept = %.4f ± %.4f\n  R²        = %.4f\n", cal.intercept, cal.intercept_se, cal.r2)
    outliers = [names[i] for i in eachindex(names) if cal.used[i] && iszero(cal.weights[i])]
    isempty(outliers) || (s *= "  outliers (not used): " * join(outliers, ", ") * "\n")
    cal.converged || (s *= "  warning: the robust fit did not converge\n")
    return s
end

function fit_curve(m::NonLinearLeastSquare, name::String, fit_noise::Bool)
    fit!(m; fit_noise = fit_noise)
    if isfitted(m)
        hasconverged(m) || return CurveResult(name, "Failed to converge", NaN, NaN, NaN, NaN)
        return CurveResult(name, "Converged", get_x0(m), get_logx0(m), AICc(m), reducedchisq(m))
    end
    status = hasamplification(m) ? "Not fitted" : "No amplification (Noise)"
    return CurveResult(name, status, NaN, NaN, NaN, NaN)
end

# A branch rather than a Dict of model types, so that each call has a concrete type
function fit_curve(model::String, y::Vector{Float64}, name::String, fit_noise::Bool)
    model == "birch" && return fit_curve(BirchModel(y), name, fit_noise)
    model == "logistic" && return fit_curve(LogisticModel(y), name, fit_noise)
    model == "richards" && return fit_curve(RichardsModel(y), name, fit_noise)
    return fit_curve(MAK2Model(y), name, fit_noise)
end

function _csv(x::Float64, fmt::Symbol)
    isnan(x) && return "N/A"
    fmt === :e && return @sprintf("%.8e", x)
    fmt === :g && return @sprintf("%.4g", x)
    return @sprintf("%.6f", x)
end

_row(r::CurveResult, quantity::String) =
    r.status == "Converged" ? @sprintf("%-18s %-25s %-12.4e %-10.2f %-10.2f %-12.4e%s\n", r.curve, r.status, r.x0, r.log_x0, r.aicc, r.red_chisq, quantity) :
    @sprintf("%-18s %-25s %-12s %-10s %-10s %-12s%s\n", r.curve, r.status, "N/A", "N/A", "N/A", "N/A", quantity)

function main(args::Vector{String})
    opts = try
        parse_args(args)
    catch e
        err("Error: " * errmsg(e) * "\n\n" * HELP)
        return 1
    end

    if opts.help || isempty(opts.input)
        out(HELP)
        return 0
    end

    if !(opts.model in MODEL_NAMES)
        err("Error: Invalid model '" * opts.model * "'. Allowed: " * join(MODEL_NAMES, ", ") * "\n")
        return 1
    end

    if haskey(DatasetDetails, opts.input)
        ds = try
            qpcRDataset(opts.input)
        catch e
            err("Error loading built-in dataset '" * opts.input * "': " * errmsg(e) * "\n")
            return 1
        end
        colnames, data = ds.colnames, ds.data
        header = "Built-in Dataset: " * opts.input * "\nSetup: " * ds.details.Setup * "\nDetails: " * ds.details.Details * "\n"
    else
        colnames, data = try
            load_custom_file(opts.input)
        catch e
            err("Error loading custom file '" * opts.input * "': " * errmsg(e) * "\n")
            return 1
        end
        header = "Input File: " * opts.input * "\n"
    end

    num_rows, num_cols = size(data)
    if num_cols < 1
        err("Error: No columns found in data.\n")
        return 1
    end
    start_col = num_cols == 1 ? 1 : 2
    curves = num_cols == 1 ? ["Curve1"] : colnames[2:end]

    if isempty(opts.standards)
        std_names, std_quantities = String[], Float64[]
    else
        std_names, std_quantities = try
            load_standards(opts.standards)
        catch e
            err("Error loading standards file '" * opts.standards * "': " * errmsg(e) * "\n")
            return 1
        end
    end
    has_cal = !isempty(std_names)

    out(LINE * "\n" * header)
    out("Model: " * uppercase(opts.model) * "\nCycles detected: " * string(num_rows) * "\nCurves to fit: " * string(length(curves)) * "\n")
    has_cal && out("Standards: " * string(length(std_names)) * " (" * opts.standards * ")\n")
    out(LINE * "\n")

    results = CurveResult[fit_curve(opts.model, data[:, start_col - 1 + i], name, opts.fit_noise) for (i, name) in enumerate(curves)]

    quantity = fill(NaN, length(results))
    if has_cal
        cal = try
            fit_calibration(results, curves, std_names, std_quantities, opts.fixed_slope, opts.robust)
        catch e
            err("Error: " * errmsg(e) * "\n")
            return 1
        end
        quantity .= [quantify(cal, r.x0) for r in results]
        out(calibration_summary(cal, std_names) * LINE * "\n")
    end

    out(@sprintf("%-18s %-25s %-12s %-10s %-10s %-12s%s\n", "Curve", "Status", "x0", "log(x0)", "AICc", "Red.Chisq", has_cal ? " Quantity" : ""))
    out(DASHES * "\n")
    for (r, q) in zip(results, quantity)
        out(_row(r, has_cal ? " " * _csv(q, :g) : ""))
    end
    out(DASHES * "\n")

    if !isempty(opts.output)
        known = Dict{String, Float64}(zip(std_names, std_quantities))
        csv = IOBuffer()
        print(csv, "Curve,Status,x0,log_x0,AICc,Reduced_Chisq", has_cal ? ",Known_quantity,Quantity\n" : "\n")
        for (r, q) in zip(results, quantity)
            print(csv, r.curve, ",", r.status, ",", _csv(r.x0, :e), ",", _csv(r.log_x0, :f), ",", _csv(r.aicc, :f), ",", _csv(r.red_chisq, :e))
            has_cal && print(csv, ",", _csv(get(known, r.curve, NaN), :e), ",", _csv(q, :e))
            print(csv, "\n")
        end
        try
            write(opts.output, take!(csv))
        catch e
            err("Error saving results to '" * opts.output * "': " * errmsg(e) * "\n")
            return 1
        end
        out("Results successfully saved to: " * opts.output * "\n")
    end
    return 0
end

# Entry point: `@main` (Julia ≥ 1.11, needed by juliac); on older versions, run when the script is the program
@static if isdefined(Base, Symbol("@main"))
    @main
else
    abspath(PROGRAM_FILE) == (@__FILE__) && exit(main(ARGS))
end
