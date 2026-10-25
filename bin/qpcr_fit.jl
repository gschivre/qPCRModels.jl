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
  --fit-noise             Force fitting even on curves detected as noise (no amplification)
  -h, --help              Show this help message and exit

Examples:
  qpcr_fit batsch1
  qpcr_fit qpcR_datasets/boggy.txt --model mak2 -o boggy_results.csv
"""

struct Options
    input::String # empty if not given
    model::String
    output::String # empty if not given
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
    input, model, output, fit_noise, help = "", "birch", "", false, false
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "-h" || arg == "--help"
            help = true
        elseif arg == "-m" || arg == "--model"
            i + 1 <= length(args) || throw(CLIError("Missing value for option --model"))
            model = lowercase(args[i+1])
            i += 1
        elseif arg == "-o" || arg == "--output"
            i + 1 <= length(args) || throw(CLIError("Missing value for option --output"))
            output = args[i+1]
            i += 1
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
    return Options(input, model, output, fit_noise, help)
end

function load_custom_file(filepath::String)
    isfile(filepath) || throw(CLIError("File not found: " * abspath(filepath)))

    lines = filter(!isempty, strip.(split(read(filepath, String), '\n')))
    isempty(lines) && throw(CLIError("File is empty: " * filepath))

    first_line = lines[1]
    delim = ','
    if occursin('\t', first_line)
        delim = '\t'
    elseif occursin(';', first_line)
        delim = ';'
    elseif occursin(' ', first_line) && !occursin(',', first_line)
        delim = ' '
    end

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
            data[i-1, j] = v
        end
    end

    return colnames, data
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

_csv(x::Float64, fmt::Symbol) = isnan(x) ? "N/A" : fmt === :e ? @sprintf("%.8e", x) : @sprintf("%.6f", x)

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

    out(LINE * "\n" * header)
    out("Model: " * uppercase(opts.model) * "\nCycles detected: " * string(num_rows) * "\nCurves to fit: " * string(length(curves)) * "\n")
    out(LINE * "\n")
    out(@sprintf("%-18s %-25s %-12s %-10s %-10s %-12s\n", "Curve", "Status", "x0", "log(x0)", "AICc", "Red.Chisq"))
    out(DASHES * "\n")

    results = CurveResult[]
    for (i, curve_name) in enumerate(curves)
        r = fit_curve(opts.model, data[:, start_col-1+i], curve_name, opts.fit_noise)
        push!(results, r)
        if r.status == "Converged"
            out(@sprintf("%-18s %-25s %-12.4e %-10.2f %-10.2f %-12.4e\n", r.curve, r.status, r.x0, r.log_x0, r.aicc, r.red_chisq))
        else
            out(@sprintf("%-18s %-25s %-12s %-10s %-10s %-12s\n", r.curve, r.status, "N/A", "N/A", "N/A", "N/A"))
        end
    end
    out(DASHES * "\n")

    if !isempty(opts.output)
        csv = IOBuffer()
        print(csv, "Curve,Status,x0,log_x0,AICc,Reduced_Chisq\n")
        for r in results
            print(csv, r.curve, ",", r.status, ",", _csv(r.x0, :e), ",", _csv(r.log_x0, :f), ",", _csv(r.aicc, :f), ",", _csv(r.red_chisq, :e), "\n")
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
