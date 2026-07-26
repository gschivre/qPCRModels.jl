#!/usr/bin/env julia
# qpcr_fit.jl - CLI tool for qPCR model fitting.

using Printf
using qPCRModels

# Constants
const VALID_MODELS = Dict(
    "birch" => BirchModel,
    "logistic" => LogisticModel,
    "richards" => RichardsModel,
    "mak2" => MAK2Model
)

function print_help()
    println("""
qPCR Model Fitting CLI Tool
===========================
Usage:
  julia --project=@__DIR__ qpcr_fit.jl <dataset_or_filepath> [options]

Arguments:
  <dataset_or_filepath>   Path to a TSV/CSV file containing qPCR data, OR
                          the name of a built-in dataset (e.g., batsch1, boggy, testdat).

Options:
  -m, --model <name>      Model to fit: "birch" (default), "logistic", "richards", or "mak2"
  -o, --output <path>     Path to save results as a CSV file
  --fit-noise             Force fitting even on curves detected as noise (no amplification)
  -h, --help              Show this help message and exit

Examples:
  julia --project=. qpcr_fit.jl batsch1
  julia --project=. qpcr_fit.jl qpcR_datasets/boggy.txt --model mak2 -o boggy_results.csv
""")
end

function parse_args(args)
    options = Dict{Symbol, Any}(
        :input => nothing,
        :model => "birch",
        :output => nothing,
        :fit_noise => false,
        :help => false
    )
    
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "-h" || arg == "--help"
            options[:help] = true
        elseif arg == "-m" || arg == "--model"
            if i + 1 <= length(args)
                options[:model] = lowercase(args[i+1])
                i += 1
            else
                error("Missing value for option --model")
            end
        elseif arg == "-o" || arg == "--output"
            if i + 1 <= length(args)
                options[:output] = args[i+1]
                i += 1
            else
                error("Missing value for option --output")
            end
        elseif arg == "--fit-noise"
            options[:fit_noise] = true
        elseif startswith(arg, "-")
            error("Unknown option: $arg")
        else
            if isnothing(options[:input])
                options[:input] = arg
            else
                error("Multiple input arguments specified")
            end
        end
        i += 1
    end
    return options
end

function load_custom_file(filepath)
    if !isfile(filepath)
        error("File not found: $filepath")
    end
    
    lines = filter(!isempty, strip.(readlines(filepath)))
    if isempty(lines)
        error("File is empty: $filepath")
    end
    
    # Detect delimiter
    first_line = lines[1]
    delim = ','
    if occursin('\t', first_line)
        delim = '\t'
    elseif occursin(';', first_line)
        delim = ';'
    elseif occursin(' ', first_line) && !occursin(',', first_line) && !occursin('\t', first_line)
        delim = ' '
    end
    
    # Parse header
    header = split(first_line, delim; keepempty = false)
    colnames = String.(strip.(header))
    
    # Parse numeric data
    num_rows = length(lines) - 1
    num_cols = length(colnames)
    data = zeros(Float64, num_rows, num_cols)
    
    for (i, line) in enumerate(lines[2:end])
        tokens = split(line, delim; keepempty = false)
        if length(tokens) < num_cols
            error("Row $(i+1) has fewer columns than the header ($(length(tokens)) vs $num_cols)")
        end
        for j in 1:num_cols
            data[i, j] = parse(Float64, tokens[j])
        end
    end
    
    return colnames, data
end

function main()
    options = try
        parse_args(ARGS)
    catch e
        println(stderr, "Error: ", e.msg)
        println(stderr)
        print_help()
        exit(1)
    end
    
    if options[:help] || isnothing(options[:input])
        print_help()
        exit(0)
    end
    
    input = options[:input]
    model_name = options[:model]
    output_path = options[:output]
    fit_noise = options[:fit_noise]
    
    if !(model_name in keys(VALID_MODELS))
        println(stderr, "Error: Invalid model '$model_name'. Allowed: ", join(keys(VALID_MODELS), ", "))
        exit(1)
    end
    
    model_type = VALID_MODELS[model_name]
    
    # Load dataset
    colnames = String[]
    data = Matrix{Float64}(undef, 0, 0)
    dataset_details = ""
    
    # Check if built-in
    if input in keys(DatasetDetails)
        try
            ds = qpcRDataset(input)
            colnames = ds.colnames
            data = ds.data
            setup = ds.details.Setup
            details = ds.details.Details
            dataset_details = "Built-in Dataset: $input\nSetup: $setup\nDetails: $details\n"
        catch e
            println(stderr, "Error loading built-in dataset: ", e)
            exit(1)
        end
    else
        try
            colnames, data = load_custom_file(input)
        catch e
            println(stderr, "Error loading custom file '$input': ", e)
            exit(1)
        end
    end
    
    num_rows, num_cols = size(data)
    if num_cols < 1
        println(stderr, "Error: No columns found in data.")
        exit(1)
    end
    
    start_col = 2
    if num_cols == 1
        start_col = 1
        curves = ["Curve1"]
    else
        curves = colnames[start_col:end]
    end
    
    println("="^88)
    if !isempty(dataset_details)
        print(dataset_details)
    else
        println("Input File: $input")
    end
    println("Model: ", uppercase(model_name))
    println("Cycles detected: $num_rows")
    println("Curves to fit: ", length(curves))
    println("="^88)
    
    @printf("%-18s %-25s %-12s %-10s %-10s %-12s\n", "Curve", "Status", "x0", "log(x0)", "AICc", "Red.Chisq")
    println("-"^88)
    
    results = []
    
    for (i, curve_name) in enumerate(curves)
        col_idx = start_col - 1 + i
        y = data[:, col_idx]
        
        m = model_type(y)
        fit!(m; fit_noise = fit_noise)
        
        status = "Failed"
        x0_val = "N/A"
        log_x0_val = "N/A"
        aicc_val = "N/A"
        red_chisq_val = "N/A"
        
        if isfitted(m)
            if hasconverged(m)
                status = "Converged"
                x0_val = get_x0(m)
                log_x0_val = get_logx0(m)
                aicc_val = AICc(m)
                red_chisq_val = reducedchisq(m)
            else
                status = "Failed to converge"
            end
        else
            if !hasamplification(m)
                status = "No amplification (Noise)"
            else
                status = "Not fitted"
            end
        end
        
        push!(results, (
            curve = curve_name,
            status = status,
            x0 = x0_val,
            log_x0 = log_x0_val,
            aicc = aicc_val,
            red_chisq = red_chisq_val
        ))
        
        if status == "Converged"
            @printf("%-18s %-25s %-12.4e %-10.2f %-10.2f %-12.4e\n", curve_name, status, x0_val, log_x0_val, aicc_val, red_chisq_val)
        else
            @printf("%-18s %-25s %-12s %-10s %-10s %-12s\n", curve_name, status, "N/A", "N/A", "N/A", "N/A")
        end
    end
    println("-"^88)
    
    if !isnothing(output_path)
        try
            open(output_path, "w") do f
                write(f, "Curve,Status,x0,log_x0,AICc,Reduced_Chisq\n")
                for r in results
                    x0_str = r.x0 isa Number ? @sprintf("%.8e", r.x0) : "N/A"
                    log_x0_str = r.log_x0 isa Number ? @sprintf("%.6f", r.log_x0) : "N/A"
                    aicc_str = r.aicc isa Number ? @sprintf("%.6f", r.aicc) : "N/A"
                    red_chisq_str = r.red_chisq isa Number ? @sprintf("%.8e", r.red_chisq) : "N/A"
                    
                    write(f, "$(r.curve),$(r.status),$(x0_str),$(log_x0_str),$(aicc_str),$(red_chisq_str)\n")
                end
            end
            println("Results successfully saved to: $output_path")
        catch e
            println(stderr, "Error saving results to '$output_path': ", e)
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
