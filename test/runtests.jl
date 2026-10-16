using Test
using qPCRModels
using qPCRModels: RichardsLogistic
using LogExpFunctions: logistic, log1mexp
using Random
using Statistics

# Exact Birch curve: bisection on the implicit solution ln x - (1 + b) ln(1 - x) = ln(2) (c - tₘ + b).
# Independent of the package root finder, so it can serve as ground truth.
function birch_x(c, b, tₘ)
    g(x) = log(x) - (1 + b) * log1p(-x) - log(2) * (c - tₘ + b)
    lo, hi = 0.0, 1.0
    for _ in 1:2000
        mid = (lo + hi) / 2
        (mid == lo || mid == hi) && break
        g(mid) < 0 ? (lo = mid) : (hi = mid)
    end
    return (lo + hi) / 2
end

# Cy₀ from its geometric definition: tangent at the steepest point of a normalized curve meets x = 0
function cy0_numeric(xfun)
    cs = range(-20, 80; length = 200_001)
    xs = xfun.(cs)
    s = diff(xs) ./ step(cs)
    k = argmax(s)
    return (cs[k] + cs[k + 1]) / 2 - (xs[k] + xs[k + 1]) / 2 / s[k]
end

col(name, c) = (d = qpcRDataset(name); d[:, findfirst(==(c), d.colnames)])

const ALL_MODELS = (BirchModel, LogisticModel, RichardsModel, MAK2Model)

@testset "qPCRModels.jl Package Tests" begin
    @testset "Datasets" begin
        ds = qpcRDataset("batsch1")
        @test size(ds, 2) == length(ds.colnames) == 16 # cycles + 5 dilutions × 3 replicates
        @test ds.colnames[1] == "Cycles"
        @test ds[:, 1] == 1:size(ds, 1)
        @test ds.details == DatasetDetails["batsch1"]

        files = Set(first(splitext(f)) for f in readdir(joinpath(pkgdir(qPCRModels), "qpcR_datasets")) if endswith(f, ".txt"))
        @test files == Set(keys(DatasetDetails))

        all32 = load_qpcR_dataset(Float32)
        @test Set(keys(all32)) == Set(keys(DatasetDetails))
        @test eltype(all32["boggy"]) == Float32
        @test all32["boggy"] ≈ qpcRDataset("boggy")

        @test_throws Exception qpcRDataset("not_a_dataset")
    end

    @testset "Utility functions" begin
        y = [0.0, 1.0, 2.0, 3.0, 4.0]
        @test cycletreshold(y, 2.5) == 4
        @test cycletreshold(y, 2.5; linear_interp = true) ≈ 3.5
        @test cycletreshold(y, 0.0; linear_interp = true) == 1 # reached at the first cycle
        @test isnan(cycletreshold(y, 10.0)) # never reached

        @test hac(collect(0.0:10.0); linear_interp = true) ≈ 6.0 # 50% of the range is reached at cycle 6

        # second-derivative maximum of a logistic curve is at tₘ - ln(2 + √3) / a (smoothing bias < 0.3 cycle)
        for (a, tₘ) in ((0.5, 22.0), (0.7, 25.3), (0.4, 30.0))
            @test sdm(logistic.(a .* ((1:45) .- tₘ)); parabolic_interp = true) ≈ tₘ - log(2 + sqrt(3)) / a atol = 0.3
        end
        @test sdm(logistic.(0.5 .* ((1:45) .- 22.0))) == 19 # whole cycle without interpolation

        @test mad([1.0, 2.0, 3.0, 4.0, 100.0]) ≈ 1.4826022185056018 # robust to the outlier

        flat = 1.0 .+ 0.01 .* sin.(7.3 .* (1:40)) # deterministic baseline wobble
        @test isnoise(flat)
        @test !isnoise(qpcRDataset("batsch1")[:, 2])
        @test isnoise(qpcRDataset("batsch1")[:, 2]; detection_threshold = 1e6)
    end

    @testset "Solvers" begin
        x, _ = newtsafe_rootfinding(x -> x^2 - 2, x -> 2x, 0.0, 2.0)
        @test x ≈ sqrt(2)
        x, _ = newtsafe_rootfinding(x -> x^2, x -> 2x, 0.0, 2.0; offset = -3.0) # solves x² = 3
        @test x ≈ sqrt(3)
        @test_throws AssertionError newtsafe_rootfinding(x -> x^2 + 1, x -> 2x, 0.0, 2.0)
        # root within one ulp of a large bracket end (MAK2 trial step): used to spin until RFmaxiter, as eps(z) > RFtol
        z₀ = 38.617540143048885
        z, converged = @test_logs newtsafe_rootfinding(z -> qPCRModels.expinti(z) - qPCRModels.expinti(z₀), z -> exp(z) / z, z₀, 42.8; x₀ = z₀, offset = -5 * log(2))
        @test converged && (z ≈ z₀)
        x, converged = @test_logs newtsafe_rootfinding(x -> x^2 - 2, x -> 2x, 0.0, 2.0; RFtol = 0.0) # to machine precision
        @test converged && (x ≈ sqrt(2))

        m = BirchModel(qpcRDataset("batsch1")[:, 2])
        chisq_lm, converged = @test_logs (:warn, "Reached maximum number of iterations") match_mode = :any LevenbergMarquardt!(m; LMmaxiter = 1)
        @test !converged
        @test isfinite(chisq_lm)
    end

    @testset "Model fitting on real data" begin
        y = qpcRDataset("batsch1")[:, 2]
        for ModelType in ALL_MODELS
            m = ModelType(y)
            @test !isfitted(m)
            @test hasamplification(m)
            fit!(m)
            @test isfitted(m)
            @test hasconverged(m)
            @test isfinite(get_x0(m))
            @test get_x0(m) ≈ exp(get_logx0(m))
            @test isfinite(AICc(m))
            @test isfinite(reducedchisq(m))
            ModelType === MAK2Model || @test isfinite(get_cy0(m)) # no Cy₀ for MAK2
            @test length(get_predicted(m)) == size(m.xmat, 1)
        end
    end

    @testset "Parameter recovery on exact curves" begin
        # Birch: logistic (b = 0), asymmetric (b > 0), sharp (b < 0) and late curves (tₘ ≳ 48 used to collapse to a flat line)
        for (b, tₘ, N) in ((0.0, 22.0, 45), (1.5, 24.0, 45), (4.0, 20.0, 45), (-0.6, 26.0, 45), (1.5, 52.0, 65), (-0.9, 52.0, 80))
            m = BirchModel(2.0 .+ 10.0 .* birch_x.(1:N, b, tₘ))
            fit!(m)
            @test hasconverged(m)
            @test collect(get_params(m)) ≈ [b, tₘ] atol = 1e-8
            @test get_logx0(m) ≈ log(birch_x(0, b, tₘ)) atol = 1e-8
            @test get_linearparams(m) .* m.σ ≈ [2.0, 10.0] # baseline and scale
        end

        a, tₘ = 0.55, 21.0
        m = LogisticModel(2.0 .+ 10.0 .* logistic.(a .* ((1:45) .- tₘ)))
        fit!(m)
        @test collect(get_params(m)) ≈ [a, tₘ] atol = 1e-8
        @test get_x0(m) ≈ logistic(-a * tₘ)

        a, tₘ, ν = 0.6, 23.0, 0.4
        m = RichardsModel(2.0 .+ 10.0 .* RichardsLogistic.(a .* ((1:45) .- tₘ), ν))
        fit!(m)
        @test collect(get_params(m)) ≈ [a, tₘ, ν] atol = 1e-5
    end

    @testset "Birch with b = 0 is the logistic with a = ln 2" begin
        tₘ = 23.0
        y = 2.0 .+ 10.0 .* logistic.(log(2) .* ((1:45) .- tₘ))
        mb, ml = BirchModel(y), LogisticModel(y)
        fit!(mb)
        fit!(ml)
        @test get_params(mb)[1] ≈ 0 atol = 1e-8
        @test get_logx0(mb) ≈ get_logx0(ml) atol = 1e-8
        @test get_cy0(mb) ≈ get_cy0(ml) atol = 1e-8
        @test get_predicted(mb) ≈ get_predicted(ml)
    end

    @testset "Birch x₀ across parameter space" begin
        # covers every branch of the starting-value heuristic; x₀ must satisfy the implicit solution at c = 0
        y = qpcRDataset("batsch1")[:, 2]
        for b in (-0.9, -0.6, 0.0, 2.0, 60.0, 200.0), tₘ in (-3.0, 0.5, 20.0, 40.0, 50.0)
            m = BirchModel(y)
            m.params .= [log1p(b), tₘ]
            z₀ = get_logx0(m)
            @test z₀ - (1 + b) * log1mexp(z₀) ≈ log(2) * (b - tₘ) atol = 1e-8
            @test get_x0(m) ≈ exp(z₀) # used to lose accuracy for late curves (absolute tolerance in x)
        end
    end

    @testset "Birch saturating curves freeze on the plateau" begin
        # b → -1 (no re-annealing) reaches the plateau abruptly; values used to turn NaN in the VarPro matrix
        y = qpcRDataset("batsch1")[:, 2]
        for b in (-0.999, -1 + 1e-12), tₘ in (5.0, 20.0)
            m = BirchModel(y)
            m.params .= [log1p(b), tₘ]
            r, J = residualsJacobian(m)
            x = m.xmat[:, 2]
            @test all(isfinite, r) && all(isfinite, J)
            @test issorted(x)
            @test x[end] ≈ 1
        end
    end

    @testset "Cy₀ matches the tangent-line definition" begin
        y = qpcRDataset("batsch1")[:, 2]
        m = BirchModel(y)
        m.params .= [log1p(2.0), 24.0]
        @test get_cy0(m) ≈ cy0_numeric(c -> birch_x(c, 2.0, 24.0)) atol = 1e-6
        m = LogisticModel(y)
        m.params .= [log(0.5), 22.0]
        @test get_cy0(m) ≈ cy0_numeric(c -> logistic(0.5 * (c - 22.0))) atol = 1e-6
        m = RichardsModel(y)
        m.params .= [log(0.6), 23.0, log(0.4)]
        @test get_cy0(m) ≈ cy0_numeric(c -> RichardsLogistic(0.6 * (c - 23.0), 0.4)) atol = 1e-6
    end

    @testset "Residuals, χ² and information criteria" begin
        y = qpcRDataset("batsch1")[:, 2]
        for ModelType in ALL_MODELS
            m = ModelType(y)
            r0 = residuals(m; scaled = false) # before fitting: initial parameters
            @test get_linearparams(m) == m.linearparams # computed on demand before fitting
            @test length(r0) == size(m.xmat, 1)
            @test isfinite(chisq(m))
            ModelType in (BirchModel, MAK2Model) && @test isfinite(chisq(m; RFtol = 1e-12)) # keywords used to throw
            fit!(m)
            n, k = size(m.xmat, 1), length(m.params) + length(m.linearparams)
            r = residuals(m)
            @test sum(abs2, r) ≈ chisq(m) rtol = 1e-6
            @test residuals(m; scaled = false) ≈ y[1:n] .- get_predicted(m; scaled = false)
            @test residuals(m; scaled = false) ≈ m.σ .* r
            @test get_predicted(m; scaled = false) ≈ m.σ .* get_predicted(m)
            @test reducedchisq(m) ≈ chisq(m) / (n - k)
            @test AICc(m) ≈ n * log(chisq(m) / n) + 2k + n * (1 + log(2π)) + 2 * (k^2 + k) / (n - k - 1)
        end
    end

    @testset "Noise curves are not fitted unless asked" begin
        flat = 1.0 .+ 0.01 .* sin.(7.3 .* (1:40))
        for ModelType in (BirchModel, LogisticModel, RichardsModel)
            m = ModelType(flat)
            @test !hasamplification(m)
            fit!(m)
            @test !isfitted(m)
            fit!(m; fit_noise = true)
            @test isfitted(m)
        end
    end

    @testset "Parameters outside Float64 range" begin
        y = qpcRDataset("batsch1")[:, 2]
        m = BirchModel(y)
        m.params .= [1e6, 20.0] # b = expm1(1e6) overflows
        @test all(isnan, get_linearparams!(m))
        r, _ = residualsJacobian(m)
        @test all(isinf, r) # LevenbergMarquardt! rejects such a step

        m = MAK2Model(y)
        m.params .= [1e4] # x₀ = e¹⁰⁰⁰⁰
        @test all(isnan, get_linearparams!(m))
        r, _ = residualsJacobian(m)
        @test all(isinf, r)
    end

    @testset "MAK2 x₀ proportional to template on recursion data" begin
        # noise-free curves from the MAK2 recursion xₙ = xₙ₋₁ + ln(1 + xₙ₋₁) for a 4-fold dilution series
        x0s = 1e-3 ./ 4.0 .^ (0:4)
        logx0 = map(x0s) do x0
            m = MAK2Model(1.0 .+ 10.0 .* accumulate((x, _) -> x + log1p(x), 1:40; init = x0))
            fit!(m)
            get_logx0(m)
        end
        slope = ([ones(5) log.(x0s)] \ logx0)[2]
        @test slope ≈ 1 atol = 0.01
    end

    @testset "MAK2 fitting window" begin
        y = qpcRDataset("boggy")[:, 2]
        m = MAK2Model(y)
        @test_throws AssertionError update_n!(m, 4) # 3 parameters and AICc need n ≥ 5
        update_n!(m, 12)
        @test m.n == 12
        @test size(m.xmat, 1) == 12
        @test !isfitted(m)
    end

    @testset "MAK2 refit! over all windows" begin
        # windows that are pure baseline or reach the plateau used to push x₀ out of Float64 range
        y = qpcRDataset("boggy")[:, 2]
        m = MAK2Model(y)
        aicc = refit!(m)
        @test isfitted(m)
        @test any(isfinite, aicc)
        @test isfinite(get_logx0(m))
    end

    @testset "Missing cycles" begin
        y = qpcRDataset("boggy")[:, 2]
        ytrunc = [y; fill(NaN, 5)] # run stopped early
        for ModelType in ALL_MODELS
            m = ModelType(ytrunc)
            @test length(m.rawdata) == length(y)
            fit!(m)
            @test hasconverged(m)
            ref = ModelType(y)
            fit!(ref)
            @test get_logx0(m) ≈ get_logx0(ref) # same fit as the curve without the missing cycles
        end
        ygap = copy(y)
        ygap[10] = NaN
        @test_throws ArgumentError BirchModel(ygap)
        @test_throws ArgumentError MAK2Model(fill(NaN, 40))
    end

    @testset "Standard errors" begin
        y = qpcRDataset("dil4reps94")[:, 2]
        for ModelType in ALL_MODELS
            m = ModelType(y)
            @test all(isnan, get_stderrors(m)) # not fitted
            fit!(m)
            se = get_stderrors(m)
            @test keys(se) == (qPCRModels._param_names(m)..., :F_bg, :F_max, :logx0)
            @test all(s -> isfinite(s) && s > 0, se)
            # analytic gradient of ln x₀ against central differences
            p = copy(m.params)
            g = map(eachindex(p)) do k
                h = 1e-6 * max(1, abs(p[k]))
                m.params .= p; m.params[k] += h; up = get_logx0(m)
                m.params .= p; m.params[k] -= h; dn = get_logx0(m)
                (up - dn) / 2h
            end
            m.params .= p
            @test qPCRModels._logx0_gradient(m) ≈ g rtol = 1e-6
        end
        # standard errors match the spread of refits on curves simulated around the fit
        m = BirchModel(y)
        fit!(m)
        ŷ, s = get_predicted(m; scaled = false), m.σ * sqrt(reducedchisq(m))
        rng = MersenneTwister(1)
        sims = map(1:300) do _
            mm = BirchModel(ŷ .+ s .* randn(rng, length(ŷ)))
            fit!(mm)
            (get_params(mm)..., (mm.σ .* mm.linearparams)..., get_logx0(mm))
        end
        se = get_stderrors(m)
        for (k, name) in enumerate(keys(se))
            @test std(getindex.(sims, k)) ≈ se[name] rtol = 0.2
        end
    end

    @testset "Calibration curve" begin
        q = repeat(10.0 .^ (1:5); inner = 3)
        y = -6 .+ log10.(q) .+ 0.01 .* sin.(1:15) # log₁₀ x₀ with small deterministic scatter
        y[4] += 1.5 # outlier
        cal = calibration(q, exp10.(y))
        @test cal.slope ≈ 1 atol = 0.01
        @test cal.intercept ≈ -6 atol = 0.03
        @test cal.weights[4] == 0 && all(>(0.5), cal.weights[[1:3; 5:15]])
        @test cal.slope_se > 0 && cal.intercept_se > 0 && cal.converged
        ols = calibration(q, exp10.(y); robust = false)
        @test abs(ols.slope - 1) > abs(cal.slope - 1) # pulled by the outlier
        @test all(==(1), ols.weights)
        fixed = calibration(q, exp10.(y); fixed_slope = true)
        @test fixed.slope == 1 && fixed.slope_se == 0
        @test fixed.intercept ≈ -6 atol = 0.03
        @test quantify(cal, exp10(cal.intercept + cal.slope * log10(250.0))) ≈ 250
        # perfect line: zero MAD must not give NaN weights
        exact = calibration(q, exp10.(-6 .+ log10.(q)))
        @test exact.slope ≈ 1 && all(==(1), exact.weights)

        x0 = exp10.(y)
        x0[2] = NaN
        x0[3] = 0.0
        @test calibration(q, x0).used == [i ∉ (2, 3) for i in 1:15]
        @test_throws DimensionMismatch calibration(q[1:3], x0)
        @test_throws ArgumentError calibration(-q, x0)
        @test_throws ArgumentError calibration(fill(10.0, 15), x0) # a single quantity
        @test_throws ArgumentError calibration(q[1:2], x0[1:2])

        # from fitted models: unconverged models are left out
        ds = qpcRDataset("dil4reps94")
        cols = [findfirst(==("$(g)_$(r)"), ds.colnames) for g in ("F15", "F150", "F1500", "F15000") for r in 1:5]
        models = [BirchModel(ds[:, j]) for j in cols]
        foreach(fit!, models[2:end])
        rel = [1 / parse(Int, split(ds.colnames[j], "_")[1][2:end]) for j in cols] # F15 holds the most template
        cal = calibration(rel, models)
        @test !cal.used[1] && all(cal.used[2:end])
        @test cal.slope ≈ 1 atol = 0.05
        @test quantify(cal, models[end]) ≈ rel[end] rtol = 0.3
    end

    @testset "Display" begin
        y = qpcRDataset("dil4reps94")[:, 2]
        plain(x) = sprint(show, MIME"text/plain"(), x)
        for ModelType in ALL_MODELS
            m = ModelType(y)
            @test occursin("not fitted", plain(m))
            @test occursin("not fitted", sprint(show, m))
            fit!(m)
            txt = plain(m)
            @test startswith(txt, string(nameof(ModelType), " · ")) && occursin("converged", txt)
            @test occursin("log₁₀ x₀", txt) && occursin("RSE = ", txt)
            # the first digits of the values line up
            lines = filter(l -> occursin(" ± ", l), split(txt, "\n"))
            @test length(unique(textwidth(l[1:prevind(l, findfirst(isdigit, l[3:end]) + 2)]) for l in lines)) == 1
            @test sprint(show, m) == string(nameof(ModelType), "(", sprint(show, m)[length(string(nameof(ModelType))) + 2:end])
            @test endswith(sprint(show, m), "converged)")
        end
        @test occursin("no amplification", plain(BirchModel(1.0 .+ 0.01 .* sin.(7.3 .* (1:40)))))
        q = repeat(10.0 .^ (1:5); inner = 3)
        cal = calibration(q, exp10.(-6 .+ log10.(q) .+ 0.01 .* sin.(1:15)))
        txt = plain(cal)
        @test occursin("robust", txt) && occursin("15 points (0 outliers)", txt) && occursin("slope", txt)
        @test !any(startswith("  slope"), split(plain(calibration(q, exp10.(-6 .+ log10.(q)); fixed_slope = true)), "\n"))
        @test startswith(sprint(show, cal), "CalibrationCurve(slope = ")
    end

    @testset "Drifting curves that pass isnoise do not crash" begin
        # smooth baseline drift without amplification: LM used to step out of Float64 range (Birch)
        # or hit an exactly singular damped Hessian (Richards)
        for (name, c, ModelType) in (("htPCR", "P48.W12", BirchModel), ("vermeulen1", "CPSG3.1495", BirchModel),
                                     ("htPCR", "P67.W91", RichardsModel))
            m = ModelType(col(name, c))
            fit!(m)
            @test isfitted(m)
        end
    end

    @testset "Command-line interface" begin
        script = joinpath(pkgdir(qPCRModels), "bin", "qpcr_fit.jl")
        cli = Module(:CLI)
        Base.include(cli, script) # main() only runs when the script is the program
        @test cli.parse_args(["boggy"])[:input] == "boggy"
        opts = cli.parse_args(["data.csv", "-m", "MAK2", "--output", "out.csv", "--fit-noise"])
        @test opts[:model] == "mak2" && opts[:output] == "out.csv" && opts[:fit_noise]
        @test cli.parse_args(["-h"])[:help]
        @test_throws ErrorException cli.parse_args(["a", "b"])
        @test_throws ErrorException cli.parse_args(["a", "--bogus"])
        @test_throws ErrorException cli.parse_args(["a", "-m"])
        @test_throws ErrorException cli.parse_args(["a", "-o"])

        dir = mktempdir()
        amp = 1.0 .+ 10.0 .* logistic.(0.6 .* ((1:40) .- 22.0)) .+ 0.01 .* sin.(3.1 .* (1:40))
        flat = 1.0 .+ 0.01 .* sin.(7.3 .* (1:40))
        for (delim, ext) in (('\t', "tsv"), (',', "csv"), (';', "txt"), (' ', "dat"))
            path = joinpath(dir, "curves.$ext")
            write(path, join(["Cycle", "A", "B"], delim) * "\n" * join((join((c, amp[c], flat[c]), delim) for c in 1:40), "\n") * "\n\n")
            colnames, data = cli.load_custom_file(path)
            @test colnames == ["Cycle", "A", "B"]
            @test size(data) == (40, 3) && data[:, 2] ≈ amp # blank trailing line ignored
        end
        short = joinpath(dir, "short.csv")
        write(short, "Cycle,A\n1,2.0\n2\n")
        @test_throws ErrorException cli.load_custom_file(short)
        @test_throws ErrorException cli.load_custom_file(joinpath(dir, "missing.csv"))
        empty = joinpath(dir, "empty.csv")
        write(empty, "\n")
        @test_throws ErrorException cli.load_custom_file(empty)

        # end-to-end runs in a separate process, with the test environment
        function run_cli(args...)
            out, err = IOBuffer(), IOBuffer()
            cmd = `$(Base.julia_cmd()) --startup-file=no --project=$(Base.active_project()) $script $args`
            p = run(pipeline(ignorestatus(cmd); stdout = out, stderr = err))
            return (code = p.exitcode, out = String(take!(out)), err = String(take!(err)))
        end
        res = run_cli("--help")
        @test res.code == 0 && occursin("Usage:", res.out)
        res = run_cli("boggy", "--model", "nope")
        @test res.code == 1 && occursin("Invalid model 'nope'", res.err)
        res = run_cli("boggy", "--bogus")
        @test res.code == 1 && occursin("Unknown option: --bogus", res.err)
        res = run_cli(joinpath(dir, "missing.csv"))
        @test res.code == 1 && occursin("File not found", res.err)

        csv = joinpath(dir, "boggy.csv")
        res = run_cli("boggy", "-o", csv)
        @test res.code == 0 && occursin("Built-in Dataset: boggy", res.out) && occursin("Results successfully saved", res.out)
        rows = readlines(csv)
        @test rows[1] == "Curve,Status,x0,log_x0,AICc,Reduced_Chisq"
        @test length(rows) == size(qpcRDataset("boggy"), 2) # header + one row per curve
        fields = split(rows[2], ',')
        @test fields[1] == "F1.1" && fields[2] == "Converged"
        ref = BirchModel(qpcRDataset("boggy")[:, 2])
        fit!(ref)
        @test parse(Float64, fields[4]) ≈ get_logx0(ref) atol = 1e-5 # same fit as the package

        csv = joinpath(dir, "custom.csv")
        res = run_cli(joinpath(dir, "curves.tsv"), "-m", "logistic", "-o", csv)
        @test res.code == 0 && occursin("Model: LOGISTIC", res.out)
        @test startswith.(readlines(csv)[2:3], ("A,Converged,", "B,No amplification (Noise),N/A")) == [true, true]
        res = run_cli(joinpath(dir, "curves.tsv"), "-m", "logistic", "--fit-noise", "-o", csv)
        @test !occursin("No amplification", read(csv, String))
    end
end
