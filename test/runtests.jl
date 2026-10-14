using Test
using qPCRModels

@testset "qPCRModels.jl Package Tests" begin
    @testset "Dataset Loading" begin
        @test "batsch1" in keys(DatasetDetails)
        ds = qpcRDataset("batsch1")
        @test size(ds, 1) > 0
        @test size(ds, 2) > 1
        @test !isempty(ds.colnames)
    end

    @testset "Utility Functions" begin
        y = [0.1, 0.2, 0.5, 1.2, 3.5, 8.0, 15.0, 20.0, 22.0, 22.5]
        ct = cycletreshold(y, 10.0; linear_interp=true)
        @test !isnan(ct)
        @test ct > 1.0
        
        h = hac(y)
        @test !isnan(h)
        
        s = sdm(y)
        @test s > 0
        
        @test isnoise(rand(10) .* 0.01) == true
    end

    @testset "Model Fitting" begin
        ds = qpcRDataset("batsch1")
        y = ds[:, 2] # first curve

        for ModelType in (BirchModel, LogisticModel, RichardsModel, MAK2Model)
            m = ModelType(y)
            @test !isfitted(m)
            fit!(m)
            @test isfitted(m)
            @test hasconverged(m)
            @test !isnan(get_x0(m))
            @test !isnan(get_logx0(m))
            @test !isnan(AICc(m))
            @test !isnan(reducedchisq(m))
            @test length(get_predicted(m)) > 0
        end
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

    @testset "Missing cycles" begin
        y = qpcRDataset("boggy")[:, 2]
        ytrunc = [y; fill(NaN, 5)] # run stopped early
        for ModelType in (BirchModel, LogisticModel, RichardsModel, MAK2Model)
            m = ModelType(ytrunc)
            @test length(m.rawdata) == length(y)
            fit!(m)
            @test hasconverged(m)
            ref = ModelType(y)
            fit!(ref)
            @test get_logx0(m) ≈ get_logx0(ref) # same fit as the curve without the missing cycles
        end
        ygap = copy(y); ygap[10] = NaN
        @test_throws ArgumentError BirchModel(ygap)
        @test_throws ArgumentError MAK2Model(fill(NaN, 40))
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
end
