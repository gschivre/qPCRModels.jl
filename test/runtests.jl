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
end
