using HypothesisTests
using Test
using DelimitedFiles
using StatsBase

@testset "Correlation" begin
    # Columns are line number, calcium, iron
    nutrient = readdlm(joinpath(@__DIR__, "data", "nutrient.txt"))[:, 1:3]
    w = CorrelationTest(nutrient[:,2], nutrient[:,3])
    let out = sprint(show, w)
        @test occursin("reject h_0", out) && !occursin("fail to", out)
    end
    let ci = confint(w)
        @test first(ci) ≈ 0.3327138 atol=1e-6
        @test last(ci) ≈ 0.4546639 atol=1e-6
    end
    @test nobs(w) == 737
    @test dof(w) == 735
    @test pvalue(w) < 1e-25

    x = CorrelationTest(nutrient[:,1], nutrient[:,2])
    @test occursin("fail to reject", sprint(show, x))
    let ci = confint(x)
        @test first(ci) ≈ -0.1105478 atol=1e-6
        @test last(ci) ≈ 0.0336730 atol=1e-6
    end
    @test pvalue(x) ≈ 0.2948405 atol=1e-6
end

@testset "Partial correlation" begin
    # Columns are information, similarities, arithmetic, picture completion
    wechsler = readdlm(joinpath(@__DIR__, "data", "wechsler.txt"))[:,2:end]
    w = CorrelationTest(wechsler[:,1], wechsler[:,2], wechsler[:,3:4])
    let out = sprint(show, w)
        @test occursin("reject h_0", out) && !occursin("fail to", out)
    end
    let ci = confint(w)
        @test first(ci) ≈ 0.4963917 atol=1e-6
        @test last(ci) ≈ 0.8447292 atol=1e-6
    end
    @test nobs(w) == 37
    @test dof(w) == 33
    @test pvalue(w) < 0.00001

    X = [ 2 1 0
          4 2 0
         15 3 1
         20 4 1]
    x = CorrelationTest(view(X,:,1), view(X,:,2), view(X,:,3))
    @test occursin("fail to reject", sprint(show, x))
    @test confint(x) == (-1.0, 1.0)
    @test nobs(x) == 4
    @test dof(x) == 1
    @test pvalue(x) ≈ 0.25776212 atol=1e-6
end
