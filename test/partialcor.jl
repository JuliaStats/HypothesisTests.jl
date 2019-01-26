using HypothesisTests
using Test
using DelimitedFiles
using StatsBase

@testset "Partial correlation" begin
    # Columns are information, similarities, arithmetic, picture completion
    wechsler = readdlm(joinpath(@__DIR__, "data", "wechsler.txt"))[:,2:end]
    w = PartialCorTest(wechsler[:,1], wechsler[:,2], wechsler[:,3:4])
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
    x = PartialCorTest(view(X,:,1), view(X,:,2), view(X,:,3))
    @test occursin("fail to reject", sprint(show, x))
    @test confint(x) == (-1.0, 1.0)
    @test nobs(x) == 4
    @test dof(x) == 1
    @test pvalue(x) ≈ 0.25776212 atol=1e-6
end
