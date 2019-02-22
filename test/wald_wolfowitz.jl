using HypothesisTests, Test
using Distributions

@testset "WW test" begin

@testset "Independent Observations" begin
    # Get set of independent observations
    x = randn(100)
    tst = WaldWolfowitzTest(x)

    # Should not be significant dependence
    @test pvalue(tst) > 0.05

    # Test consistency of z-statistics
    @test pvalue(tst) == pvalue(Normal(tst.μ, tst.σ), tst.N_runs)
    @test pvalue(tst, tail=:left) == pvalue(Normal(tst.μ, tst.σ), tst.N_runs, tail=:left)
    @test pvalue(tst, tail=:right) == pvalue(Normal(tst.μ, tst.σ), tst.N_runs, tail=:right)
    show(IOBuffer(), tst)
end

@testset "Dependent Observations" begin
    # Get set of independent observations
    x = randn(1)
    for i in 2:100
        append!(x, x[i-1] + randn())
    end
    tst = WaldWolfowitzTest(x)
    
    # Should have significant dependence
    @test pvalue(tst) < 0.05

    # Test consistency of z-statistics
    @test pvalue(tst) == pvalue(Normal(tst.μ, tst.σ), tst.N_runs)
    @test pvalue(tst, tail=:left) == pvalue(Normal(tst.μ, tst.σ), tst.N_runs, tail=:left)
    @test pvalue(tst, tail=:right) == pvalue(Normal(tst.μ, tst.σ), tst.N_runs, tail=:right)
    show(IOBuffer(), tst)
end


end