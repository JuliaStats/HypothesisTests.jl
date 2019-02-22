using .HypothesisTests, Test
using Distributions

@testset "Wald-Wolfowitz test" begin

@testset "Independent Observations" begin
    # Get set of independent observations
    x = randn(1000)
    tst = WaldWolfowitzTest(x)

    # Should not be significant dependence
    @test pvalue(tst) > 0.05

    # Test consistency of z-statistics
    @test pvalue(tst) == pvalue(Normal(tst.μ, tst.σ), tst.nruns)
    @test pvalue(tst, tail=:left) == pvalue(Normal(tst.μ, tst.σ), tst.nruns, tail=:left)
    @test pvalue(tst, tail=:right) == pvalue(Normal(tst.μ, tst.σ), tst.nruns, tail=:right)
    show(IOBuffer(), tst)
end

@testset "Dependent Observations" begin
    # Get set of dependent observations
    x = 1:1000
    tst = WaldWolfowitzTest(x)
    
    # Should have significant dependence
    @test pvalue(tst) < 0.05

    # Test consistency of z-statistics
    @test pvalue(tst) == pvalue(Normal(tst.μ, tst.σ), tst.nruns)
    @test pvalue(tst, tail=:left) == pvalue(Normal(tst.μ, tst.σ), tst.nruns, tail=:left)
    @test pvalue(tst, tail=:right) == pvalue(Normal(tst.μ, tst.σ), tst.nruns, tail=:right)
    show(IOBuffer(), tst)
end


end