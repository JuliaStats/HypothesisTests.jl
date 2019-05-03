using HypothesisTests, Test
using Distributions

@testset "Wald-Wolfowitz test" begin

@testset "Many-Valued Observations" begin
    # Get set of dependent observations
    x = 1:1000
    tst = WaldWolfowitzTest(x)
    
    @test tst.z ≈ -31.575 atol=1e-3
    @test pvalue(tst) ≈ 8.052e-219 atol=1e-222

    # Test consistency of z-statistics
    @test pvalue(tst) == pvalue(Normal(tst.μ, tst.σ), tst.nruns)
    @test pvalue(tst, tail=:left) == pvalue(Normal(tst.μ, tst.σ), tst.nruns, tail=:left)
    @test pvalue(tst, tail=:right) == pvalue(Normal(tst.μ, tst.σ), tst.nruns, tail=:right)
    show(IOBuffer(), tst)
end

@testset "Two-Valued Observations" begin
    # equivalent data as above (half under median half over)
    x = [falses(500); trues(500)]
    tst = WaldWolfowitzTest(x)
    
    @test tst.z ≈ -31.575 atol=1e-3
    @test pvalue(tst) ≈ 8.052e-219 atol=1e-222

    # Test consistency of z-statistics
    @test pvalue(tst) == pvalue(Normal(tst.μ, tst.σ), tst.nruns)
    @test pvalue(tst, tail=:left) == pvalue(Normal(tst.μ, tst.σ), tst.nruns, tail=:left)
    @test pvalue(tst, tail=:right) == pvalue(Normal(tst.μ, tst.σ), tst.nruns, tail=:right)
    show(IOBuffer(), tst)
end


end