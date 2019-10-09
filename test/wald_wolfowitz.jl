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
    expected_output = """
    Wald-Wolfowitz Test
    -------------------
    Population details:
        parameter of interest:   Number of runs
        value under h_0:         501.0
        point estimate:          2

    Test summary:
        outcome with 95% confidence: reject h_0
        two-sided p-value:           <1e-99

    Details:
        number of runs:  2
        z-statistic:     -31.575338477995764
    """
    output = sprint(show, tst)
    @test output == expected_output
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
    expected_output = """
    Wald-Wolfowitz Test
    -------------------
    Population details:
        parameter of interest:   Number of runs
        value under h_0:         501.0
        point estimate:          2

    Test summary:
        outcome with 95% confidence: reject h_0
        two-sided p-value:           <1e-99

    Details:
        number of runs:  2
        z-statistic:     -31.575338477995764
    """
    output = sprint(show, tst)
    @test output == expected_output
end

end
