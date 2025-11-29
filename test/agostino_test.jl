using HypothesisTests, Test

@testset "D'Agostino Test" begin
    # Check errors
    @test_throws ArgumentError AgostinoTest(rand(5))     # too few samples
    @test_throws ArgumentError AgostinoTest(rand(50000)) # too many samples

    # Check correctness of results (benchmarked against `moments` package in R)
    y = [
        2.44816, 1.086045, 0.154618, 0.778503, 0.923376,
        1.18671, -0.495009, -0.151178, 1.662281, 2.615421,
    ]

    result = @inferred AgostinoTest(y)
    @test result.z_statistic ≈ 0.28356 atol = 1.0e-6
    @test result.skewness ≈ 0.158581 atol = 1.0e-6
    @test nobs(result) == 10
    @test pvalue(result) ≈ 0.776748 atol = 1.0e-6
    @test pvalue(result, tail = :left) ≈ 0.611626 atol = 1.0e-6
    @test pvalue(result, tail = :right) ≈ 0.388374 atol = 1.0e-6
    @test AgostinoTest(exp.(y)).z_statistic ≈ 2.10965 atol = 1.0e-5
    @test AgostinoTest(-1 .* exp.(y)).z_statistic ≈ -2.10965 atol = 1.0e-5
end
