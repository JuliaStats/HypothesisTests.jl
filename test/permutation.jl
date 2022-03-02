using HypothesisTests, Test, Statistics, StableRNGs

@testset "Permutation" begin
    @testset "ExactPermutationTest" begin
        @test pvalue(ExactPermutationTest([1,2,3], [4,5,6], mean), tail=:both) ≈ 0.1
        @test pvalue(ExactPermutationTest([4,5,6], [1,2,3], mean), tail=:both) ≈ 0.1
        @test pvalue(ExactPermutationTest([1,2,3], [4,5,6], mean), tail=:left) ≈ 0.05
        @test pvalue(ExactPermutationTest([4,5,6], [1,2,3], mean), tail=:right) ≈ 0.05
    end

    @testset "ApproximatePermutationTest" begin
        rng = StableRNG(12345)
        @test pvalue(ApproximatePermutationTest(rng, [1,2,3], [4,5,6], mean, 100), tail=:both) ≈ 0.14
        @test pvalue(ApproximatePermutationTest(rng, [4,5,6], [1,2,3], mean, 100), tail=:both) ≈ 0.08
        @test pvalue(ApproximatePermutationTest(rng, [1,2,3], [4,5,6], mean, 100), tail=:left) ≈ 0.05
        @test pvalue(ApproximatePermutationTest(rng, [4,5,6], [1,2,3], mean, 100), tail=:right) ≈ 0.05
        # Test that the non-rng method works
        @test ApproximatePermutationTest([4,5,6], [1,2,3], mean, 100) isa HypothesisTests.PermutationTest{Float64}
    end
end
