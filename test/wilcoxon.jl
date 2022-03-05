using HypothesisTests, Test
using HypothesisTests: default_tail

@testset "Wilcoxon" begin
@testset "Basic exact test" begin
    @test abs(@inferred(pvalue(ExactSignedRankTest([1:10;], [2:2:20;]))) - 0.0020) <= 1e-4
    @test abs(@inferred(pvalue(ExactSignedRankTest([2:2:20;], [1:10;]))) - 0.0020) <= 1e-4
    @test abs(@inferred(pvalue(ExactSignedRankTest([1:10;], [2:2:16; -1; 1]))) - 0.4316) <= 1e-4
    @test abs(@inferred(pvalue(ExactSignedRankTest([2:2:16; -1; 1], [1:10;]))) - 0.4316) <= 1e-4
	@test default_tail(ExactSignedRankTest([1:10;], [2:2:20;])) == :both
	show(IOBuffer(), ExactSignedRankTest([1:10;], [2:2:20;]))
end

@testset "Exact with ties" begin
    @test abs(@inferred(pvalue(ExactSignedRankTest([1:10;], [1:10;]))) - 1) <= 1e-4
    @test abs(@inferred(pvalue(ExactSignedRankTest([1:10;], [2:11;]))) - 0.0020) <= 1e-4
    @test abs(@inferred(pvalue(ExactSignedRankTest([2:11;], [1:10;]))) - 0.0020) <= 1e-4
    @test abs(@inferred(pvalue(ExactSignedRankTest(1:10, [1:5; ones(5)]))) - 0.0625) <= 1e-4
    @test abs(@inferred(pvalue(ExactSignedRankTest([1:5; ones(5)], [1:10;]))) - 0.0625) <= 1e-4
	show(IOBuffer(), ExactSignedRankTest([1:10;], [1:10;]))
end

@testset "Approximate test" begin
    @test abs(@inferred(pvalue(ApproximateSignedRankTest([1:10;], [2:2:20;]))) - 0.005922) <= 1e-6
    @test abs(@inferred(pvalue(ApproximateSignedRankTest([2:2:20;], [1:10;]))) - 0.005922) <= 1e-6
    @test abs(@inferred(pvalue(ApproximateSignedRankTest([1:10;], [2:2:16; -1; 1]))) - 0.4148) <= 1e-4
    @test abs(@inferred(pvalue(ApproximateSignedRankTest([2:2:16; -1; 1], [1:10;]))) - 0.4148) <= 1e-4
	@test default_tail(ApproximateSignedRankTest([1:10;], [2:2:20;])) == :both
	show(IOBuffer(), ApproximateSignedRankTest([1:10;], [2:2:20;]))
end

@testset "Approximate with ties" begin
    @test abs(@inferred(pvalue(ApproximateSignedRankTest([1:10;], [1:10;]))) - 1) <= 1e-4
    @test abs(@inferred(pvalue(ApproximateSignedRankTest([1:10;], [2:11;]))) - 0.001904) <= 1e-6
    @test abs(@inferred(pvalue(ApproximateSignedRankTest([2:11;], [1:10;]))) - 0.001904) <= 1e-6
    @test abs(@inferred(pvalue(ApproximateSignedRankTest([1:10;], [1:5; ones(5)]))) - 0.05906) <= 1e-5
    @test abs(@inferred(pvalue(ApproximateSignedRankTest([1:5; ones(5)], 1:10))) - 0.05906) <= 1e-5
	show(IOBuffer(), ApproximateSignedRankTest([1:10;], [1:10;]))
end

@testset "Tests for automatic selection" begin
    @test abs(@inferred(pvalue(SignedRankTest([1:10;], [2:2:20;]))) - 0.0020) <= 1e-4
    @test abs(@inferred(pvalue(SignedRankTest([1:10;], [2:11;]))) - 0.0020) <= 1e-4
	@test default_tail(SignedRankTest([1:10;], [2:2:20;])) == :both
	show(IOBuffer(), SignedRankTest([1:10;], [2:2:20;]))
end

@testset "One Sample tests" begin
	# P-value computed using R wilcox.test
    @test abs(@inferred(pvalue(SignedRankTest([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15] .- 10.1))) - 0.09460449) <= 1e-4
	# P-value computed using R wilcox.test
    @test abs(@inferred(pvalue(SignedRankTest([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16] .- 10.1))) - 0.1928101) <= 1e-4
end

@testset "One Sample tests with ties" begin
	# P-value computed using R package exactRankTests wilcox.exact
    @test abs(@inferred(pvalue(SignedRankTest([1,2,3,4,5,6,7,10,10,10,10,10,13,14,15] .- 10.1))) - 0.04052734) <= 1e-4
	# P-value computed using R wilcox.test
    @test abs(@inferred(pvalue(SignedRankTest([1,2,3,4,5,6,7,10,10,10,10,10,13,14,15,16] .- 10.1))) - 0.1021964) <= 1e-4
end

@testset "Issue 128" begin
    @test @inferred(pvalue(SignedRankTest([54.5, 54.5, 95.0, 51.5]), tail=:left))  == 1
    @test @inferred(pvalue(SignedRankTest([54.5, 54.5, 95.0, 51.5]), tail=:right)) == 0.0625
end

@testset "Test confidence interval" begin
    x = [-7.8, -6.9, -4.7, 3.7, 6.5, 8.7, 9.1, 10.1, 10.8, 13.6, 14.4, 16.6, 20.2, 22.4, 23.5]
    @test isapprox(@inferred(confint(ExactSignedRankTest(x)))[1], 3.3, atol=1e-4)
    @test isapprox(@inferred(confint(ExactSignedRankTest(x)))[2], 15.5, atol=1e-4)
    @test isapprox(@inferred(confint(ApproximateSignedRankTest(x)))[1], 3.3, atol=1e-4)
    @test isapprox(@inferred(confint(ApproximateSignedRankTest(x)))[2], 15.5, atol=1e-4)
    @test isapprox(@inferred(confint(SignedRankTest(x); tail=:left))[1], 4.45, atol=1e-4)
    @test isapprox(@inferred(confint(SignedRankTest(x); tail=:right))[2], 14.45, atol=1e-4)
end
end
