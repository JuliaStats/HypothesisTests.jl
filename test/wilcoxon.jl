using HypothesisTests, Test
using HypothesisTests: default_tail, population_param_of_interest
using Statistics: median

@testset "Wilcoxon" begin
@testset "Basic exact test" begin
    @test default_tail(ExactSignedRankTest([1:10;], [2:2:20;])) == :both
	show(IOBuffer(), ExactSignedRankTest([1:10;], [2:2:20;]))

    # Two-sided
    for kwargs in ((), (; tail = :both))
        @test abs(@inferred(pvalue(ExactSignedRankTest([1:10;], [2:2:20;]); kwargs...)) - 0.0020) <= 1e-4
        @test abs(@inferred(pvalue(ExactSignedRankTest([2:2:20;], [1:10;]); kwargs...)) - 0.0020) <= 1e-4
        @test abs(@inferred(pvalue(ExactSignedRankTest([1:10;], [2:2:16; -1; 1]); kwargs...)) - 0.4316) <= 1e-4
        @test abs(@inferred(pvalue(ExactSignedRankTest([2:2:16; -1; 1], [1:10;]); kwargs...)) - 0.4316) <= 1e-4
    end

    # Left tail
    @test abs(@inferred(pvalue(ExactSignedRankTest([1:10;], [2:2:20;]); tail = :left)) - 0.0009) <= 1e-4
    @test abs(@inferred(pvalue(ExactSignedRankTest([2:2:20;], [1:10;]); tail = :left)) - 1) <= 1e-4
    @test abs(@inferred(pvalue(ExactSignedRankTest([1:10;], [2:2:16; -1; 1]); tail = :left)) - 0.2158) <= 1e-4
    @test abs(@inferred(pvalue(ExactSignedRankTest([2:2:16; -1; 1], [1:10;]); tail = :left)) - 0.8125) <= 1e-4

    # Right tail
    @test abs(@inferred(pvalue(ExactSignedRankTest([1:10;], [2:2:20;]); tail = :right)) - 1) <= 1e-4
    @test abs(@inferred(pvalue(ExactSignedRankTest([2:2:20;], [1:10;]); tail = :right)) - 0.0009) <= 1e-4
    @test abs(@inferred(pvalue(ExactSignedRankTest([1:10;], [2:2:16; -1; 1]); tail = :right)) - 0.8125) <= 1e-4
    @test abs(@inferred(pvalue(ExactSignedRankTest([2:2:16; -1; 1], [1:10;]); tail = :right)) - 0.2158) <= 1e-4
end

@testset "Exact with ties" begin
    show(IOBuffer(), ExactSignedRankTest([1:10;], [1:10;]))

    # Two-sided
    for kwargs in ((), (; tail = :both))
        @test abs(@inferred(pvalue(ExactSignedRankTest([1:10;], [1:10;]); kwargs...)) - 1) <= 1e-4
        @test abs(@inferred(pvalue(ExactSignedRankTest([1:10;], [2:11;]); kwargs...)) - 0.0020) <= 1e-4
        @test abs(@inferred(pvalue(ExactSignedRankTest([2:11;], [1:10;]); kwargs...)) - 0.0020) <= 1e-4
        @test abs(@inferred(pvalue(ExactSignedRankTest(1:10, [1:5; ones(5)]); kwargs...)) - 0.0625) <= 1e-4
        @test abs(@inferred(pvalue(ExactSignedRankTest([1:5; ones(5)], [1:10;]); kwargs...)) - 0.0625) <= 1e-4
    end

    # Left tail
    @test abs(@inferred(pvalue(ExactSignedRankTest([1:10;], [1:10;]); tail = :left)) - 1) <= 1e-4
    @test abs(@inferred(pvalue(ExactSignedRankTest([1:10;], [2:11;]); tail = :left)) - 0.0009) <= 1e-4
    @test abs(@inferred(pvalue(ExactSignedRankTest([2:11;], [1:10;]); tail = :left)) - 1) <= 1e-4
    @test abs(@inferred(pvalue(ExactSignedRankTest(1:10, [1:5; ones(5)]); tail = :left)) - 1) <= 1e-4
    @test abs(@inferred(pvalue(ExactSignedRankTest([1:5; ones(5)], [1:10;]); tail = :left)) - 0.0312) <= 1e-4

    # Right tail
    @test abs(@inferred(pvalue(ExactSignedRankTest([1:10;], [1:10;]); tail = :right)) - 1) <= 1e-4
    @test abs(@inferred(pvalue(ExactSignedRankTest([1:10;], [2:11;]); tail = :right)) - 1) <= 1e-4
    @test abs(@inferred(pvalue(ExactSignedRankTest([2:11;], [1:10;]); tail = :right)) - 0.0009) <= 1e-4
    @test abs(@inferred(pvalue(ExactSignedRankTest(1:10, [1:5; ones(5)]); tail = :right)) - 0.0312) <= 1e-4
    @test abs(@inferred(pvalue(ExactSignedRankTest([1:5; ones(5)], [1:10;]); tail = :right)) - 1) <= 1e-4
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
    # R: wilcox.test(x, conf.int = TRUE, exact = TRUE)  ->  (3.3, 15.5)
    @test isapprox(@inferred(confint(ExactSignedRankTest(x)))[1], 3.3, atol=1e-4)
    @test isapprox(@inferred(confint(ExactSignedRankTest(x)))[2], 15.5, atol=1e-4)
    # R: wilcox.test(x, conf.int = TRUE, exact = FALSE)  ->  (3.05004, 15.50001).
    # The approximate test gets the approximate interval; it used to report the exact one.
    @test isapprox(@inferred(confint(ApproximateSignedRankTest(x)))[1], 3.05, atol=1e-4)
    @test isapprox(@inferred(confint(ApproximateSignedRankTest(x)))[2], 15.5, atol=1e-4)
    @test isapprox(@inferred(confint(SignedRankTest(x); tail=:left))[1], 4.45, atol=1e-4)
    @test isapprox(@inferred(confint(SignedRankTest(x); tail=:right))[2], 14.45, atol=1e-4)
end

@testset "Hodges-Lehmann estimate" begin
    x = [-7.8, -6.9, -4.7, 3.7, 6.5, 8.7, 9.1, 10.1, 10.8, 13.6, 14.4, 16.6, 20.2, 22.4, 23.5]
    walsh = [(x[i] + x[j]) / 2 for i in eachindex(x) for j in i:length(x)]
    # R reports 9.675 as the estimate for this sample
    @test @inferred(hodgeslehmann(SignedRankTest(x))) ≈ median(walsh) ≈ 9.675
    @test hodgeslehmann(ExactSignedRankTest(x)) == hodgeslehmann(ApproximateSignedRankTest(x))
    # the estimate that is reported is the one the interval is built around
    lo, hi = confint(ExactSignedRankTest(x))
    @test lo <= hodgeslehmann(ExactSignedRankTest(x)) <= hi
    @test population_param_of_interest(ExactSignedRankTest(x))[3] == hodgeslehmann(ExactSignedRankTest(x))

    # ties: R reports -2.1
    h = [1, 2, 3, 4, 5, 6, 7, 10, 10, 10, 10, 10, 13, 14, 15] .- 10.1
    @test hodgeslehmann(SignedRankTest(h)) ≈ -2.1
end

@testset "Zero differences are dropped consistently" begin
    # the statistic already ignores zeros; the interval and the estimate now do too
    d = [0.0, 0, 0, 0.5, 0.5, 1, -0.5, -1, 1.5, -1.5, 0.5, 0, 1, -0.5, 2, 0, 0.5, -1, 1, 0.5]
    nz = d[d .!= 0]
    @test SignedRankTest(d).W == SignedRankTest(nz).W
    @test pvalue(SignedRankTest(d)) == pvalue(SignedRankTest(nz))
    @test confint(SignedRankTest(d); level=0.9) == confint(SignedRankTest(nz); level=0.9)
    @test hodgeslehmann(SignedRankTest(d)) == hodgeslehmann(SignedRankTest(nz))
    # R: wilcox.test(d, conf.int = TRUE, conf.level = 0.9) -> (-0.24999, 0.75005), est 0.5
    @test confint(SignedRankTest(d); level=0.9) == (-0.25, 0.75)
    @test hodgeslehmann(SignedRankTest(d)) == 0.5
    # all-zero input still yields the degenerate interval rather than an error
    @test confint(SignedRankTest(zeros(6))) == (0.0, 0.0)
end

@testset "method keyword" begin
    x = [-7.8, -6.9, -4.7, 3.7, 6.5, 8.7, 9.1, 10.1, 10.8, 13.6, 14.4, 16.6, 20.2, 22.4, 23.5]
    big = collect(1.0:60)                      # n > 50, untied: auto picks approximate
    @test SignedRankTest(x) isa ExactSignedRankTest
    @test SignedRankTest(big) isa ApproximateSignedRankTest
    # :auto reproduces the automatic rule
    @test SignedRankTest(x; method=:auto) isa ExactSignedRankTest
    @test SignedRankTest(big; method=:auto) isa ApproximateSignedRankTest
    # and both may be forced
    @test SignedRankTest(x; method=:approximate) isa ApproximateSignedRankTest
    @test SignedRankTest(big; method=:exact) isa ExactSignedRankTest
    # a callable sees the decision inputs
    @test SignedRankTest(big; method = s -> s.n_nonzero <= 100 ? :exact : :approximate) isa
        ExactSignedRankTest
    @test SignedRankTest(x; method = s -> s.ties ? :exact : :approximate) isa
        ApproximateSignedRankTest
    @test SignedRankTest(x, x .+ 1; method=:approximate) isa ApproximateSignedRankTest

    @test_throws ArgumentError SignedRankTest(x; method=:fastest)
    @test_throws ArgumentError SignedRankTest(x; method = s -> :fastest)
end

@testset "Exact branch refuses where it is not computable" begin
    # StatsFuns.signrankcdf is silently invalid past MAX_EXACT_SIGNRANK_N (StatsFuns#219)
    @test HypothesisTests.MAX_EXACT_SIGNRANK_N == 71
    ok = collect(1.0:71)
    @test pvalue(SignedRankTest(ok; method=:exact)) isa Float64
    @test confint(SignedRankTest(ok; method=:exact)) isa Tuple
    bad = collect(1.0:72)
    @test_throws ArgumentError pvalue(SignedRankTest(bad; method=:exact))
    @test_throws ArgumentError confint(SignedRankTest(bad; method=:exact))
    # ... and the tied branch will not enumerate without bound
    tied = repeat([1.0, 2.0], 20)
    @test_throws ArgumentError pvalue(SignedRankTest(tied; method=:exact))
end
end
