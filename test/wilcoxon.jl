using HypothesisTests, Test
using HypothesisTests: default_tail, population_param_of_interest
using Statistics: median, quantile
using Distributions: Normal

# A parameterised route rule written the natural way, as a struct you can store, reuse
# and print. It is callable, and it is not a subtype of `Function`: defining a call
# method does not place a type under `Function` unless it says so.
struct SignedRankExactBelow
    n::Int
end
(r::SignedRankExactBelow)(s) = s.n_nonzero <= r.n ? :exact : :approximate

# A vector whose indices do not start at 1, enough to check that the pairwise-estimate
# loop refuses one rather than reading out of bounds. Avoids a dependency on OffsetArrays
# for a single assertion.
struct ZeroBasedVector{T} <: AbstractVector{T}
    data::Vector{T}
end
Base.size(v::ZeroBasedVector) = size(v.data)
Base.axes(v::ZeroBasedVector) = (Base.IdentityUnitRange(0:(length(v.data) - 1)),)
Base.getindex(v::ZeroBasedVector, i::Int) = v.data[i + 1]

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

@testset "Two-sided p-values are probabilities" begin
    # Doubling a discrete tail overshoots at `W == n(n+1)/4`, which is attainable
    # whenever `n(n+1)/2` is even. Sweep every sign pattern of every untied sample
    # up to n = 12; 220 of the 8190 used to come back above 1.
    signed_ranks(n, bits) =
        [(bits >> (i - 1)) & 1 == 1 ? float(i) : -float(i) for i in 1:n]
    overshoots = [(n, bits) for n in 1:12 for bits in 0:(2^n - 1)
                  if pvalue(ExactSignedRankTest(signed_ranks(n, bits))) > 1]
    @test isempty(overshoots)
    # the worst of them: W = 3 = 3*4/4, both tails 0.625
    @test pvalue(ExactSignedRankTest([1.0, 2.0, -3.0])) == 1
    # and unchanged away from the null mean
    @test pvalue(ExactSignedRankTest([1.0, 2.0, 3.0])) ≈ 0.25
end

@testset "Exact with ties" begin
    # every difference is zero here, so the interval shown cannot have 95% coverage and
    # `show` says so
    @test_logs (:warn, r"not attainable") show(IOBuffer(), ExactSignedRankTest([1:10;], [1:10;]))

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
    # 2/2^10 exactly, from the tied enumeration; the approximate route gives
    # 0.0019042, which an atol of 1e-4 against 0.0020 would also have accepted
    @test @inferred(pvalue(SignedRankTest([1:10;], [2:11;]))) == 0.001953125
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
    # The approximate test now gets the approximate interval; it used to report the exact one.
    @test isapprox(@inferred(confint(ApproximateSignedRankTest(x)))[1], 3.05, atol=1e-4)
    @test isapprox(@inferred(confint(ApproximateSignedRankTest(x)))[2], 15.5, atol=1e-4)
    # one-sided, at level 0.95, is the corresponding endpoint of the two-sided 0.90
    # interval, and the tail names the alternative: :left is location below the null, so it
    # keeps the upper bound. R: wilcox.test(x, alternative = "less", conf.int = TRUE) ->
    # (-Inf, 14.45), and "greater" -> (4.45, Inf). See #368.
    @test @inferred(confint(SignedRankTest(x); tail=:left))[1] == -Inf
    @test isapprox(@inferred(confint(SignedRankTest(x); tail=:left))[2], 14.45, atol=1e-4)
    @test isapprox(@inferred(confint(SignedRankTest(x); tail=:right))[1], 4.45, atol=1e-4)
    @test @inferred(confint(SignedRankTest(x); tail=:right))[2] == Inf
    # and the pair matches what pvalue means by the same tail: the :left bound is the
    # acceptance limit of the left-tailed test
    u = confint(SignedRankTest(x); tail=:left)[2]
    @test pvalue(SignedRankTest(x .- (u - 0.02)); tail=:left) > 0.05
    @test pvalue(SignedRankTest(x .- (u + 0.02)); tail=:left) < 0.05
    l = confint(SignedRankTest(x); tail=:right)[1]
    @test pvalue(SignedRankTest(x .- (l + 0.02)); tail=:right) > 0.05
    @test pvalue(SignedRankTest(x .- (l - 0.02)); tail=:right) < 0.05

    # one-sided p-values on the same sample, against R:
    # wilcox.test(x, alternative = "less")    -> 0.9986877441
    # wilcox.test(x, alternative = "greater") -> 0.001678466797
    @test isapprox(pvalue(ExactSignedRankTest(x); tail=:left), 0.9986877441, atol=1e-9)
    @test isapprox(pvalue(ExactSignedRankTest(x); tail=:right), 0.001678466797, atol=1e-9)
    # and the approximate route, R with exact = FALSE, correct = TRUE:
    # "less" -> 0.9975337639, "greater" -> 0.002938062707
    @test isapprox(pvalue(ApproximateSignedRankTest(x); tail=:left), 0.9975337639, atol=1e-9)
    @test isapprox(pvalue(ApproximateSignedRankTest(x); tail=:right), 0.002938062707, atol=1e-9)
    # approximate one-sided intervals: R solves numerically for (-Inf, 14.45001612) and
    # (4.450022698, Inf); these are the order statistics it lands beside
    @test confint(ApproximateSignedRankTest(x); tail=:left)[1] == -Inf
    @test isapprox(confint(ApproximateSignedRankTest(x); tail=:left)[2], 14.45, atol=1e-4)
    @test isapprox(confint(ApproximateSignedRankTest(x); tail=:right)[1], 4.45, atol=1e-4)
    @test confint(ApproximateSignedRankTest(x); tail=:right)[2] == Inf

    # the sample from #368, all three tails against R:
    # wilcox.test(z, conf.int = TRUE, alternative = "less")    -> (-Inf, 0.95)
    # wilcox.test(z, conf.int = TRUE, alternative = "greater") -> (-0.45, Inf)
    # wilcox.test(z, conf.int = TRUE)                          -> (-0.7, 1.1)
    z = [1.2, 2.3, 3.1, 4.4, 2.8, 3.9, 5.1, 2.2, 3.3, 4.0] .- 3
    @test all(isapprox.(confint(ExactSignedRankTest(z); tail=:left), (-Inf, 0.95); atol=1e-9))
    @test all(isapprox.(confint(ExactSignedRankTest(z); tail=:right), (-0.45, Inf); atol=1e-9))
    @test all(isapprox.(confint(ExactSignedRankTest(z)), (-0.7, 1.1); atol=1e-9))
end

@testset "Tied one-sided p-values against exactRankTests" begin
    # R's exactRankTests::wilcox.exact computes the same conditional distribution as the
    # tied enumeration here. Its one-sided p-values are convention-free and pin ours
    # exactly; its two-sided can differ from doubling only where the conditional
    # distribution is asymmetric, which the sign-flip symmetry rules out one-sample.
    h = [1, 2, 3, 4, 5, 6, 7, 10, 10, 10, 10, 10, 13, 14, 15] .- 10.1
    # wilcox.exact(h): 0.04052734375; "less": 0.02026367188; "greater": 0.9828491211
    @test pvalue(ExactSignedRankTest(h)) == 0.04052734375
    @test isapprox(pvalue(ExactSignedRankTest(h); tail=:left), 0.02026367188, atol=1e-9)
    @test isapprox(pvalue(ExactSignedRankTest(h); tail=:right), 0.9828491211, atol=1e-9)
    # with zeros as well as ties; wilcox.exact drops the zeros likewise
    # wilcox.exact(d): 0.3071899414; "less": 0.8592834473; "greater": 0.1535949707
    d = [0, 0, 0, 0.5, 0.5, 1, -0.5, -1, 1.5, -1.5, 0.5, 0, 1, -0.5, 2, 0, 0.5, -1, 1, 0.5]
    @test isapprox(pvalue(SignedRankTest(d); tail=:left), 0.8592834473, atol=1e-9)
    @test isapprox(pvalue(SignedRankTest(d); tail=:right), 0.1535949707, atol=1e-9)
end

@testset "Hodges-Lehmann estimate" begin
    x = [-7.8, -6.9, -4.7, 3.7, 6.5, 8.7, 9.1, 10.1, 10.8, 13.6, 14.4, 16.6, 20.2, 22.4, 23.5]
    walsh = [(x[i] + x[j]) / 2 for i in eachindex(x) for j in i:length(x)]
    # R reports 9.675 as the estimate for this sample
    @test @inferred(hodgeslehmann(SignedRankTest(x))) ≈ median(walsh) ≈ 9.675
    @test hodgeslehmann(ExactSignedRankTest(x)) == hodgeslehmann(ApproximateSignedRankTest(x))
    lo, hi = confint(ExactSignedRankTest(x))
    @test lo <= hodgeslehmann(ExactSignedRankTest(x)) <= hi
    # the estimate that is reported is the one the interval is built around
    @test population_param_of_interest(ExactSignedRankTest(x))[3] == hodgeslehmann(ExactSignedRankTest(x))
    @test population_param_of_interest(ExactSignedRankTest(x))[3] != median(x)

    # ties: R's wilcox.test reports -2.1; exactRankTests::wilcox.exact reports -2.35,
    # which is its conditional estimator, not this one
    h = [1, 2, 3, 4, 5, 6, 7, 10, 10, 10, 10, 10, 13, 14, 15] .- 10.1
    @test hodgeslehmann(SignedRankTest(h)) ≈ -2.1
end

@testset "One observation" begin
    # One difference gives one Walsh average, so the only interval is that point, and
    # the index rule has nothing to choose between: the scan this bisection replaced was
    # handed the empty range 1:0 and `argmin` raised. R's `wilcox.test` returns (1.5, 1.5).
    #
    # Every interval here warns, and should: one observation cannot cover 95% by either
    # route. The exact route names the coverage it does reach, 0 two-sided and 0.5
    # one-sided; the approximate route reports its endpoint falling outside the set
    # instead, since unattainability is not what its rule establishes.
    t = SignedRankTest([1.5])
    @test (@test_logs (:warn, r"not attainable") confint(t)) == (1.5, 1.5)
    @test (@test_logs (:warn, r"not attainable") confint(t; level=0.99)) == (1.5, 1.5)
    # `tail = :left` keeps the endpoint inverting the test of the same name (#368)
    @test (@test_logs (:warn, r"not attainable") confint(t; tail=:left)) == (-Inf, 1.5)
    @test (@test_logs (:warn, r"not attainable") confint(t; tail=:right)) == (1.5, Inf)
    @test hodgeslehmann(t) == 1.5
    @test pvalue(t) == 1.0
    @test (@test_logs (:warn, r"outside the set of pairwise estimates") confint(ApproximateSignedRankTest([1.5]))) == (1.5, 1.5)
    # and through the two-sample form
    @test (@test_logs (:warn, r"not attainable") confint(SignedRankTest([4.0], [2.5]))) == (1.5, 1.5)
    # the two-sample tests already agreed: one against one is the one difference
    @test (@test_logs (:warn, r"not attainable") confint(MannWhitneyUTest([4.0], [2.5]))) == (1.5, 1.5)
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
    # the interval is read off the 120 Walsh averages of the 15 non-zero differences,
    # not the 210 of all 20
    @test length(HypothesisTests.walsh_averages(d)) == 210
    @test length(HypothesisTests.signedrank_pairwise_estimates(d)) == 120
    # all-zero input still yields the degenerate interval rather than an error, and says
    # that the level was not met
    @test (@test_logs (:warn, r"not attainable") confint(SignedRankTest(zeros(6)))) == (0.0, 0.0)
end

@testset "Unattainable levels warn" begin
    # the widest interval a sample of n has is (V_(1), V_(m)), covering 1 - 2^(1-n), so
    # 0.95 is out of reach below n = 6 and 0.99 below n = 8
    small = [1.4, -0.6, 2.3, 0.8, -1.9]                       # n = 5, widest is 0.9375
    ci = @test_logs (:warn, r"not attainable") confint(ExactSignedRankTest(small))
    walsh = sort([(small[i] + small[j]) / 2 for i in eachindex(small) for j in i:length(small)])
    @test ci == (first(walsh), last(walsh))                   # the widest available
    # R returns the same interval, (-1.9, 2.3), but silently: its conf.level attribute
    # still claims 0.95, and it warns only once the shortfall exceeds alpha/2, as at
    # the 2 v 2 case below
    @test ci == (-1.9, 2.3)
    @test_logs (:warn, r"not attainable") confint(ExactSignedRankTest(small); level=0.99)
    # R returns (-3, -1) with an achieved conf.level attribute of 2/3, which is the
    # attainable coverage the warning here names: 1 - 2 P(U <= 0) = 1 - 2/6
    ci22 = @test_logs (:warn, r"not attainable") confint(ExactMannWhitneyUTest([1.0, 2.0], [3.0, 4.0]))
    @test ci22 == (-3.0, -1.0)
    # the approximate route warns for its own reason, and does not claim unattainability:
    # at n = 8, level = 0.99 it asks for k = -1 while the exact route reaches 0.9922
    @test_logs (:warn, r"outside the set of pairwise estimates") confint(ApproximateSignedRankTest(small))
    @test_logs (:warn, r"outside the set of pairwise estimates") confint(
        ApproximateSignedRankTest([1.4, -0.6, 2.3, 0.8, -1.9, 3.1, 2.7, -0.3]); level=0.99)
    # and no warning once the level is attainable
    @test_logs confint(ExactSignedRankTest([1.4, -0.6, 2.3, 0.8, -1.9, 3.1]))
    @test_logs confint(ExactSignedRankTest(collect(1.0:15)); level=0.99)
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

    # a callable need not be a `Function`; the keyword documents "a callable"
    @test !(SignedRankExactBelow(100) isa Function)
    @test SignedRankTest(big; method = SignedRankExactBelow(100)) isa ExactSignedRankTest
    @test SignedRankTest(big; method = SignedRankExactBelow(10)) isa ApproximateSignedRankTest

    @test_throws ArgumentError SignedRankTest(x; method=:fastest)
    @test_throws ArgumentError SignedRankTest(x; method = s -> :fastest)
    @test_throws ArgumentError SignedRankTest(x; method = 3)
    # a callable this cannot be called with is rejected with the same message rather
    # than raising a MethodError from inside the constructor
    @test_throws ArgumentError SignedRankTest(x; method = (a, b) -> :exact)

    # the choice reaches the interval, not just the type
    @test confint(SignedRankTest(x; method=:exact)) !=
          confint(SignedRankTest(x; method=:approximate))

    # the callable is handed the documented fields
    seen = Ref{Any}(nothing)
    SignedRankTest([1.0, 2, 3, 3, 4]; method = s -> (seen[] = s; :exact))
    @test issetequal(keys(seen[]), (:n, :n_nonzero, :ties, :tie_adjustment))
    @test seen[].n == 5
    @test seen[].n_nonzero == 5
    @test seen[].ties
    @test seen[].tie_adjustment == SignedRankTest([1.0, 2, 3, 3, 4]).tie_adjustment
    # and zeros are excluded from n_nonzero but not from n
    SignedRankTest([0.0, 1, 2, 3]; method = s -> (seen[] = s; :exact))
    @test (seen[].n, seen[].n_nonzero, seen[].ties) == (4, 3, false)
end

@testset "Interval invariants" begin
    for n in (8, 13, 21, 34), seed in (1, 2)
        v = [sin(seed * 1.7 + i * 0.9) * 3 + cos(i * 0.31) for i in 1:n]
        for t in (SignedRankTest(v; method=:exact), SignedRankTest(v; method=:approximate))
            # the estimate the interval is built around lies inside it
            lo, hi = confint(t)
            @test lo <= hodgeslehmann(t) <= hi
            # raising the level cannot narrow the interval (the approximate route warns at
            # n = 8, where its rule wants an endpoint outside the pairwise estimates)
            wide = n == 8 && t isa ApproximateSignedRankTest ?
                (@test_logs (:warn, r"outside the set of pairwise estimates") confint(t; level=0.99)) :
                confint(t; level=0.99)
            @test wide[1] <= lo && wide[2] >= hi
            # one-sided bounds are the corresponding two-sided endpoints, and open on the
            # side the tail's alternative allows (#368)
            @test confint(t; tail=:left)[1] == -Inf
            @test confint(t; tail=:right)[2] == Inf
            @test confint(t; tail=:left)[2] <= hi
            @test confint(t; tail=:right)[1] >= lo
        end
    end
end

@testset "Element types other than Float64" begin
    ints = [3, -1, 4, -1, 5, 9, -2, 6, 5, 3]
    ref = confint(SignedRankTest(Float64.(ints)))
    for T in (Int, Int32, Float64, Float32, Float16, BigFloat, Rational{Int})
        x = T.(ints)
        t = SignedRankTest(x)
        @test t.median == median(x)                       # value survives the conversion
        @test all(isapprox.(confint(t), ref; atol = 1e-2)) # Float16 is the loose one
        @test hodgeslehmann(t) ≈ hodgeslehmann(SignedRankTest(Float64.(ints))) atol = 1e-2
    end
    # every element type above is accepted, and every one is stored as Float64: the
    # fields are converted rather than parametrised
    @test SignedRankTest(Float32.(ints)).median isa Float64
    @test SignedRankTest(BigFloat.(ints)).median isa Float64
    @test SignedRankTest(BigFloat.(ints)).median == median(BigFloat.(ints))
    @test SignedRankTest(ints).median isa Float64
    @test SignedRankTest(Float32.(ints)).vals isa Vector{Float64}
    # and neither type carries a parameter, so these spellings stay concrete
    @test isconcretetype(ExactSignedRankTest)
    @test isconcretetype(ApproximateSignedRankTest)
    # both routes, and the two-sample form
    @test ExactSignedRankTest(Float32.(ints)) isa ExactSignedRankTest
    @test ApproximateSignedRankTest(Float32.(ints)) isa ApproximateSignedRankTest
    @test SignedRankTest(Float32[1, 2, 3, 4, 5, 6], Float32[2, 3, 4, 5, 6, 8]) isa
        HypothesisTests.HypothesisTest
    # `level` is declared `Real` too, so it accepts any real rather than only Float64
    for t in (ExactSignedRankTest(Float64.(ints)), ApproximateSignedRankTest(Float64.(ints)))
        lref = confint(t; level=0.90)
        @test confint(t; level=0.90f0) == lref
        @test confint(t; level=9//10) == lref
        @test_throws ArgumentError confint(t; level=1//2)
    end
end

@testset "walsh_averages states its indexing assumption" begin
    # `walsh_averages` indexes over `1:n` with `@inbounds`, so a vector indexed otherwise
    # is refused rather than read out of bounds
    @test_throws ArgumentError HypothesisTests.walsh_averages(ZeroBasedVector([1.0, 2.0, 3.0]))
    @test HypothesisTests.walsh_averages([1.0, 2.0, 3.0]) == [1.0, 1.5, 2.0, 2.0, 2.5, 3.0]
    # `cross_differences` iterates values, so it needs no such guard
    @test HypothesisTests.cross_differences(ZeroBasedVector([3.0]), ZeroBasedVector([1.0])) == [2.0]
end

@testset "Enumeration and pairwise estimates are bounded" begin
    @test HypothesisTests.MAX_PAIRWISE_ESTIMATES == 1_000_000
    @test HypothesisTests.check_estimate_count(10, "Walsh averages") === nothing
    @test_throws HypothesisTests.ComputationTooLarge HypothesisTests.check_estimate_count(10^7, "Walsh averages")

    # the enumeration bound bites at exactly MAX_EXACT_ENUMERATION_N non-zero values
    N = HypothesisTests.MAX_EXACT_ENUMERATION_N
    @test N == 30
    @test HypothesisTests.check_exact_enumeration(N) === nothing
    @test_throws HypothesisTests.ComputationTooLarge HypothesisTests.check_exact_enumeration(N + 1)
    # and on the two-sample side at binomial(nx + ny, min(nx, ny)) > 2^N
    @test HypothesisTests.check_exact_enumeration(2, 2) === nothing
    @test_throws HypothesisTests.ComputationTooLarge HypothesisTests.check_exact_enumeration(30, 30)
    # a binomial too large to represent is bounded rather than thrown from
    @test HypothesisTests.rank_binomial(10^6, 10^6) === nothing
    @test_throws HypothesisTests.ComputationTooLarge HypothesisTests.check_exact_enumeration(10^6, 10^6)

    # the signed rank exact interval carries no time bound of its own: the bisection is
    # cheap enough that the memory bound above is the one it meets. n = 300 forms
    # m = 45,150 Walsh averages and inverts in under a tenth of a second, where the
    # m/2-recursion scan this branch replaced took 39 s and the branch below refused
    # past n = 223 (MAX_EXACT_CI_ESTIMATES bounds the two-sample route alone now; its
    # boundary is asserted in test/mann_whitney.jl)
    @test confint(ExactSignedRankTest(collect(1.0:300))) isa Tuple
    # n = 2000, the size #97 reports hanging, is stopped by the memory bound on
    # materialising the set at all
    big = SignedRankTest(collect(1.0:2000))
    @test_throws HypothesisTests.ComputationTooLarge confint(big)
    # the test is still printable, without its interval line
    out = sprint(show, big)
    @test occursin("Wilcoxon signed rank", out) || occursin("Wilcoxon", out)
    @test !occursin("confidence interval", out)
    # and a sample inside the bound still gets one
    @test confint(SignedRankTest(collect(1.0:60))) isa Tuple
end

@testset "Exact route past the automatic thresholds" begin
    # StatsFuns 2.2.1 fixed the Int overflow in signrankcdf (StatsFuns#219, #221), so
    # the exact route stays valid well past where the automatic rule abandons it
    big = collect(1.0:100)
    @test SignedRankTest(big) isa ApproximateSignedRankTest
    @test 0 <= pvalue(SignedRankTest(big; method=:exact)) <= 1
    # the null distribution is symmetric about n(n+1)/4; this identity is what failed
    # on StatsFuns <= 2.2.0, where the counts overflowed and the cdf became garbage
    let n = 100, M = div(n * (n + 1), 2)
        for w in (1000, 2525, 4000)
            @test HypothesisTests.signrankcdf(n, w) +
                  HypothesisTests.signrankcdf(n, M - w - 1) ≈ 1
        end
    end
    # the tied branch will not enumerate without bound, though
    tied = repeat([1.0, 2.0], 20)
    @test_throws HypothesisTests.ComputationTooLarge pvalue(SignedRankTest(tied; method=:exact))
end

@testset "Reflection under negation" begin
    # The two-sample reflection is swapping the samples; the one-sample one is negating
    # the sample, which negates the location the test is about. Exact for the same
    # reason: negation is exact, and the index rules see only n and the tie adjustment.
    for x in ([1.4, -2.6, 3.1, -0.7, 5.2, -4.8, 2.9],        # no ties, no zeros
              [1.0, -1.0, 2.0, 2.0, -2.0, 3.0, 3.0, -3.0],   # tied absolute values
              [0.0, 1.5, -2.5, 0.0, 3.5, -0.5, 4.0, 0.0])    # zeros, which are dropped
        for T in (ExactSignedRankTest, ApproximateSignedRankTest)
            a, b = T(x), T(-x)
            @test hodgeslehmann(a) == -hodgeslehmann(b)
            @test population_param_of_interest(a)[3] ==
                  -population_param_of_interest(b)[3]
            # W sums the positive signed ranks, so the pair exhausts the rank total
            # over the non-zero observations rather than negating
            nz = length(a.ranks)
            @test a.W + b.W == nz * (nz + 1) / 2
            @test pvalue(a) == pvalue(b)
            @test pvalue(a; tail=:left) == pvalue(b; tail=:right)
            for level in (0.90, 0.95, 0.99)
                # these samples are small enough that the higher levels are out of
                # reach, so both routes warn here. That behaviour is pinned in
                # "Unattainable levels warn"; in this testset it is noise, and the
                # identity being checked is the reflection.
                (lo_a, hi_a), (lo_b, hi_b) =
                    Base.CoreLogging.with_logger(Base.CoreLogging.NullLogger()) do
                        confint(a; level=level), confint(b; level=level)
                    end
                @test lo_a == -hi_b
                @test hi_a == -lo_b
            end
            # compare the finite ends: under the #368 convention `:left` is an upper
            # bound and `:right` a lower one, so `[1]` on both would be -Inf either way
            @test confint(a; tail=:left)[2] == -confint(b; tail=:right)[1]
        end
    end
end

@testset "AbstractVector inputs" begin
    # ranges and views used to fail with a MethodError: the structs store Vector,
    # and nothing converted on the way in
    @test pvalue(SignedRankTest(-2:7)) == pvalue(SignedRankTest([-2:7;]))
    @test pvalue(ExactSignedRankTest(-2:7)) == pvalue(ExactSignedRankTest([-2:7;]))
    @test pvalue(ApproximateSignedRankTest(1.0:20.0)) ==
          pvalue(ApproximateSignedRankTest([1.0:20.0;]))
    v = [1.0, -2.0, 3.0, -4.0, 5.0, 6.0]
    @test pvalue(SignedRankTest(view(v, 2:5))) == pvalue(SignedRankTest(v[2:5]))
end

@testset "No Int64 overflow at two million observations" begin
    # nz(nz+1)(2nz+1) wrapped in Int64 from nz = 2^21, making sigma three orders of
    # magnitude too small with no error raised; t^3 in the tie adjustment wrapped the
    # same way for a single group of 2^21 equal values
    nz = 2^21
    t = ApproximateSignedRankTest(collect(1.0:nz))
    @test t.sigma ≈ Float64(sqrt(big(nz) * (nz + 1) * (2 * nz + 1) // 24)) rtol = 1e-12
    @test HypothesisTests.tiedrank_adj(fill(1.0, nz))[2] ≈ Float64(big(nz)^3 - nz) rtol = 1e-12
end

@testset "Attainable levels do not warn" begin
    small = [1.4, -0.6, 2.3, 0.8, -1.9]              # n = 5: P(W <= 0) = 1/32
    # alpha/2 = 1/32 exactly: k = 0 attains the request, so there is nothing to warn
    # about. R returns the same (-1.9, 2.3) at conf.level = 0.9375, silently.
    @test (@test_logs confint(ExactSignedRankTest(small); level = 1 - 2/32)) == (-1.9, 2.3)
end

@testset "One-sided warnings name the one-sided request" begin
    small = [1.4, -0.6, 2.3, 0.8, -1.9]
    logs, ci = Test.collect_test_logs() do
        confint(ExactSignedRankTest(small); level = 0.99, tail = :left)
    end
    @test ci == (-Inf, 2.3)
    @test length(logs) == 1
    # the request was 0.99, not the two-sided 0.98 this used to report, and the bound
    # misses by its one overshooting tail only: 1 - 1/32, not 1 - 2/32
    @test logs[1].kwargs[:requested] == 0.99
    @test logs[1].kwargs[:attainable] == 1 - 1/32
end

@testset "One-sided intervals on the tied routes against R" begin
    h = [1, 2, 3, 4, 5, 6, 7, 10, 10, 10, 10, 10, 13, 14, 15] .- 10.1
    # base R falls back to the normal approximation under ties and lands on the same
    # order statistics as this route, which inverts the untied lattice:
    # wilcox.test(h, conf.int = TRUE, alternative = "less") -> (-Inf, -0.09997),
    # "greater" -> (-4.10001, Inf). exactRankTests::wilcox.exact, which inverts the
    # tied conditional distribution instead, gives (-4.6, Inf) for "greater".
    te = ExactSignedRankTest(h)
    @test confint(te; tail = :left)[1] == -Inf
    @test isapprox(confint(te; tail = :left)[2], -0.1, atol = 1e-9)
    @test isapprox(confint(te; tail = :right)[1], -4.1, atol = 1e-9)
    @test confint(te; tail = :right)[2] == Inf

    # the approximate test's tie-corrected sigma reaches one order statistic further
    # on the left: R gives (-4.60008, -0.09997) with exact = FALSE
    ta = ApproximateSignedRankTest(h)
    @test all(isapprox.(confint(ta), (-4.6, -0.1); atol = 1e-9))
end

@testset "Auto boundaries at n = 50 and 51, against R" begin
    s50 = quantile.(Normal(), (1:50) ./ 51) .+ 0.3
    s51 = quantile.(Normal(), (1:51) ./ 52) .+ 0.3
    @test SignedRankTest(s50) isa ExactSignedRankTest
    @test SignedRankTest(s51) isa ApproximateSignedRankTest
    # R switches at n < 50, so its exact side must be forced to compare:
    # wilcox.test(qnorm((1:50)/51) + 0.3, exact = TRUE)  -> 0.0371666050798
    # wilcox.test(qnorm((1:51)/52) + 0.3, exact = FALSE) -> 0.0345394645593
    @test isapprox(pvalue(SignedRankTest(s50)), 0.0371666050798, atol = 1e-10)
    @test isapprox(pvalue(SignedRankTest(s51)), 0.0345394645593, atol = 1e-10)
end

@testset "Tied auto boundary counts non-zero observations" begin
    tied16 = repeat([1.0, -1.0, 2.5, 3.5], 4)        # 16 non-zero, tied
    @test SignedRankTest(tied16) isa ApproximateSignedRankTest
    @test SignedRankTest(tied16[1:15]) isa ExactSignedRankTest
end

@testset "Levels at the edges of (0.5, 1)" begin
    x = [-7.8, -6.9, -4.7, 3.7, 6.5, 8.7, 9.1, 10.1, 10.8, 13.6, 14.4, 16.6, 20.2, 22.4, 23.5]
    t = ExactSignedRankTest(x)
    @test_throws ArgumentError confint(t; level = 0.5)
    @test_throws ArgumentError confint(t; level = 1.0)
    # R: wilcox.test(x, conf.int = TRUE, conf.level = 0.51) -> (7.85, 11.75)
    @test all(isapprox.(confint(t; level = 0.51), (7.85, 11.75); atol = 1e-9))
    # R: alternative = "less", conf.level = 0.99 -> (-Inf, 16.3)
    @test confint(t; level = 0.99, tail = :left)[1] == -Inf
    @test isapprox(confint(t; level = 0.99, tail = :left)[2], 16.3, atol = 1e-9)
end
end
