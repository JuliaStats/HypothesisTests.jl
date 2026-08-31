using HypothesisTests, Test
using HypothesisTests: default_tail, population_param_of_interest
using Statistics: median

# A parameterised route rule written the natural way, as a struct you can store, reuse
# and print. It is callable, and it is not a subtype of `Function`: defining a call
# method does not place a type under `Function` unless it says so.
struct SignedRankExactBelow
    n::Int
end
(r::SignedRankExactBelow)(s) = s.n_nonzero <= r.n ? :exact : :approximate

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
    # PENDING #361: the approximate test still reports the exact interval, because
    # `calculate_ci` inverts the exact null distribution whichever test it is called on.
    # R with exact = FALSE gives (3.05004, 15.50001).
    @test isapprox(@inferred(confint(ApproximateSignedRankTest(x)))[1], 3.3, atol=1e-4)
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
    lo, hi = confint(ExactSignedRankTest(x))
    @test lo <= hodgeslehmann(ExactSignedRankTest(x)) <= hi
    # PENDING #363: the reported estimate is still the sample median, not the
    # Hodges-Lehmann estimate the "pseudomedian" label names.
    @test population_param_of_interest(ExactSignedRankTest(x))[3] == median(x) == 10.1
    @test population_param_of_interest(ExactSignedRankTest(x))[3] != hodgeslehmann(ExactSignedRankTest(x))

    # ties: R reports -2.1
    h = [1, 2, 3, 4, 5, 6, 7, 10, 10, 10, 10, 10, 13, 14, 15] .- 10.1
    @test hodgeslehmann(SignedRankTest(h)) ≈ -2.1
end

@testset "One observation" begin
    # One difference gives one Walsh average, so the only interval is that point, and
    # the scan that picks an order statistic has nothing to pick between: it was handed
    # the empty range 1:0 and `argmin` raised. R's `wilcox.test` returns (1.5, 1.5).
    t = SignedRankTest([1.5])
    @test confint(t) == (1.5, 1.5)
    @test confint(t; level=0.99) == (1.5, 1.5)
    @test confint(t; tail=:left) == (1.5, Inf)
    @test confint(t; tail=:right) == (-Inf, 1.5)
    @test hodgeslehmann(t) == 1.5
    @test pvalue(t) == 1.0
    @test confint(ApproximateSignedRankTest([1.5])) == (1.5, 1.5)
    # and through the two-sample form
    @test confint(SignedRankTest([4.0], [2.5])) == (1.5, 1.5)
    # the two-sample tests already agreed: one against one is the one difference
    @test confint(MannWhitneyUTest([4.0], [2.5])) == (1.5, 1.5)
end

@testset "Zero differences" begin
    d = [0.0, 0, 0, 0.5, 0.5, 1, -0.5, -1, 1.5, -1.5, 0.5, 0, 1, -0.5, 2, 0, 0.5, -1, 1, 0.5]
    nz = d[d .!= 0]
    # the statistic and the p-value already ignore zeros
    @test SignedRankTest(d).W == SignedRankTest(nz).W
    @test pvalue(SignedRankTest(d)) == pvalue(SignedRankTest(nz))
    # and so does the estimate, which is formed from the same pairwise estimates
    @test hodgeslehmann(SignedRankTest(d)) == hodgeslehmann(SignedRankTest(nz)) == 0.5

    # PENDING #362: the interval does not. It is built from the Walsh averages of all
    # 20 differences rather than the 15 non-zero ones, so it describes a different
    # sample from the p-value beside it. R gives (-0.24999, 0.75005) at this level.
    @test confint(SignedRankTest(d); level=0.9) != confint(SignedRankTest(nz); level=0.9)
    @test confint(SignedRankTest(d); level=0.9) == (0.0, 0.5)
    @test length(HypothesisTests.walsh_averages(d)) == 210
    @test length(HypothesisTests.signedrank_pairwise_estimates(d)) == 120

    # PENDING #362: and so the estimate can fall outside the interval printed beside it.
    # Here the p-value and the estimate describe the seven non-zero differences while the
    # interval is built from all fourteen, so it lands entirely below the estimate. The
    # interval of the sample the estimate does describe contains it. This happens on
    # roughly a tenth of samples that contain zeros, so it is worth a sample of its own.
    e = [0.0, 3, 2, 0, 0, 3, 0, -1, 0, 0, 3, 0, 3, 3]
    enz = e[e .!= 0]
    @test pvalue(SignedRankTest(e)) == pvalue(SignedRankTest(enz))
    @test hodgeslehmann(SignedRankTest(e)) == hodgeslehmann(SignedRankTest(enz)) == 3.0
    lo, hi = confint(SignedRankTest(e))
    @test (lo, hi) == (0.0, 1.5)
    @test !(lo <= hodgeslehmann(SignedRankTest(e)) <= hi)
    @test confint(SignedRankTest(enz)) == (1.0, 3.0)
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

    # PENDING #361: the choice reaches the type but not yet the interval
    @test confint(SignedRankTest(x; method=:exact)) ==
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

@testset "Element types other than Float64" begin
    ints = [3, -1, 4, -1, 5, 9, -2, 6, 5, 3]
    ref = confint(SignedRankTest(Float64.(ints)))
    for T in (Int, Int32, Float64, Float32, Float16, BigFloat, Rational{Int})
        x = T.(ints)
        t = SignedRankTest(x)
        @test t.median == median(x)                       # and not silently widened
        @test all(isapprox.(confint(t), ref; atol = 1e-2)) # Float16 is the loose one
        @test hodgeslehmann(t) ≈ hodgeslehmann(SignedRankTest(Float64.(ints))) atol = 1e-2
    end
    # the median keeps its own type rather than being truncated to Float64
    @test SignedRankTest(Float32.(ints)).median isa Float32
    @test SignedRankTest(BigFloat.(ints)).median isa BigFloat
    @test SignedRankTest(BigFloat.(ints)).median == median(BigFloat.(ints))
    # integers still give a Float64 median, as median() does
    @test SignedRankTest(ints).median isa Float64
    # both concrete types, and the two-sample form
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

@testset "Construction snapshots the sample" begin
    # `vals` is what `confint` and `hodgeslehmann` read, while every other field is
    # computed at construction, so keeping a view of the caller's array would let the
    # interval and the estimate describe a different sample from the p-value. This was
    # already so on master, where `confint` read `vals`. See
    # JuliaStats/HypothesisTests.jl#375.
    for f in (SignedRankTest,
              x -> SignedRankTest(x; method = :approximate),
              ExactSignedRankTest,
              ApproximateSignedRankTest)
        xs = [1.0, 2.0, 3.0, 4.0, 9.0, -2.0, 0.5]
        t = f(xs)
        @test t.vals !== xs
        p, ci, hl = pvalue(t), confint(t), hodgeslehmann(t)
        xs[1] += 100.0
        @test pvalue(t) == p
        @test confint(t) == ci
        @test hodgeslehmann(t) == hl
    end
    # the seven-argument approximate constructor can be called directly, so it takes
    # the snapshot itself rather than relying on its callers to have done so
    xs = [1.0, 2.0, 3.0, 4.0, 9.0, -2.0, 0.5]
    t = ApproximateSignedRankTest(xs, HypothesisTests.signedrankstats(xs)...)
    @test t.vals !== xs
    ci = confint(t)
    xs[1] += 100.0
    @test confint(t) == ci
end

@testset "Enumeration and pairwise estimates are bounded" begin
    @test HypothesisTests.MAX_PAIRWISE_ESTIMATES == 100_000_000
    @test HypothesisTests.check_estimate_count(10, "Walsh averages") === nothing
    @test_throws ArgumentError HypothesisTests.check_estimate_count(10^9, "Walsh averages")

    # the enumeration bound bites at exactly MAX_EXACT_ENUMERATION_N non-zero values
    N = HypothesisTests.MAX_EXACT_ENUMERATION_N
    @test N == 25
    @test HypothesisTests.check_exact_enumeration(N) === nothing
    @test_throws ArgumentError HypothesisTests.check_exact_enumeration(N + 1)
    # and on the two-sample side at binomial(nx + ny, min(nx, ny)) > 2^N
    @test HypothesisTests.check_exact_enumeration(2, 2) === nothing
    @test_throws ArgumentError HypothesisTests.check_exact_enumeration(30, 30)
    # a binomial too large to represent is bounded rather than thrown from
    @test HypothesisTests.rank_binomial(10^6, 10^6) === nothing
    @test_throws ArgumentError HypothesisTests.check_exact_enumeration(10^6, 10^6)
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
    @test_throws ArgumentError pvalue(SignedRankTest(tied; method=:exact))
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
                lo_a, hi_a = confint(a; level=level)
                lo_b, hi_b = confint(b; level=level)
                @test lo_a == -hi_b
                @test hi_a == -lo_b
            end
            @test confint(a; tail=:left)[1] == -confint(b; tail=:right)[2]
        end
    end
end
end
