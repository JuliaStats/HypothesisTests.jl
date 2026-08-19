using HypothesisTests, Test
using HypothesisTests: default_tail, population_param_of_interest
using Statistics: median

# Callable, and deliberately not a subtype of `Function`. See test/wilcoxon.jl.
struct MannWhitneyExactBelow
    n::Int
end
(r::MannWhitneyExactBelow)(s) = s.nx + s.ny <= r.n ? :exact : :approximate

@testset "Mann-Whitney" begin
@testset "Basic exact test" begin
    test = ExactMannWhitneyUTest([1:10;], [2.1:2:21;])
    @test test.nx == 10
    @test test.ny == 10
    @test test.ranks == [1, 2, 4, 5, 7, 8, 10, 11, 13, 14, 3, 6, 9, 12, 15, 16, 17, 18, 19, 20]
    @test iszero(test.tie_adjustment)
    @test test.U == 20
    @test test.median == -5.6

    @test default_tail(test) == :both
    @test repr(test) == """
    Exact Mann-Whitney U test
    -------------------------
    Population details:
        parameter of interest:   Location parameter (pseudomedian)
        value under h_0:         0
        point estimate:          -5.6
        95% confidence interval: (-11.1, -0.1)

    Test summary:
        outcome with 95% confidence: reject h_0
        two-sided p-value:           0.0232
    
    Details:
        number of observations in each group: [10, 10]
        Mann-Whitney-U statistic:             20.0
        rank sums:                            [75.0, 135.0]
        adjustment for ties:                  0.0
    """

    # Two-sided
    for kwargs in ((), (; tail = :both))
        @test abs(@inferred(pvalue(ExactMannWhitneyUTest([1:10;], [2.1:2:21;]); kwargs...)) - 0.0232) <= 1e-4
        @test abs(@inferred(pvalue(ExactMannWhitneyUTest([2.1:2:21;], [1:10;]); kwargs...)) - 0.0232) <= 1e-4
        @test abs(@inferred(pvalue(ExactMannWhitneyUTest([1.5:10:100;], [2.1:2:21;]); kwargs...)) - 0.0068) <= 1e-4
        @test abs(@inferred(pvalue(ExactMannWhitneyUTest([2.1:2:21;], [1.5:10:100;]); kwargs...)) - 0.0068) <= 1e-4
    end

    # Left tail
    @test abs(@inferred(pvalue(ExactMannWhitneyUTest([1:10;], [2.1:2:21;]); tail = :left)) - 0.0116) <= 1e-4
    @test abs(@inferred(pvalue(ExactMannWhitneyUTest([2.1:2:21;], [1:10;]); tail = :left)) - 0.9907) <= 1e-4
    @test abs(@inferred(pvalue(ExactMannWhitneyUTest([1.5:10:100;], [2.1:2:21;]); tail = :left)) - 0.9974) <= 1e-4
    @test abs(@inferred(pvalue(ExactMannWhitneyUTest([2.1:2:21;], [1.5:10:100;]); tail = :left)) - 0.0034) <= 1e-4


    # Right tail
    @test abs(@inferred(pvalue(ExactMannWhitneyUTest([1:10;], [2.1:2:21;]); tail = :right)) - 0.9907) <= 1e-4
    @test abs(@inferred(pvalue(ExactMannWhitneyUTest([2.1:2:21;], [1:10;]); tail = :right)) - 0.0116) <= 1e-4
    @test abs(@inferred(pvalue(ExactMannWhitneyUTest([1.5:10:100;], [2.1:2:21;]); tail = :right)) - 0.0034) <= 1e-4
    @test abs(@inferred(pvalue(ExactMannWhitneyUTest([2.1:2:21;], [1.5:10:100;]); tail = :right)) - 0.9974) <= 1e-4
end

@testset "Exact with ties" begin
    show(IOBuffer(), ExactMannWhitneyUTest([1:10;], [1:10;]))

    # Two-sided
    for kwargs in ((), (; tail = :both))
        @test abs(@inferred(pvalue(ExactMannWhitneyUTest([1:10;], [1:10;]); kwargs...)) - 1) <= 1e-4
        @test abs(@inferred(pvalue(ExactMannWhitneyUTest([1:10;], [2:11;]); kwargs...)) - 0.5096) <= 1e-4
        @test abs(@inferred(pvalue(ExactMannWhitneyUTest([2:11;], [1:10;]); kwargs...)) - 0.5096) <= 1e-4
        @test abs(@inferred(pvalue(ExactMannWhitneyUTest([1:10;], [1:5; ones(5)]); kwargs...)) - 0.0057) <= 1e-4
        @test abs(@inferred(pvalue(ExactMannWhitneyUTest([1:5; ones(5)], [1:10;]); kwargs...)) - 0.0057) <= 1e-4
    end

    # Left tail
    @test abs(@inferred(pvalue(ExactMannWhitneyUTest([1:10;], [1:10;]); tail = :left)) - 0.5296) <= 1e-4
    @test abs(@inferred(pvalue(ExactMannWhitneyUTest([1:10;], [2:11;]); tail = :left)) - 0.2548) <= 1e-4
    @test abs(@inferred(pvalue(ExactMannWhitneyUTest([2:11;], [1:10;]); tail = :left)) - 0.7634) <= 1e-4
    @test abs(@inferred(pvalue(ExactMannWhitneyUTest([1:10;], [1:5; ones(5)]); tail = :left)) - 0.9979) <= 1e-4
    @test abs(@inferred(pvalue(ExactMannWhitneyUTest([1:5; ones(5)], [1:10;]); tail = :left)) - 0.0028) <= 1e-4

    # Right tail
    @test abs(@inferred(pvalue(ExactMannWhitneyUTest([1:10;], [1:10;]); tail = :right)) - 0.5296) <= 1e-4
    @test abs(@inferred(pvalue(ExactMannWhitneyUTest([1:10;], [2:11;]); tail = :right)) - 0.7634) <= 1e-4
    @test abs(@inferred(pvalue(ExactMannWhitneyUTest([2:11;], [1:10;]); tail = :right)) - 0.2548) <= 1e-4
    @test abs(@inferred(pvalue(ExactMannWhitneyUTest([1:10;], [1:5; ones(5)]); tail = :right)) - 0.0028) <= 1e-4
    @test abs(@inferred(pvalue(ExactMannWhitneyUTest([1:5; ones(5)], [1:10;]); tail = :right)) - 0.9979) <= 1e-4
end

@testset "Exact with ties and unequal lengths" begin
    test1 = ExactMannWhitneyUTest([1:10;], [2:2:24;])
    test2 = ExactMannWhitneyUTest([2:2:24;], [1:10;])
    @test test1.nx == test2.ny == 10
    @test test2.nx == test1.ny == 12
    @test test1.ranks == test2.ranks[[13:22; 1:12]]
    @test test1.tie_adjustment == test2.tie_adjustment
    @test test1.U == 22.5
    @test test1.U + test2.U == 120
    @test test1.median == -test2.median

    @test repr(test1) == """
    Exact Mann-Whitney U test
    -------------------------
    Population details:
        parameter of interest:   Location parameter (pseudomedian)
        value under h_0:         0
        point estimate:          -7.5
        95% confidence interval: (-14.0, -1.0)

    Test summary:
        outcome with 95% confidence: reject h_0
        two-sided p-value:           0.0120

    Details:
        number of observations in each group: [10, 12]
        Mann-Whitney-U statistic:             22.5
        rank sums:                            [77.5, 175.5]
        adjustment for ties:                  30.0
    """
    @test repr(test2) == """
    Exact Mann-Whitney U test
    -------------------------
    Population details:
        parameter of interest:   Location parameter (pseudomedian)
        value under h_0:         0
        point estimate:          7.5
        95% confidence interval: (1.0, 14.0)

    Test summary:
        outcome with 95% confidence: reject h_0
        two-sided p-value:           0.0120

    Details:
        number of observations in each group: [12, 10]
        Mann-Whitney-U statistic:             97.5
        rank sums:                            [175.5, 77.5]
        adjustment for ties:                  30.0
    """

    # Two-sided
    for kwargs in ((), (; tail = :both))
        @test abs(@inferred(pvalue(ExactMannWhitneyUTest([1:10;], [2:2:24;]); kwargs...)) - 0.0120) <= 1e-4
        @test abs(@inferred(pvalue(ExactMannWhitneyUTest([2:2:24;], [1:10;]); kwargs...)) - 0.0120) <= 1e-4
    end

    # Left tail
    @test abs(@inferred(pvalue(ExactMannWhitneyUTest([1:10;], [2:2:24;]); tail = :left)) - 0.0060) <= 1e-4
    @test abs(@inferred(pvalue(ExactMannWhitneyUTest([2:2:24;], [1:10;]); tail = :left)) - 0.9949) <= 1e-4

    # Right tail
    @test abs(@inferred(pvalue(ExactMannWhitneyUTest([1:10;], [2:2:24;]); tail = :right)) - 0.9949) <= 1e-4
    @test abs(@inferred(pvalue(ExactMannWhitneyUTest([2:2:24;], [1:10;]); tail = :right)) - 0.0060) <= 1e-4
end

@testset "Approximate test" begin
    @test abs(@inferred(pvalue(ApproximateMannWhitneyUTest([1:10;], [2.1:2:21;]))) - 0.0257) <= 1e-4
    @test abs(@inferred(pvalue(ApproximateMannWhitneyUTest([2.1:2:21;], [1:10;]))) - 0.0257) <= 1e-4
    @test abs(@inferred(pvalue(ApproximateMannWhitneyUTest([1.5:10:100;], [2.1:2:21;]))) - 0.0091) <= 1e-4
    @test abs(@inferred(pvalue(ApproximateMannWhitneyUTest([2.1:2:21;], [1.5:10:100;]))) - 0.0091) <= 1e-4
    @test default_tail(ApproximateMannWhitneyUTest([1:10;], [2.1:2:21;])) == :both
	show(IOBuffer(), ApproximateMannWhitneyUTest([1:10;], [2.1:2:21;]))
end

@testset "Approximate with ties" begin
    @test abs(@inferred(pvalue(ApproximateMannWhitneyUTest([1:10;], [1:10;]))) - 1) <= 1e-4
    @test abs(@inferred(pvalue(ApproximateMannWhitneyUTest([1:10;], [2:11;]))) - 0.4948) <= 1e-4
    @test abs(@inferred(pvalue(ApproximateMannWhitneyUTest([2:11;], [1:10;]))) - 0.4948) <= 1e-4
    @test abs(@inferred(pvalue(ApproximateMannWhitneyUTest([1:10;], [1:5; ones(5)]))) - 0.0076) <= 1e-4
    @test abs(@inferred(pvalue(ApproximateMannWhitneyUTest([1:5; ones(5)], [1:10;]))) - 0.0076) <= 1e-4
	show(IOBuffer(), ApproximateMannWhitneyUTest([1:10;], [1:10;]))
end

@testset "Tests for automatic selection" begin
    @test abs(@inferred(pvalue(MannWhitneyUTest([1:10;], [2.1:2:21;]))) - 0.0232) <= 1e-4
    @test abs(@inferred(pvalue(MannWhitneyUTest([1:10;], [2:11;]))) - 0.4948) <= 1e-4
	show(IOBuffer(), MannWhitneyUTest([1:10;], [2.1:2:21;]))
end

@testset "Issue #39" begin
    @test abs(@inferred(pvalue(ExactMannWhitneyUTest(Float32[1:10;], Float32[2:11;]))) - 0.5096) <= 1e-4
end

@testset "Issue #113" begin
    @test abs(pvalue(ApproximateMannWhitneyUTest(Float32[1:10;], Float32[2:11;])) - 0.4948) <= 1e-4
end

@testset "Issue #126 pvalues above 1" begin
	A = [1.34937, 1.75722,0.276514, 1.04546, 1.69085, 0.738085, 2.36313]
	B = [2.62325, 1.16533, 1.1327, 0.728714]

	m = ExactMannWhitneyUTest(A,B)
	p = @inferred(pvalue(m; tail = :both))
	@test p == 1

	A = [12,10,7,6,3,1]
	B = [11,9,8,5,4,2]

	m = MannWhitneyUTest(A,B)
	p = @inferred(pvalue(m; tail = :both))

	@test p == 1

	m = ApproximateMannWhitneyUTest(A, B)
	p = @inferred(pvalue(m; tail = :both))
end

@testset "Confidence interval and Hodges-Lehmann estimate" begin
    # R: wilcox.test(x, y, conf.int = TRUE) -> (-11.1, -0.1), estimate -5.6
    x = [1:10;]
    y = [2.1:2:21;]
    @test all(isapprox.(@inferred(confint(ExactMannWhitneyUTest(x, y))), (-11.1, -0.1); atol=1e-8))
    @test @inferred(hodgeslehmann(ExactMannWhitneyUTest(x, y))) ≈ -5.6
    @test hodgeslehmann(ExactMannWhitneyUTest(x, y)) ≈ median([xi - yj for xi in x, yj in y])

    # the samples are kept so that the cross-differences can be formed
    t = MannWhitneyUTest(x, y)
    @test t.x == x
    @test t.y == y

    # ties: R -> (-13.99997, -1.00003), i.e. the order statistics (-14, -1)
    @test confint(ExactMannWhitneyUTest([1:10;], [2:2:24;])) == (-14.0, -1.0)
    @test confint(ExactMannWhitneyUTest([2:2:24;], [1:10;])) == (1.0, 14.0)

    # approximate, with ties. R (correct = TRUE) -> (-0.36998, 1.18297), estimate 0.56
    a = [1.83, 0.50, 1.62, 2.48, 1.68, 1.88, 1.55, 3.06, 1.30]
    b = [0.878, 0.647, 0.598, 2.05, 1.06, 1.29, 1.06, 3.14, 1.29]
    ta = MannWhitneyUTest(a, b)
    @test ta isa ApproximateMannWhitneyUTest
    @test all(isapprox.(confint(ta), (-0.37, 1.183); atol = 1e-4))
    @test hodgeslehmann(ta) ≈ 0.56

    # PENDING #363: the reported estimate is still the difference of sample medians,
    # not the Hodges-Lehmann estimate the "pseudomedian" label names and the interval
    # above is built around. R reports 0.56005616 for this sample.
    @test population_param_of_interest(ta)[3] ≈ 0.62
    @test population_param_of_interest(ta)[3] == median(a) - median(b)
    @test population_param_of_interest(ta)[3] != hodgeslehmann(ta)

    # one-sided bounds
    @test confint(ExactMannWhitneyUTest(x, y); tail=:left)[2] == Inf
    @test confint(ExactMannWhitneyUTest(x, y); tail=:right)[1] == -Inf

    # a narrower interval for a lower confidence level
    lo95, hi95 = confint(ExactMannWhitneyUTest(x, y); level=0.95)
    lo80, hi80 = confint(ExactMannWhitneyUTest(x, y); level=0.80)
    @test lo95 <= lo80 <= hi80 <= hi95
end

@testset "method keyword" begin
    x = [1:10;]
    y = [2.1:2:21;]
    big1 = collect(1.0:40)
    big2 = collect(1.5:1:45)               # nx + ny > 50, untied: auto picks approximate
    @test MannWhitneyUTest(x, y) isa ExactMannWhitneyUTest
    @test MannWhitneyUTest(big1, big2) isa ApproximateMannWhitneyUTest
    @test MannWhitneyUTest(x, y; method=:auto) isa ExactMannWhitneyUTest
    @test MannWhitneyUTest(x, y; method=:approximate) isa ApproximateMannWhitneyUTest
    @test MannWhitneyUTest(big1, big2; method=:exact) isa ExactMannWhitneyUTest
    @test MannWhitneyUTest(big1, big2; method = s -> s.nx + s.ny <= 100 ? :exact : :approximate) isa
        ExactMannWhitneyUTest
    # a callable need not be a `Function`; the keyword documents "a callable"
    @test !(MannWhitneyExactBelow(100) isa Function)
    @test MannWhitneyUTest(big1, big2; method = MannWhitneyExactBelow(100)) isa
        ExactMannWhitneyUTest
    @test MannWhitneyUTest(big1, big2; method = MannWhitneyExactBelow(10)) isa
        ApproximateMannWhitneyUTest

    @test_throws ArgumentError MannWhitneyUTest(x, y; method=:fastest)
    @test_throws ArgumentError MannWhitneyUTest(x, y; method = s -> :fastest)
    @test_throws ArgumentError MannWhitneyUTest(x, y; method = 3)
    # a callable this cannot be called with is rejected with the same message rather
    # than raising a MethodError from inside the constructor
    @test_throws ArgumentError MannWhitneyUTest(x, y; method = (a, b) -> :exact)

    # the callable is handed the documented fields
    seen = Ref{Any}(nothing)
    MannWhitneyUTest([1.0, 2, 3], [3.0, 4, 5, 6]; method = s -> (seen[] = s; :exact))
    @test issetequal(keys(seen[]), (:nx, :ny, :ties, :tie_adjustment))
    @test (seen[].nx, seen[].ny) == (3, 4)
    @test seen[].ties
    MannWhitneyUTest([1.0, 2, 3], [4.0, 5, 6]; method = s -> (seen[] = s; :exact))
    @test !seen[].ties
end

@testset "Samples are carried through faithfully" begin
    # the interval and the estimator need the values, not just the ranks
    x = [3.0, 1.0, 4.0, 1.5, 5.0]
    y = [2.0, 7.0, 1.0, 8.0]
    for t in (MannWhitneyUTest(x, y),
              MannWhitneyUTest(x, y; method=:approximate))
        @test t.x == x
        @test t.y == y
        @test hodgeslehmann(t) ≈ median([xi - yj for xi in x, yj in y])
    end
    # eltype is preserved rather than widened to Float64
    t32 = MannWhitneyUTest(Float32[1:10;], Float32[2.5f0:2:21;])
    @test eltype(confint(t32)) === Float32
    @test hodgeslehmann(t32) isa Float32
    # mixed eltypes promote instead of erroring
    tm = MannWhitneyUTest([1, 2, 3, 4, 5, 6], [2.5, 3.5, 4.5, 5.5, 6.5, 7.5])
    @test eltype(confint(tm)) === Float64
end

@testset "Exact route past the automatic thresholds" begin
    # StatsFuns 2.2.1 fixed the Int overflow in wilcoxcdf (StatsFuns#220, #221). The
    # frontier used to be binomial(nx + ny, nx): (30, 35) was the largest valid split
    # and (31, 34) the smallest invalid one. Both are fine now, as is (40, 45), where
    # binomial(85, 40) is about 2.9e24 and the old code could not even form the count.
    for (nx, ny) in ((30, 35), (31, 34), (40, 45))
        t = MannWhitneyUTest(collect(1.0:nx), collect(1.5:1:(ny + 0.5)); method=:exact)
        @test 0 <= pvalue(t) <= 1
        lo, hi = confint(t)
        @test lo <= hodgeslehmann(t) <= hi
    end
    # the null distribution is symmetric about nx*ny/2
    for (nx, ny, u) in ((40, 45, 500), (40, 45, 900), (31, 34, 300))
        @test HypothesisTests.wilcoxcdf(nx, ny, u) +
              HypothesisTests.wilcoxcdf(nx, ny, nx * ny - u - 1) ≈ 1
    end
    # the tied branch will not enumerate without bound, though
    tied = MannWhitneyUTest(repeat([1.0, 2.0], 15), repeat([1.0, 3.0], 15); method=:exact)
    @test_throws ArgumentError pvalue(tied)
end

@testset "show past the pairwise bound" begin
    # `show` reads the interval off nx * ny cross-group differences, which `confint`
    # refuses to form past MAX_PAIRWISE_ESTIMATES. The name, the statistic and the
    # p-value are still worth printing, so it is the line that goes, not the display.
    nx = ny = 10_001  # nx * ny = 100_020_001
    t = MannWhitneyUTest(collect(1.0:nx), collect(1.5:1:(ny + 0.5)))
    @test nx * ny > HypothesisTests.MAX_PAIRWISE_ESTIMATES
    @test_throws ArgumentError confint(t)
    out = sprint(show, t)
    @test occursin("Approximate Mann-Whitney U test", out)
    @test occursin("two-sided p-value", out)
    @test !occursin("confidence interval", out)

    # and where it can afford one it still prints it
    @test occursin("95% confidence interval",
                   sprint(show, MannWhitneyUTest([1.0, 3, 5, 7], [2.0, 4, 6, 8])))
end

@testset "Reflection under argument swap" begin
    # Swapping the samples negates the shift the test is about, so every quantity
    # carrying that sign flips and every quantity that does not is unchanged. The
    # arithmetic is exact rather than approximate: IEEE subtraction gives
    # xi - yj == -(yj - xi), negation is exact, and the index rules see only nx, ny
    # and the tie adjustment, all three symmetric in the two samples. So `==` here,
    # not `≈`: an off-by-one in either index rule breaks it.
    #
    # The last two samples are separated on purpose. Where U sits near its null mean
    # the normal index rule cannot tell that mean from the centred statistic, and the
    # reflection holds either way; on a clearly shifted sample it does not, which is
    # what makes this testset detect the confusion its `confint` comment warns about.
    for (x, y) in (([1.0, 4, 6, 9, 11, 14], [2.0, 5, 7, 8, 12, 13]),   # no ties
                   ([1.0, 2, 2, 5, 5, 5, 8], [2.0, 3, 5, 5, 9, 9]),    # tied within and across
                   ([1.0, 2, 3], [10.0, 11, 12, 13]),                  # disjoint, unequal sizes
                   ([10.0, 11, 12, 13, 14, 15], [1.0, 2, 3, 4, 5, 6]), # shifted well clear
                   ([5.0, 5, 6, 8, 9, 9, 12], [1.0, 2, 2, 4, 6, 7]))   # shifted, and tied
        for T in (ExactMannWhitneyUTest, ApproximateMannWhitneyUTest)
            a, b = T(x, y), T(y, x)
            @test hodgeslehmann(a) == -hodgeslehmann(b)
            @test population_param_of_interest(a)[3] ==
                  -population_param_of_interest(b)[3]
            # U is complementary rather than antisymmetric: it counts pairs one way
            @test a.U + b.U == length(x) * length(y)
            @test pvalue(a) == pvalue(b)
            @test pvalue(a; tail=:left) == pvalue(b; tail=:right)
            # the two-sided interval reflects about zero, ends exchanged
            for level in (0.90, 0.95, 0.99)
                lo_a, hi_a = confint(a; level=level)
                lo_b, hi_b = confint(b; level=level)
                @test lo_a == -hi_b
                @test hi_a == -lo_b
            end
            # and a lower bound becomes the negated upper bound
            @test confint(a; tail=:left)[1] == -confint(b; tail=:right)[2]
        end
    end
end
end
