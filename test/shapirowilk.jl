using HypothesisTests, LinearAlgebra, Test, Random
using StableRNGs

@testset "ShapiroWilk" begin
    @testset "ShapiroWilkCoefs" begin
        @test HypothesisTests.ShapiroWilkCoefs(3).A == [sqrt(2.0) / 2.0]
        @test length(HypothesisTests.ShapiroWilkCoefs(3)) == 3

        swc = HypothesisTests.ShapiroWilkCoefs(10)
        @test length(swc) == 10
        @test lastindex(swc) == 10
        @test length(swc.A) == 5

        @test swc.A[4] == swc[4] == -swc[7]
        @test swc.A[2] == swc[2] == -swc[9]

        swc = HypothesisTests.ShapiroWilkCoefs(11)
        @test length(swc) == 11
        @test lastindex(swc) == 11
        @test length(swc.A) == 5

        @test swc[6] == 0.0
        @test swc.A[5] == swc[5] == -swc[7]
        @test swc.A[3] == swc[3] == -swc[9]
        @test swc.A[1] == swc[1] == -swc[11]

        #anti-symmetry
        swc = HypothesisTests.ShapiroWilkCoefs(20)
        @test all([swc[i] == -swc[end-i+1] for i in eachindex(swc)])

        # Values obtained by calling `_swilkfort` fortran subroutine directly.

        swc10 = HypothesisTests.ShapiroWilkCoefs(10)
        res = swc10.A .- [0.573715, 0.32897, 0.214349, 0.122791, 0.0400887]
        @test norm(res, 1) ≈ 0.0 atol = length(swc10) * eps(Float32)

        swc20 = HypothesisTests.ShapiroWilkCoefs(20)
        res = swc20.A .- [0.473371, 0.32174, 0.255663, 0.208297, 0.16864, 0.133584, 0.101474, 0.0712893, 0.0423232, 0.0140351]
        @test norm(res, 1) ≈ 0.0 atol = length(swc20) * eps(Float32)

        rng = StableRNG(0x5bca7c69b794f8ce)
        X = sort(rand(rng, Float64, 10))
        # W, pval, _ = swilkfort(X)
        W = 0.9088434774710951
        # pval = 0.2731410626084226
        @test HypothesisTests.swstat(X, swc10) ≈ W atol = eps(Float32)
    end

    @testset "ShapiroWilk" begin
        # syntactic tests
        @test_throws ArgumentError ShapiroWilkTest([1, 2])
        @test_throws ArgumentError ShapiroWilkTest([1, 2, 3], lower_uncensored=4)
        @test_throws DimensionMismatch ShapiroWilkTest([1, 2, 3], HypothesisTests.ShapiroWilkCoefs(4))

        t = ShapiroWilkTest([1, 2, 3])
        @test t.W == 1.0
        @test HypothesisTests.pvalue(t) == 1.0

        str = sprint(show, t)
        @test occursin("parameter of interest:   Squared correlation of sorted data and the uncorrelated expected order statistics of the normal distribution (W)", str)
        @test occursin("fail to reject h_0", str)
        @test occursin("number of observations: 3", str)
        @test occursin("censored ratio:         0.0", str)
        @test occursin("W-statistic:            1.0", str)

        # testing different cases of N
        for N in (3, 5, 11, 12)
            rng = StableRNG(0x5bca7c69b794f8ce)
            X = sort(randn(rng, N))
            t = ShapiroWilkTest(X, sorted=true)

            # analytic properties from Shapiro-Wilk 1965:
            # Lemma 1: Scale and origin invariance:
            scale, shift = rand(rng, 2)
            @test t.W ≈ ShapiroWilkTest(X .+ shift).W
            @test t.W ≈ ShapiroWilkTest(scale .* X).W
            @test t.W ≈ ShapiroWilkTest(scale .* X .+ shift).W
            # Lemma 2, 3: upper and lower bounds
            @test N * t.coefs[1]^2 / (N - 1) ≤ t.W ≤ 1.0
        end

        @testset "Worked Example" begin
            # **Worked Example** (Section 4) from
            # PATRICK ROYSTON Approximating the Shapiro-Wilk W-test for non-normality
            # *Statistics and Computing* (1992) **2**, 117-119

            X = [48.4, 49.0, 59.5, 59.6, 60.7, 88.8, 98.2, 109.4, 169.1, 227.1]
            swc = HypothesisTests.ShapiroWilkCoefs(length(X))
            @test norm(swc.A .- [0.5737, 0.3290, 0.2143, 0.1228, 0.0401], Inf) < 5.0e-5
            W = HypothesisTests.swstat(X, swc)
            @test W ≈ 0.8078 atol = 2.9e-5

            t = ShapiroWilkTest(X)
            @test t.W == W
            @test HypothesisTests.pvalue(t) ≈ 0.018 atol = 4.7e-5

            @test iszero(HypothesisTests.censored_ratio(t))
            @test length(t.coefs) == length(X)
        end
    end
end

