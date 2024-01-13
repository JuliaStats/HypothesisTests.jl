using HypothesisTests, LinearAlgebra, Test, Random
using StableRNGs

@testset "Shapiro-Wilk" begin
    @testset "shapiro_wilk_coefs" begin
        @test HypothesisTests.shapiro_wilk_coefs(3) == [sqrt(2.0) / 2.0, 0.0, -sqrt(2.0) / 2.0]

        swc = HypothesisTests.shapiro_wilk_coefs(10)
        @test swc[4] == -swc[7]
        @test swc[2] == -swc[9]

        swc = HypothesisTests.shapiro_wilk_coefs(11)

        @test swc[6] == 0.0
        @test swc[5] == -swc[7]
        @test swc[3] == -swc[9]
        @test swc[1] == -swc[11]

        #anti-symmetry
        swc = HypothesisTests.shapiro_wilk_coefs(20)
        @test all([swc[i] == -swc[end - i + 1] for i in eachindex(swc)])

        # Values obtained by calling `_swilkfort` fortran subroutine directly.

        swc10 = HypothesisTests.shapiro_wilk_coefs(10)
        res = swc10[1:5] .- [0.573715, 0.32897, 0.214349, 0.122791, 0.0400887]
        @test norm(res, 1) ≈ 0.0 atol = length(swc10) * eps(Float32)

        swc20 = HypothesisTests.shapiro_wilk_coefs(20)
        res = swc20[1:10] .-
              [0.473371, 0.32174, 0.255663, 0.208297, 0.16864, 0.133584, 0.101474,
               0.0712893, 0.0423232, 0.0140351]
        @test norm(res, 1) ≈ 0.0 atol = length(swc20) * eps(Float32)

        rng = StableRNG(0x5bca7c69b794f8ce)
        X = sort(rand(rng, Float64, 10))
        # W, pval, _ = swilkfort(X)
        W = 0.9088434774710951
        # pval = 0.2731410626084226
        @test HypothesisTests.unsafe_swstat(X, swc10) ≈ W atol = eps(Float32)
    end

    @testset "Shapiro-Wilk" begin
        # syntactic tests

        @test_throws ArgumentError ShapiroWilkTest([1, 2])
        @test_throws ArgumentError("at least 3 samples are required, got 2") ShapiroWilkTest([1, 2], HypothesisTests.shapiro_wilk_coefs(3))
        @test_throws ArgumentError ShapiroWilkTest([1, 2, 3], censored=4)
        @test_throws DimensionMismatch ShapiroWilkTest([1, 2, 3],
                                                       HypothesisTests.shapiro_wilk_coefs(4))

        @test_throws ArgumentError("sample doesn't seem to be sorted or is constant (up to numerical accuracy)") ShapiroWilkTest([1,1,1])
        @test_throws ArgumentError("sample is constant (up to numerical accuracy)") ShapiroWilkTest([1,1,1], sorted=false)

        t = ShapiroWilkTest([1, 2, 3])
        @test t.W == 1.0
        @test pvalue(t) == 1.0

        @test_throws "censored samples not implemented yet" pvalue(ShapiroWilkTest(1:4, censored=1))

        str = sprint(show, t)
        @test str ==
              """Shapiro-Wilk normality test
              ---------------------------
              Population details:
                  parameter of interest:   Squared correlation of data and expected order statistics of N(0,1) (W)
                  value under h_0:         1.0
                  point estimate:          1.0

              Test summary:
                  outcome with 95% confidence: fail to reject h_0
                  one-sided p-value:           1.0000

              Details:
                  number of observations: 3
                  censored ratio:         0.0
                  W-statistic:            1.0
              """

        # testing different cases of N
        for N in (3, 5, 11, 12)
            rng = StableRNG(0x5bca7c69b794f8ce)
            X = sort(randn(rng, N))
            t = ShapiroWilkTest(X; sorted=true)

            # analytic properties from Shapiro-Wilk 1965:
            # Lemma 1: Scale and origin invariance:
            scale, shift = rand(rng, 2)
            @test t.W ≈ ShapiroWilkTest(X .+ shift).W
            @test t.W ≈ ShapiroWilkTest(scale .* X).W
            @test t.W ≈ ShapiroWilkTest(scale .* X .+ shift).W
            # Lemma 2, 3: upper and lower bounds
            @test N * t.coefs[1]^2 / (N - 1) ≤ t.W ≤ 1.0

            # test the computation of pvalue in those cases
            @test pvalue(t) ≥ 0.05
        end

        @testset "Worked Example" begin
            # **Worked Example** (Section 4) from
            # PATRICK ROYSTON Approximating the Shapiro-Wilk W-test for non-normality
            # *Statistics and Computing* (1992) **2**, 117-119

            X = [48.4, 49.0, 59.5, 59.6, 60.7, 88.8, 98.2, 109.4, 169.1, 227.1]
            swc = HypothesisTests.shapiro_wilk_coefs(length(X))
            @test norm(swc[1:5] .- [0.5737, 0.3290, 0.2143, 0.1228, 0.0401], Inf) < 5.0e-5
            W = HypothesisTests.unsafe_swstat(X, swc)
            @test W ≈ 0.8078 atol = 2.9e-5

            t = ShapiroWilkTest(X)
            @test t.W == W
            @test pvalue(t) ≈ 0.018 atol = 4.7e-5
            @test pvalue(t) ≈ pvalue(ShapiroWilkTest(X, sorted=true))
            @test iszero(HypothesisTests.censored_ratio(t))
            @test length(t.coefs) == length(X)

            # test for un-sorted sample
            X2 = X[[9,8,2,3,4,5,7,10,1,6]]
            t2 = ShapiroWilkTest(X2)
            @test_throws ArgumentError ShapiroWilkTest(X2, sorted=true)
            @test t2.W ≈ t.W
            @test pvalue(t2) ≈ pvalue(t)
            X3 = X[[2,8,9,3,4,5,7,10,1,6]]
            t3 = ShapiroWilkTest(X3)
            @test t3.W ≈ t.W
            @test pvalue(t3) ≈ pvalue(t)
            @test pvalue(t3) ≠ pvalue(ShapiroWilkTest(X3, sorted=true))
        end
    end
end
