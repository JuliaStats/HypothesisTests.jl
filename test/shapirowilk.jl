using HypothesisTests, LinearAlgebra, Test, Random
using StableRNGs
@testset "ShapiroWilk" begin
    @testset "SWcoeffs" begin
        @test HypothesisTests.SWCoeffs(3).A == [sqrt(2.0) / 2.0]
        @test length(HypothesisTests.SWCoeffs(3)) == 3

        S = HypothesisTests.SWCoeffs(10)
        @test length(S) == 10
        @test lastindex(S) == 10
        @test length(S.A) == 5

        @test S.A[4] == S[4] == -S[7]
        @test S.A[2] == S[2] == -S[9]

        S = HypothesisTests.SWCoeffs(11)
        @test length(S) == 11
        @test lastindex(S) == 11
        @test length(S.A) == 5

        @test S[6] == 0.0
        @test S.A[5] == S[5] == -S[7]
        @test S.A[3] == S[3] == -S[9]
        @test S.A[1] == S[1] == -S[11]

        #anti-symmetry
        S = HypothesisTests.SWCoeffs(20)
        @test all([S[i] == -S[end-i+1] for i in eachindex(S)])

        # Values obtained by calling `_swilkfort` fortran subroutine directly.

        SWc10 = HypothesisTests.SWCoeffs(10)
        res = SWc10.A .- [0.573715, 0.32897, 0.214349, 0.122791, 0.0400887]
        @test norm(res, 1) ≈ 0.0 atol = length(SWc10) * eps(Float32)

        SWc20 = HypothesisTests.SWCoeffs(20)
        res = SWc20.A .- [0.473371, 0.32174, 0.255663, 0.208297, 0.16864, 0.133584, 0.101474, 0.0712893, 0.0423232, 0.0140351]
        @test norm(res, 1) ≈ 0.0 atol = length(SWc20) * eps(Float32)

        rng = StableRNG(0x5bca7c69b794f8ce)
        X = sort(rand(rng, Float64, 10))
        # W, pval, _ = swilkfort(X)
        W = 0.9088434774710951
        # pval = 0.2731410626084226
        @test HypothesisTests.swstat(X, SWc10) ≈ W atol = eps(Float32)
    end

    @testset "ShapiroWilk" begin
        # syntactic tests
        @test_throws ArgumentError ShapiroWilkTest([1, 2])
        @test_throws ArgumentError ShapiroWilkTest([1, 2, 3], N1=4)
        @test_throws DimensionMismatch ShapiroWilkTest([1, 2, 3], HypothesisTests.SWCoeffs(4))

        t = ShapiroWilkTest([1, 2, 3])
        @test t.W == 1.0 # W is at most one but numerics makes it nextfloat(1.0)
        @test pvalue(t) == 1.0
        @test sprint(show, t) isa String

        # testing different cases of N
        for N in (3, 11, 1000)
            rng = StableRNG(0x5bca7c69b794f8ce)
            X = sort(randn(rng, N))
            t = ShapiroWilkTest(X, sample_sorted=true)

            # analytic properties from Shapiro-Wilk 1965:
            # Lemma 1: Scale and origin invariance:
            scale, shift = rand(rng, 2)
            @test t.W ≈ ShapiroWilkTest(X .+ shift).W
            @test t.W ≈ ShapiroWilkTest(scale .* X).W
            @test t.W ≈ ShapiroWilkTest(scale .* X .+ shift).W
            # Lemma 2, 3: upper and lower bounds
            @test N * t.SWc[1]^2 / (N - 1) ≤ t.W ≤ 1.0
        end

        @testset "Worked Example" begin
            # **Worked Example** (Section 4) from
            # PATRICK ROYSTON Approximating the Shapiro-Wilk W-test for non-normality
            # *Statistics and Computing* (1992) **2**, 117-119

            X = [48.4, 49.0, 59.5, 59.6, 60.7, 88.8, 98.2, 109.4, 169.1, 227.1]
            SWc = HypothesisTests.SWCoeffs(length(X))
            @test norm(SWc.A .- [0.5737, 0.3290, 0.2143, 0.1228, 0.0401], Inf) < 5.0e-5
            W = HypothesisTests.swstat(X, SWc)
            @test W ≈ 0.8078 atol = 2.9e-5

            t = ShapiroWilkTest(X)
            @test t.W == W
            @test pvalue(t) ≈ 0.018 atol = 4.7e-5

            @test t.N1 == length(X)
            @test length(t.SWc) == length(X)
        end
    end
end

