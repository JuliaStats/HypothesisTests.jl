using HypothesisTests, Distributions, Test, Random
using HypothesisTests: default_tail
using StableRNGs

@testset "Anderson-Darling" begin
@testset "One sample test" begin
    n = 1000
    rng = StableRNG(1984948)

    d = Normal(1000, 100)
    x = rand(rng, d, n)
    t = OneSampleADTest(x, d)
    @test t.A² ≈ 1.2960 atol=0.1^4
    @test pvalue(t) ≈ 0.2336 atol=0.1^4
    @test default_tail(t) == :right

    d = DoubleExponential()
    x = rand(rng, d, n)
    t = OneSampleADTest(x, Normal(mean(d), std(d)))
    @test t.A² ≈ 11.0704 atol=0.1^4
    @test pvalue(t) ≈ 0.0 atol=0.1^4
    t = OneSampleADTest(x, d)
    @test t.A² ≈ 0.8968 atol=0.1^4
    @test pvalue(t) ≈ 0.4162 atol=0.1^4

    d = Cauchy()
    x = rand(rng, Cauchy(), n)
    t = OneSampleADTest(x, Normal())
    @test pvalue(t) ≈ 0.0 atol=0.1^4
    t = OneSampleADTest(x, d)
    @test pvalue(t) ≈ 0.9640 atol=0.1^4

    d = LogNormal()
    x = rand(rng, d, n)
    t = OneSampleADTest(x, Normal(mean(d), std(d)))
    @test pvalue(t) ≈ 0.0 atol=0.1^4
    t = OneSampleADTest(x, d)
    @test pvalue(t) ≈ 0.9123 atol=0.1^4

    d = Uniform(-pi, 2pi)
    x = rand(rng, d, n)
    t = OneSampleADTest(x, d)
    @test pvalue(t) ≈ 0.2337 atol=0.1^4

    x = rand(rng, Uniform(0, 1.8), n)
    t = OneSampleADTest(x, Uniform())
    @test pvalue(t) ≈ 0.0 atol=0.1^4

    x = rand(rng, Exponential(), n)
    t = OneSampleADTest(x, Exponential())
    @test pvalue(t) ≈ 0.8579 atol=0.1^4
end

@testset "k-sample test" begin
    samples = Any[
        [38.7, 41.5, 43.8, 44.5, 45.5, 46.0, 47.7, 58.0],
        [39.2, 39.3, 39.7, 41.4, 41.8, 42.9, 43.3, 45.8],
        [34.0, 35.0, 39.0, 40.0, 43.0, 43.0, 44.0, 45.0],
        [34.0, 34.8, 34.8, 35.4, 37.2, 37.8, 41.2, 42.8]
    ]

    t = KSampleADTest(samples...)
    @test t.A²k ≈ 8.3926 atol=0.1^4
    @test t.σ ≈ 1.2038 atol=0.1^4
    @test pvalue(t) ≈ 0.0022 atol=0.1^4
    @test default_tail(t) == :right

    ts = KSampleADTest(samples..., nsim = 20000);
    @test pvalue(ts) ≈ 0.00155 atol=0.1^3

    t = KSampleADTest(samples..., modified = false)
    @test t.A²k ≈ 8.3559 atol=0.1^4
    @test t.σ ≈ 1.2038 atol=0.1^4
    @test pvalue(t) ≈ 0.0023 atol=0.1^4

    ts = KSampleADTest(samples..., modified = false, nsim = 20000);
    @test pvalue(ts) ≈ 0.00150 atol=0.1^3
end

@testset "more tests" begin
    rng =  StableRNG(31412455)
    samples = Any[rand(rng, Normal(), 50), rand(rng, Normal(0.5), 30)]
    t = KSampleADTest(samples...)
    @test pvalue(t) < 0.05

    samples = Any[rand(rng, Normal(), 50), rand(rng, Normal(), 30), rand(rng, Normal(), 20)]
    t = KSampleADTest(samples...)
    @test pvalue(t) > 0.05

    @test pvalue(OneSampleADTest(vcat(rand(rng, Normal(),500), rand(rng, Beta(2,2),500)), Beta(2,2))) ≈ 0.0 atol=0.1^4

    n = 1000
    x = rand(rng, Exponential(), n)
    @test pvalue(KSampleADTest(rand(rng, 1000), randn(rng, 1000))) ≈ 0.0 atol=eps()
    @test pvalue(KSampleADTest(x,x,x,x,x,x)) ≈ 1.0
end

end
