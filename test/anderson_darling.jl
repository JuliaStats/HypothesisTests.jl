using HypothesisTests, Distributions, Test, Random
using HypothesisTests: default_tail

@testset "Anderson-Darling" begin
@testset "One sample test" begin
    n = 1000
    Random.seed!(1984948)

    d = Normal(1000, 100)
    x = rand(d, n)
    t = OneSampleADTest(x, d)
    @test t.A² ≈ 0.3078 atol=0.1^4
    @test pvalue(t) ≈ 0.9321 atol=0.1^4
    @test default_tail(t) == :right

    d = DoubleExponential()
    x = rand(d, n)
    t = OneSampleADTest(x, Normal(mean(d), std(d)))
    @test t.A² ≈ 10.9678 atol=0.1^4
    @test pvalue(t) ≈ 0.0 atol=0.1^4
    t = OneSampleADTest(x, d)
    @test t.A² ≈ 0.21942 atol=0.1^4
    @test pvalue(t) ≈ 0.9842 atol=0.1^4

    d = Cauchy()
    x = rand(Cauchy(), n)
    t = OneSampleADTest(x, Normal())
    @test pvalue(t) ≈ 0.0 atol=0.1^4
    t = OneSampleADTest(x, d)
    @test pvalue(t) ≈ 0.4228 atol=0.1^4

    d = LogNormal()
    x = rand(d, n)
    t = OneSampleADTest(x, Normal(mean(d), std(d)))
    @test pvalue(t) ≈ 0.0 atol=0.1^4
    t = OneSampleADTest(x, d)
    @test pvalue(t) ≈ 0.1569 atol=0.1^4

    d = Uniform(-pi, 2pi)
    x = rand(d, n)
    t = OneSampleADTest(x, d)
    @test pvalue(t) ≈ 0.5489 atol=0.1^4

    x = rand(Uniform(0, 1.8), n)
    t = OneSampleADTest(x, Uniform())
    @test pvalue(t) ≈ 0.0 atol=0.1^4

    x = rand(Exponential(), n)
    t = OneSampleADTest(x, Exponential())
    @test pvalue(t) ≈ 0.1153 atol=0.1^4
end

@testset "k-sample test" begin
    Random.seed!(948574875)
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
    Random.seed!(31412455)
    samples = Any[rand(Normal(), 50), rand(Normal(0.5), 30)]
    t = KSampleADTest(samples...)
    @test pvalue(t) < 0.05

    samples = Any[rand(Normal(), 50), rand(Normal(), 30), rand(Normal(), 20)]
    t = KSampleADTest(samples...)
    @test pvalue(t) > 0.05

    @test pvalue(OneSampleADTest(vcat(rand(Normal(),500), rand(Beta(2,2),500)), Beta(2,2))) ≈ 0.0 atol=0.1^4

    n = 1000
    x = rand(Exponential(), n)
    @test pvalue(KSampleADTest(rand(1000),randn(1000))) ≈ 0.0 atol=eps()
    @test pvalue(KSampleADTest(x,x,x,x,x,x)) ≈ 1.0
end

end
