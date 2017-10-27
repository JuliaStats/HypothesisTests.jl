using HypothesisTests, Distributions, Base.Test
using HypothesisTests: default_tail

# One sample test
n = 1000
srand(1984948)

x = rand(Normal(), n)
t = OneSampleADTest(x, Normal())
@test isapprox(t.A², 0.2013, atol=0.1^4)
@test isapprox(pvalue(t), 0.99, atol=0.1^4)
@test default_tail(t) == :right

x = rand(DoubleExponential(), n)
t = OneSampleADTest(x, Normal())
@test isapprox(t.A², 10.7439, atol=0.1^4)
@test isapprox(pvalue(t), 0.0, atol=0.1^4)

x = rand(Cauchy(), n)
t = OneSampleADTest(x, Normal())
@test isapprox(pvalue(t), 0.0, atol=0.1^4)

x = rand(LogNormal(), n)
t = OneSampleADTest(x, Normal())
@test isapprox(pvalue(t), 0.0, atol=0.1^4)

x = rand(Uniform(), n)
t = OneSampleADTest(x, Uniform())
@test isapprox(pvalue(t), 0.145, atol=0.1^4)

x = rand(Uniform(0, 1.8), n)
t = OneSampleADTest(x, Uniform())
@test isapprox(pvalue(t), 0.5112, atol=0.1^4)

x = rand(Exponential(), n)
t = OneSampleADTest(x, Exponential())
@test isapprox(pvalue(t), 0.9427, atol=0.1^4)

# k-sample test
srand(948574875)
samples = Any[
    [38.7, 41.5, 43.8, 44.5, 45.5, 46.0, 47.7, 58.0],
    [39.2, 39.3, 39.7, 41.4, 41.8, 42.9, 43.3, 45.8],
    [34.0, 35.0, 39.0, 40.0, 43.0, 43.0, 44.0, 45.0],
    [34.0, 34.8, 34.8, 35.4, 37.2, 37.8, 41.2, 42.8]
]

t = KSampleADTest(samples...)
@test isapprox(t.A²k, 8.3926, atol=0.1^4)
@test isapprox(t.σ, 1.2038, atol=0.1^4)
@test isapprox(pvalue(t), 0.0022, atol=0.1^4)
@test default_tail(t) == :right

ts = KSampleADTest(samples..., nsim = 20000);
@test isapprox(pvalue(ts), 0.00155, atol=0.1^3)

t = KSampleADTest(samples..., modified = false)
@test isapprox(t.A²k, 8.3559, atol=0.1^4)
@test isapprox(t.σ, 1.2038, atol=0.1^4)
@test isapprox(pvalue(t), 0.0023, atol=0.1^4)

ts = KSampleADTest(samples..., modified = false, nsim = 20000);
@test isapprox(pvalue(ts), 0.00150, atol=0.1^3)

# more tests
srand(31412455)
samples = Any[rand(Normal(), 50), rand(Normal(0.5), 30)]
t = KSampleADTest(samples...)
@test pvalue(t) < 0.05

samples = Any[rand(Normal(), 50), rand(Normal(), 30), rand(Normal(), 20)]
t = KSampleADTest(samples...)
@test pvalue(t) > 0.05

@test isapprox(pvalue(OneSampleADTest(vcat(rand(Normal(),500), rand(Beta(2,2),500)), Beta(2,2))), 0.0, atol=0.1^4)
