using HypothesisTests, Distributions, Base.Test

# One sample test
n = 1000
srand(1984948)

x = rand(Normal(), n)
t = OneSampleADTest(x, Normal())
@test_approx_eq_eps t.A² 0.2013 0.1^4
@test_approx_eq_eps pvalue(t) 0.8811 0.1^4

x = rand(DoubleExponential(), n)
t = OneSampleADTest(x, Normal())
@test_approx_eq_eps t.A² 10.7439 0.1^4
@test_approx_eq_eps pvalue(t) 0.0 0.1^4

x = rand(Cauchy(), n)
t = OneSampleADTest(x, Normal())
@test_approx_eq_eps t.A² 278.0190 0.1^4

x = rand(LogNormal(), n)
t = OneSampleADTest(x, Normal())
@test_approx_eq_eps t.A² 85.5217 0.1^4

# k-sample test
samples = Any[
    [38.7, 41.5, 43.8, 44.5, 45.5, 46.0, 47.7, 58.0],
    [39.2, 39.3, 39.7, 41.4, 41.8, 42.9, 43.3, 45.8],
    [34.0, 35.0, 39.0, 40.0, 43.0, 43.0, 44.0, 45.0],
    [34.0, 34.8, 34.8, 35.4, 37.2, 37.8, 41.2, 42.8]
]

t = KSampleADTest(samples...)
@test_approx_eq_eps t.A²k 8.3926 0.1^4
@test_approx_eq_eps t.σ 1.2038 0.1^4
@test_approx_eq_eps pvalue(t) 0.0020 0.1^4

t = KSampleADTest(samples..., modified = false)
@test_approx_eq_eps t.A²k 8.3559 0.1^4
@test_approx_eq_eps t.σ 1.2038 0.1^4
@test_approx_eq_eps pvalue(t) 0.0021 0.1^4

srand(31412455)
samples = Any[rand(Normal(), 50), rand(Normal(0.5), 30)]
t = KSampleADTest(samples...)
@test pvalue(t) < 0.05

samples = Any[rand(Normal(), 50), rand(Normal(), 30), rand(Normal(), 20)]
t = KSampleADTest(samples...)
@test pvalue(t) > 0.05
