using HypothesisTests, Distributions, Base.Test

n = 1000
srand(1984948)

x = rand(Normal(), n)
t = OneSampleADTest(x, Normal())
@test_approx_eq_eps t.A² 0.2013 0.1^4

x = rand(DoubleExponential(), n)
t = OneSampleADTest(x, Normal())
@test_approx_eq_eps t.A² 10.7439 0.1^4

x = rand(Cauchy(), n)
t = OneSampleADTest(x, Normal())
@test_approx_eq_eps t.A² 278.0190 0.1^4

x = rand(LogNormal(), n)
t = OneSampleADTest(x, Normal())
@test_approx_eq_eps t.A² 85.5217 0.1^4
