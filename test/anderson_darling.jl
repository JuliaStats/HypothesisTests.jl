using HypothesisTests, Distributions, Base.Test

n = 1000
srand(1984948)

x = rand(Normal(), n)
t = OneSampleADTest(x, Normal())
@test_approx_eq t.A² 0.2012884979634464

x = rand(DoubleExponential(), n)
t = OneSampleADTest(x, Normal())
@test_approx_eq t.A² 10.743875548302867

x = rand(Cauchy(), n)
t = OneSampleADTest(x, Normal())
@test_approx_eq t.A² 278.0189969882192

x = rand(LogNormal(), n)
t = OneSampleADTest(x, Normal())
@test_approx_eq t.A² 85.52170553455244
