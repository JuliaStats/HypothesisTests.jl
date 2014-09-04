using HypothesisTests, Base.Test

# http://en.wikipedia.org/wiki/Fisher%27s_exact_test
t = HypothesisTests.FisherExactTest(1, 9, 11, 3)
@test_approx_eq HypothesisTests.pvalue(t; tail=:both) 0.002759456185220088
@test_approx_eq HypothesisTests.pvalue(t; tail=:left) 0.001379728092610043
@test_approx_eq HypothesisTests.pvalue(t; tail=:right) 0.9999663480953023
show(IOBuffer(), t)

# http://www.physics.csbsju.edu/stats/exact.html
t = HypothesisTests.FisherExactTest(7, 12, 0, 5)
@test_approx_eq HypothesisTests.pvalue(t; tail=:both) 0.2720685111989459
@test_approx_eq HypothesisTests.pvalue(t; tail=:left) 1.0
@test_approx_eq HypothesisTests.pvalue(t; tail=:right) 0.1455862977602108
show(IOBuffer(), t)
t = HypothesisTests.FisherExactTest(12, 7, 5, 0)
@test_approx_eq HypothesisTests.pvalue(t; tail=:both) 0.2720685111989459
@test_approx_eq HypothesisTests.pvalue(t; tail=:left) 0.1455862977602108
@test_approx_eq HypothesisTests.pvalue(t; tail=:right) 1.0
show(IOBuffer(), t)
t = HypothesisTests.FisherExactTest(0, 5, 7, 12)
@test_approx_eq HypothesisTests.pvalue(t; tail=:both) 0.2720685111989459
@test_approx_eq HypothesisTests.pvalue(t; tail=:left) 0.1455862977602108
@test_approx_eq HypothesisTests.pvalue(t; tail=:right) 1.0
show(IOBuffer(), t)
t = HypothesisTests.FisherExactTest(5, 0, 12, 7)
@test_approx_eq HypothesisTests.pvalue(t; tail=:both) 0.2720685111989459
@test_approx_eq HypothesisTests.pvalue(t; tail=:left) 1.0
@test_approx_eq HypothesisTests.pvalue(t; tail=:right) 0.1455862977602108
show(IOBuffer(), t)

# http://www.stata.com/support/faqs/statistics/fishers-exact-test/
t = HypothesisTests.FisherExactTest(2, 31, 136, 15532)
@test_approx_eq HypothesisTests.pvalue(t; tail=:both) 0.03390271476034175
@test_approx_eq HypothesisTests.pvalue(t; tail=:left) 0.9970112864705307
@test_approx_eq HypothesisTests.pvalue(t; tail=:right) 0.03390271476034175
show(IOBuffer(), t)

t = HypothesisTests.FisherExactTest(1, 1, 1, 1)
@test HypothesisTests.pvalue(t, tail=:both) <= 1
show(IOBuffer(), t)
