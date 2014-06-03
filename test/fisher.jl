using HypothesisTests, Base.Test

# http://en.wikipedia.org/wiki/Fisher%27s_exact_test
t = HypothesisTests.FisherExactTest(1, 9, 11, 3)
@test_approx_eq HypothesisTests.pvalue(t; method=:minlike, tail=:both) 0.002759456185220088
@test_approx_eq HypothesisTests.pvalue(t; method=:central, tail=:both) 0.0027594561852200853
@test_approx_eq HypothesisTests.pvalue(t; tail=:left) 0.001379728092610043
@test_approx_eq HypothesisTests.pvalue(t; tail=:right) 0.9999663480953023
show(IOBuffer(), t)

# http://www.physics.csbsju.edu/stats/exact.html
t = HypothesisTests.FisherExactTest(7, 12, 0, 5)
@test_approx_eq HypothesisTests.pvalue(t; method=:minlike, tail=:both) 0.2720685111989459
@test_approx_eq HypothesisTests.pvalue(t; tail=:left) 1.0
@test_approx_eq HypothesisTests.pvalue(t; tail=:right) 0.1455862977602108
show(IOBuffer(), t)
t = HypothesisTests.FisherExactTest(12, 7, 5, 0)
@test_approx_eq HypothesisTests.pvalue(t; method=:minlike, tail=:both) 0.2720685111989459
@test_approx_eq HypothesisTests.pvalue(t; tail=:left) 0.1455862977602108
@test_approx_eq HypothesisTests.pvalue(t; tail=:right) 1.0
show(IOBuffer(), t)
t = HypothesisTests.FisherExactTest(0, 5, 7, 12)
@test_approx_eq HypothesisTests.pvalue(t; method=:minlike, tail=:both) 0.2720685111989459
@test_approx_eq HypothesisTests.pvalue(t; tail=:left) 0.1455862977602108
@test_approx_eq HypothesisTests.pvalue(t; tail=:right) 1.0
show(IOBuffer(), t)
t = HypothesisTests.FisherExactTest(5, 0, 12, 7)
@test_approx_eq HypothesisTests.pvalue(t; method=:minlike, tail=:both) 0.2720685111989459
@test_approx_eq HypothesisTests.pvalue(t; tail=:left) 1.0
@test_approx_eq HypothesisTests.pvalue(t; tail=:right) 0.1455862977602108
show(IOBuffer(), t)

# http://www.stata.com/support/faqs/statistics/fishers-exact-test/
t = HypothesisTests.FisherExactTest(2, 31, 136, 15532)
@test_approx_eq HypothesisTests.pvalue(t; method=:minlike, tail=:both) 0.03390271476034175
@test_approx_eq HypothesisTests.pvalue(t; tail=:left) 0.9970112864705307
@test_approx_eq HypothesisTests.pvalue(t; tail=:right) 0.03390271476034175
show(IOBuffer(), t)

# http://www.utstat.toronto.edu/~brunner/oldclass/312f12/lectures/312f12FisherWithR.pdf
# http://vassarstats.net/odds2x2.html
t = HypothesisTests.FisherExactTest(4, 1, 20, 1)
#@test_approx_eq [ci(t)...] [0.0102,3.9081]
#@test_approx_eq [ci(t)...] [0.002439905, 19.594803004]
@test_approx_eq pvalue(t; method=:minlike) 0.353846153846154
@test_approx_eq pvalue(t; tail=:left) 0.353846153846154
@test_approx_eq pvalue(t; tail=:right) 0.9692307692307692
