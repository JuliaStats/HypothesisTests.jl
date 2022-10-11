using HypothesisTests, Test
using HypothesisTests: default_tail

@testset "Fisher" begin
t = @inferred(HypothesisTests.FisherExactTest(1, 1, 1, 1))
@test t.ω ≈ 1.0
@test pvalue(t; tail=:left) ≈ 0.8333333333333337
@test pvalue(t; tail=:right) ≈ 0.8333333333333337
@test pvalue(t; method=:central) ≈ 1.0
@test pvalue(t; method=:minlike) ≈ 1.0
@test default_tail(t) == :both
@test_ci_approx confint(t; tail=:left) (0.0, 76.24918299781056)
@test_ci_approx confint(t; tail=:right) (0.013114894621608135, Inf)
@test_ci_approx confint(t; method=:central) (0.006400016357911029, 156.2496006379585)
#@test_approx_eq [confint(t; method=:minlike)...] [0.0131, 76.2492]
show(IOBuffer(), t)

# http://en.wikipedia.org/wiki/Fisher%27s_exact_test
t = HypothesisTests.FisherExactTest(1, 9, 11, 3)
@test t.ω ≈ 0.03720908483238119
@test pvalue(t; tail=:left) ≈ 0.0013797280926100427
@test pvalue(t; tail=:right) ≈ 0.9999663480953023
@test pvalue(t; method=:central) ≈ 0.0027594561852200853
@test pvalue(t; method=:minlike) ≈ 0.002759456185220088
@test_ci_approx confint(t; tail=:left) (0.0, 0.32600296913224913)
@test_ci_approx confint(t; tail=:right) (0.0012948958389639856, Inf)
@test_ci_approx confint(t; method=:central) (0.0006360029488751071, 0.42586647569637387)
#@test_approx_eq [confint(t; method=:minlike)...] [0.0013, 0.3567]
show(IOBuffer(), t)

# http://www.physics.csbsju.edu/stats/exact.html
t = HypothesisTests.FisherExactTest(7, 12, 0, 5)
@test t.ω ≈ Inf
@test pvalue(t; tail=:left) ≈ 1.0
@test pvalue(t; tail=:right) ≈ 0.1455862977602108
@test pvalue(t; method=:central, tail=:both) ≈ 0.2911725955204216
@test pvalue(t; method=:minlike, tail=:both) ≈ 0.2720685111989459
@test_ci_approx confint(t; tail=:left) (0.0, Inf)
@test_ci_approx confint(t; tail=:right) (0.539859556284207, Inf)
@test_ci_approx confint(t; method=:central) (0.39239937500428096, Inf)
#@test_approx_eq [confint(t; method=:minlike)...] [0.5171, Inf]
show(IOBuffer(), t)

t = HypothesisTests.FisherExactTest(12, 7, 5, 0)
@test t.ω ≈ 0.0
@test pvalue(t; tail=:left) ≈ 0.1455862977602108
@test pvalue(t; tail=:right) ≈ 1.0
@test pvalue(t; method=:central, tail=:both) ≈ 0.29117259552042146
@test pvalue(t; method=:minlike, tail=:both) ≈ 0.2720685111989459
@test_ci_approx confint(t; tail=:left) (0.0, 1.852333608546057)
@test_ci_approx confint(t; tail=:right) (0.0, Inf)
@test_ci_approx confint(t; method=:central) (0.0, 2.5484240386190433)
#@test_approx_eq [confint(t; method=:minlike)...] [0.0, 1.9338]
show(IOBuffer(), t)

t = HypothesisTests.FisherExactTest(0, 5, 7, 12)
@test t.ω ≈ 0.0
@test pvalue(t; tail=:left) ≈ 0.1455862977602108
@test pvalue(t; tail=:right) ≈ 1.0
@test pvalue(t; method=:central, tail=:both) ≈ 0.29117259552042146
@test pvalue(t; method=:minlike, tail=:both) ≈ 0.2720685111989459
@test_ci_approx confint(t; tail=:left) (0.0, 1.8523336085460567)
@test_ci_approx confint(t; tail=:right) (0.0, Inf)
@test_ci_approx confint(t; method=:central) (0.0, 2.5484240386190433)
#@test_approx_eq [confint(t; method=:minlike)...] [0.0, 1.9338]
show(IOBuffer(), t)

t = HypothesisTests.FisherExactTest(5, 0, 12, 7)
@test t.ω ≈ Inf
@test pvalue(t; tail=:left) ≈ 1.0
@test pvalue(t; tail=:right) ≈ 0.1455862977602108
@test pvalue(t; method=:central, tail=:both) ≈ 0.29117259552042146
@test pvalue(t; method=:minlike, tail=:both) ≈ 0.2720685111989459
@test_ci_approx confint(t; tail=:left) (0.0, Inf)
@test_ci_approx confint(t; tail=:right) (0.5398595562842079, Inf)
@test_ci_approx confint(t; method=:central) (0.39239937500428096, Inf)
#@test_approx_eq [confint(t; method=:minlike)...] [0.5171, Inf]
show(IOBuffer(), t)

# http://www.stata.com/support/faqs/statistics/fishers-exact-test/
t = HypothesisTests.FisherExactTest(2, 31, 136, 15532)
@test t.ω ≈ 7.3653226779913386
@test pvalue(t; tail=:left) ≈ 0.9970112864705307
@test pvalue(t; tail=:right) ≈ 0.03390271476034175
@test pvalue(t; method=:central, tail=:both) ≈ 0.06780542952068347
@test pvalue(t; method=:minlike, tail=:both) ≈ 0.03390271476034175
@test_ci_approx confint(t; tail=:left) (0.0,25.20252454804777)
@test_ci_approx confint(t; tail=:right) (1.2436262312601785, Inf)
@test_ci_approx confint(t; method=:central) (0.8458141614657836, 29.44434308524672)
#@test_approx_eq [confint(t; method=:minlike)...] [1.2436, 28.3557]
show(IOBuffer(), t)

# http://www.utstat.toronto.edu/~brunner/oldclass/312f12/lectures/312f12FisherWithR.pdf
# http://vassarstats.net/odds2x2.html
t = HypothesisTests.FisherExactTest(4, 1, 20, 1)
@test t.ω ≈ 0.21821789023599236
@test pvalue(t; tail=:left) ≈ 0.353846153846154
@test pvalue(t; tail=:right) ≈ 0.9692307692307692
@test pvalue(t; method=:central) ≈ 0.7076923076923074
@test pvalue(t; method=:minlike) ≈ 0.353846153846154
@test_ci_approx confint(t; tail=:left) (0.0, 9.594302003876502)
@test_ci_approx confint(t; tail=:right) (0.004963263361921223, Inf)
@test_ci_approx confint(t; method=:central) (0.002430190787475382, 19.59477744071154)
#@test_approx_eq [confint(t; method=:minlike)...] [0.005, 9.5943]
show(IOBuffer(), t)

# Corner cases gh #276
t = HypothesisTests.FisherExactTest(5, 0, 5, 0)
@test pvalue(t; tail=:left) ≈ 1
@test pvalue(t; tail=:right) ≈ 1
@test pvalue(t; method=:central) ≈ 1
@test pvalue(t; method=:minlike) ≈ 1
@test_ci_approx confint(t; tail=:left) (0.0, Inf)
@test_ci_approx confint(t; tail=:right) (0.0, Inf)
@test_ci_approx confint(t; method=:central) (0.0, Inf)

t = HypothesisTests.FisherExactTest(0, 5, 0, 5)
@test pvalue(t; tail=:left) ≈ 1
@test pvalue(t; tail=:right) ≈ 1
@test pvalue(t; method=:central) ≈ 1
@test pvalue(t; method=:minlike) ≈ 1
@test_ci_approx confint(t; tail=:left) (0.0, Inf)
@test_ci_approx confint(t; tail=:right) (0.0, Inf)
@test_ci_approx confint(t; method=:central) (0.0, Inf)

t = HypothesisTests.FisherExactTest(1, 1, 1, 1)
@test HypothesisTests.pvalue(t, tail=:both) <= 1
show(IOBuffer(), t)
end
