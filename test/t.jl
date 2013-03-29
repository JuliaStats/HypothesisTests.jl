using HypothesisTests, Test

## ONE SAMPLE T-TEST

# One sample
@test pvalue(OneSampleTTest([-5:5])) == 1
@test abs(pvalue(OneSampleTTest([-5:10])) - 0.0530) <= 1e-4
@test abs(pvalue(OneSampleTTest([-10:5])) - 0.0530) <= 1e-4

# Paired samples
@test abs(pvalue(OneSampleTTest([1, 1, 2, 1, 0], [0, 1, 1, 1, 0])) - 0.1778) <= 1e-4

## TWO SAMPLE T-TESTS

# From http://en.wikipedia.org/w/index.php?title=Student%27s_t-test&oldid=526762741

a1 = [30.02, 29.99, 30.11, 29.97, 30.01, 29.99]
a2 = [29.89, 29.93, 29.72, 29.98, 30.02, 29.98]

tst = EqualVarianceTTest(a1, a2)
@test tst.df == 10
@test abs(tst.t - 1.959) <= 1e-3
@test abs(pvalue(tst) - 0.078) <= 1e-3

tst = UnequalVarianceTTest(a1, a2)
@test abs(tst.df - 7.03) <= 0.01
@test abs(tst.t - 1.959) <= 1e-3
@test abs(pvalue(tst) - 0.091) <= 1e-3