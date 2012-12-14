load("test")
load("src/HypothesisTests")
using HypothesisTests, Test

## ONE SAMPLE T-TEST

# One sample
@test p_value(OneSampleTTest, [-5:5]) == 1
@test abs(p_value(OneSampleTTest, [-5:10]) - 0.0530) <= 1e-4
@test abs(p_value(OneSampleTTest, [-10:5]) - 0.0530) <= 1e-4

# Paired samples
@test abs(p_value(OneSampleTTest, [1, 1, 2, 1, 0], [0, 1, 1, 1, 0]) - 0.1778) <= 1e-4

## TWO SAMPLE T-TESTS

# From http://en.wikipedia.org/w/index.php?title=Student%27s_t-test&oldid=526762741

a1 = [30.02, 29.99, 30.11, 29.97, 30.01, 29.99]
a2 = [29.89, 29.93, 29.72, 29.98, 30.02, 29.98]

@test df(EqualVarianceTTest, a1, a2) == 10
@test abs(test_statistic(EqualVarianceTTest, a1, a2) - 1.959) <= 1e-3
@test abs(p_value(EqualVarianceTTest, a1, a2) - 0.078) <= 1e-3
test = EqualVarianceTTest(a1, a2)
@test abs(test.sd - 0.084) <= 1e-3

@test abs(df(UnequalVarianceTTest, a1, a2) - 7.03) <= 0.01
@test abs(test_statistic(UnequalVarianceTTest, a1, a2) - 1.959) <= 1e-3
@test abs(p_value(UnequalVarianceTTest, a1, a2) - 0.091) <= 1e-3
test = UnequalVarianceTTest(a1, a2)
@test abs(test.sd - 0.0485) <= 1e-3